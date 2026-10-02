import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
import urllib.request

import numpy as np
import pandas as pd
import scipy

from .core import checksum, holm_two, load_expression, load_metadata, load_signature, rank_effect, rank_test, score_expression

GEO_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253564/suppl/GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz"


def run(args):
    out = Path(args.output)
    if out.exists() and any(out.iterdir()):
        raise ValueError("Output directory must be empty to preserve previous runs.")
    expression = load_expression(args.expression)
    specs = [load_signature(path) for path in args.signatures]
    if len(specs) != 2 or {s.get("role") for s in specs} != {"primary", "secondary"}:
        raise ValueError("Provide exactly one primary RSS and one secondary immune specification.")
    if len({s["name"] for s in specs}) != 2:
        raise ValueError("Signature names must be unique.")
    primary = next(s for s in specs if s["role"] == "primary")
    if len(primary["genes"]) != 34:
        raise ValueError("The primary published RSS requires exactly 34 genes.")
    # Calculate both signatures before reading clinical responses.
    scores = pd.DataFrame({s["name"]: score_expression(expression, s) for s in specs})
    adapted = pd.DataFrame({s["name"]: score_expression(expression, s, adapted_rank=True) for s in specs})
    included, flow = load_metadata(args.metadata, set(expression.columns))
    combination = included[included.arm == "durvalumab_sbrt"]
    results, sensitivities = [], []
    for spec in specs:
        name = spec["name"]
        y = combination.set_index("sample_id").mpr
        values = scores.loc[y.index, name]
        result = rank_test(values[y == 1], values[y == 0])
        results.append({"signature": name, "role": spec["role"], **result})
        ranked = adapted.loc[y.index, name]
        sensitivities.append({"signature": name, "analysis": "adapted_rank", "omitted_sample": "",
                              "rank_probability": rank_effect(ranked[y == 1], ranked[y == 0])})
        for sample in y.index:
            keep = y.index != sample
            remainder_y, remainder = y[keep], values[keep]
            if remainder_y.nunique() == 2:
                sensitivities.append({"signature": name, "analysis": "leave_one_out", "omitted_sample": sample,
                                      "rank_probability": rank_effect(remainder[remainder_y == 1], remainder[remainder_y == 0])})
    adjusted = holm_two([r["p_exact_two_sided"] for r in results])
    for result, p in zip(results, adjusted): result["p_holm"] = p
    out.mkdir(parents=True, exist_ok=True)
    scores.rename_axis("sample_id").to_csv(out / "scores.csv")
    adapted.rename_axis("sample_id").to_csv(out / "adapted_rank_scores.csv")
    flow.to_csv(out / "sample_flow.csv", index=False)
    pd.DataFrame(results).to_csv(out / "statistics.csv", index=False)
    pd.DataFrame(sensitivities).to_csv(out / "sensitivity.csv", index=False)
    qc = pd.DataFrame({"gene": expression.index, "mean_fpkm": expression.mean(axis=1),
                       "constant": expression.nunique(axis=1) == 1})
    qc.to_csv(out / "gene_qc.csv", index=False)
    inputs = [args.expression, args.metadata, *args.signatures]
    manifest = {"created_utc": datetime.now(timezone.utc).isoformat(),
                "inputs": [{"path": str(path), "sha256": checksum(path)} for path in inputs],
                "versions": {"python": platform.python_version(), "numpy": np.__version__,
                             "pandas": pd.__version__, "scipy": scipy.__version__},
                "seed": 4370, "bootstrap_replicates": 2000,
                "bootstrap_method": "stratified percentile; descriptive uncertainty in a small cohort",
                "expression_samples": expression.shape[1], "eligible_samples": len(included),
                "signature_specs": specs,
                "interaction": "not implemented in initial milestone; no treatment-benefit claim"}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Saved analysis to {out}")


def download(args):
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    target = out / GEO_URL.rsplit("/", 1)[1]
    if target.exists():
        raise ValueError("Download already exists; preserve the archived input.")
    temporary = target.with_suffix(".partial")
    try:
        with urllib.request.urlopen(GEO_URL, timeout=60) as response, temporary.open("wb") as stream:
            while block := response.read(1024 * 1024): stream.write(block)
        load_expression(temporary, compression="gzip")
        temporary.rename(target)
    finally:
        temporary.unlink(missing_ok=True)
    (out / "download_manifest.json").write_text(json.dumps({"url": GEO_URL,
        "downloaded_utc": datetime.now(timezone.utc).isoformat(), "sha256": checksum(target)}, indent=2) + "\n")
    print(target)


def main():
    parser = argparse.ArgumentParser(description="Fixed radiotherapy-signature transfer analysis")
    sub = parser.add_subparsers(dest="command", required=True)
    fetch = sub.add_parser("download", help="Archive the GSE253564 processed FPKM matrix")
    fetch.add_argument("--output", default="data/raw")
    fetch.set_defaults(func=download)
    analysis = sub.add_parser("analyze")
    analysis.add_argument("--expression", required=True)
    analysis.add_argument("--metadata", required=True)
    analysis.add_argument("--signatures", nargs=2, required=True)
    analysis.add_argument("--output", default="results/run")
    analysis.set_defaults(func=run)
    args = parser.parse_args()
    try:
        args.func(args)
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(2, f"Error: {exc}\n")


if __name__ == "__main__":
    main()
