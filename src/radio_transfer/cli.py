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
from .explore import run as explore
from .preclinical import CCLE_SOURCE, RADIATION_SOURCE, load_radiation_response, run as preclinical
from .report import write_analysis_report
from .visualize import run as visualize

GEO_DATASETS = {
    "baseline": {
        "accession": "GSE253564",
        "filename": "GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253564/suppl/GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz",
    },
    "post-treatment": {
        "accession": "GSE248378",
        "filename": "GSE248378_Durva_Post_FPKMs.txt.gz",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248378/suppl/GSE248378_Durva_Post_FPKMs.txt.gz",
    },
}

ADDITIONAL_DATASETS = {
    "radiation-response": {
        "accession": RADIATION_SOURCE["accession"],
        "filename": "41467_2016_BFncomms11428_MOESM708_ESM.xlsx",
        "url": RADIATION_SOURCE["supplement"],
        "format": "excel",
    },
    "ccle-expression": {
        "accession": CCLE_SOURCE["accession"],
        "filename": CCLE_SOURCE["filename"],
        "url": CCLE_SOURCE["url"],
        "format": "gct",
    },
}

DOWNLOAD_DATASETS = {**GEO_DATASETS, **ADDITIONAL_DATASETS}


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
    write_analysis_report(out, scores, included, results)
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
    selected = list(DOWNLOAD_DATASETS) if args.dataset == "all" else [args.dataset]
    targets = [out / DOWNLOAD_DATASETS[name]["filename"] for name in selected]
    if any(target.exists() for target in targets):
        raise ValueError("At least one requested download already exists; preserve the archived input.")

    downloaded = []
    for name, target in zip(selected, targets):
        details = DOWNLOAD_DATASETS[name]
        temporary = target.with_suffix(".partial")
        try:
            with urllib.request.urlopen(details["url"], timeout=60) as response, temporary.open("wb") as stream:
                while block := response.read(1024 * 1024):
                    stream.write(block)
            if name in GEO_DATASETS:
                load_expression(temporary, compression="gzip")
            elif details["format"] == "excel":
                load_radiation_response(temporary)
            elif details["format"] == "gct":
                header = pd.read_csv(temporary, sep="\t", compression="gzip", skiprows=2, nrows=0)
                if list(header.columns[:2]) != ["Name", "Description"] or len(header.columns) < 4:
                    raise ValueError("Downloaded CCLE file is not the expected GCT gene-expression matrix.")
            temporary.rename(target)
        finally:
            temporary.unlink(missing_ok=True)
        downloaded.append({"dataset": name, "accession": details["accession"], "filename": target.name,
                           "url": details["url"], "sha256": checksum(target)})

    manifest_path = out / "download_manifest.json"
    previous = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    previous_records = previous.get("downloads", [])
    if not previous_records and previous.get("url"):
        legacy = next((name for name, record in DOWNLOAD_DATASETS.items() if record["url"] == previous["url"]), None)
        if legacy:
            previous_records = [{"dataset": legacy, "accession": DOWNLOAD_DATASETS[legacy]["accession"],
                                 "filename": DOWNLOAD_DATASETS[legacy]["filename"], "url": previous["url"],
                                 "sha256": previous.get("sha256", "")}]
    previous_records = [record for record in previous_records if record.get("dataset") not in selected]
    manifest = {"downloaded_utc": datetime.now(timezone.utc).isoformat(),
                "downloads": previous_records + downloaded}
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    for target in targets:
        print(target)


def main():
    parser = argparse.ArgumentParser(description="Fixed radiotherapy-signature transfer analysis")
    sub = parser.add_subparsers(dest="command", required=True)
    fetch = sub.add_parser("download", help="Archive the public expression and radiation-response datasets")
    fetch.add_argument("--output", default="data/raw")
    fetch.add_argument("--dataset", choices=[*DOWNLOAD_DATASETS, "all"], default="baseline",
                       help="Choose one source dataset or all GEO and preclinical files")
    fetch.set_defaults(func=download)
    analysis = sub.add_parser("analyze")
    analysis.add_argument("--expression", required=True)
    analysis.add_argument("--metadata", required=True)
    analysis.add_argument("--signatures", nargs=2, required=True)
    analysis.add_argument("--output", default="results/run")
    analysis.set_defaults(func=run)
    preclinical_analysis = sub.add_parser("preclinical", help="Test the locked RSS against cell-line radiation survival")
    preclinical_analysis.add_argument("--radiation-response", required=True, help="Yard et al. Supplementary Data 1 workbook")
    preclinical_analysis.add_argument("--expression", required=True, help="CCLE RPKM GCT.gz matrix")
    preclinical_analysis.add_argument("--signature", required=True, help="Locked RSS JSON specification")
    preclinical_analysis.add_argument("--output", default="results/preclinical")
    preclinical_analysis.set_defaults(func=preclinical)
    figures = sub.add_parser("visualize", help="Create exploratory plots and explanatory notes")
    figures.add_argument("--expression", required=True)
    figures.add_argument("--metadata", help="Optional verified clinical metadata CSV")
    figures.add_argument("--signatures", nargs=2, required=True)
    figures.add_argument("--output", default="results/visualization")
    figures.add_argument("--dataset-label", default="Pretreatment baseline")
    figures.add_argument("--accession", default="GSE253564")
    figures.set_defaults(func=visualize)

    exploration = sub.add_parser("explore", help="Inventory and visualize baseline and post-treatment data")
    exploration.add_argument("--baseline-expression", required=True)
    exploration.add_argument("--post-treatment-expression", required=True)
    exploration.add_argument("--metadata", help="Optional verified baseline clinical metadata CSV")
    exploration.add_argument("--signatures", nargs=2, required=True)
    exploration.add_argument("--output", default="results/data_exploration")
    exploration.set_defaults(func=explore)
    args = parser.parse_args()
    try:
        args.func(args)
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(2, f"Error: {exc}\n")


if __name__ == "__main__":
    main()
