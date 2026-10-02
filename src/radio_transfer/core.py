import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata


ARMS = {"durvalumab", "durvalumab_sbrt"}


def checksum(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_expression(path, compression="infer"):
    frame = pd.read_csv(path, sep="\t", index_col=0, compression=compression)
    frame.index = frame.index.astype(str)
    if frame.empty or frame.index.has_duplicates or frame.columns.has_duplicates:
        raise ValueError("Expression matrix must have unique gene and sample identifiers.")
    frame = frame.apply(pd.to_numeric, errors="raise")
    if not np.isfinite(frame.to_numpy()).all() or (frame < 0).any().any():
        raise ValueError("FPKM values must be finite, nonnegative and complete.")
    return frame


def load_signature(path):
    spec = json.loads(Path(path).read_text())
    for key in ("name", "source", "source_verified", "locked_at", "direction", "transform", "genes"):
        if key not in spec:
            raise ValueError(f"Signature missing {key}.")
    if spec["source_verified"] is not True or not spec["source"] or not spec["locked_at"]:
        raise ValueError("Signature source and scoring specification must be verified and locked before analysis.")
    if spec["direction"] not in (-1, 1) or spec["transform"] not in ("identity", "log2p1", "log2p1_gene_zscore"):
        raise ValueError("Unsupported score direction or expression transformation.")
    genes = spec["genes"]
    if not genes or len({g["gene"] for g in genes}) != len(genes):
        raise ValueError("Signature genes must be nonempty and unique.")
    if spec.get("zscore_ddof", 0) not in (0, 1):
        raise ValueError("zscore_ddof must be 0 or 1.")
    weights = np.array([g["weight"] for g in genes], dtype=float)
    if not np.isfinite(weights).all() or not np.any(weights):
        raise ValueError("Signature weights must be finite and not all zero.")
    return spec


def score_expression(expression, spec, adapted_rank=False):
    genes = [g["gene"] for g in spec["genes"]]
    missing = sorted(set(genes) - set(expression.index))
    if missing:
        raise ValueError(f"Missing required signature genes: {', '.join(missing)}")
    weights = np.array([g["weight"] for g in spec["genes"]], dtype=float)
    values = expression.loc[genes].to_numpy(dtype=float)
    if adapted_rank:
        # Rank over the full measured gene universe within each sample.
        universe = expression.rank(axis=0, method="average", pct=True)
        values = universe.loc[genes].to_numpy()
        weights = np.sign(weights) / np.abs(np.sign(weights)).sum()
    else:
        if spec["transform"].startswith("log2p1"):
            values = np.log2(values + 1)
        if spec["transform"] == "log2p1_gene_zscore":
            if values.shape[1] <= spec.get("zscore_ddof", 0):
                raise ValueError("Insufficient samples for the chosen z-score SD divisor.")
            sd = values.std(axis=1, ddof=spec.get("zscore_ddof", 0))
            if np.any(sd == 0):
                raise ValueError("Cannot standardize constant signature genes.")
            values = (values - values.mean(axis=1, keepdims=True)) / sd[:, None]
    return pd.Series(spec["direction"] * (weights @ values), index=expression.columns, name=spec["name"])


def load_metadata(path, sample_ids):
    data = pd.read_csv(path, dtype=str, keep_default_na=False)
    required = {"sample_id", "patient_id", "baseline", "arm", "mpr", "label_source", "verified"}
    if not required.issubset(data.columns):
        raise ValueError(f"Metadata requires columns: {sorted(required)}")
    if data.sample_id.duplicated().any() or (data.sample_id == "").any():
        raise ValueError("Metadata sample IDs must be unique and nonempty.")
    if not set(data.sample_id).issubset(sample_ids):
        raise ValueError("Metadata contains samples absent from expression matrix.")
    reasons = []
    for row in data.itertuples():
        problems = []
        if row.baseline != "1": problems.append("baseline_not_verified")
        if row.arm not in ARMS: problems.append("unresolved_arm")
        if row.mpr not in {"0", "1"}: problems.append("unresolved_mpr")
        if row.verified != "1" or not row.label_source: problems.append("unverified_label")
        if not row.patient_id: problems.append("unresolved_patient")
        reasons.append(";".join(problems))
    data["exclusion_reason"] = reasons
    absent = sorted(set(sample_ids) - set(data.sample_id))
    if absent:
        data = pd.concat([data, pd.DataFrame({"sample_id": absent, "exclusion_reason": "missing_metadata"})], ignore_index=True).fillna("")
    included = data[data.exclusion_reason == ""].copy()
    if included.patient_id.duplicated().any():
        raise ValueError("Multiple eligible samples per patient; reconcile before analysis.")
    included["mpr"] = included.mpr.astype(int)
    return included, data


def rank_effect(positive, negative):
    differences = np.asarray(positive)[:, None] - np.asarray(negative)[None, :]
    return float(np.mean((differences > 0) + 0.5 * (differences == 0)))


def rank_test(positive, negative, seed=4370, bootstraps=2000):
    positive, negative = np.asarray(positive, dtype=float), np.asarray(negative, dtype=float)
    if not len(positive) or not len(negative):
        raise ValueError("Both MPR and non-MPR groups are required.")
    if not np.isfinite(np.r_[positive, negative]).all():
        raise ValueError("Scores must be finite.")
    pooled = np.r_[positive, negative]
    n = len(positive)
    assignments = math.comb(len(pooled), n)
    if assignments > 200000:
        raise ValueError("Exact test exceeds 200,000 assignments; revise the prespecified plan.")
    ranks = rankdata(pooled)
    center = n * (len(pooled) + 1) / 2
    observed = abs(ranks[:n].sum() - center)
    extreme = sum(abs(ranks[list(indices)].sum() - center) >= observed - 1e-12
                  for indices in itertools.combinations(range(len(pooled)), n))
    rng = np.random.default_rng(seed)
    effects = [rank_effect(rng.choice(positive, n, replace=True),
                           rng.choice(negative, len(negative), replace=True))
               for _ in range(bootstraps)]
    low, high = np.quantile(effects, [0.025, 0.975])
    return {"n_mpr": n, "n_non_mpr": len(negative), "median_mpr": float(np.median(positive)),
            "median_non_mpr": float(np.median(negative)), "rank_probability": rank_effect(positive, negative),
            "bootstrap_ci_low": float(low), "bootstrap_ci_high": float(high),
            "p_exact_two_sided": extreme / assignments, "assignments": assignments}


def holm_two(p_values):
    if len(p_values) != 2:
        raise ValueError("The prespecified Holm family contains RSS and immune tests.")
    order = np.argsort(p_values)
    adjusted = np.empty(2)
    adjusted[order[0]] = min(1, 2 * p_values[order[0]])
    adjusted[order[1]] = max(adjusted[order[0]], p_values[order[1]])
    return adjusted.tolist()
