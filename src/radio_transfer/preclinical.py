"""Cross-platform evaluation of the locked RSS in a public irradiated cell-line panel."""

from html import escape
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

from .core import checksum, holm_two, load_signature, score_expression


RADIATION_SOURCE = {
    "article": "https://doi.org/10.1038/ncomms11428",
    "supplement": "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fncomms11428/MediaObjects/41467_2016_BFncomms11428_MOESM708_ESM.xlsx",
    "accession": "Supplementary Data 1",
    "reported_total_cell_lines": 533,
    "reported_nsclc": 89,
    "reported_luad": 39,
}

CCLE_SOURCE = {
    "accession": "CCLE RNA-seq RPKM 2018-09-29",
    "url": "https://data.broadinstitute.org/ccle/CCLE_RNAseq_genes_rpkm_20180929.gct.gz",
    "filename": "CCLE_RNAseq_genes_rpkm_20180929.gct.gz",
}

# This explicit mapping reproduces the paper's 89-line NSCLC grouping from the
# supplement's Subhistology field while excluding SCLC and unspecified labels.
NSCLC_SUBHISTOLOGIES = {
    "adenocarcinoma",
    "squamous_cell_carcinoma",
    "non_small_cell_carcinoma",
    "large_cell_carcinoma",
    "bronchioloalveolar_adenocarcinoma",
    "mixed_adenosquamous_carcinoma",
    "mucoepidermoid_carcinoma",
}
LUAD_SUBHISTOLOGIES = {"adenocarcinoma", "bronchioloalveolar_adenocarcinoma"}

RADIATION_PARAMETERS = [
    ("Master_ccl_id", "Unique identifier assigned to the radiation-response record.", "Retain as a stable record key in the complete source audit."),
    ("Cell Line", "Cell-line name used in the radiation-survival experiment.", "Match to CCLE after normalization, together with Site."),
    ("Site", "Primary tissue site for the source cell line.", "Restrict the prespecified cohort to lung and support the exact CCLE join."),
    ("Histology", "Broad tumor histology label from the source panel.", "Audit histology inclusion and characterize the source cohort."),
    ("Subhistology", "Specific tumor subtype label from the source panel.", "Define NSCLC/LUAD groups and adjust the sensitivity analysis."),
    ("Culture_media", "Culture medium used for the cell line in the source experiment.", "Document experimental conditions; not a patient covariate."),
    ("Snp_fp_status", "SNP fingerprint comparison status for cell-line identity.", "Review identity provenance and flag lines without a matched reference."),
    ("AUC", "Integral survival after irradiation, scaled from 0 (sensitive) to 7 (resistant).", "Continuous primary cell-intrinsic radiation-response endpoint."),
    ("CCLE Name", "Ensembl gene identifier in the CCLE GCT Name field.", "Preserve in the RSS gene-to-expression mapping audit."),
    ("CCLE Description", "Gene-symbol annotation in the CCLE GCT Description field.", "Join the locked RSS gene symbols to RPKM expression rows."),
    ("RPKM", "CCLE RNA-seq expression measurement per cell line.", "Transform as log2(RPKM + 1), standardize by gene, then apply the fixed RSS."),
    ("ccle_sample_id", "CCLE sample column matched to one source cell-line/site pair.", "Retain exact join decisions; exclude missing or ambiguous matches."),
]


def _compact_identifier(value):
    return re.sub(r"[^A-Z0-9]", "", str(value).strip().upper())


def load_radiation_response(path):
    """Read the publisher's supplementary sheet and apply the prespecified histology labels."""
    frame = pd.read_excel(path, sheet_name=0, header=1, engine="openpyxl")
    required = {"Master_ccl_id", "Cell Line", "Site", "Histology", "Subhistology", "AUC"}
    if not required.issubset(frame.columns):
        raise ValueError(f"Radiation response workbook requires columns: {sorted(required)}")
    frame = frame.dropna(subset=["Cell Line", "Site", "AUC"]).copy()
    frame["Cell Line"] = frame["Cell Line"].astype(str).str.strip()
    for field in ("Site", "Histology", "Subhistology"):
        frame[field] = frame[field].fillna("").astype(str).str.strip()
    frame["AUC"] = pd.to_numeric(frame["AUC"], errors="raise")
    if not np.isfinite(frame["AUC"].to_numpy()).all():
        raise ValueError("Radiation AUC values must be finite.")
    frame["ccle_join_key"] = frame["Cell Line"].map(_compact_identifier) + frame["Site"].map(_compact_identifier)
    frame["subhistology_key"] = frame["Subhistology"].str.casefold().str.strip()
    frame["is_nsclc"] = frame["Site"].str.casefold().eq("lung") & frame["subhistology_key"].isin(NSCLC_SUBHISTOLOGIES)
    frame["is_luad"] = frame["Site"].str.casefold().eq("lung") & frame["subhistology_key"].isin(LUAD_SUBHISTOLOGIES)
    frame["reported_endpoint"] = frame["AUC"]
    return frame


def match_ccle_columns(radiation, columns):
    """Create an exact normalized cell-line/site join and retain every join decision."""
    candidates = {}
    for column in columns:
        candidates.setdefault(_compact_identifier(column), []).append(column)
    matched = []
    statuses = []
    for key in radiation["ccle_join_key"]:
        options = candidates.get(key, [])
        if len(options) == 1:
            matched.append(options[0])
            statuses.append("matched_exact_normalized_name_and_site")
        elif len(options) > 1:
            matched.append("")
            statuses.append("ambiguous_ccle_identifier")
        else:
            matched.append("")
            statuses.append("no_ccle_column_match")
    result = radiation.copy()
    result["ccle_sample_id"] = matched
    result["expression_match_status"] = statuses
    duplicate_pairs = result.loc[result.ccle_sample_id.ne(""), "ccle_sample_id"].duplicated(keep=False)
    if duplicate_pairs.any():
        duplicate_ids = set(result.loc[result.ccle_sample_id.ne(""), "ccle_sample_id"][duplicate_pairs])
        mask = result.ccle_sample_id.isin(duplicate_ids)
        result.loc[mask, "ccle_sample_id"] = ""
        result.loc[mask, "expression_match_status"] = "ambiguous_multiple_radiation_records"
    return result


def load_ccle_signature(path, signature_genes, chunk_rows=4000):
    """Load only signature genes from the large CCLE GCT, using its Description symbol field."""
    header = pd.read_csv(path, sep="\t", compression="gzip", skiprows=2, nrows=0)
    columns = list(header.columns)
    if len(columns) < 3 or columns[:2] != ["Name", "Description"]:
        raise ValueError("CCLE GCT must have Name and Description columns after its two metadata rows.")
    sample_columns = columns[2:]
    selected = []
    for chunk in pd.read_csv(path, sep="\t", compression="gzip", skiprows=2,
                             usecols=["Name", "Description", *sample_columns], chunksize=chunk_rows,
                             low_memory=False):
        chunk["Description"] = chunk["Description"].fillna("").astype(str).str.strip()
        hit = chunk[chunk["Description"].isin(signature_genes)]
        if not hit.empty:
            selected.append(hit)
    if not selected:
        raise ValueError("No locked RSS genes were found in the CCLE Description field.")
    rows = pd.concat(selected, ignore_index=True)
    for sample in sample_columns:
        rows[sample] = pd.to_numeric(rows[sample], errors="coerce")
    # Sum transcript/gene rows that map to the same HGNC symbol before log transform.
    expression = rows.groupby("Description", sort=False)[sample_columns].sum(min_count=1)
    expression.index.name = "gene"
    missing = sorted(set(signature_genes) - set(expression.index))
    if missing:
        raise ValueError(f"CCLE expression is missing locked RSS symbols: {', '.join(missing)}")
    if (expression < 0).any().any():
        raise ValueError("Selected CCLE RPKM values must be nonnegative.")
    # Do not impute expression. A line with any missing RSS value cannot be scored.
    expression = expression.dropna(axis=1, how="any")

    mapping = rows.groupby("Description", sort=False)["Name"].apply(
        lambda series: ";".join(sorted(set(series.astype(str))))
    )
    mapping = pd.DataFrame({"gene": mapping.index, "ensembl_ids": mapping.values})
    return expression, mapping, sample_columns


def association_statistics(score, auc, seed=4370, permutations=10000, bootstraps=2000):
    score = np.asarray(score, dtype=float)
    auc = np.asarray(auc, dtype=float)
    if len(score) != len(auc) or len(score) < 3 or not np.isfinite(score).all() or not np.isfinite(auc).all():
        raise ValueError("Association analysis needs at least three complete, finite score/AUC pairs.")
    observed = float(spearmanr(score, auc).statistic)
    if not np.isfinite(observed):
        raise ValueError("RSS or AUC is constant in the selected cell-line subset.")
    rng = np.random.default_rng(seed)
    permuted = np.empty(permutations)
    for index in range(permutations):
        permuted[index] = spearmanr(score, rng.permutation(auc)).statistic
    p_value = float((1 + np.count_nonzero(np.abs(permuted) >= abs(observed))) / (permutations + 1))
    boot = []
    for _ in range(bootstraps):
        indices = rng.integers(0, len(score), size=len(score))
        value = spearmanr(score[indices], auc[indices]).statistic
        if np.isfinite(value):
            boot.append(value)
    low, high = np.quantile(boot, [0.025, 0.975]) if boot else (np.nan, np.nan)
    return {"n": len(score), "spearman_rho": observed, "bootstrap_ci_low": float(low),
            "bootstrap_ci_high": float(high), "p_permutation_two_sided": p_value,
            "permutations": permutations, "valid_bootstrap_replicates": len(boot)}


def _partial_rank_correlation(score, auc, histology):
    score_rank = rankdata(score)
    auc_rank = rankdata(auc)
    labels = np.asarray(histology, dtype=str)
    levels = sorted(set(labels))
    design = np.column_stack([np.ones(len(labels)), *[(labels == level).astype(float) for level in levels[1:]]])
    score_residual = score_rank - design @ np.linalg.lstsq(design, score_rank, rcond=None)[0]
    auc_residual = auc_rank - design @ np.linalg.lstsq(design, auc_rank, rcond=None)[0]
    if np.std(score_residual) == 0 or np.std(auc_residual) == 0:
        return float("nan")
    return float(np.corrcoef(score_residual, auc_residual)[0, 1])


def histology_adjusted_statistics(score, auc, histology, seed=4372, permutations=10000, bootstraps=2000):
    """Partial Spearman with permutation and bootstrap resampling within histology groups."""
    score = np.asarray(score, dtype=float)
    auc = np.asarray(auc, dtype=float)
    histology = np.asarray(histology, dtype=str)
    if len(score) != len(auc) or len(score) != len(histology) or len(score) < 10:
        raise ValueError("Histology-adjusted sensitivity analysis needs at least ten complete cell lines.")
    if not np.isfinite(score).all() or not np.isfinite(auc).all():
        raise ValueError("Histology-adjusted sensitivity analysis requires finite values.")
    observed = _partial_rank_correlation(score, auc, histology)
    if not np.isfinite(observed):
        raise ValueError("The histology-adjusted rank association is not estimable.")
    rng = np.random.default_rng(seed)
    auc_ranks = rankdata(auc)
    groups = [np.flatnonzero(histology == level) for level in sorted(set(histology))]
    null = np.empty(permutations)
    for index in range(permutations):
        permuted = auc_ranks.copy()
        for group in groups:
            permuted[group] = rng.permutation(permuted[group])
        null[index] = _partial_rank_correlation(score, permuted, histology)
    p_value = float((1 + np.count_nonzero(np.abs(null) >= abs(observed))) / (permutations + 1))
    boot = []
    for _ in range(bootstraps):
        indices = np.concatenate([rng.choice(group, len(group), replace=True) for group in groups])
        value = _partial_rank_correlation(score[indices], auc[indices], histology[indices])
        if np.isfinite(value):
            boot.append(value)
    low, high = np.quantile(boot, [0.025, 0.975]) if boot else (np.nan, np.nan)
    return {"n": len(score), "spearman_rho": observed, "bootstrap_ci_low": float(low),
            "bootstrap_ci_high": float(high), "p_permutation_two_sided": p_value,
            "permutations": permutations, "valid_bootstrap_replicates": len(boot),
            "adjustment": "partial Spearman; ranked score and AUC residualized against source subhistology"}


def gene_component_associations(expression, signature, samples, auc):
    """Explore which weighted RSS components align with AUC; q-values use BH over 34 genes."""
    genes = [entry["gene"] for entry in signature["genes"]]
    weights = np.asarray([entry["weight"] for entry in signature["genes"]], dtype=float)
    values = np.log2(expression.loc[genes, samples].to_numpy(dtype=float) + 1)
    sd = values.std(axis=1, ddof=signature.get("zscore_ddof", 1))
    if np.any(sd == 0):
        raise ValueError("A constant RSS gene cannot be standardized in the CCLE panel.")
    z = (values - values.mean(axis=1, keepdims=True)) / sd[:, None]
    components = (signature["direction"] * weights[:, None]) * z
    auc = np.asarray(auc, dtype=float)
    rows = []
    for index, entry in enumerate(signature["genes"]):
        rho, p_value = spearmanr(components[index], auc)
        rows.append({"gene": entry["gene"], "source_weight": entry["weight"],
                     "oriented_component_weight": signature["direction"] * entry["weight"],
                     "component_auc_spearman_rho": float(rho), "p_approximate": float(p_value)})
    order = np.argsort([row["p_approximate"] for row in rows])
    adjusted = np.empty(len(rows), dtype=float)
    running = 1.0
    for rank in range(len(rows) - 1, -1, -1):
        original = order[rank]
        running = min(running, rows[original]["p_approximate"] * len(rows) / (rank + 1))
        adjusted[original] = min(1.0, running)
    for row, q_value in zip(rows, adjusted):
        row["q_bh_34_genes"] = float(q_value)
    return pd.DataFrame(rows).sort_values("component_auc_spearman_rho", key=lambda series: series.abs(), ascending=False)


def _write_figure(data, results, output):
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.4), squeeze=False)
    axes = axes[0]
    panels = [("NSCLC", data[data.is_nsclc], results[0]), ("LUAD", data[data.is_luad], results[1])]
    for axis, (label, subset, result) in zip(axes, panels):
        colors = plt.get_cmap("tab10")
        for index, histology in enumerate(sorted(subset.Subhistology.unique())):
            group = subset[subset.Subhistology == histology]
            axis.scatter(group.rss_score, group.AUC, s=47, alpha=0.85,
                         color=colors(index % 10), edgecolor="white", linewidth=0.55,
                         label=histology.replace("_", " "))
        axis.set_title(label, loc="left", weight="bold", color="#173047")
        axis.set_xlabel("Locked RSS score (higher = predicted more radiosensitive)")
        axis.set_ylabel("Integral radiation survival (0 sensitive → 7 resistant)")
        axis.grid(color="#e6ecf0", linewidth=0.7)
        axis.spines[["top", "right"]].set_visible(False)
        axis.text(0.03, 0.97,
                  f"n = {result['n']}\nSpearman ρ = {result['spearman_rho']:.3f}\n"
                  f"95% bootstrap CI: {result['bootstrap_ci_low']:.3f} to {result['bootstrap_ci_high']:.3f}\n"
                  f"Permutation p = {result['p_permutation_two_sided']:.4g}\nHolm p = {result['p_holm']:.4g}",
                  transform=axis.transAxes, va="top", ha="left", fontsize=8.8,
                  bbox={"boxstyle": "round,pad=.45", "facecolor": "white", "edgecolor": "#dce5ec"})
        if subset.Subhistology.nunique() <= 7:
            axis.legend(frameon=False, fontsize=7.5, loc="best", title="Source subtype")
    fig.suptitle("Does a fixed breast RSS transfer to lung cancer radiation survival?", x=0.05, ha="left",
                 fontsize=14, color="#173047", weight="bold")
    fig.text(0.05, 0.005, "Each point is one cell line. Higher AUC means greater in-vitro resistance; no patient-level prediction is implied.",
             fontsize=8.8, color="#607186")
    fig.tight_layout(rect=(0, 0.05, 1, 0.91))
    fig.savefig(output / "rss_vs_radiation_survival.png", dpi=190, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _write_panel_overview(radiation, output):
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.2))
    axes[0].hist(radiation.AUC, bins=24, color="#95a9b7", alpha=0.75,
                 label=f"All cell lines (n={len(radiation)})")
    lung = radiation[radiation.Site.str.casefold().eq("lung")]
    nsclc = radiation[radiation.is_nsclc]
    luad = radiation[radiation.is_luad]
    axes[0].hist(nsclc.AUC, bins=18, color="#187d78", alpha=0.62,
                 label=f"NSCLC source set (n={len(nsclc)})")
    axes[0].hist(luad.AUC, bins=14, color="#d87846", alpha=0.58,
                 label=f"LUAD source set (n={len(luad)})")
    axes[0].set_title("Integral radiation survival", loc="left", weight="bold")
    axes[0].set_xlabel("AUC (0 sensitive → 7 resistant)")
    axes[0].set_ylabel("Cell lines")
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].grid(axis="y", color="#e6ecf0", linewidth=0.7)
    counts = lung.Subhistology.replace("", "unspecified").value_counts().sort_values()
    axes[1].barh([label.replace("_", " ") for label in counts.index], counts.values, color="#598293")
    axes[1].set_title(f"All lung-site source records (n={len(lung)})", loc="left", weight="bold")
    axes[1].set_xlabel("Cell lines")
    axes[1].grid(axis="x", color="#e6ecf0", linewidth=0.7)
    axes[1].spines[["top", "right"]].set_visible(False)
    fig.suptitle("Radiotherapy data inventory: phenotype distribution and lung histologies", x=0.05,
                 ha="left", fontsize=14, color="#173047", weight="bold")
    fig.text(0.05, 0.005, "Lung site includes small-cell and unspecified histologies; the prespecified NSCLC set excludes them.",
             fontsize=8.7, color="#607186")
    fig.tight_layout(rect=(0, 0.05, 1, 0.91))
    fig.savefig(output / "radiation_data_overview.png", dpi=190, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _write_report(output, results, sensitivity, gene_results, flow, mapping, signature, source_counts):
    stats_rows = "".join(
        f"<tr><th>{escape(row['analysis'])}</th><td>{row['n']}</td><td>{row['spearman_rho']:.3f}</td>"
        f"<td>{row['bootstrap_ci_low']:.3f}–{row['bootstrap_ci_high']:.3f}</td>"
        f"<td>{row['p_permutation_two_sided']:.4g}</td><td>{row['p_holm']:.4g}</td></tr>"
        for row in results
    )
    flow_rows = "".join(f"<tr><td>{escape(row['stage'])}</td><td>{row['n']}</td><td>{escape(row['note'])}</td></tr>" for row in flow)
    mapping_rows = "".join(
        f"<tr><td>{escape(row.gene)}</td><td>{escape(row.ensembl_ids)}</td><td>{escape('present' if row.gene in set(mapping.gene) else 'missing')}</td></tr>"
        for row in mapping.itertuples()
    )
    n_nsclc_source = source_counts["nsclc_source"]
    n_luad_source = source_counts["luad_source"]
    genes_preview = gene_results.head(10)
    gene_rows = "".join(
        f"<tr><th>{escape(row.gene)}</th><td>{row.oriented_component_weight:.3g}</td>"
        f"<td>{row.component_auc_spearman_rho:.3f}</td><td>{row.q_bh_34_genes:.4g}</td></tr>"
        for row in genes_preview.itertuples()
    )
    html = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<meta name="description" content="External preclinical evaluation of a fixed breast-cancer radiation-response gene signature in NSCLC cell lines.">
<title>Preclinical radiation-response validation</title>
<style>
:root{{--ink:#173047;--muted:#607186;--paper:#f3f6f8;--line:#dce5ec;--teal:#187d78;--white:#fff}}
*{{box-sizing:border-box}}body{{margin:0;background:var(--paper);color:var(--ink);font:15px/1.55 system-ui,-apple-system,"Segoe UI",sans-serif}}
header{{padding:42px max(22px,calc((100vw - 1060px)/2));background:linear-gradient(120deg,#10263d,#1d5364);color:#fff}}header p{{color:#d8e8ef;margin:8px 0 0}}
main{{max-width:1060px;margin:26px auto;padding:0 20px}}section{{background:#fff;border:1px solid var(--line);border-radius:14px;padding:22px;margin-bottom:16px}}
h1{{font-size:clamp(26px,4vw,39px);line-height:1.1;margin:7px 0}}h2{{font-size:20px;margin:0 0 12px}}.eyebrow{{text-transform:uppercase;letter-spacing:.12em;font-size:11px;font-weight:750;color:#8ed0c9;margin:0}}
.muted{{color:var(--muted)}}.callout{{background:#fff7e8;border-left:4px solid #d78b45;padding:13px 16px;border-radius:5px}}
img{{max-width:100%;height:auto;border:1px solid var(--line);border-radius:10px}}.table-wrap{{overflow:auto}}table{{border-collapse:collapse;width:100%;min-width:700px}}th,td{{text-align:left;border-bottom:1px solid var(--line);padding:10px;vertical-align:top}}
a{{color:var(--teal)}}footer{{color:var(--muted);padding:8px 0 30px}}
</style></head><body><header><p class="eyebrow">Aim 1 · external functional evaluation</p><h1>Locked RSS and radiation survival in lung cancer cell lines</h1>
<p>Published breast-cancer signature, evaluated without refitting against an independent in-vitro radiation survival phenotype.</p></header>
<main><section><h2>Study question</h2><p>Is the fixed 34-gene RSS associated with integral survival after irradiation in NSCLC cell lines, and does the association persist in lung adenocarcinoma?</p>
<p class="muted">The signature is locked at <code>{escape(signature['locked_at'])}</code>. Expression was log2(RPKM + 1), then standardized per gene across the measured CCLE cell-line panel; the source weights and configured direction were applied unchanged.</p></section>
<section><h2>Radiotherapy data and available parameters</h2><img src="radiation_data_overview.png" alt="Integral radiation survival distribution and histology counts for the source cell-line data">
<p>The workbook contains cell-line IDs, tissue site, histology/subhistology, culture medium, SNP fingerprinting status, and integral-survival AUC. The expression file adds Ensembl IDs, gene-symbol descriptions, and RPKM. <a href="radiation_parameter_inventory.csv">Open the parameter dictionary</a> or <a href="cell_line_join_audit.csv">inspect all source rows and joins</a>.</p></section>
<section><h2>Association results</h2><img src="rss_vs_radiation_survival.png" alt="Locked radiosensitivity score plotted against radiation survival in NSCLC and LUAD cell lines">
<div class="table-wrap"><table><thead><tr><th>Prespecified group</th><th>Matched cell lines</th><th>Spearman ρ</th><th>95% bootstrap CI</th><th>Permutation p</th><th>Holm p</th></tr></thead><tbody>{stats_rows}</tbody></table></div>
<p>Higher AUC means greater integral survival and relative resistance (the publisher scales the endpoint from 0, completely sensitive, to 7, completely resistant). The directional hypothesis is therefore a negative RSS–AUC association. P-values come from two-sided 10,000-draw permutations; the bootstrap interval resamples cell lines. Holm adjustment covers the NSCLC and LUAD tests.</p></section>
<section><h2>Histology-adjusted sensitivity</h2><p>Partial Spearman association after residualizing ranked RSS and AUC against the source subhistology: <strong>ρ = {sensitivity['spearman_rho']:.3f}</strong> (95% stratified bootstrap CI {sensitivity['bootstrap_ci_low']:.3f}–{sensitivity['bootstrap_ci_high']:.3f}; within-subhistology permutation p = {sensitivity['p_permutation_two_sided']:.4g}). This checks whether the pooled NSCLC result is explained only by differences between subtypes. It is a sensitivity analysis, not an additional primary endpoint.</p></section>
<section><h2>RSS component exploration</h2><p>These 34 exploratory associations ask whether individual weighted score components align with the phenotype. Approximate Spearman p-values are Benjamini–Hochberg adjusted across all RSS genes; they do not replace the fixed-score primary test.</p>
<div class="table-wrap"><table><thead><tr><th>Gene</th><th>Oriented source weight</th><th>Component–AUC ρ</th><th>BH q</th></tr></thead><tbody>{gene_rows}</tbody></table></div><p><a href="gene_component_associations.csv">Download all component results</a></p></section>
<section><h2>Sample flow and identifier audit</h2><p>The source article reports {RADIATION_SOURCE['reported_total_cell_lines']} profiled lines, including {RADIATION_SOURCE['reported_nsclc']} NSCLC and {RADIATION_SOURCE['reported_luad']} LUAD. This implementation selected {n_nsclc_source} NSCLC-labelled records and {n_luad_source} LUAD-labelled records from the supplementary histology fields.</p>
<div class="table-wrap"><table><thead><tr><th>Stage</th><th>Records</th><th>Definition</th></tr></thead><tbody>{flow_rows}</tbody></table></div>
<p>Only normalized exact cell-line plus site matches are included. Unmatched and ambiguous identifiers remain in <a href="cell_line_join_audit.csv">the full join audit</a>; no expression value is imputed.</p></section>
<section><h2>Gene identifier mapping</h2><p>CCLE GCT <code>Description</code> symbols were joined to the locked RSS. If multiple Ensembl rows share a symbol, their RPKM values were summed before transformation.</p>
<div class="table-wrap"><table><thead><tr><th>HGNC symbol</th><th>Ensembl row(s)</th><th>Status</th></tr></thead><tbody>{mapping_rows}</tbody></table></div></section>
<section class="callout"><strong>Interpretation limits.</strong> This is an observational cell-line association. The phenotype is in-vitro integral survival from a multi-dose assay, not a clinical response, prescribed dose, or patient radiosensitivity measurement. The clinical MPR cohort is analyzed separately as a small translation aim; no causal treatment-benefit or clinical utility claim follows from either dataset.</section>
<section><h2>Reproducibility files</h2><ul><li><a href="statistics.csv">Statistics and uncertainty</a></li><li><a href="cell_line_scores.csv">Matched scores and radiation phenotypes</a></li><li><a href="cell_line_join_audit.csv">All source records and matching decisions</a></li><li><a href="gene_mapping.csv">Signature gene crosswalk</a></li><li><a href="cohort_flow.csv">Cohort flow</a></li><li><a href="manifest.json">Checksums, software versions, and settings</a></li></ul></section>
<footer>Exploratory academic analysis. Not a clinical diagnostic or radiotherapy decision tool.</footer></main></body></html>'''
    (output / "preclinical_report.html").write_text(html)


def run(args):
    output = Path(args.output)
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output directory must be empty to preserve previous runs.")
    signature = load_signature(args.signature)
    if signature.get("role") != "primary" or len(signature["genes"]) != 34:
        raise ValueError("The preclinical analysis requires the locked 34-gene primary RSS.")
    response = load_radiation_response(args.radiation_response)
    expression, mapping, ccle_columns = load_ccle_signature(
        args.expression, [gene["gene"] for gene in signature["genes"]]
    )
    response = match_ccle_columns(response, ccle_columns)
    exact_join = response.expression_match_status.eq("matched_exact_normalized_name_and_site")
    n_exact_joins = int(exact_join.sum())
    complete_expression = response.ccle_sample_id.isin(expression.columns)
    response.loc[exact_join & ~complete_expression, "expression_match_status"] = "incomplete_rss_expression"
    matched = response[response.expression_match_status == "matched_exact_normalized_name_and_site"].copy()
    eligible = matched[matched.is_nsclc | matched.is_luad].copy()
    score = score_expression(expression, signature)
    eligible["rss_score"] = eligible.ccle_sample_id.map(score)
    eligible = eligible.dropna(subset=["rss_score", "AUC"]).copy()
    if eligible.ccle_sample_id.duplicated().any():
        raise ValueError("A matched CCLE expression profile appears more than once in the radiation-response panel.")

    nsclc = eligible[eligible.is_nsclc]
    luad = eligible[eligible.is_luad]
    if len(nsclc) < 10 or len(luad) < 10:
        raise ValueError(f"Too few exact matches for analysis: NSCLC n={len(nsclc)}, LUAD n={len(luad)}.")
    raw_results = []
    for label, subset in (("NSCLC", nsclc), ("LUAD", luad)):
        raw_results.append({"analysis": label, **association_statistics(
            subset.rss_score.to_numpy(), subset.AUC.to_numpy(), seed=4370 if label == "NSCLC" else 4371
        )})
    adjusted = holm_two([result["p_permutation_two_sided"] for result in raw_results])
    for result, adjusted_p in zip(raw_results, adjusted):
        result["p_holm"] = adjusted_p
    sensitivity = {"analysis": "NSCLC subhistology-adjusted sensitivity", **histology_adjusted_statistics(
        nsclc.rss_score.to_numpy(), nsclc.AUC.to_numpy(), nsclc.Subhistology.to_numpy()
    )}
    gene_results = gene_component_associations(expression, signature, list(nsclc.ccle_sample_id), nsclc.AUC.to_numpy())

    output.mkdir(parents=True, exist_ok=True)
    flow = [
        {"stage": "Radiation workbook records", "n": len(response), "note": "Nonempty cell-line, site, and AUC rows in Supplementary Data 1."},
        {"stage": "Lung site", "n": int(response.Site.str.casefold().eq("lung").sum()), "note": "Site field equals lung."},
        {"stage": "NSCLC histologies", "n": int(response.is_nsclc.sum()), "note": "Explicit carcinoma subhistologies; small-cell and unspecified NS labels are excluded."},
        {"stage": "LUAD histologies", "n": int(response.is_luad.sum()), "note": "Adenocarcinoma plus bronchioloalveolar adenocarcinoma."},
        {"stage": "Exact normalized CCLE matches", "n": n_exact_joins, "note": "Cell line and tissue site both match after case/punctuation normalization."},
        {"stage": "Complete RSS expression profiles", "n": int((exact_join & complete_expression).sum()), "note": "All 34 required genes are present; incomplete cell lines are not imputed."},
        {"stage": "Matched NSCLC", "n": len(nsclc), "note": "Complete 34-gene RSS and radiation AUC."},
        {"stage": "Matched LUAD", "n": len(luad), "note": "Complete 34-gene RSS and radiation AUC."},
    ]
    response.rename(columns={"AUC": "radiation_auc"}).to_csv(output / "cell_line_join_audit.csv", index=False)
    eligible.to_csv(output / "cell_line_scores.csv", index=False)
    mapping.to_csv(output / "gene_mapping.csv", index=False)
    pd.DataFrame(flow).to_csv(output / "cohort_flow.csv", index=False)
    gene_results.to_csv(output / "gene_component_associations.csv", index=False)
    pd.DataFrame(RADIATION_PARAMETERS, columns=["parameter", "meaning", "project_use"]).to_csv(
        output / "radiation_parameter_inventory.csv", index=False
    )
    pd.DataFrame([*raw_results, sensitivity]).to_csv(output / "statistics.csv", index=False)
    _write_panel_overview(response, output)
    _write_figure(eligible, raw_results, output)
    _write_report(output, raw_results, sensitivity, gene_results, flow, mapping, signature,
                  {"nsclc_source": int(response.is_nsclc.sum()), "luad_source": int(response.is_luad.sum())})
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": [{"path": str(path), "sha256": checksum(path)} for path in (args.radiation_response, args.expression, args.signature)],
        "radiation_source": RADIATION_SOURCE,
        "expression_source": CCLE_SOURCE,
        "radiation_rows": len(response),
        "ccle_sample_columns": len(ccle_columns),
        "rss_gene_count": len(signature["genes"]),
        "signature": signature,
        "histology_definitions": {"NSCLC": sorted(NSCLC_SUBHISTOLOGIES), "LUAD": sorted(LUAD_SUBHISTOLOGIES)},
        "join": "Exact cell-line plus site after case and punctuation normalization; duplicate matches excluded; no imputation.",
        "transformation": "log2(RPKM + 1), then per-gene z-score across CCLE sample columns complete for all 34 RSS genes (sample SD); source weights/direction unchanged.",
        "statistics": {"effect": "Spearman rho", "p_value": "two-sided permutation, 10000 draws", "interval": "95% percentile bootstrap, 2000 cell-line resamples", "seed": 4370, "multiplicity": "Holm adjustment for NSCLC and LUAD tests"},
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Saved preclinical radiation-response analysis to {output}")

