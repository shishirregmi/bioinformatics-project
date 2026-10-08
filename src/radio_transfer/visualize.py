"""Create exploratory, publication-ready views of the locked transfer dataset."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from .core import checksum, load_expression, load_metadata, load_signature, score_expression


ARM_COLORS = {
    "durvalumab": "#5577aa",
    "durvalumab_sbrt": "#d97745",
    "unlabeled": "#9298a3",
}
ARM_LABELS = {
    "durvalumab": "Durvalumab",
    "durvalumab_sbrt": "Durvalumab + SBRT",
    "unlabeled": "No verified arm label",
}


def _load_specs(paths):
    specs = [load_signature(path) for path in paths]
    if len(specs) != 2 or {spec.get("role") for spec in specs} != {"primary", "secondary"}:
        raise ValueError("Provide exactly one primary RSS and one secondary immune specification.")
    if len({spec["name"] for spec in specs}) != 2:
        raise ValueError("Signature names must be unique.")
    primary = next(spec for spec in specs if spec["role"] == "primary")
    if len(primary["genes"]) != 34:
        raise ValueError("The primary published RSS requires exactly 34 genes.")
    return specs


def _sample_arm(sample_id, annotations):
    if annotations is None or sample_id not in annotations.index:
        return "unlabeled"
    return annotations.loc[sample_id, "arm"]


def _expression_distributions(expression, annotations, out, dataset_label):
    values = np.log2(expression.to_numpy(dtype=float) + 1)
    fig, ax = plt.subplots(figsize=(max(12, expression.shape[1] * 0.38), 6.5))
    boxes = ax.boxplot(
        [values[:, index] for index in range(values.shape[1])],
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#172033", "linewidth": 1.3},
        whiskerprops={"color": "#667085"},
        capprops={"color": "#667085"},
    )
    for index, box in enumerate(boxes["boxes"]):
        arm = _sample_arm(expression.columns[index], annotations)
        box.set_facecolor(ARM_COLORS.get(arm, ARM_COLORS["unlabeled"]))
        box.set_alpha(0.8)
    ax.set_xticks(np.arange(1, expression.shape[1] + 1))
    ax.set_xticklabels(expression.columns, rotation=90, fontsize=8)
    ax.set_ylabel("log2(FPKM + 1)")
    ax.set_xlabel(f"{dataset_label} tumor sample")
    ax.set_title(f"Expression distribution across {dataset_label.lower()} samples")
    ax.grid(axis="y", color="#d9dee8", linewidth=0.7, alpha=0.8)
    present_arms = {_sample_arm(sample, annotations) for sample in expression.columns}
    handles = [
        Line2D([0], [0], color=color, marker="s", linestyle="", markersize=9, label=label)
        for arm, label in ARM_LABELS.items()
        if arm in present_arms
        for color in [ARM_COLORS[arm]]
    ]
    if handles:
        ax.legend(handles=handles, frameon=False, loc="best")
    fig.tight_layout()
    fig.savefig(out / "sample_expression_distributions.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def _pca(expression, annotations, out, dataset_label):
    transformed = np.log2(expression.to_numpy(dtype=float) + 1)
    variances = transformed.var(axis=1, ddof=1)
    variable = np.flatnonzero(np.isfinite(variances) & (variances > 0))
    if len(variable) < 2 or expression.shape[1] < 3:
        raise ValueError("PCA needs at least two variable genes and three samples.")
    selected = variable[np.argsort(variances[variable])[-min(500, len(variable)):]]
    matrix = transformed[selected]
    centered = matrix - matrix.mean(axis=1, keepdims=True)
    standardized = centered / matrix.std(axis=1, ddof=1, keepdims=True)
    standardized = standardized.T
    left, singular_values, _ = np.linalg.svd(standardized, full_matrices=False)
    coordinates = left[:, :2] * singular_values[:2]
    proportions = singular_values**2 / np.sum(singular_values**2)

    fig, ax = plt.subplots(figsize=(9.5, 7))
    arms = [_sample_arm(sample, annotations) for sample in expression.columns]
    for arm in dict.fromkeys(arms):
        indices = [index for index, value in enumerate(arms) if value == arm]
        ax.scatter(
            coordinates[indices, 0], coordinates[indices, 1],
            s=58, color=ARM_COLORS.get(arm, ARM_COLORS["unlabeled"]),
            edgecolor="white", linewidth=0.7, label=ARM_LABELS.get(arm, arm), zorder=3,
        )
        for index in indices:
            ax.annotate(expression.columns[index], coordinates[index], xytext=(4, 4),
                        textcoords="offset points", fontsize=7, color="#344054")
    ax.axhline(0, color="#e5e7eb", linewidth=0.8, zorder=0)
    ax.axvline(0, color="#e5e7eb", linewidth=0.8, zorder=0)
    ax.set_xlabel(f"PC1 ({proportions[0] * 100:.1f}% of selected-gene variance)")
    ax.set_ylabel(f"PC2 ({proportions[1] * 100:.1f}% of selected-gene variance)")
    ax.set_title(f"Unsupervised PCA of {dataset_label.lower()} tumor expression")
    ax.legend(frameon=False, loc="best")
    ax.text(
        0.01, -0.15,
        f"Top {len(selected)} variable genes; log2(FPKM + 1), then per-gene z-score. "
        "PCA is descriptive and does not use response labels.",
        transform=ax.transAxes, fontsize=8, color="#667085",
    )
    fig.tight_layout()
    fig.savefig(out / "expression_pca.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    return {"genes_used": int(len(selected)), "pc1_variance_fraction": float(proportions[0]),
            "pc2_variance_fraction": float(proportions[1])}


def _signature_coverage(expression, specs, out):
    available = set(map(str, expression.index))
    rows = []
    for spec in specs:
        expected = [gene["gene"] for gene in spec["genes"]]
        missing = [gene for gene in expected if gene not in available]
        rows.append({
            "signature": spec["name"], "role": spec["role"],
            "genes_expected": len(expected), "genes_matched_exact_symbol": len(expected) - len(missing),
            "genes_missing_exact_symbol": len(missing), "missing_symbols": ";".join(missing),
        })
    report = pd.DataFrame(rows)
    report.to_csv(out / "signature_gene_coverage.csv", index=False)

    fig, ax = plt.subplots(figsize=(9, 3.8))
    y_positions = np.arange(len(report))
    for index, row in report.iterrows():
        ax.barh(index, row.genes_matched_exact_symbol, color="#327a73", height=0.55)
        ax.barh(index, row.genes_missing_exact_symbol, left=row.genes_matched_exact_symbol,
                color="#e7a35a", height=0.55)
        ax.text(row.genes_expected + 0.25, index,
                f"{row.genes_matched_exact_symbol} / {row.genes_expected} symbols found",
                va="center", fontsize=9, color="#344054")
    ax.set_yticks(y_positions)
    ax.set_yticklabels([f"{row.signature} ({row.role})" for row in report.itertuples()])
    ax.invert_yaxis()
    ax.set_xlim(0, max(report.genes_expected) * 1.55)
    ax.set_xlabel("Signature genes")
    ax.set_title("Exact gene-symbol coverage in the downloaded expression matrix")
    ax.text(0, -0.42,
            "No aliases or Entrez-ID substitutions are guessed; see signature_gene_coverage.csv for unmatched symbols.",
            transform=ax.transAxes, fontsize=8, color="#667085")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(axis="x", color="#e5e7eb", linewidth=0.7)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(out / "signature_gene_coverage.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    return report


def _signature_weights(specs, out):
    fig, axes = plt.subplots(2, 1, figsize=(12, 11), squeeze=False)
    for ax, spec in zip(axes[:, 0], sorted(specs, key=lambda item: item["role"])):
        genes = spec["genes"]
        weights = np.asarray([gene["weight"] for gene in genes], dtype=float)
        order = np.argsort(weights)
        labels = [genes[index]["gene"] for index in order]
        values = weights[order]
        colors = ["#5079a5" if value >= 0 else "#d37a4a" for value in values]
        ax.barh(np.arange(len(values)), values, color=colors, height=0.72)
        ax.set_yticks(np.arange(len(values)))
        ax.set_yticklabels(labels, fontsize=7)
        ax.axvline(0, color="#344054", linewidth=0.8)
        ax.set_xlabel("Published raw coefficient")
        ax.set_title(f"{spec['name']}: {len(genes)} fixed genes ({spec['role']} signature)", loc="left")
        ax.grid(axis="x", color="#e5e7eb", linewidth=0.7)
        ax.set_axisbelow(True)
    fig.suptitle("What contributes to each locked signature score", y=1.01, fontsize=14)
    fig.text(0.12, 0.005,
             "Bars are source coefficients, not patient expression. The configured direction is −1, "
             "so exported scores reverse each raw weighted sum.", fontsize=9, color="#667085")
    fig.tight_layout()
    fig.savefig(out / "signature_coefficients.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def _signature_heatmap(expression, specs, annotations, out, dataset_label):
    genes = [gene["gene"] for spec in specs for gene in spec["genes"]]
    available = set(map(str, expression.index))
    missing = [gene for gene in genes if gene not in available]
    if missing:
        return False
    transformed = np.log2(expression.loc[genes].to_numpy(dtype=float) + 1)
    means = transformed.mean(axis=1, keepdims=True)
    standard_deviations = transformed.std(axis=1, ddof=1, keepdims=True)
    if np.any(standard_deviations == 0):
        return False
    standardized = (transformed - means) / standard_deviations
    order = list(range(len(expression.columns)))
    if annotations is not None:
        order.sort(key=lambda i: (
            annotations.loc[expression.columns[i], "arm"] if expression.columns[i] in annotations.index else "z_unlabeled",
            int(annotations.loc[expression.columns[i], "mpr"]) if expression.columns[i] in annotations.index else 9,
            expression.columns[i],
        ))
    fig, ax = plt.subplots(figsize=(max(11, expression.shape[1] * 0.38), 11.5))
    image = ax.imshow(standardized[:, order], aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5,
                      interpolation="nearest")
    ax.set_yticks(np.arange(len(genes)))
    ax.set_yticklabels(genes, fontsize=7)
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels([expression.columns[index] for index in order], rotation=90, fontsize=7)
    ax.set_xlabel(f"{dataset_label} tumor sample")
    ax.set_title(f"Expression of fixed signature genes in {dataset_label.lower()} samples")
    rss_count = len(specs[next(i for i, spec in enumerate(specs) if spec["role"] == "primary")]["genes"])
    ax.axhline(rss_count - 0.5, color="#172033", linewidth=1.5)
    colorbar = fig.colorbar(image, ax=ax, fraction=0.025, pad=0.02)
    colorbar.set_label("Per-gene z-score across all baseline samples")
    fig.text(0.12, 0.012,
        "Blue = lower and red = higher expression relative to that gene's dataset mean. "
        "This is an exploratory view, not an outcome-based clustering result.",
             fontsize=8, color="#667085")
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(out / "signature_expression_heatmap.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    return True


def _response_scores(expression, specs, coverage, included, out):
    if included is None:
        return None
    complete = set(coverage.loc[coverage.genes_missing_exact_symbol == 0, "signature"])
    combo = included[included.arm == "durvalumab_sbrt"]
    if combo.empty or combo.mpr.nunique() != 2:
        return None
    selected = [spec for spec in specs if spec["name"] in complete]
    if not selected:
        return None
    fig, axes = plt.subplots(1, len(selected), figsize=(6 * len(selected), 5.4), squeeze=False)
    y = combo.set_index("sample_id").mpr
    for ax, spec in zip(axes[0], selected):
        scores = score_expression(expression, spec)
        values_by_group = [scores.loc[y.index[y == value]].to_numpy() for value in (0, 1)]
        ax.boxplot(values_by_group, positions=[0, 1], widths=0.48, showfliers=False,
                   medianprops={"color": "#172033", "linewidth": 1.4},
                   boxprops={"color": "#667085"}, whiskerprops={"color": "#667085"},
                   capprops={"color": "#667085"})
        rng = np.random.default_rng(4370)
        for group, values in enumerate(values_by_group):
            jitter = rng.uniform(-0.08, 0.08, len(values))
            ax.scatter(group + jitter, values, color=["#5577aa", "#d97745"][group],
                       edgecolor="white", linewidth=0.6, s=52, zorder=3)
        ax.set_xticks([0, 1], [f"No MPR\n(n={sum(y == 0)})", f"MPR\n(n={sum(y == 1)})"])
        ax.set_ylabel("Configured oriented score")
        ax.set_title(spec["name"])
        ax.grid(axis="y", color="#e5e7eb", linewidth=0.7)
        ax.text(0.02, 0.98, "Higher score follows configured favorable direction",
                transform=ax.transAxes, va="top", fontsize=7, color="#667085")
    fig.suptitle("Exploratory signature scores in the verified durvalumab + SBRT group", y=1.02)
    fig.text(0.02, -0.02,
             "Descriptive sample-level comparison only. The cohort is small; this does not establish "
             "prediction, treatment benefit, or a clinical cutoff.", fontsize=8, color="#667085")
    fig.tight_layout()
    fig.savefig(out / "signature_scores_by_mpr.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    return {"combination_arm_samples": int(len(combo)), "mpr": int(sum(y == 1)),
            "no_mpr": int(sum(y == 0)), "sample_ids": list(y.index)}


def _write_notes(out, summary, coverage):
    lines = [
        "# Visualization proof of concept",
        "",
        f"This run contains **{summary['sample_count']} samples** and **{summary['measured_gene_count']} measured gene rows**.",
        f"The source is the {summary['dataset_label'].lower()} {summary['accession']} FPKM matrix. Figures are descriptive; they are not a clinical validation.",
        "",
        "## Figures",
        "",
        "- `sample_expression_distributions.png`: compares log2(FPKM+1) distributions sample by sample. It can expose unusual sample-wide expression patterns, but is not a complete RNA-seq quality-control assessment.",
        "- `expression_pca.png`: projects the 500 most variable measured genes after log transformation and per-gene standardization. It helps reveal broad sample structure or possible outliers; it does not use response labels.",
        "- `signature_gene_coverage.png` and `signature_gene_coverage.csv`: show whether the fixed source genes are present under exact gene-symbol matching. Missing genes need explicit identifier reconciliation before scoring.",
        "- `signature_coefficients.png`: displays the fixed source weights. It explains how a score is assembled; these coefficients are not expression values.",
    ]
    if summary["heatmap_generated"]:
        lines.append("- `signature_expression_heatmap.png`: displays per-gene standardized expression for the locked signature genes across the whole baseline cohort.")
    else:
        lines.append("- Signature heatmap was skipped because at least one required gene symbol was unmatched or constant. Do not silently drop unmatched genes.")
    if summary["response_plot_generated"]:
        lines.append("- `signature_scores_by_mpr.png`: shows oriented signature scores for verified MPR and non-MPR samples in the durvalumab + SBRT arm. Points are individual samples; no cutoff is applied.")
    else:
        lines.append("- No MPR comparison was made because complete verified clinical labels for both response groups were not available or signature genes were unmatched.")
    lines.extend([
        "",
        "## What this can support",
        "",
        "These plots can help explain the study design, inspect the shape of the expression matrix, check whether the published signature genes can be mapped to the target data, and decide whether the fixed scores merit a cautious exploratory comparison.",
        "",
        "The response comparison is a small-cohort association check. It cannot show that a signature predicts individual outcomes, proves benefit from adding radiation, or is ready for treatment decisions. The configured log2(FPKM+1) and per-gene z-score steps are an explicit RNA-seq adaptation of the source scoring method.",
        "",
        "## Gene coverage",
        "",
    ])
    for row in coverage.itertuples():
        lines.append(f"- **{row.signature} ({row.role}):** {row.genes_matched_exact_symbol}/{row.genes_expected} exact symbols matched.")
    if summary["clinical_labels_present"]:
        lines.extend(["", "Clinical labels were read only from the supplied metadata file after the pipeline's verification checks."])
    else:
        lines.extend(["", "No clinical outcome labels were supplied for this run. No response group was inferred from sample IDs or published counts."])
    (out / "visualization_notes.md").write_text("\n".join(lines) + "\n")


def run(args):
    out = Path(args.output)
    if out.exists() and any(out.iterdir()):
        raise ValueError("Output directory must be empty to preserve previous visualizations.")
    expression = load_expression(args.expression)
    specs = _load_specs(args.signatures)
    dataset_label = getattr(args, "dataset_label", "Pretreatment baseline")
    accession = getattr(args, "accession", "GSE253564")
    annotations = None
    included = None
    if args.metadata:
        included, flow = load_metadata(args.metadata, set(expression.columns))
        annotations = included.set_index("sample_id")
    out.mkdir(parents=True, exist_ok=True)
    if args.metadata:
        flow.to_csv(out / "sample_flow.csv", index=False)

    _expression_distributions(expression, annotations, out, dataset_label)
    pca_summary = _pca(expression, annotations, out, dataset_label)
    coverage = _signature_coverage(expression, specs, out)
    _signature_weights(specs, out)
    heatmap_generated = _signature_heatmap(expression, specs, annotations, out, dataset_label)
    response_summary = _response_scores(expression, specs, coverage, included, out)
    summary = {
        "source": f"NCBI GEO {accession} {dataset_label.lower()} FPKM matrix",
        "accession": accession,
        "dataset_label": dataset_label,
        "expression_path": str(args.expression), "expression_sha256": checksum(args.expression),
        "clinical_metadata_path": str(args.metadata) if args.metadata else None,
        "clinical_metadata_sha256": checksum(args.metadata) if args.metadata else None,
        "sample_count": int(expression.shape[1]), "measured_gene_count": int(expression.shape[0]),
        "pca": pca_summary, "clinical_labels_present": annotations is not None,
        "verified_clinical_sample_count": int(len(included)) if included is not None else 0,
        "response_group_summary": response_summary,
        "heatmap_generated": bool(heatmap_generated),
        "response_plot_generated": response_summary is not None,
        "signature_coverage": coverage.to_dict(orient="records"),
        "figures": [
            "sample_expression_distributions.png", "expression_pca.png",
            "signature_gene_coverage.png", "signature_coefficients.png",
            *(["signature_expression_heatmap.png"] if heatmap_generated else []),
            *(["signature_scores_by_mpr.png"] if response_summary is not None else []),
        ],
    }
    (out / "visualization_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    _write_notes(out, summary, coverage)
    print(f"Saved exploratory figures and notes to {out}")

