"""Profile both public expression matrices and create an interpretable data inventory."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .core import checksum, load_expression, load_signature
from .visualize import run as visualize


DATASETS = {
    "baseline": {
        "label": "Pretreatment baseline",
        "accession": "GSE253564",
        "reported_samples": 32,
        "timepoint": "Before neoadjuvant durvalumab, with or without SBRT",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253564/suppl/GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz",
    },
    "post_treatment": {
        "label": "Post-treatment resection",
        "accession": "GSE248378",
        "reported_samples": 46,
        "timepoint": "At surgery after neoadjuvant therapy",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248378/suppl/GSE248378_Durva_Post_FPKMs.txt.gz",
    },
}


PARAMETERS = [
    ("gene_id", "Expression matrix", "Gene identifier in the first column; the project expects exact symbols for signature matching.", "Map matrix rows to fixed signature genes; reconcile unmatched IDs explicitly."),
    ("sample_id", "Expression matrix / GEO", "Identifier for one RNA-seq profile; expression-matrix columns are samples.", "Link expression to GEO and verified clinical records; never infer outcomes from the ID."),
    ("FPKM", "Expression matrix", "Processed fragments-per-kilobase-per-million expression value for a gene in a sample.", "Describe expression distributions, PCA, gene coverage, and fixed signature scores."),
    ("RNA-seq processing", "GEO sample records", "Illumina NovaSeq 6000; reads aligned to hg19 with STAR 2.4.0f1; gene FPKMs estimated with Cufflinks 2.0.2 and GENCODE v19.", "Check assay and reference compatibility before comparing with other cohorts or remapping gene identifiers."),
    ("timepoint", "GEO series / clinical table", "Pretreatment biopsy or post-treatment resection, depending on series.", "Keep the two timepoints separate unless patient pairing is verified."),
    ("patient_id", "Verified clinical table", "De-identified study participant identifier.", "Establish patient-level linkage and detect repeated samples; do not publish direct identifiers."),
    ("baseline", "Verified clinical table", "Whether the profile is a pretreatment baseline sample (1/0).", "Restrict the prespecified baseline signature-transfer analysis."),
    ("arm", "Verified clinical table", "Durvalumab alone or durvalumab plus SBRT.", "Describe randomized treatment groups; an arm label is not a patient-specific delivered-dose record."),
    ("mpr", "Verified clinical table", "Major pathologic response indicator (1/0).", "Exploratory outcome grouping only after sample-to-patient linkage and source verification."),
    ("label_source", "Verified clinical table", "Source document or table supporting an individual clinical label.", "Audit clinical mapping and avoid guessed labels."),
    ("verified", "Verified clinical table", "Whether the clinical mapping has been checked (1/0).", "Exclude unresolved labels from outcome plots and analysis."),
    ("gene", "Locked signature configuration", "Gene included in the fixed RSS or immune signature.", "Check exact-symbol coverage and explain the score components."),
    ("weight", "Locked signature configuration", "Published coefficient attached to a signature gene.", "Compute the prespecified weighted score; not a patient measurement."),
    ("direction", "Locked signature configuration", "Orientation applied to the weighted sum (+1 or -1).", "Ensure larger exported scores point in the configured favorable direction."),
    ("transform", "Locked signature configuration", "Expression scale used by the score (identity, log2p1, or log2p1_gene_zscore).", "Make platform adaptation explicit and reproducible."),
    ("RT regimen", "Trial protocol / article", "Combination arm used 8 Gy per fraction for 3 daily fractions (24 Gy total), initiated with the first durvalumab cycle.", "Describe the trial's radiotherapy exposure at protocol level; do not treat it as individual dosimetry."),
]


def _profile_matrix(key, path, output, metadata_path, signatures):
    expression = load_expression(path)
    details = DATASETS[key]
    label = details["label"]
    child = output / key
    visualize_args = type("VisualizationArgs", (), {
        "expression": str(path),
        "metadata": str(metadata_path) if key == "baseline" and metadata_path else None,
        "signatures": [str(item) for item in signatures],
        "output": str(child),
        "dataset_label": label,
        "accession": details["accession"],
    })()
    visualize(visualize_args)

    values = expression.to_numpy(dtype=float)
    sample_rows = []
    for sample, column in expression.items():
        sample_values = column.to_numpy(dtype=float)
        sample_rows.append({
            "dataset": key,
            "accession": details["accession"],
            "sample_id": sample,
            "mean_fpkm": float(sample_values.mean()),
            "median_fpkm": float(np.median(sample_values)),
            "p95_fpkm": float(np.quantile(sample_values, 0.95)),
            "genes_nonzero": int(np.count_nonzero(sample_values)),
            "genes_nonzero_fraction": float(np.count_nonzero(sample_values) / len(sample_values)),
        })
    gene_rows = pd.DataFrame({
        "dataset": key,
        "accession": details["accession"],
        "gene_id": expression.index,
        "mean_fpkm": expression.mean(axis=1).to_numpy(),
        "median_fpkm": expression.median(axis=1).to_numpy(),
        "sd_fpkm": expression.std(axis=1, ddof=1).to_numpy(),
        "samples_nonzero": (expression > 0).sum(axis=1).to_numpy(),
        "samples_nonzero_fraction": (expression > 0).mean(axis=1).to_numpy(),
    })

    summary = {
        "dataset": key,
        "accession": details["accession"],
        "label": label,
        "timepoint": details["timepoint"],
        "expression_file": str(path),
        "sha256": checksum(path),
        "sample_count": int(expression.shape[1]),
        "gene_count": int(expression.shape[0]),
        "reported_study_samples": details["reported_samples"],
        "fpkm_min": float(values.min()),
        "fpkm_median": float(np.median(values)),
        "fpkm_max": float(values.max()),
        "zero_fraction": float(np.mean(values == 0)),
        "median_genes_detected_per_sample": float(np.median([row["genes_nonzero"] for row in sample_rows])),
    }
    return summary, pd.DataFrame(sample_rows), gene_rows


def _metadata_fields(path):
    if not path:
        return []
    data = pd.read_csv(path, dtype=str, keep_default_na=False)
    fields = []
    for column in data.columns:
        values = [str(value) for value in data[column] if str(value).strip()]
        fields.append({
            "dataset": "baseline_clinical_metadata",
            "field": column,
            "nonempty_rows": len(values),
            "missing_or_empty_rows": int(len(data) - len(values)),
            "unique_nonempty_values": len(set(values)),
            "project_use": next((use for name, _, _, use in PARAMETERS if name == column), "Inspect as a candidate study annotation; verify source and meaning before analysis."),
        })
    return fields


def _write_overview_figure(summaries, sample_tables, output):
    fig, axes = plt.subplots(1, 3, figsize=(14, 5.2))
    labels = [item["label"] for item in summaries]
    colors = ["#397c8c", "#d7824b"]

    axes[0].bar(labels, [item["sample_count"] for item in summaries], color=colors)
    axes[0].set_title("Expression profiles available")
    axes[0].set_ylabel("Samples")
    axes[0].tick_params(axis="x", rotation=18)
    for index, item in enumerate(summaries):
        axes[0].text(index, item["sample_count"], f"n={item['sample_count']}", ha="center", va="bottom")

    axes[1].bar(labels, [item["gene_count"] for item in summaries], color=colors)
    axes[1].set_title("Measured gene rows")
    axes[1].set_ylabel("Genes")
    axes[1].tick_params(axis="x", rotation=18)

    distributions = [table["median_fpkm"].to_numpy() for table in sample_tables]
    axes[2].boxplot(distributions, showfliers=False)
    axes[2].set_xticks([1, 2], labels)
    for index, values in enumerate(distributions, start=1):
        jitter = np.linspace(-0.07, 0.07, len(values)) if len(values) > 1 else np.zeros(len(values))
        axes[2].scatter(index + jitter, values, color=colors[index - 1], edgecolor="white", linewidth=0.5, s=35, zorder=3)
    axes[2].set_title("Per-sample median FPKM")
    axes[2].set_ylabel("Median over measured genes")
    axes[2].tick_params(axis="x", rotation=18)
    axes[2].text(0.02, -0.27, "QC summary only; not differential expression.", transform=axes[2].transAxes, fontsize=8, color="#667085")

    fig.suptitle("Public RNA-seq matrices: baseline and post-treatment", fontsize=14)
    fig.tight_layout()
    fig.savefig(output / "data_overview.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def _write_inventory_doc(output, summaries, metadata_fields):
    lines = [
        "# Project data inventory",
        "",
        "This report profiles the two public processed RNA-seq matrices separately. Baseline and post-treatment samples are not pooled or treated as matched pairs unless patient-level linkage has been verified.",
        "",
        "## Expression datasets",
        "",
        "| Dataset | Timepoint | Samples in matrix | Gene rows | Expression values |",
        "| --- | --- | ---: | ---: | --- |",
    ]
    for item in summaries:
        lines.append(f"| [{item['accession']}](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={item['accession']}) | {item['label']} | {item['sample_count']} | {item['gene_count']} | FPKM |")
    lines.extend([
        "",
        "The study paper reports 32 pretreatment profiles and 46 post-treatment tissues. The current GEO record for GSE248378 lists 29 samples; check Table S1 and the GEO sample list before claiming that the deposited post-treatment matrix contains every study tissue.",
        "",
        "## Parameters and how the project can use them",
        "",
        "| Parameter | Source | Meaning | Project use |",
        "| --- | --- | --- | --- |",
    ])
    for name, source, meaning, use in PARAMETERS:
        lines.append(f"| `{name}` | {source} | {meaning} | {use} |")
    lines.extend([
        "",
        "The matrix itself contains gene identifiers, sample identifiers, and FPKM values. Clinical fields are supplied separately through the verified metadata CSV template; GEO/clinical labels must be reconciled to expression sample IDs before outcome comparisons.",
        "",
        "## Radiotherapy data available here",
        "",
        "The public expression data support treatment-arm annotation and the protocol-level SBRT regimen: 8 Gy × 3 consecutive daily fractions (24 Gy total) with the first durvalumab cycle. They do not include individual RT treatment plans, RTDOSE/RTSTRUCT DICOM objects, dose-volume histograms, delivered-dose deviations, or target/organ-at-risk dosimetry in these FPKM matrices. Therefore the project can compare a fixed gene signature against a verified arm or response label, but cannot model per-patient dose-response from these matrices alone.",
        "",
        "## Useful source links",
        "",
        "- [GSE253564 baseline expression matrix and sample records](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253564)",
        "- [Direct baseline FPKM matrix](" + DATASETS["baseline"]["url"] + ")",
        "- [GSE248378 post-treatment expression matrix and sample records](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248378)",
        "- [Direct post-treatment FPKM matrix](" + DATASETS["post_treatment"]["url"] + ")",
        "- [Baseline raw reads: BioProject PRJNA1066291](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1066291)",
        "- [Post-treatment raw reads: BioProject PRJNA1043632](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1043632)",
        "- [Study paper and supplementary files, including Table S1](https://pmc.ncbi.nlm.nih.gov/articles/PMC10982989/)",
        "- [RSS source paper](https://doi.org/10.1158/1078-0432.CCR-18-0825); locked coefficients are recorded in `config/rss.json` and `config/immune.json`.",
        "",
        "## Generated files",
        "",
        "- `data_overview.png`: counts of samples and genes plus sample-level median-FPKM QC summaries.",
        "- `parameter_inventory.csv`: parameter definitions and project uses.",
        "- `signature_parameters.csv`: source, role, direction, transformation, lock date, and gene count for each signature.",
        "- `sample_parameter_profile.csv`: per-sample mean/median/95th percentile FPKM and number/fraction of detected genes.",
        "- `gene_parameter_profile.csv`: per-gene expression summaries across each dataset.",
        "- `clinical_metadata_fields.csv`: column names, completeness, and unique-value counts from the supplied clinical metadata, if provided. Values are not copied into this inventory.",
        "- `baseline/` and `post_treatment/`: individual dataset figures and run notes.",
        "",
        "The FPKM summaries are exploratory QC views. A between-timepoint difference is not a paired treatment effect and can reflect cohort composition, tissue collection, or processing differences.",
    ])
    if metadata_fields:
        lines.extend(["", "## Clinical metadata fields found in this run", ""])
        lines.append("| Field | Nonempty rows | Unique values |")
        lines.append("| --- | ---: | ---: |")
        for field in metadata_fields:
            lines.append(f"| `{field['field']}` | {field['nonempty_rows']} | {field['unique_nonempty_values']} |")
    (output / "data_inventory.md").write_text("\n".join(lines) + "\n")


def run(args):
    output = Path(args.output)
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output directory must be empty to preserve previous exploration reports.")
    output.mkdir(parents=True, exist_ok=True)
    signatures = [Path(path) for path in args.signatures]
    specs = [load_signature(path) for path in signatures]
    if len(specs) != 2 or {spec.get("role") for spec in specs} != {"primary", "secondary"}:
        raise ValueError("Provide exactly one primary RSS and one secondary immune specification.")

    summaries, sample_tables, gene_tables = [], [], []
    for key, path in (("baseline", args.baseline_expression), ("post_treatment", args.post_treatment_expression)):
        summary, sample_table, gene_table = _profile_matrix(key, Path(path), output, args.metadata, signatures)
        summaries.append(summary)
        sample_tables.append(sample_table)
        gene_tables.append(gene_table)

    pd.concat(sample_tables, ignore_index=True).to_csv(output / "sample_parameter_profile.csv", index=False)
    pd.concat(gene_tables, ignore_index=True).to_csv(output / "gene_parameter_profile.csv", index=False)
    pd.DataFrame(PARAMETERS, columns=["parameter", "source", "meaning", "project_use"]).to_csv(output / "parameter_inventory.csv", index=False)
    signature_rows = [{
        "name": spec["name"],
        "role": spec["role"],
        "gene_count": len(spec["genes"]),
        "source": spec["source"],
        "source_verified": spec["source_verified"],
        "locked_at": spec["locked_at"],
        "direction": spec["direction"],
        "transform": spec["transform"],
        "zscore_ddof": spec.get("zscore_ddof", ""),
    } for spec in specs]
    pd.DataFrame(signature_rows).to_csv(output / "signature_parameters.csv", index=False)
    metadata_fields = _metadata_fields(args.metadata)
    if metadata_fields:
        pd.DataFrame(metadata_fields).to_csv(output / "clinical_metadata_fields.csv", index=False)
    _write_overview_figure(summaries, sample_tables, output)
    _write_inventory_doc(output, summaries, metadata_fields)

    manifest = {
        "datasets": summaries,
        "clinical_metadata": ({"path": str(args.metadata), "sha256": checksum(args.metadata)} if args.metadata else None),
        "signature_files": [{"path": str(path), "sha256": checksum(path)} for path in signatures],
        "notes": [
            "Baseline and post-treatment datasets are summarized separately.",
            "GSE248378 currently lists fewer GEO samples than the study's reported post-treatment tissue count; reconcile before claiming completeness.",
            "RT regimen is protocol-level. No individual dose plan or dose-volume data are included in the FPKM matrices.",
        ],
    }
    (output / "data_inventory.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Saved data inventory and per-dataset figures to {output}")
