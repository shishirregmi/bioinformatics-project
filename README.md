# Cross-scale radiotherapy-signature transfer in lung cancer

This repository implements a reproducible master’s-project framework to test whether a **locked 34-gene breast-cancer radiosensitivity signature** transfers to lung cancer radiation response. It links a large preclinical radiation-survival panel to CCLE RNA-seq, then examines a separate, much smaller NSCLC patient cohort as a translational check.

The project does not train a classifier or claim clinical utility. It is designed to test one biological question across two settings: **does baseline expression captured by the published signature track intrinsic radiation survival in lung cancer cells, and is its direction compatible with pathologic response in a small clinical cohort?**

## Project aims

### Aim 1 — External preclinical radiation-response evaluation

Use the fixed RSS, without selecting genes or tuning a cutoff, to test its association with the continuous integral-survival phenotype in the published NSCLC cell-line panel. The paper reports 89 NSCLC and 39 lung adenocarcinoma (LUAD) lines among 533 profiled cancer cell lines. LUAD is a prespecified subgroup analysis. A histology-adjusted rank-correlation sensitivity analysis checks whether a pooled NSCLC result is explained by differences between cell-line subtypes. Exploratory gene-component associations show whether the fixed score’s individual weighted components point in a consistent direction.

The primary effect is Spearman’s ρ between the oriented RSS and integral radiation survival. The source endpoint ranges from 0 (completely sensitive) to 7 (completely resistant), so the directional hypothesis is a negative association. Two-sided permutation p-values and bootstrap intervals quantify uncertainty; Holm adjustment covers the NSCLC and LUAD association tests.

### Aim 2 — Clinical translation check

Use the pretreatment GSE253564 RNA-seq profiles and source-verified clinical mapping to compare the locked RSS with major pathologic response (MPR) in the durvalumab + SBRT group. The initial cohort contains 16 eligible patients in that arm (10 MPR, 6 without MPR). The published 4-gene immune score is secondary and is used only in the patient cohort. This analysis is exploratory because of its small sample size; it does not estimate a treatment interaction, causal radiation effect, or patient-level predictive accuracy.

The GSE248378 post-treatment matrix and complete data explorer are included for context and teaching. Post-treatment samples are not treated as baseline samples or matched pairs unless patient linkage is confirmed.

See [`docs/master_project_proposal.md`](docs/master_project_proposal.md) for the full research rationale, aims, analysis design, feasibility, milestones, and limitations. The statistical contract is in [`docs/analysis_plan.md`](docs/analysis_plan.md).
The [data and visualization guide](docs/visualization_poc.md) explains the expression and radiotherapy fields and how to read the generated views. The exact locked-score provenance and RNA-seq adaptation are described in [`docs/cui_scoring.md`](docs/cui_scoring.md).

## Public datasets and parameters

| Source | Parameters available | Use in this project |
| --- | --- | --- |
| [Yard et al. radiation-response study](https://doi.org/10.1038/ncomms11428) and [Supplementary Data 1 workbook](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fncomms11428/MediaObjects/41467_2016_BFncomms11428_MOESM708_ESM.xlsx) | Cell-line ID, name, cancer site, histology/subhistology, culture medium, SNP fingerprinting status, integral-survival AUC | Continuous in-vitro radiation-response endpoint; prespecified NSCLC and LUAD subsets |
| [CCLE RNA-seq RPKM matrix](https://data.broadinstitute.org/ccle/CCLE_RNAseq_genes_rpkm_20180929.gct.gz) | Ensembl `Name`, gene-symbol `Description`, RPKM per cell-line sample | Locked 34-gene RSS expression score, joined by exact normalized cell-line name and site |
| [GSE253564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253564) | Pretreatment tumor RNA-seq FPKM per sample and gene | Clinical translation analysis after exact sample-to-patient label verification |
| [GSE248378](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248378) | Post-treatment resection RNA-seq FPKM per sample and gene | Separate expression exploration; not combined with baseline outcomes |
| [Trial article and supplements](https://pmc.ncbi.nlm.nih.gov/articles/PMC10982989/) | Sample identifiers, treatment arm, pathology response and MPR labels; protocol regimen | Source crosswalk for the clinical translation analysis; trial SBRT was 8 Gy × 3 daily fractions |
| [`config/rss.json`](config/rss.json) and [`config/immune.json`](config/immune.json) | Published genes, coefficients, direction, transform, and source provenance | Locked primary RSS and secondary patient-only immune score |

**How the radiation endpoint works:** Yard et al. integrated survival measurements across 1, 2, 3, 4, 5, 6, 8, and 10 Gy and rescaled the area from 0 (completely sensitive) to 7 (completely resistant). It is a cell-culture phenotype, not a patient dose, tumor dose-volume histogram, or clinical response. The trial RNA-seq data do not include patient-specific RT plans, RTDOSE/RTSTRUCT DICOM, delivered dose, or organ-at-risk dosimetry.

The radiation workbook includes lung-site records beyond NSCLC, including small-cell lung cancer and unspecified histologies. The pipeline selects a documented set of explicit non-small-cell subhistologies, audits all source rows, and keeps only exact expression joins. It does not infer labels from sample counts or use fuzzy name matching. CCLE Ensembl rows are mapped through the GCT `Description` symbol field; repeated rows for a symbol are summed before log transformation. Lines with incomplete RSS expression are excluded without imputation.

### First source-checked readout

The downloaded radiation workbook contained 533 records, including 89 explicit NSCLC and 39 LUAD lines. The normalized exact CCLE join matched 518/533 panel records overall; 86 NSCLC and 37 LUAD records had both a matched expression profile and complete 34-gene RSS values. The NSCLC RSS–AUC correlation was **ρ = −0.027** (95% bootstrap CI −0.246 to 0.204; permutation *p* = 0.806; Holm *p* = 1.000). In LUAD it was **ρ = −0.024** (95% CI −0.378 to 0.316; *p* = 0.889; Holm *p* = 1.000). The histology-adjusted sensitivity was also near zero (ρ = −0.047; 95% CI −0.273 to 0.172; within-subhistology permutation *p* = 0.669).

This first real-data run shows no evidence that the fixed breast RSS transfers to the in-vitro lung radiation-survival phenotype under this RNA-seq adaptation. The small clinical pilot in the merged proof of concept showed an inverse RSS–MPR association in the combination arm (rank probability 0.15; exact *p* = 0.0225; Holm *p* = 0.045, with 10 MPR and 6 non-MPR patients). That estimate is too small and context-dependent to support prediction. The preclinical/clinical difference motivates a careful thesis on cross-cancer and cross-scale portability; it does not establish a mechanism or a clinical biomarker.

## Run everything in one Docker run

From the repository root:

```bash
docker build -t radio-transfer .
mkdir -p data results
docker run --rm -p 8000:8000 \
  -v "$PWD/data:/app/data" \
  -v "$PWD/results:/app/results" \
  radio-transfer
```

Open [http://localhost:8000](http://localhost:8000). The first run downloads both GEO matrices, the radiation-response supplement, and the CCLE RPKM GCT file (about 137 MB compressed). It then generates the data explorer, preclinical and clinical reports, and the project landing page. Downloaded files are cached in `data/raw/`; timestamped reports are written to `results/`. Stop the server with Ctrl+C.

To use a corrected clinical mapping, place a verified CSV at `data/clinical.csv` and run:

```bash
docker run --rm -p 8000:8000 \
  -v "$PWD/data:/app/data" \
  -v "$PWD/results:/app/results" \
  -e METADATA=/app/data/clinical.csv \
  radio-transfer
```

The default clinical mapping is crosswalked to the published trial table; its source and rules are in [`docs/clinical_metadata.md`](docs/clinical_metadata.md). The bundled labels are never inferred from expression sample names.

## Run selected steps with Python

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .

# Download all four public sources
radio-transfer download --dataset all --output data/raw

# Primary preclinical analysis
radio-transfer preclinical \
  --radiation-response data/raw/41467_2016_BFncomms11428_MOESM708_ESM.xlsx \
  --expression data/raw/CCLE_RNAseq_genes_rpkm_20180929.gct.gz \
  --signature config/rss.json \
  --output results/preclinical

# Clinical translation analysis
radio-transfer analyze \
  --expression data/raw/GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz \
  --metadata config/clinical.csv \
  --signatures config/rss.json config/immune.json \
  --output results/clinical

# Full baseline/post-treatment data explorer
radio-transfer explore \
  --baseline-expression data/raw/GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz \
  --post-treatment-expression data/raw/GSE248378_Durva_Post_FPKMs.txt.gz \
  --signatures config/rss.json config/immune.json \
  --output results/exploration
```

The default workflow is also configurable with `DATA_DIR`, `OUTPUT_DIR`, `PORT`, and `METADATA` environment variables.

## Outputs

| Output | Contents |
| --- | --- |
| `index.html` | Project landing page with links to all views |
| `exploration/index.html` | Baseline and post-treatment data inventory, parameter descriptions, plots, and downloads |
| `preclinical/preclinical_report.html` | NSCLC and LUAD score–survival associations, subtype adjustment, gene-component analysis, and limitations |
| `preclinical/radiation_data_overview.png` | AUC distribution and lung-site histology counts for the full irradiation panel |
| `preclinical/radiation_parameter_inventory.csv` | Meaning and intended use of each radiation or expression parameter |
| `preclinical/cell_line_join_audit.csv` | Every radiation-panel record with histology selection and CCLE join status |
| `preclinical/cell_line_scores.csv` | Matched cell-line subtype, AUC, CCLE ID, and fixed RSS |
| `preclinical/statistics.csv` | Correlations, permutation p-values, bootstrap intervals, and Holm adjustment |
| `preclinical/gene_component_associations.csv` | Exploratory per-gene component correlations and 34-gene BH q-values |
| `preclinical/manifest.json` | Input SHA-256 checksums, source definitions, transforms, and statistical settings |
| `analysis/analysis_report.html` | Clinical MPR comparisons and small-cohort interpretation limits |
| `analysis/sample_flow.csv` | Verified clinical inclusion and exclusion audit |

## Reproducibility and tests

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

The test suite checks the data crosswalk, transformation, histology filters, exact identifier matching, permutation/bootstrap calculations, report generation, and Docker dashboard using synthetic fixtures. Tests validate the software, not biological conclusions. The actual-data report records checksums and the numbers that pass each matching step.

## Interpretation boundaries

- The primary evidence is an association with a preclinical cell-culture radiation phenotype.
- The signature is fixed from a breast-cancer source; no genes, coefficients, or cutoffs are refit in NSCLC.
- A cell-line association is not a patient-level radiation-response estimate.
- MPR in the clinical cohort reflects the assigned regimen and tumor/immune context; it does not isolate radiation’s causal effect.
- No patient-specific dosimetry is present in these public RNA-seq matrices.
- No classifier, ROC analysis, tuned threshold, clinical utility claim, or treatment-benefit interaction is part of this project.
