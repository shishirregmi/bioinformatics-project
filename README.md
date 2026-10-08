# Cross-cancer radiotherapy signature transfer

An exploratory, reproducible analysis of whether a **fixed breast-cancer radiosensitivity signature** is associated with major pathologic response (MPR) in pretreatment non-small-cell lung cancer (NSCLC) samples from a neoadjuvant durvalumab ± stereotactic body radiotherapy (SBRT) trial.

The project is intended to make the data, scoring choices, and limits clear for review with a professor. It is not a clinical prediction tool.

## Research question

- **Primary:** Is the published 34-gene radiosensitivity signature (RSS) associated with MPR within the durvalumab + SBRT arm?
- **Secondary:** Is the published 4-gene immune signature associated with MPR in that same arm?
- **Exploratory only:** A treatment-by-score interaction may be considered after all sample labels are verified and an appropriate small-sample method is prespecified. No treatment-benefit claim is made by this project.

The signatures, source links, coefficients, direction, preprocessing adaptation, and lock date are recorded in [`config/rss.json`](config/rss.json) and [`config/immune.json`](config/immune.json). The source scoring process was developed for breast-cancer data. This project uses an explicitly documented RNA-seq adaptation and does not reuse the source cohort's cutoffs as NSCLC thresholds.

## Data

| Data | What it contains | Project use |
| --- | --- | --- |
| [GSE253564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253564) | Processed pretreatment tumor RNA-seq FPKM matrix; the paper reports 32 samples, 16 per arm | Primary signature scoring and response analysis after sample-to-patient labels are verified |
| [GSE248378](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248378) | Processed post-treatment resection RNA-seq FPKM matrix | Separate post-treatment expression exploration; do not treat as matched pairs until patient linkage is verified |
| [GEO sample records](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253564) | Sample identifiers and assay/provenance details | Check expression identifiers and specimen provenance |
| Clinical metadata CSV | Sample ID, patient ID, baseline status, treatment arm, MPR, source, verification flag | Select eligible baseline samples and define MPR groups; labels are never inferred from filenames or counts |
| Trial protocol | Combination arm used durvalumab with SBRT at 8 Gy × 3 daily fractions (24 Gy total) | Describe the study-level RT regimen; this is not a record of an individual's delivered dose |

The study paper reports **46 post-treatment tissues**, while GSE248378 has historically listed **29 GEO samples**. Reconcile the series and the paper's Table S1 before calling the deposited post-treatment matrix complete. This project does not contain patient-specific RT plans, RTDOSE/RTSTRUCT DICOM, dose-volume histograms, or organ-at-risk dosimetry.

### Source links

- [Direct GSE253564 pretreatment FPKM matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253564/suppl/GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz)
- [Direct GSE248378 post-treatment FPKM matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248378/suppl/GSE248378_Durva_Post_FPKMs.txt.gz)
- [Study article and supplementary files (including Table S1)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10982989/)
- [Cui et al. radiosensitivity signature source paper](https://doi.org/10.1158/1078-0432.CCR-18-0825)
- [Cui signature supplement](https://aacr.figshare.com/articles/journal_contribution/22470557)
- [Visualization and parameter guide](docs/visualization_poc.md)
- [Analysis contract](docs/analysis_plan.md)
- [Scoring details](docs/cui_scoring.md)

## Run the full workflow with Docker

From the repository root, build the image and run the container once:

```bash
docker build -t radio-transfer .
mkdir -p data results
docker run --rm -p 8000:8000 \
  -v "$PWD/data:/app/data" \
  -v "$PWD/results:/app/results" \
  radio-transfer
```

Open [http://localhost:8000](http://localhost:8000). On the first run, the container downloads both public GEO matrices, builds the data inventory and figures, runs the fixed-signature MPR analysis using the reviewed public mapping in `config/clinical.csv`, and starts the local web page. It saves downloads under `data/raw/` and reports under a timestamped folder in `results/`. Stop the server with Ctrl+C.

The bundled mapping is built from Table S1 and documented in [`docs/clinical_metadata.md`](docs/clinical_metadata.md). To supply a different or corrected mapping, place it at `data/clinical.csv` and run with `-e METADATA=/app/data/clinical.csv`.

### Run with Python instead

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .

radio-transfer download --dataset all --output data/raw
radio-transfer explore \
  --baseline-expression data/raw/GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz \
  --post-treatment-expression data/raw/GSE248378_Durva_Post_FPKMs.txt.gz \
  --signatures config/rss.json config/immune.json \
  --output results/exploration
```

To run the outcome analysis directly:

```bash
radio-transfer analyze \
  --expression data/raw/GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz \
  --metadata data/clinical.csv \
  --signatures config/rss.json config/immune.json \
  --output results/analysis
```

To generate figures for one matrix without the combined inventory, use `radio-transfer visualize --expression FILE --signatures config/rss.json config/immune.json --output results/figures`.

## Clinical metadata requirements

The repository includes [`config/clinical.csv`](config/clinical.csv), a verified minimal mapping for the 32 baseline expression profiles. It is crosswalked to the study's Table S1 and uses the published MPR definition; see [`docs/clinical_metadata.md`](docs/clinical_metadata.md) for the rules and validation counts. For a corrected or alternative source, start from [`config/metadata_template.csv`](config/metadata_template.csv) and record the evidence for each label in `label_source`.

Required fields:

| Field | Required value or meaning |
| --- | --- |
| `sample_id` | Exact column name from the pretreatment expression matrix |
| `patient_id` | De-identified participant ID used to verify one eligible baseline sample per patient |
| `baseline` | `1` for verified pretreatment samples, otherwise `0` |
| `arm` | `durvalumab` or `durvalumab_sbrt` |
| `mpr` | `1` for MPR, `0` for no MPR; leave unresolved values blank until verified |
| `pathology_response_signed_pct` | Optional source audit field; the bundled mapping retains the signed Pathology Response value from Table S1 |
| `label_source` | Paper/table/record that supports the sample and clinical labels |
| `verified` | `1` only after the sample-to-patient, baseline, arm, and MPR mapping has been checked |

All verified baseline records are reported in `sample_flow.csv`; unresolved, excluded, or expression-missing records remain visible with exclusion reasons. Duplicate eligible patient samples stop the analysis and require reconciliation. The paper's published group counts are checks on the crosswalk, **not** a way to assign labels to sample IDs.

## Score and statistical methods

1. Calculate each locked signature score for all pretreatment expression samples before joining MPR labels.
2. Apply `log2(FPKM + 1)`, then standardize each gene across the pretreatment samples using the sample standard deviation. The signatures' configured direction is applied so the exported score direction is consistent. This is an explicit platform adaptation; the original work used breast-cancer data and the supplement does not specify every scaling choice needed for this RNA-seq matrix.
3. Within the durvalumab + SBRT arm, compare MPR and no-MPR scores using a two-sided exact rank permutation test. Report group medians, a rank probability (MPR score greater than no-MPR score, with half weight for ties), and a descriptive stratified bootstrap interval. Holm adjustment covers the two prespecified signatures.
4. Report a signed within-sample rank adaptation and leave-one-out effects as sensitivity checks. These do not replace the locked weighted scores.

### First real-cohort run

The first source-checked run used 32 baseline profiles, including 10 MPR and 6 no-MPR samples in the durvalumab + SBRT arm. For the RSS, the rank probability was **0.15** (95% bootstrap interval 0.00–0.383), with median oriented scores of −3.61 for MPR and 8.75 for no MPR (exact *p* = 0.0225; Holm-adjusted *p* = 0.0450). Since larger oriented scores point in the published favorable direction, this is an inverse association in this RNA-seq adaptation. The secondary immune signature had rank probability **0.267** (95% interval 0.050–0.550; exact and Holm-adjusted *p* = 0.147), which does not show clear evidence of association.

The cohort is small. Bootstrap intervals are unstable, a non-significant result is not evidence of no biological relationship, and an association does not establish radiation-specific causation. No feature selection, classifier training, cutoff tuning, ROC analysis, predictive accuracy, or clinical validation is performed.

## What the reports contain

The combined explorer is available at `exploration/index.html`; the Docker landing page is at the run folder's `index.html`. When the clinical labels validate against the baseline samples, the analysis report is at `analysis/analysis_report.html`.

| Output | Contents |
| --- | --- |
| `data_overview.png` | Sample counts, gene-row counts, and per-sample median FPKM overview |
| `baseline/`, `post_treatment/` | Separate expression distributions, PCA, signature coverage, coefficients, heatmaps, and figure notes |
| `parameter_inventory.csv` | Meaning and intended use of each project parameter |
| `sample_parameter_profile.csv` | Per-sample mean/median/95th percentile FPKM and detected-gene counts |
| `gene_parameter_profile.csv` | Per-gene expression summaries for each dataset |
| `data_inventory.md`, `data_inventory.json` | Dataset descriptions, source links, counts, checksums, and limitations |
| `analysis_report.html`, `score_by_mpr.png` | Outcome score plot and accessible summary when analysis runs |
| `statistics.csv` | Medians, rank probability, interval, exact p-value, and Holm-adjusted p-value |
| `scores.csv`, `adapted_rank_scores.csv` | Per-sample locked scores and the separate rank sensitivity scores |
| `sample_flow.csv` | Metadata inclusion/exclusion audit |
| `sensitivity.csv` | Adapted-rank and leave-one-out sensitivity values |
| `gene_qc.csv`, `manifest.json` | Expression QC plus input hashes, versions, signature settings, and seed |

## Tests

Run the unit and synthetic end-to-end checks:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

The tests use synthetic cohort data and verify scoring, metadata exclusions, exact tests, and report generation. Passing tests validate software behavior; they do not validate the clinical sample mapping or the scientific hypothesis.

## Current analysis status

The pipeline, visual explorer, score calculations, source-based baseline clinical mapping, and synthetic tests are implemented. Each Docker run generates the real-cohort reports from the public expression matrices and the checked-in minimal clinical mapping. Review the input checksums, mapping source, counts, and limitations in each report before presenting outcome results.
