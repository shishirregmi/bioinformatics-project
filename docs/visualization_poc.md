# Data inventory and visualization proof of concept

## What data the project uses

The study has two processed RNA-seq FPKM matrices. Each matrix is organized with gene identifiers as rows and tumor sample IDs as columns. The baseline series contains 32 pretreatment samples. The post-treatment GEO series currently lists 29 samples. The paper describes 46 post-treatment tissue samples, so reconcile the deposited GEO sample set with the paper's Table S1 before calling it complete.

The project also uses two locked gene signatures: a 34-gene breast-cancer radiosensitivity signature (RSS) and a 4-gene immune signature. The JSON configurations contain gene symbols, published weights, score direction, expression transformation, source, and lock/verification fields. A minimal public clinical crosswalk for all 32 baseline profiles is included in `config/clinical.csv`, derived from the study's Table S1 and documented in `docs/clinical_metadata.md`. Labels are not inferred from expression values, sample names, or published cohort counts.

| Data block | Parameters available | How it can help this project |
| --- | --- | --- |
| Baseline expression, GSE253564 | Gene ID, sample ID, FPKM | Describe pretreatment NSCLC expression, check signature-gene coverage, calculate the fixed signature scores, and inspect broad sample structure. |
| Post-treatment expression, GSE248378 | Gene ID, sample ID, FPKM | Explore expression after neoadjuvant therapy and signature-gene expression at resection. Keep it separate from baseline unless patient pairing is verified. |
| GEO sample records | Sample title/source, tissue or sample characteristics, RNA/library/sequencing details, platform, and processing notes | Verify specimen/timepoint and technical provenance; these records describe the deposited RNA-seq samples, not the full patient-level clinical table. |
| Clinical metadata CSV | `sample_id`, `patient_id`, `baseline`, `arm`, `mpr`, `label_source`, `verified` | Verify sample-to-patient links, define treatment groups, and make a cautious exploratory MPR comparison. The project currently requires these fields for the outcome analysis. |
| Locked signature configurations | `gene`, `weight`, `direction`, `transform`, `source`, `source_verified`, `locked_at`, `role`, optional `zscore_ddof` | Document and reproduce the score definition; exact-symbol coverage shows if the target matrix contains all source genes. |
| Trial radiotherapy protocol | SBRT treatment arm; 8 Gy × 3 daily fractions (24 Gy total), delivered with the first durvalumab treatment | Describe the intended RT exposure for this trial and explain the combination-arm label. This is a study-level regimen, not an individual patient's delivered dose record. |

Both RNA-seq series report NovaSeq 6000 data processed against hg19 using STAR 2.4.0f1, Cufflinks 2.0.2 FPKM quantification, and GENCODE v19. Check that annotation and identifier background when reconciling signature genes or comparing another cohort.

## Visualize and profile both expression datasets

Download both processed matrices and run the combined explorer with the source-verified baseline mapping:

```bash
radio-transfer download --dataset all --output data/raw
radio-transfer explore \
  --baseline-expression data/raw/GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz \
  --post-treatment-expression data/raw/GSE248378_Durva_Post_FPKMs.txt.gz \
  --metadata config/clinical.csv \
  --signatures config/rss.json config/immune.json \
  --output results/data_exploration
```

To use a different or corrected clinical mapping, start from `config/metadata_template.csv` and provide the completed CSV:

```bash
radio-transfer explore \
  --baseline-expression data/raw/GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz \
  --post-treatment-expression data/raw/GSE248378_Durva_Post_FPKMs.txt.gz \
  --metadata data/clinical.csv \
  --signatures config/rss.json config/immune.json \
  --output results/data_exploration
```

A single matrix can still be visualized with `radio-transfer visualize --expression FILE --signatures config/rss.json config/immune.json`. It supports `--dataset-label` and `--accession` so figures identify their timepoint and GEO series.

## Launch the whole workflow with Docker

Build the image once, then run one container command to download both matrices if missing, generate all figures and inventories, run the MPR analysis with `config/clinical.csv`, and serve the web report:

```bash
docker build -t radio-transfer .
mkdir -p data results
docker run --rm -p 8000:8000 \
  -v "$PWD/data:/app/data" \
  -v "$PWD/results:/app/results" \
  radio-transfer
```

Open [http://localhost:8000](http://localhost:8000). The container prints the report folder, persists downloaded matrices under `data/`, and saves each run under a timestamped directory in `results/`. It uses the bundled mapping by default; override it with `-e METADATA=/app/data/clinical.csv` if needed. The landing page links to the explorer and statistical report. Press Ctrl+C to stop serving the page.

## Generated views and files

The combined command creates a run-level inventory, separate figures for each timepoint, and a response-analysis report. The `Entrez.ID` field in the GEO files is gene annotation, not a sample, and is removed before expression summaries. Baseline and post-treatment samples are not pooled or assumed to be matched.

| Output | What it represents | How it may help |
| --- | --- | --- |
| `data_overview.png` | Sample and measured-gene counts plus each sample's median FPKM | Gives a quick inventory of both matrices and a simple expression QC summary. It is not a differential-expression result. |
| `baseline/sample_expression_distributions.png`, `post_treatment/sample_expression_distributions.png` | Per-sample gene expression distributions after `log2(FPKM + 1)` | Highlights unusually broad or shifted distributions in either set. This is not a full RNA-seq QC report. |
| `baseline/expression_pca.png`, `post_treatment/expression_pca.png` | Separate PCA views of each timepoint's 500 most variable genes | Shows within-dataset structure or possible outliers. PCA does not use MPR labels. The two PCAs are not on a shared scale. |
| `*/signature_gene_coverage.png` and `*/signature_gene_coverage.csv` | Exact gene-symbol matches for the RSS and immune signature | Shows whether each fixed score can be calculated without dropping or guessing genes. |
| `*/signature_coefficients.png` | Source coefficient for each signature gene | Explains the fixed score definition. Bars are model weights, not patient expression. |
| `*/signature_expression_heatmap.png` | Per-gene standardized expression for signature genes | Shows patterns across each dataset. It is exploratory and not outcome-trained clustering. |
| `baseline/signature_scores_by_mpr.png` | Scores in verified MPR and non-MPR samples of the durvalumab + SBRT arm | Provides a small-cohort descriptive comparison only when complete verified labels and genes are present. |
| `sample_parameter_profile.csv` | Per sample: mean, median and 95th percentile FPKM, plus number/fraction of nonzero gene values | Helps inspect every sample's measured expression range and gene detection. |
| `gene_parameter_profile.csv` | Per gene and dataset: mean, median, SD, and nonzero sample count/fraction | Lets the project inspect all measured gene rows and identify variable or rarely detected genes. |
| `parameter_inventory.csv` | Field descriptions and project uses | Gives the professor a data dictionary of the expression, clinical, signature, and radiotherapy fields. |
| `clinical_metadata_fields.csv` | Field names, completeness, and unique-value counts from supplied metadata | Shows what clinical parameters were supplied without copying any values into the inventory. |
| `signature_parameters.csv`, `data_inventory.json`, `data_inventory.md` | Signature definitions, input hashes, matrix dimensions, interpretation notes, and source links | Makes the exploration reviewable and reproducible. |
| `index.html` | Browser page with the figures, data dictionary, and inventory file links | View all generated plots and parameter definitions in one place. |
| `analysis/analysis_report.html`, `analysis/statistics.csv` | Fixed-signature MPR score plot and exact-test results | Review exploratory association estimates and their limits. |

## What the radiotherapy information does and does not contain

The public RNA-seq matrices can be linked to randomized study arms and outcomes after clinical sample mapping. The trial paper describes the combination arm as durvalumab plus SBRT at 8 Gy for 3 daily fractions. These expression matrices do not contain patient-specific treatment plans, RTDOSE or RTSTRUCT DICOM files, dose-volume histograms, delivered-dose deviations, or target/organ-at-risk dosimetry. A patient-level dose-response analysis would require those separate RT records.

## Interpretation limits

The visualizations describe the data and test whether the published fixed signatures can be mapped to NSCLC expression. An MPR score plot is only a hypothesis-generating association in a small cohort; it cannot establish individual prediction, prove radiation caused the response, or guide treatment. The project uses log2(FPKM + 1) and per-gene z-scores as an explicit RNA-seq adaptation; it does not claim to exactly reproduce the breast-cancer training preprocessing or reuse breast-specific cutoffs.

The baseline study reports 16 pretreatment profiles per arm. Within the 16 baseline combination-arm profiles, the paper reports 10 MPR and 6 non-MPR. Do not assign these labels to sample IDs by count or filename; verify each mapping against the clinical source table.

## Links to the data and sources

- [GSE253564: baseline series, sample records, and processed data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253564)
- [Direct baseline FPKM matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253564/suppl/GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz)
- [GSE248378: post-treatment series, sample records, and processed data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248378)
- [Direct post-treatment FPKM matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE248nnn/GSE248378/suppl/GSE248378_Durva_Post_FPKMs.txt.gz)
- [Baseline raw sequencing reads: BioProject PRJNA1066291](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1066291)
- [Post-treatment raw sequencing reads: BioProject PRJNA1043632](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1043632)
- [Study article and supplementary files](https://pmc.ncbi.nlm.nih.gov/articles/PMC10982989/); [direct Table S1 workbook](https://pmc.ncbi.nlm.nih.gov/articles/instance/10982989/bin/mmc2.xlsx)
- [Cui et al. source RSS paper](https://doi.org/10.1158/1078-0432.CCR-18-0825) and [supplemental source deposit](https://aacr.figshare.com/articles/journal_contribution/22470557)
