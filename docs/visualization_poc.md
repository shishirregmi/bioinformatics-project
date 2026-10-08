# Data and visualization guide

The project connects two kinds of radiotherapy data: an in-vitro radiation-survival panel with matched cell-line expression, and a small clinical NSCLC study with pretreatment RNA-seq and pathologic response. The preclinical cell-line panel is the primary biological test; the patient data provide a separate translation check.

## Data blocks and parameters

| Data block | Parameters available | How the project uses them |
| --- | --- | --- |
| [Yard et al. Supplementary Data 1](https://doi.org/10.1038/ncomms11428) | `Master_ccl_id`, `Cell Line`, `Site`, `Histology`, `Subhistology`, `Culture_media`, `Snp_fp_status`, `AUC` | Filter explicit NSCLC/LUAD source labels; treat integral-survival AUC as the continuous radiation-response endpoint; audit source line and histology information. |
| [CCLE RNA-seq RPKM GCT](https://data.broadinstitute.org/ccle/CCLE_RNAseq_genes_rpkm_20180929.gct.gz) | `Name` (Ensembl ID), `Description` (gene symbol), cell-line sample columns, RPKM | Match by normalized exact cell-line name and site; map the fixed RSS genes and calculate the expression score. |
| [GSE253564](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253564) | Gene symbol, sample ID, FPKM; GEO sample records | Describe baseline tumor expression, check signature coverage, and run the verified MPR translation analysis. |
| [GSE248378](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248378) | Gene symbol, sample ID, FPKM; GEO sample records | Explore post-treatment resection expression separately. Do not assume sample pairing or treat it as baseline. |
| Clinical metadata and trial source table | `sample_id`, `patient_id`, `baseline`, `arm`, `mpr`, label source, verification status; trial SBRT regimen | Verify sample-to-patient labels and treatment context; the regimen was 8 Gy × 3 daily fractions. It is not individual delivered-dose information. |
| Locked signature JSON files | Gene symbols, published weights, score direction, transform, source and lock date | Apply the published 34-gene RSS as the primary score and the 4-gene immune score as secondary in patients only. |

### What the radiation-response AUC represents

Yard et al. integrated cell survival after 1, 2, 3, 4, 5, 6, 8, and 10 Gy. The endpoint is scaled from 0 (completely sensitive) to 7 (completely resistant). It summarizes a multi-dose in-vitro assay; it is not a patient dose, tumor dose-volume histogram, or clinical outcome.

The workbook contains 533 cell-line records. The article reports 89 NSCLC and 39 LUAD lines. Lung-site rows also include SCLC and unspecified histologies; the project uses explicit histology labels rather than treating all lung-site records as NSCLC. The full audit retains unmatched and ambiguous rows. Only one-to-one normalized exact joins to CCLE are used; there is no fuzzy matching or expression imputation.

CCLE’s GCT uses Ensembl identifiers in `Name` and gene-symbol annotations in `Description`. The project records this crosswalk for the 34 signature genes. If more than one Ensembl row has the same symbol, the RPKM values are summed before log transformation. Cell lines missing any required RSS expression value are excluded from scoring.

### Clinical RT information available

The trial paper describes a protocol-level regimen for the combination arm: durvalumab plus SBRT at 8 Gy for 3 daily fractions (24 Gy total). The public RNA-seq matrices and clinical mapping do not contain individual RT plans, RTDOSE/RTSTRUCT DICOM, dose-volume histograms, delivered-dose deviations, or organ-at-risk dosimetry. A patient dose-response analysis would require those separate records.

## Visualize all data with Docker

```bash
docker build -t radio-transfer .
mkdir -p data results
docker run --rm -p 8000:8000 \
  -v "$PWD/data:/app/data" \
  -v "$PWD/results:/app/results" \
  radio-transfer
```

Open [http://localhost:8000](http://localhost:8000). The first run downloads the two GEO matrices, the radiation-response workbook, and the CCLE expression file, then creates the patient RNA-seq explorer, radiation-data overview, preclinical association report, clinical MPR report, and dashboard. CCLE is about 137 MB compressed. Files are cached under `data/raw/`; timestamped outputs are saved under `results/`.

To add a corrected clinical mapping, place a verified CSV at `data/clinical.csv` and set `-e METADATA=/app/data/clinical.csv` on `docker run`. The default mapping is source-checked and explained in [`clinical_metadata.md`](clinical_metadata.md).

## Run an individual visualization or analysis

Download all public sources once:

```bash
radio-transfer download --dataset all --output data/raw
```

Run the primary cell-line analysis:

```bash
radio-transfer preclinical \
  --radiation-response data/raw/41467_2016_BFncomms11428_MOESM708_ESM.xlsx \
  --expression data/raw/CCLE_RNAseq_genes_rpkm_20180929.gct.gz \
  --signature config/rss.json \
  --output results/preclinical
```

Run the two expression matrices through the patient-data explorer:

```bash
radio-transfer explore \
  --baseline-expression data/raw/GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz \
  --post-treatment-expression data/raw/GSE248378_Durva_Post_FPKMs.txt.gz \
  --signatures config/rss.json config/immune.json \
  --output results/exploration
```

## Views and files

| View or file | What it represents | How it helps |
| --- | --- | --- |
| `index.html` | Project dashboard with links to each report | Presents both study aims and their interpretation limits together. |
| `exploration/index.html` | Baseline and post-treatment expression views, data dictionary, and file inventory | Inspect every available RNA-seq parameter without mixing timepoints. |
| `exploration/data_overview.png` | GEO sample counts, measured-gene counts, and per-sample median FPKM | Quick inventory and expression QC summary. |
| `exploration/baseline/` and `post_treatment/` | Expression distributions, separate PCA, gene coverage, coefficients, heatmaps, and optional MPR plot | Visualize each GEO dataset and its fixed signatures. |
| `preclinical/radiation_data_overview.png` | AUC distribution and lung-site subhistology counts across the radiation panel | See the radiotherapy measurements and why lung-site count is not the NSCLC analysis count. |
| `preclinical/radiation_parameter_inventory.csv` | Definitions and intended uses for AUC, cell-line metadata, CCLE IDs, Ensembl IDs, symbols, and RPKM | Give a professor a compact radiation-data dictionary. |
| `preclinical/preclinical_report.html` | Fixed RSS vs. AUC in NSCLC and LUAD, histology-adjusted sensitivity, and gene-component associations | Review the primary external test and exploratory interpretation. |
| `preclinical/cell_line_join_audit.csv` | All source cell lines, histology flags, CCLE matching status, and identifiers | Inspect included, unmatched, ambiguous, and incomplete records. |
| `preclinical/cell_line_scores.csv` | Matched cell-line subtype, AUC, expression column, and RSS | Review the analysis cohort and per-line measurements. |
| `preclinical/statistics.csv` | Spearman estimates, permutation p-values, bootstrap intervals, and Holm adjustment | Understand effect size and uncertainty for the planned tests. |
| `preclinical/gene_component_associations.csv` | All 34 weighted-component associations and BH q-values | Explore whether individual score components align with the phenotype. |
| `analysis/analysis_report.html` | Small clinical MPR comparison and limitations | Treat patient results as a separate exploratory translation check. |
| `parameter_inventory.csv`, `sample_parameter_profile.csv`, `gene_parameter_profile.csv` | Expression and clinical parameter descriptions and summaries | Inspect individual RNA-seq samples and gene measurements. |

## Interpret the results carefully

The primary association is preclinical and cell-intrinsic. The clinical arm has 16 patients (10 MPR, 6 no MPR), so it cannot support prediction, treatment interaction, or causal separation of SBRT from durvalumab and immune effects. Neither dataset includes patient-level dosimetry. The pipeline does not train a classifier, tune cutoffs, or claim clinical utility. A null or inverse association is reported as observed; it is not corrected by changing the signature.

## Direct source links

- [Yard et al. study](https://doi.org/10.1038/ncomms11428)
- [Yard et al. Supplementary Data 1 workbook](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fncomms11428/MediaObjects/41467_2016_BFncomms11428_MOESM708_ESM.xlsx)
- [CCLE RNA-seq RPKM matrix, compressed](https://data.broadinstitute.org/ccle/CCLE_RNAseq_genes_rpkm_20180929.gct.gz)
- [CCLE data directory](https://data.broadinstitute.org/ccle/)
- [GSE253564 pretreatment series and matrix](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253564)
- [GSE248378 post-treatment series and matrix](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248378)
- [Clinical NSCLC study and supplementary files](https://pmc.ncbi.nlm.nih.gov/articles/PMC10982989/)
- [Cui et al. 34-gene RSS source](https://doi.org/10.1158/1078-0432.CCR-18-0825)
- [Cui RSS supplement](https://aacr.figshare.com/articles/journal_contribution/22470557)
