# Prespecified analysis plan

## Version and scope

This plan separates the primary external preclinical test from the exploratory clinical translation check. The published RSS is treated as fixed. No genes, weights, cutoffs, or model settings are selected from the outcome data.

## Research question

Does a breast-cancer radiosensitivity score transfer to NSCLC cell-line radiation survival, and is the score direction compatible with MPR in a separate small NSCLC cohort?

## Aim 1: preclinical radiation-survival association

### Population

Use Yard et al. Supplementary Data 1 records whose `Site` is `lung` and whose source `Subhistology` is one of the explicit non-small-cell labels in `preclinical.py`. Small-cell lung cancer and unspecified labels are excluded. LUAD is defined as `adenocarcinoma` plus `bronchioloalveolar_adenocarcinoma`. The report shows the source count and the count remaining after exact CCLE matching; it does not force the expression-linked count to equal the source count.

### Joining and expression scoring

Normalize case and punctuation from `Cell Line` plus `Site` and match against CCLE sample columns. Keep one-to-one exact normalized matches only. Do not fuzzy match or impute. Read gene symbols from the CCLE GCT `Description` column, retain the Ensembl `Name` identifiers for the audit, and sum repeated Ensembl rows sharing an RSS symbol before transformation. Drop cell-line profiles missing any of the 34 RSS genes.

For each signature gene, transform RPKM as `log2(RPKM + 1)` and standardize across CCLE cell-line profiles complete for the fixed RSS using the sample standard deviation. Apply the locked coefficient and direction from `config/rss.json` unchanged. No outcome-informed normalization, gene selection, or threshold is used.

### Endpoints and tests

The primary endpoint is the source study’s integral-survival AUC, which integrates survival across 1, 2, 3, 4, 5, 6, 8, and 10 Gy and is scaled from 0 (completely sensitive) to 7 (completely resistant). The primary effect is Spearman’s rank correlation between RSS and AUC in matched NSCLC lines. The prespecified LUAD subgroup repeats the same estimate.

For each estimate, use 10,000 two-sided Monte Carlo permutations of AUC among cell lines and 2,000 percentile bootstrap resamples of cell lines. The random seeds are fixed and recorded. Holm adjustment covers the two planned tests (NSCLC and LUAD). An estimate with an uncomputable or constant score is an explicit analysis failure, not grounds to change the locked scoring rule.

### Prespecified sensitivity and interpretation analyses

1. **Histology-adjusted sensitivity:** rank RSS and AUC, residualize both rank variables against source subhistology indicators, and correlate residuals. Permute AUC ranks within subhistology and bootstrap within subhistology to preserve observed subtype composition. This checks whether a pooled association is driven only by between-subtype differences. It is not a third confirmatory hypothesis test.
2. **RSS component exploration:** correlate each oriented weighted gene contribution with AUC in NSCLC. Report all 34 correlations and Benjamini–Hochberg q-values; label them exploratory and do not use them to change the signature.
3. **Join sensitivity:** report all source lines, exact matches, ambiguous identifiers, and complete-expression exclusions. The primary analysis includes only one-to-one exact normalized joins.

## Aim 2: clinical translation check

Use pretreatment GSE253564 expression and [`config/clinical.csv`](../config/clinical.csv), which is crosswalked to the source trial table. Require verified sample ID, patient ID, baseline status, treatment arm, and MPR label. Calculate locked RSS and immune scores before joining labels. The primary clinical comparison remains RSS versus MPR within the durvalumab + SBRT arm; the locked immune score is secondary.

Use the existing two-sided exact label-permutation comparison and rank probability. Bootstrap intervals are descriptive because the combination arm has only 16 eligible patients (10 MPR and 6 no MPR). Holm adjustment covers the two patient signatures. Do not estimate treatment-by-score interaction, classifier accuracy, an ROC curve, clinical utility, or a radiation-specific causal effect.

## Missingness and exclusions

- Missing required RSS gene coverage stops the preclinical score rather than triggering gene substitution.
- A CCLE row missing one or more required expression values is excluded from scoring without imputation and is visible in join or sample-flow outputs.
- Unmatched and ambiguous cell-line records remain in `cell_line_join_audit.csv`.
- Clinical rows with an unverified label, unresolved patient, absent expression profile, unresolved outcome, or non-baseline status remain in `sample_flow.csv` with an exclusion reason.
- Repeated eligible patient samples stop the clinical analysis pending reconciliation.

## Reproducibility records

Each report includes source paths, SHA-256 checksums, signature specification, cell-line/sample counts, transform choices, random seeds, software versions, and analysis parameters. Unit and synthetic end-to-end tests verify implementation behavior. They do not validate biological association or source labels.

## Interpretation boundaries

The cell-line AUC measures in-vitro survival after a multi-dose irradiation assay. It does not represent patient-level tumor response or delivered dose. The MPR label comes from a combined neoadjuvant treatment setting and cannot isolate SBRT’s causal effect. The project has no patient-specific RT plan, RTDOSE, RTSTRUCT, dose-volume histogram, or organ-at-risk dose data. Findings are exploratory and are not for clinical decision-making.
