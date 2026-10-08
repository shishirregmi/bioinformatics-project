# Exploratory analysis contract

Primary source: Cui et al. (2018), https://doi.org/10.1158/1078-0432.CCR-18-0825.
Target: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253564.
Trial and molecular report: https://doi.org/10.1038/s41467-023-44195-x and https://doi.org/10.1016/j.xcrm.2024.101438.

1. Archive the processed pretreatment FPKM matrix and checksum. Keep GSE248378 post-treatment specimens separate from baseline.
2. Link each pretreatment expression column to the exact Study number in Table S1. Map Study Arm 1 to durvalumab and Study Arm 2 to durvalumab plus SBRT. Use the published Pathology Response value and MPR definition to classify response; retain source and exclusions in config/clinical.csv and docs/clinical_metadata.md. The published 32-profile and 10/6 arm-2 group counts are cross-checks, never a rule for assigning labels.
3. The baseline matrix contains an Entrez.ID annotation field in addition to gene symbols and sample profiles. Drop that annotation from the numeric expression matrix; do not count it as a sample. Do not silently alias or aggregate gene symbols.
4. Use the source-verified, locked gene lists, coefficients, formula, direction, and expression transformation in config/rss.json and config/immune.json. The configurations were locked before outcome analysis. Missing genes, duplicate mappings, nonfinite expression, and unsupported formulas must stop analysis.
5. Score all baseline samples without using outcome labels. The RNA-seq adaptation is log2(FPKM + 1), followed by per-gene standardization across baseline samples using sample SD. This is explicit and reproducible, but not an exact reproduction of the breast-cancer source preprocessing.
6. Compare MPR and no-MPR scores within the durvalumab + SBRT arm using a two-sided exact rank-sum label permutation test including ties. Report medians, rank probability, and descriptive 95% stratified bootstrap percentile intervals (2,000 resamples; seed 4370); adjust the two prespecified signature tests with Holm. The intervals are unstable in this small cohort.
7. Report signed within-sample gene ranks and leave-one-out effects as sensitivity analyses. These are not substitutes for the locked scores.
8. Defer treatment-by-score interaction analyses until arm labels and event counts have been independently reviewed and a suitable method is prespecified. No treatment-benefit inference is implemented.

This is a reproducible exploratory contract, not a registered clinical validation study. Synthetic tests verify software behavior; the public-data crosswalk is checked separately against its source records.