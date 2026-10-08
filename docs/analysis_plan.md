# Initial analysis contract

Primary source: Cui et al. (2018), https://doi.org/10.1158/1078-0432.CCR-18-0825.
Target: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253564.
Trial and molecular report: https://doi.org/10.1038/s41467-023-44195-x and https://doi.org/10.1016/j.xcrm.2024.101438.

1. Archive processed pretreatment FPKM and checksum. Do not combine GSE248378 post-treatment specimens with baseline profiles.
2. Reconcile expression/sample/patient identifiers. Verify baseline status, arm and trial-defined MPR against source metadata and supplements. Record each source and exclusion. Published counts (32 profiles; combination arm 10 MPR/6 non-MPR) are cross-checks, never a rule for assigning labels.
3. Verify and lock source genes, coefficients, formula, direction and transformation before outcomes are accessed. Do not silently fill missing genes. Missing genes, duplicate mappings, nonfinite expression and unsupported formulas must be resolved before analysis.
4. Score all baseline expression columns without outcome labels. Any platform adaptation must be documented. The software offers weighted-sum transforms but cannot determine which reproduces the source method.
5. Primary and secondary combination-arm comparisons use exact rank-sum label permutation, including ties, with absolute deviation from the null center. Report rank probability, group medians and descriptive 95% stratified bootstrap percentile intervals (2,000 draws; seed 4370). These intervals are unstable in tiny cohorts and do not imply clinical validation. Apply Holm to the two prespecified tests.
6. Sensitivity outputs include signed within-sample gene ranks and leave-one-out rank effects. They are not substitutes for the locked source scores.
7. Exploratory score-by-treatment interaction is deferred until labels and event counts are verified and a method handling separation is prespecified. No treatment-benefit inference is implemented.

This is a starting implementation contract, not the final time-stamped scientific preregistration. Review and freeze the completed contract and signature specifications before conducting the cohort analysis. Synthetic checks validate software behavior only.
