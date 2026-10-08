# Master’s project proposal

## Working title

**Cross-scale transfer of a breast-cancer radiosensitivity signature to NSCLC: preclinical radiation-survival evaluation and exploratory clinical translation**

## Why this is a master’s project

The original software milestone was a small clinical proof of concept with 16 patients in the treatment arm of interest. That is too small to support a master’s thesis by itself. This expanded design adds a larger, independent radiation-response phenotype for the primary biological test, a pre-registered histology subgroup, an analysis that checks subtype confounding, and component-level interpretation of the fixed signature. The patient cohort remains a separate translation check rather than the primary source of evidence.

The work requires data provenance review, cross-platform gene identifier harmonization, cell-line/sample matching, locked-score implementation, nonparametric inference, uncertainty estimation, subtype sensitivity analysis, and a careful interpretation across in-vitro and clinical contexts. It stays within one biological question and avoids model benchmarking on a cohort too small for it.

## Background and rationale

Radiotherapy response differs across cancer lineages and among tumors within a lineage. The 2016 Yard et al. study profiled radiation survival in 533 human cancer cell lines and reports 89 NSCLC lines, including 39 LUAD lines. Its integral-survival endpoint summarizes cellular survival across multiple radiation doses and was compared with clonogenic survival in a subset of lines. The published 34-gene radiosensitivity signature (RSS) was developed from breast-cancer data. Testing the fixed score against the lung-cell-line radiation phenotype asks whether its expression pattern transfers to another cancer type without retraining.

A second dataset offers a limited clinical context: baseline RNA-seq and pathologic response labels from a neoadjuvant NSCLC study of durvalumab with or without stereotactic body radiotherapy (SBRT). The combination arm is small, and MPR can reflect immune response as well as radiation. Therefore, this cohort can show whether the score’s direction is compatible with a patient endpoint; it cannot establish radiation-specific causality or predictive performance.

The first source-checked cell-line run linked 86 NSCLC and 37 LUAD lines with complete RSS expression. The NSCLC association was near zero (ρ = −0.027; 95% bootstrap CI −0.246 to 0.204; permutation *p* = 0.806), as was the LUAD estimate. The original small patient pilot showed an inverse RSS–MPR rank probability of 0.15 (exact *p* = 0.0225; 10 MPR, 6 non-MPR). These results do not support transfer or a clinical biomarker claim. Their difference makes cross-scale portability and treatment/immune context a useful research question, while remaining insufficient to identify a biological mechanism.

## Central question and hypotheses

**Central question:** Does the fixed breast-cancer RSS track intrinsic radiation survival in NSCLC cells, and is its direction compatible with pathologic response in a separate small NSCLC patient cohort?

**Primary hypothesis:** In the external cell-line panel, higher oriented RSS scores will associate with lower integral survival AUC, consistent with greater predicted radiosensitivity.

**Prespecified subgroup hypothesis:** The score–AUC association will be examined separately in LUAD to determine whether the primary result is present in a defined NSCLC subtype.

**Clinical translation question:** In the durvalumab + SBRT arm, is baseline RSS associated with MPR? This is exploratory and has no directional confirmatory claim because of the small cohort and the combined treatment/immune context.

## Specific aims

### Aim 1 — Test the locked RSS in the cell-line radiation-response panel

1. Select the explicit NSCLC subhistologies in Yard et al. Supplementary Data 1 and the defined LUAD histologies; exclude small-cell and unspecified labels.
2. Join radiation-response records to CCLE RNA-seq profiles by normalized exact cell-line name plus tissue site. Preserve unmatched and ambiguous rows in an audit file; do not use fuzzy joins.
3. Map CCLE Ensembl rows through the GCT `Description` symbol field. Sum multiple rows with the same symbol before log transformation. Exclude lines with incomplete RSS expression without imputation.
4. Apply the locked RSS weights and direction after log2(RPKM + 1) transformation and per-gene standardization across complete CCLE cell-line profiles.
5. Estimate Spearman correlation between RSS and the continuous integral-survival AUC in NSCLC, with a prespecified LUAD analysis. Use two-sided permutation tests and bootstrap intervals; adjust these two p-values with Holm’s method.
6. Perform a histology-adjusted partial Spearman sensitivity analysis with permutations restricted within source subhistology, so that the pooled result can be compared with the within-subtype association.
7. As exploratory interpretation, estimate per-gene weighted-component correlations with AUC and apply Benjamini–Hochberg correction across the 34 genes. These are not a new signature or feature-selection step.

### Aim 2 — Check clinical direction in the small patient cohort

Apply the existing source-verified clinical mapping to pretreatment GSE253564 expression. Within the durvalumab + SBRT arm, compare the primary RSS with MPR using the project’s exact rank-permutation approach. The 4-gene immune score is secondary and remains patient-only. Report sample flow, score distributions, rank probabilities, exact p-values, uncertainty, and multiplicity correction. Do not test treatment-by-score interaction or report predictive accuracy.

## Data sources and feasibility

| Dataset | Size and contents | Role |
| --- | --- | --- |
| Yard et al. Supplementary Data 1 | 533 cancer cell lines with cancer site, histology, culture medium, SNP-fingerprint status, and integral radiation-survival AUC | Primary radiation-response phenotype |
| CCLE RNA-seq RPKM GCT | 1,019 cell-line expression profiles, Ensembl gene rows, and gene-symbol descriptions | Expression source for the fixed RSS |
| GSE253564 | 32 pretreatment tumor RNA-seq profiles | Clinical translation check; 16 patients in the durvalumab + SBRT arm |
| GSE248378 | Deposited post-treatment resection RNA-seq profiles | Separate visualization and data completeness exploration; not in the primary endpoint |
| Cui et al. locked RSS | 34 genes with published weights and direction | Prespecified primary exposure |
| Cui et al. immune score | 4 genes | Secondary clinical context only |

All primary inputs are public. The container downloads and checksums them in a single run. The CCLE GCT is about 137 MB compressed; the pipeline streams it in chunks and retains only the signature genes. The project does not require patient consent, new sequencing, or new radiotherapy plans.

## Variables and interpretation

| Variable | Meaning | Use |
| --- | --- | --- |
| CCLE `Name` | Ensembl gene identifier | Preserve for identifier audit |
| CCLE `Description` | Gene-symbol annotation | Map RSS symbols to expression rows |
| RPKM | CCLE RNA-seq relative expression | Transform and standardize for RSS scoring |
| `Cell Line`, `Site` | Radiation-panel identifier and tissue origin | Exact normalized expression join and cohort eligibility |
| `Histology`, `Subhistology` | Source tumor classification | Define NSCLC/LUAD and adjust the sensitivity analysis |
| Integral-survival AUC | Integrated response to 1, 2, 3, 4, 5, 6, 8, and 10 Gy | Continuous cell-intrinsic endpoint; 0 is sensitive and 7 resistant |
| RSS score | Weighted sum using the published 34 genes, coefficients, and direction | Fixed expression-based exposure; no cutoff or refit |
| Patient sample ID, arm, MPR | Verified pretreatment sample and clinical labels | Separate small clinical translation analysis |
| Clinical immune score | Locked published 4-gene score | Secondary patient analysis; not a cell-line microenvironment score |

## Statistical analysis

- **Primary effect size:** Spearman’s ρ between oriented RSS and integral-survival AUC in NSCLC.
- **LUAD:** same prespecified effect in LUAD; Holm adjustment covers the two association tests.
- **Uncertainty:** 10,000 two-sided permutation draws and 2,000 bootstrap resamples, deterministic seeds recorded in the report.
- **Subtype sensitivity:** partial Spearman correlation after residualizing score and AUC ranks against subhistology; permute AUC ranks within subhistology and bootstrap within subhistology.
- **Component exploration:** per-gene Spearman association of the weighted RSS component with AUC; approximate p-values and BH q-values across the 34 genes.
- **Clinical translation:** exact label-permutation tests for the two locked scores in the combination arm, with Holm adjustment; descriptive bootstrap interval for rank probability.
- No expression cutoff tuning, gene selection, model fitting, cross-validation, ROC/AUC prediction metric, or treatment-interaction test is performed.

## Expected contribution

The thesis can distinguish three possible outcomes: (1) the fixed signature has a coherent association with radiation survival in lung cells and a compatible clinical direction; (2) it associates in cell lines but does not translate to the small patient endpoint, suggesting tumor-microenvironment or treatment-context differences; or (3) the signature does not transfer, which is still informative about limits of cross-cancer signature portability. Component associations and subtype adjustment help explain whether a pooled result is distributed across genes and histologies or driven by a narrow subset.

No result is promised in advance. A null or inverse association will be reported as a result, not repaired by choosing a new score, gene subset, or cutoff.

## Limitations and safeguards

1. **In-vitro to patient gap:** Cell lines omit stromal, vascular, immune, and pharmacologic context. AUC is not clinical tumor control.
2. **Small patient cohort:** Sixteen patients in the combination arm are insufficient for reliable prediction, subgroup discovery, or treatment interactions.
3. **Expression platform adaptation:** The original signature and CCLE RPKM are from different assays. The log transform and per-gene z-score are explicit adaptations, not a reproduction of the original clinical platform.
4. **Cell-line matching:** Some radiation lines may not have an exact CCLE expression match. The pipeline records all exclusions and does not infer missing expression.
5. **No patient dosimetry:** The clinical dataset gives the study-level SBRT regimen, not each patient’s delivered tumor or normal-tissue dose.
6. **Component analysis:** Individual gene correlations are exploratory and multiplicity controlled; they do not establish mechanism.
7. **No clinical validation claim:** Neither association analysis supports clinical decisions, a treatment-benefit biomarker claim, or individual risk prediction.

## Semester plan

| Weeks | Milestone |
| --- | --- |
| 1–2 | Review the two source studies, verify variable definitions, lock the analysis contract |
| 3–4 | Download/checksum source data; implement gene and cell-line identifier audits |
| 5–6 | Validate histology cohorts and reproduce expected 533/89/39 source counts |
| 7–8 | Run fixed-score primary, LUAD, and subtype-adjusted analyses |
| 9–10 | Complete component interpretation and clinical translation check |
| 11–12 | Review exclusions, sensitivity checks, visualizations, and limits with the supervisor |
| 13–14 | Draft results/discussion and prepare a reproducible presentation |

## Deliverables

- Reproducible Docker pipeline and source-checksummed manifests
- Interactive local project dashboard and complete data dictionary
- Cell-line matching and cohort-flow audit files
- Preclinical primary/subgroup/sensitivity results and component table
- Clinical translation report with verified MPR crosswalk
- Thesis methods/results figures and a limitations-first presentation

## Primary references

- Yard BD et al. [A genetic basis for the variation in the vulnerability of cancer to DNA damage](https://doi.org/10.1038/ncomms11428). *Nature Communications*. 2016.
- Cui Y et al. [A 34-gene signature predicts radiation response in breast cancer](https://doi.org/10.1158/1078-0432.CCR-18-0825). *Clinical Cancer Research*.
- [GSE253564 and GSE248378 source study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10982989/).
- [CCLE public data directory](https://data.broadinstitute.org/ccle/).
