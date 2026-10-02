# Verified Cui signatures

Cui Y et al., Clinical Cancer Research 2018;24:4754-4762. DOI: https://doi.org/10.1158/1078-0432.CCR-18-0825. Main manuscript: https://pmc.ncbi.nlm.nih.gov/articles/PMC6168425/.

The publisher's supplementary material is available at https://aacr.figshare.com/articles/journal_contribution/22470557 (file 39922043). Its metadata associates it with the paper's DOI; its 2023 repository date is a deposit date, not a new study. Download URL, checksum, size and attribution are recorded in `cui_source_provenance.json`.

## Coefficients and formula

`config/rss_source_table.csv` preserves all 34 gene symbols, Entrez IDs and coefficients from Supplement Table 1, including the original labels LINS and H2AFJ. No aliases have been substituted: reconcile these identifiers against the actual target matrix using the source Entrez IDs before analysis. `config/rss.json` uses the same rows. Coefficients are used at the precision published in the table; more precise fitted parameters were not available in the retrieved source.

The RSS is the weighted sum of standardized expression, not an exponentiated hazard or a refitted model. The IMS formula in the main article is:

`IMS = 4.7*ADRM1 + 3.6*MICB + 4.8*PSMD13 - 3.7*RFXANK`.

Lower source scores correspond to the radiation-sensitive and immune-effective directions. Both configurations use `direction: -1`, so exported oriented scores are the negatives of raw weighted sums. Thus higher exported scores are in the prespecified favorable direction. The published breast cutoffs (RSS 1.0 and IMS -3.8) are recorded as provenance only; no threshold classification is performed in NSCLC.

## Normalization and transfer

The supplement says each dataset was transformed to z-scores to address platform differences. It does not state the scaling axis or population/sample SD divisor, and it describes microarray preprocessing rather than FPKM preprocessing. The target configuration therefore explicitly chooses log2(FPKM+1), followed by per-gene standardization across all baseline samples using sample SD (`zscore_ddof: 1`). This is a documented RNA-seq adaptation of the fixed source formula, not a claim of exact numerical reproduction of the breast microarray pipeline.

Both treatment arms must be retained in the baseline input cohort used for normalization; no outcome-specific centering is allowed. Missing genes and constant genes stop scoring. The cohort composition and expression transformation must remain fixed before outcomes are examined. Changes in input cohort change z-scores. The default contract remains exploratory transfer analysis.

## Reproduce extraction

Download the publisher's DOCX using the provenance URL, then run:

```bash
python scripts/extract_cui_rss.py /path/to/supplement.docx /tmp/rss_source_table.csv
```

The script checks the complete source SHA-256 and table schema before extracting. The source file itself is not duplicated in the repository. A test checks that the source CSV and RSS JSON agree, and a hand-computed synthetic test verifies IMS weighting, sample SD and reversed direction. These verify implementation, not biological transferability.
