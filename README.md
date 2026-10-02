# Cross-cancer radiotherapy signature transfer

Initial implementation of the attached proposal: evaluate fixed Cui et al. breast-cancer radiosensitivity and immune signatures in pretreatment NSCLC data from GSE253564. The proposal's contents concern radiotherapy, despite its chemotherapy filename.

## Run

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
radio-transfer download --output data/raw
radio-transfer analyze --expression data/raw/GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz --metadata data/clinical.csv --signatures config/rss.json config/immune.json --output results/run_01
```

Docker:

```bash
docker build -t radio-transfer .
docker run --rm -v "$PWD/data:/app/data" radio-transfer download
docker run --rm -v "$PWD/data:/app/data" -v "$PWD/config:/app/config:ro" -v "$PWD/results:/app/results" radio-transfer analyze --expression data/raw/GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz --metadata data/clinical.csv --signatures config/rss.json config/immune.json
```

## Required inputs

The expression TSV has gene symbols in its first column and sample IDs in its remaining columns. Non-symbol identifiers and duplicate gene mappings must be reconciled explicitly upstream; the pipeline does not guess aliases or aggregate duplicate genes.

Create `clinical.csv` using `config/metadata_template.csv`. `baseline`, `mpr`, and `verified` use 1/0; arms are `durvalumab` or `durvalumab_sbrt`. Each label requires a source and a patient ID. Unverified/nonbaseline/unresolved samples are listed as excluded. Duplicate eligible patients stop analysis.

Create two signature JSON files using `config/signature_template.json`: one `primary` RSS with exactly 34 genes, and one `secondary` immune signature. Each gene entry is `{"gene": "GENE_SYMBOL", "weight": 0.123}`. Obtain the actual genes, coefficients, transformation, and direction from the Cui paper and supplement. `direction` is +1 or -1 so larger oriented scores correspond to the published sensitive/immune-effective direction. Set `source_verified` to true and record `locked_at` only after checking and locking the source specification **before inspecting outcomes**. These declarations are an audit record, not independent proof of verification.

Supported transforms are `identity`, `log2p1`, and `log2p1_gene_zscore` (population SD across the baseline expression cohort). Use one only if justified by the locked source method and platform adaptation. The engine is a weighted sum, not a claim that all source scoring methods are already reproduced; extend it if the supplement requires a different formula. There are deliberately no fabricated gene weights or patient response labels in this repository.

## Outputs

Scores, complete sample-flow/exclusion table, gene QC, exact two-sided rank permutation tests in the combination arm, rank probability effect sizes with stratified bootstrap percentile intervals, Holm correction across the two signatures, signed rank adaptation and leave-one-out effects, and a manifest with input SHA-256 hashes, signature specifications, versions, and seed.

Rank probability is P(MPR score > non-MPR score) + half the tie probability. It is reported as an association effect, not predictive accuracy. Signed rank adaptation ranks every measured gene within each sample and averages signed ranks; it is separate from the source score. The implementation has no label-driven gene selection, classifier fitting, threshold optimization, or ROC claims.

## Current status

Pipeline and synthetic tests are implemented. Real source coefficients, verified clinical linkage, identifier reconciliation, QC figures, and an estimability-aware exploratory treatment interaction remain to be completed. Real cohort analysis has not run. A blocked download is an error, not evidence that the dataset is unavailable. Confirm the current GEO supplementary filename before downloading if it changes.

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

See `docs/analysis_plan.md` for the initial analysis contract.
