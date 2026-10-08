# Pretreatment clinical sample mapping

`config/clinical.csv` contains the minimal public sample-to-clinical mapping needed for the prespecified baseline analysis. The mapping was checked against the pretreatment expression columns and Table S1 of the study paper.

## Crosswalk rules

1. `sample_id` is the exact `durva###` study number shared by the GSE253564 matrix and Table S1.
2. The study's `Study Arm` labels map as follows: Arm1 = durvalumab alone; Arm2 = durvalumab plus SBRT. The paper defines the treatment arms and the SBRT regimen.
3. `baseline=1` is restricted to samples marked `Pre treatment RNAseq = Yes` in Table S1.
4. Table S1's `Pathology Response` values are signed percent response values. The exact value is retained in `pathology_response_signed_pct`. MPR is set to 1 when the absolute value is at least 90, consistent with the paper's definition of MPR as no more than 10% residual viable tumor (at least 90% cancer-cell killing). No label is assigned from sample order or expected group counts.
5. The public `Study number` is used as the de-identified `patient_id`; there is one pretreatment sample per patient in this matrix.

## Checks against the paper

The crosswalk has 32 unique baseline profiles, 16 per arm. Within the 16 profiles in Arm2, the source response values classify 10 as MPR and 6 as no MPR, matching the published pretreatment counts. Arm1 contains one MPR among the 16 available baseline profiles. The validation counts check the source mapping; they were not used to assign individual labels.

## Sources

- [Study article](https://pmc.ncbi.nlm.nih.gov/articles/PMC10982989/): the Methods define MPR as residual viable tumor ≤10%; the Results report 32 pretreatment profiles and 10 MPR / 6 no-MPR in Arm2.
- [Table S1 workbook](https://pmc.ncbi.nlm.nih.gov/articles/PMC10982989/bin/mmc2.xlsx): `Study number`, `Study Arm`, `Pathology Response`, and `Pre treatment RNAseq` fields.
- [GSE253564 record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253564): pretreatment expression profiles.

The mapping contains only study sample IDs and the fields required for this analysis. It does not contain direct identifiers or patient-specific radiotherapy dosimetry.
