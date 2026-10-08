import argparse
import gzip
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from radio_transfer.preclinical import (
    association_statistics,
    load_ccle_signature,
    load_radiation_response,
    match_ccle_columns,
    run,
)


ROOT = Path(__file__).parents[1]


def write_radiation_workbook(path, rows):
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        pd.DataFrame(rows).to_excel(writer, sheet_name="Supplementary Data 1", index=False, startrow=1)


def write_gct(path, genes, samples, values):
    with gzip.open(path, "wt", encoding="utf-8") as stream:
        stream.write("#1.2\n")
        stream.write(f"{len(genes)}\t{len(samples)}\n")
        stream.write("Name\tDescription\t" + "\t".join(samples) + "\n")
        for identifier, symbol, row in zip((g[0] for g in genes), (g[1] for g in genes), values):
            stream.write(identifier + "\t" + symbol + "\t" + "\t".join(str(value) for value in row) + "\n")


class PreclinicalTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.signature = ROOT / "config/rss.json"
        self.spec = json.loads(self.signature.read_text())
        self.genes = [entry["gene"] for entry in self.spec["genes"]]
        self.rng = np.random.default_rng(2034)

    def tearDown(self):
        self.temp.cleanup()

    def test_histology_selection_and_explicit_exact_join_audit(self):
        rows = [
            {"Master_ccl_id": 1, "Cell Line": "LungA", "Site": "lung", "Histology": "carcinoma", "Subhistology": "adenocarcinoma", "AUC": 2.0},
            {"Master_ccl_id": 2, "Cell Line": "LungS", "Site": "lung", "Histology": "carcinoma", "Subhistology": "small_cell_carcinoma", "AUC": 3.0},
            {"Master_ccl_id": 3, "Cell Line": "Other", "Site": "breast", "Histology": "carcinoma", "Subhistology": "adenocarcinoma", "AUC": 4.0},
            {"Master_ccl_id": 4, "Cell Line": "Unknown", "Site": "lung", "Histology": "carcinoma", "Subhistology": "NS", "AUC": 5.0},
        ]
        path = self.root / "radiation.xlsx"
        write_radiation_workbook(path, rows)
        radiation = load_radiation_response(path)
        self.assertEqual(int(radiation.is_nsclc.sum()), 1)
        self.assertEqual(int(radiation.is_luad.sum()), 1)
        joined = match_ccle_columns(radiation, ["LUNGA_LUNG", "OTHER_BREAST", "DUPLICATE-LUNG"])
        self.assertEqual(joined.loc[0, "expression_match_status"], "matched_exact_normalized_name_and_site")
        self.assertEqual(joined.loc[1, "expression_match_status"], "no_ccle_column_match")
        self.assertTrue(joined.loc[0, "is_nsclc"])

    def test_ccle_uses_description_symbols_sums_duplicate_rows_and_drops_incomplete_lines(self):
        samples = ["CELL1_LUNG", "CELL2_LUNG"]
        gene_rows = [(f"ENSG{i:011d}.1", gene) for i, gene in enumerate(self.genes)]
        gene_rows.extend([("ENSG99999999999.1", self.genes[0])])
        values = self.rng.uniform(0.1, 20, size=(len(gene_rows), len(samples))).round(3).astype(object)
        values[-1] = [1.0, 2.0]
        values[1, 1] = "NA"
        path = self.root / "ccle.gct.gz"
        write_gct(path, gene_rows, samples, values)
        expression, mapping, all_samples = load_ccle_signature(path, self.genes)
        self.assertEqual(all_samples, samples)
        self.assertEqual(list(expression.columns), ["CELL1_LUNG"])
        self.assertEqual(len(expression.index), 34)
        self.assertIn("ENSG99999999999.1", mapping.loc[mapping.gene == self.genes[0], "ensembl_ids"].iloc[0])
        expected = float(values[0, 0]) + 1.0
        self.assertAlmostEqual(expression.loc[self.genes[0], "CELL1_LUNG"], expected)

    def test_preclinical_pipeline_emits_auditable_results_and_visualization(self):
        rows = []
        for index in range(16):
            subtype = "adenocarcinoma" if index < 12 else "squamous_cell_carcinoma"
            rows.append({
                "Master_ccl_id": index + 1,
                "Cell Line": f"LUNG{index:02d}",
                "Site": "lung",
                "Histology": "carcinoma",
                "Subhistology": subtype,
                "Culture_media": "RPMI001",
                "Snp_fp_status": "SNP-matched-reference",
                "AUC": float(0.6 + index * 0.3),
            })
        rows.extend([
            {"Master_ccl_id": 100, "Cell Line": "SCLC1", "Site": "lung", "Histology": "carcinoma", "Subhistology": "small_cell_carcinoma", "AUC": 3.1},
            {"Master_ccl_id": 101, "Cell Line": "BREAST1", "Site": "breast", "Histology": "carcinoma", "Subhistology": "ductal_carcinoma", "AUC": 2.1},
        ])
        workbook = self.root / "radiation.xlsx"
        write_radiation_workbook(workbook, rows)

        samples = [f"LUNG{index:02d}_LUNG" for index in range(16)] + ["SCLC1_LUNG", "BREAST1_BREAST", "UNRELATED_BREAST"]
        gene_rows = [(f"ENSG{i:011d}.1", gene) for i, gene in enumerate(self.genes)]
        expression_values = self.rng.lognormal(mean=1.0, sigma=0.7, size=(len(gene_rows), len(samples))).round(5)
        expression_path = self.root / "ccle.gct.gz"
        write_gct(expression_path, gene_rows, samples, expression_values)
        output = self.root / "preclinical"
        run(argparse.Namespace(radiation_response=str(workbook), expression=str(expression_path),
                               signature=str(self.signature), output=str(output)))

        statistics = pd.read_csv(output / "statistics.csv")
        flow = pd.read_csv(output / "cohort_flow.csv")
        self.assertEqual(statistics.set_index("analysis").loc["NSCLC", "n"], 16)
        self.assertEqual(statistics.set_index("analysis").loc["LUAD", "n"], 12)
        self.assertEqual(flow.loc[flow.stage == "NSCLC histologies", "n"].iloc[0], 16)
        self.assertTrue((output / "preclinical_report.html").exists())
        self.assertTrue((output / "rss_vs_radiation_survival.png").exists())
        self.assertTrue((output / "radiation_data_overview.png").exists())
        self.assertIn("AUC", set(pd.read_csv(output / "radiation_parameter_inventory.csv").parameter))
        self.assertTrue((output / "cell_line_join_audit.csv").exists())

    def test_permutation_association_is_reproducible(self):
        x = np.arange(12, dtype=float)
        y = -x
        first = association_statistics(x, y, permutations=300, bootstraps=100)
        second = association_statistics(x, y, permutations=300, bootstraps=100)
        self.assertEqual(first, second)
        self.assertEqual(first["n"], 12)
        self.assertLess(first["spearman_rho"], -0.99)
        self.assertGreaterEqual(first["bootstrap_ci_low"], -1.0)
        self.assertLessEqual(first["bootstrap_ci_high"], -0.99)


if __name__ == "__main__":
    unittest.main()
