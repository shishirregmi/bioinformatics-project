import argparse
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from radio_transfer.core import holm_two, load_expression, load_metadata, load_signature, rank_test, score_expression
from radio_transfer.cli import run


class PipelineTests(unittest.TestCase):
    def test_exact_separation_and_ties(self):
        self.assertAlmostEqual(rank_test([4, 5, 6], [1, 2, 3])["p_exact_two_sided"], .1)
        self.assertEqual(rank_test([1, 1], [1, 1])["p_exact_two_sided"], 1)
        self.assertEqual(rank_test([1, 1], [1, 1])["rank_probability"], .5)
        self.assertEqual(holm_two([.04, .01]), [.04, .02])

    def test_no_missing_gene_imputation(self):
        with self.assertRaisesRegex(ValueError, "Missing required"):
            score_expression(pd.DataFrame({"s": [1]}, index=["A"]),
                             {"genes": [{"gene": "B", "weight": 1}]})

    def test_metadata_excludes_unverified_and_missing(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "clinical.csv"
            path.write_text("sample_id,patient_id,baseline,arm,mpr,label_source,verified\ns1,p1,1,durvalumab_sbrt,1,paper,0\n")
            included, flow = load_metadata(path, {"s1", "s2"})
            self.assertTrue(included.empty)
            self.assertEqual(len(flow), 2)
            self.assertIn("unverified_label", flow.iloc[0].exclusion_reason)

    def test_bundled_public_clinical_crosswalk(self):
        path = Path(__file__).parents[1] / "config/clinical.csv"
        metadata, flow = load_metadata(path, set(pd.read_csv(path).sample_id))
        self.assertEqual(len(metadata), 32)
        combo = metadata[metadata.arm == "durvalumab_sbrt"]
        self.assertEqual((int((combo.mpr == 1).sum()), int((combo.mpr == 0).sum())), (10, 6))
        self.assertTrue(metadata.label_source.str.contains("PMC10982989 Table S1").all())
        self.assertTrue((metadata.verified == "1").all())
        source_response = metadata.pathology_response_signed_pct.astype(float).abs() >= 90
        np.testing.assert_array_equal(source_response.to_numpy(), metadata.mpr.to_numpy() == 1)

    def test_unlocked_signature_refused(self):
        path = Path(__file__).parents[1] / "config/signature_template.json"
        with self.assertRaisesRegex(ValueError, "verified and locked"):
            load_signature(path)

    def test_full_synthetic_run(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            genes = [f"G{i}" for i in range(34)]
            samples = [f"s{i}" for i in range(6)]
            expression = pd.DataFrame(np.arange(204).reshape(34, 6), index=genes, columns=samples)
            expression.to_csv(root / "expression.tsv", sep="\t")
            pd.DataFrame({"sample_id": samples, "patient_id": samples, "baseline": 1,
                          "arm": "durvalumab_sbrt", "mpr": [0, 0, 0, 1, 1, 1],
                          "label_source": "synthetic fixture", "verified": 1}).to_csv(root / "clinical.csv", index=False)
            paths = []
            for name, role, selected in [("rss", "primary", genes), ("immune", "secondary", genes[:2])]:
                spec = {"name": name, "role": role, "source": "synthetic fixture only", "source_verified": True,
                        "locked_at": "synthetic", "direction": 1, "transform": "identity",
                        "genes": [{"gene": g, "weight": 1} for g in selected]}
                path = root / f"{name}.json"
                path.write_text(json.dumps(spec))
                paths.append(str(path))
            args = argparse.Namespace(expression=str(root / "expression.tsv"), metadata=str(root / "clinical.csv"),
                                      signatures=paths, output=str(root / "out"))
            run(args)
            stats = pd.read_csv(root / "out/statistics.csv")
            self.assertEqual(stats.rank_probability.tolist(), [1., 1.])
            self.assertEqual(stats.p_holm.tolist(), [.2, .2])
            self.assertEqual(len(pd.read_csv(root / "out/sensitivity.csv")), 14)
            self.assertTrue((root / "out/analysis_report.html").exists())
            self.assertTrue((root / "out/score_by_mpr.png").exists())
            manifest = json.loads((root / "out/manifest.json").read_text())
            self.assertEqual(len(manifest["inputs"]), 4)
            with self.assertRaisesRegex(ValueError, "empty"):
                run(args)

    def test_negative_fpkm_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "expression.tsv"
            path.write_text("gene\tsample\nA\t-1\n")
            with self.assertRaisesRegex(ValueError, "nonnegative"):
                load_expression(path)

    def test_gene_annotation_column_is_not_a_sample(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "expression.tsv"
            path.write_text("Gene\tEntrez.ID\tdurva001\tdurva003\nA\t101\t2.5\t3.0\nB\t102\t4.5\t5.0\n")
            expression = load_expression(path)
            self.assertEqual(list(expression.columns), ["durva001", "durva003"])
            self.assertEqual(expression.shape, (2, 2))


if __name__ == "__main__":
    unittest.main()
