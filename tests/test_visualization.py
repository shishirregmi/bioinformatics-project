import argparse
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd

from radio_transfer.visualize import run


ROOT = Path(__file__).parents[1]


class VisualizationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.signatures = [ROOT / "config/rss.json", ROOT / "config/immune.json"]
        genes = []
        for path in self.signatures:
            genes.extend(gene["gene"] for gene in json.loads(path.read_text())["genes"])
        genes = list(dict.fromkeys(genes)) + [f"GENE_{index}" for index in range(40)]
        rng = np.random.default_rng(4370)
        expression = pd.DataFrame(
            rng.lognormal(mean=1, sigma=1, size=(len(genes), 8)),
            index=genes,
            columns=[f"sample_{index}" for index in range(8)],
        )
        self.expression_path = self.root / "expression.tsv"
        expression.to_csv(self.expression_path, sep="\t")
        self.metadata_path = self.root / "metadata.csv"
        pd.DataFrame({
            "sample_id": expression.columns,
            "patient_id": [f"patient_{index}" for index in range(8)],
            "baseline": ["1"] * 8,
            "arm": ["durvalumab"] * 4 + ["durvalumab_sbrt"] * 4,
            "mpr": ["0", "1", "0", "1", "0", "0", "1", "1"],
            "label_source": ["synthetic test fixture"] * 8,
            "verified": ["1"] * 8,
        }).to_csv(self.metadata_path, index=False)

    def tearDown(self):
        self.temp.cleanup()

    def _run(self, output, metadata=None):
        args = argparse.Namespace(
            expression=str(self.expression_path),
            metadata=str(metadata) if metadata else None,
            signatures=[str(path) for path in self.signatures],
            output=str(output),
        )
        run(args)
        return json.loads((Path(output) / "visualization_summary.json").read_text())

    def test_expression_figures_do_not_require_outcome_labels(self):
        output = self.root / "unlabeled"
        summary = self._run(output)
        self.assertTrue(summary["heatmap_generated"])
        self.assertFalse(summary["clinical_labels_present"])
        self.assertFalse(summary["response_plot_generated"])
        self.assertTrue((output / "expression_pca.png").exists())
        self.assertTrue((output / "signature_gene_coverage.png").exists())
        notes = (output / "visualization_notes.md").read_text()
        self.assertIn("No clinical outcome labels were supplied", notes)

    def test_response_plot_requires_verified_labels(self):
        output = self.root / "labeled"
        summary = self._run(output, self.metadata_path)
        self.assertTrue(summary["clinical_labels_present"])
        self.assertTrue(summary["response_plot_generated"])
        self.assertEqual(summary["response_group_summary"]["mpr"], 2)
        self.assertEqual(summary["response_group_summary"]["no_mpr"], 2)
        self.assertTrue((output / "signature_scores_by_mpr.png").exists())


if __name__ == "__main__":
    unittest.main()
