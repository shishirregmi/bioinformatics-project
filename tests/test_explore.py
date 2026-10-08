import argparse
import gzip
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from radio_transfer.explore import run
from radio_transfer.cli import download, GEO_DATASETS


ROOT = Path(__file__).parents[1]


class DataExplorationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.signatures = [ROOT / "config/rss.json", ROOT / "config/immune.json"]
        genes = []
        for path in self.signatures:
            genes.extend(gene["gene"] for gene in json.loads(path.read_text())["genes"])
        genes = list(dict.fromkeys(genes)) + [f"GENE_{index}" for index in range(40)]
        rng = np.random.default_rng(719)
        for filename, samples in (("baseline.tsv", 8), ("post.tsv", 7)):
            expression = pd.DataFrame(
                rng.lognormal(mean=1, sigma=1, size=(len(genes), samples)),
                index=genes,
                columns=[f"{filename.split('.')[0]}_{index}" for index in range(samples)],
            )
            expression.to_csv(self.root / filename, sep="\t")
        self.metadata = self.root / "metadata.csv"
        pd.DataFrame({
            "sample_id": [f"baseline_{index}" for index in range(8)],
            "patient_id": [f"patient_{index}" for index in range(8)],
            "baseline": ["1"] * 8,
            "arm": ["durvalumab", "durvalumab_sbrt"] * 4,
            "mpr": ["0", "1"] * 4,
            "label_source": ["synthetic test fixture"] * 8,
            "verified": ["1"] * 8,
        }).to_csv(self.metadata, index=False)

    def tearDown(self):
        self.temp.cleanup()

    def test_exploration_profiles_both_timepoints_and_metadata(self):
        output = self.root / "exploration"
        args = argparse.Namespace(
            baseline_expression=str(self.root / "baseline.tsv"),
            post_treatment_expression=str(self.root / "post.tsv"),
            metadata=str(self.metadata),
            signatures=[str(path) for path in self.signatures],
            output=str(output),
        )
        run(args)
        inventory = json.loads((output / "data_inventory.json").read_text())
        self.assertEqual([row["dataset"] for row in inventory["datasets"]], ["baseline", "post_treatment"])
        self.assertEqual([row["sample_count"] for row in inventory["datasets"]], [8, 7])
        self.assertTrue((output / "baseline/expression_pca.png").exists())
        self.assertTrue((output / "post_treatment/expression_pca.png").exists())
        self.assertTrue((output / "data_overview.png").exists())
        report = (output / "index.html").read_text()
        self.assertIn("Radiotherapy RNA-seq data explorer", report)
        self.assertIn("baseline/expression_pca.png", report)
        self.assertIn("gene_parameter_profile.csv", report)
        self.assertTrue((output / "gene_parameter_profile.csv").exists())
        fields = pd.read_csv(output / "clinical_metadata_fields.csv")
        self.assertIn("mpr", set(fields.field))
        notes = (output / "data_inventory.md").read_text()
        self.assertIn("29 samples", notes)
        self.assertIn("8 Gy × 3", notes)
        self.assertIn("patient-specific", notes)

    def test_post_treatment_download_uses_its_own_geo_file(self):
        output = self.root / "downloads"
        frame = pd.DataFrame({"sample_1": [1.0, 2.0], "sample_2": [3.0, 4.0]}, index=["GENE1", "GENE2"])
        source = io.BytesIO()
        with gzip.GzipFile(fileobj=source, mode="wb") as compressed:
            compressed.write(frame.to_csv(sep="\t").encode())

        class Response:
            def __init__(self, data):
                self.data = io.BytesIO(data)

            def __enter__(self):
                return self

            def __exit__(self, *args):
                return False

            def read(self, size=-1):
                return self.data.read(size)

        args = argparse.Namespace(output=str(output), dataset="post-treatment")
        with patch("radio_transfer.cli.urllib.request.urlopen", return_value=Response(source.getvalue())) as opened:
            download(args)
        opened.assert_called_once_with(GEO_DATASETS["post-treatment"]["url"], timeout=60)
        matrix_path = output / GEO_DATASETS["post-treatment"]["filename"]
        self.assertTrue(matrix_path.exists())
        manifest = json.loads((output / "download_manifest.json").read_text())
        self.assertEqual(manifest["downloads"][0]["accession"], "GSE248378")


if __name__ == "__main__":
    unittest.main()
