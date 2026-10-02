import csv
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

from radio_transfer.core import load_signature, score_expression

ROOT = Path(__file__).parents[1]


class CuiTests(unittest.TestCase):
    def test_rss_matches_source_table(self):
        spec = load_signature(ROOT / 'config/rss.json')
        with (ROOT / 'config/rss_source_table.csv').open() as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual(len(rows), 34)
        self.assertEqual(spec['direction'], -1)
        self.assertEqual(spec['genes'], [{'gene':r['gene'], 'entrez_id':r['entrez_id'], 'weight':float(r['weight'])} for r in rows])
        self.assertEqual(spec['genes'][0]['weight'], -2.03)
        self.assertEqual(next(g['weight'] for g in spec['genes'] if g['gene']=='H2AFJ'), 5.5)

    def test_immune_formula_sample_sd_and_direction(self):
        spec = load_signature(ROOT / 'config/immune.json')
        # log2(x+1) yields [0,1,2], whose sample-standardized values are [-1,0,1].
        frame = pd.DataFrame([[0,1,3]]*4, index=['ADRM1','MICB','PSMD13','RFXANK'], columns=['a','b','c'])
        scores = score_expression(frame, spec)
        np.testing.assert_allclose(scores, [9.4,0,-9.4], atol=1e-12)

    def test_constant_gene_rejected(self):
        spec = load_signature(ROOT / 'config/immune.json')
        frame = pd.DataFrame([[1,1,1]]*4, index=[g['gene'] for g in spec['genes']])
        with self.assertRaisesRegex(ValueError, 'constant'):
            score_expression(frame, spec)
