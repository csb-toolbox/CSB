import unittest
from csb.bio.io import HHpredOutputParser

class TestHHpredOutputParser(unittest.TestCase):
    def setUp(self):
        filename = 'data/test_HHpredOutputParser.hhr'
        tmp = HHpredOutputParser()
        self.hitlist = tmp.parse_file(filename)

    def test_parse_file(self):
        '''
        TODO: write tests for alignment parsing and data contained in the hits.
        '''
        self.assertEqual(len(self.hitlist), 10)

        self.assertEqual(self.hitlist.command, './hhsearch -i /ebio/abt1/kopec/projects/hhpred_on_cluster/split_fastas/d1ea0a1.hhm -d /ebio/abt1/kopec/projects/hhpred_on_cluster/astral_scop70_all_beta_only_v1.75_hhpred_database.hhm')

        self.assertEqual(self.hitlist.query_name, 'd1ea0a1 b.80.4.1 (A:1203-1472) Alpha subunit of glutamate synthase, C-terminal domain {Azospirillum brasilense [TaxId: 192]}')

        self.assertEqual(self.hitlist.neff, float('5.2'))

        self.assertEqual(self.hitlist.searched_hmms, 2936)

        self.assertRaises(IndexError, self.hitlist.__getitem__, 20)

if __name__ == '__main__':
    unittest.main()
