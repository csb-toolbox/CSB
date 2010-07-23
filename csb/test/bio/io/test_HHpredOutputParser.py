"""
@note: Naming convention for test files is slightly different. Keep test method names in camelCase
       to achieve conformance with the PEP's "respect local environment" recommendation. 
       
@todo: Write tests for alignment parsing and data contained in the hits.      
"""

import unittest
from csb.bio.io.hhpred import HHpredOutputParser

class TestHHpredOutputParser(unittest.TestCase):

    def setUp(self):
        
        filename = 'data/test_HHpredOutputParser.hhr'
        content = open(filename).read()
        tmp = HHpredOutputParser()
        
        self.hitlist = tmp.parse_file(filename)
        self.hitlist2 = tmp.parse_string(content)

    def testParseFile(self):
        """
        Test HHpredOutputParser.parse_file()
        """
        self.assertEqual(len(self.hitlist), 10)     
        self.assertEqual(self.hitlist.command, './hhsearch -i /ebio/abt1/kopec/projects/hhpred_on_cluster/split_fastas/d1ea0a1.hhm -d /ebio/abt1/kopec/projects/hhpred_on_cluster/astral_scop70_all_beta_only_v1.75_hhpred_database.hhm')
        self.assertEqual(self.hitlist.query_name, 'd1ea0a1 b.80.4.1 (A:1203-1472) Alpha subunit of glutamate synthase, C-terminal domain {Azospirillum brasilense [TaxId: 192]}')
        self.assertEqual(self.hitlist.neff, float('5.2'))
        self.assertEqual(self.hitlist.searched_hmms, 2936)   
        self.assertRaises(IndexError, self.hitlist.__getitem__, 20)

    def testParseString(self):
        """
        Test HHpredOutputParser.parse_string()
        """
        self.assertEqual(len(self.hitlist2), 10)
        self.assertEqual(self.hitlist2.command, './hhsearch -i /ebio/abt1/kopec/projects/hhpred_on_cluster/split_fastas/d1ea0a1.hhm -d /ebio/abt1/kopec/projects/hhpred_on_cluster/astral_scop70_all_beta_only_v1.75_hhpred_database.hhm')
        self.assertEqual(self.hitlist2.query_name, 'd1ea0a1 b.80.4.1 (A:1203-1472) Alpha subunit of glutamate synthase, C-terminal domain {Azospirillum brasilense [TaxId: 192]}')
        self.assertEqual(self.hitlist2.neff, float('5.2'))
        self.assertEqual(self.hitlist2.searched_hmms, 2936)
        self.assertRaises(IndexError, self.hitlist2.__getitem__, 20)

if __name__ == '__main__':
    
    unittest.main()