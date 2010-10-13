"""
@todo: Write tests for alignment parsing and data contained in the hits.
@todo: HHpredProfileParser tests.  
"""

import csb.test as test

from csb.bio.io.hhpred import HHpredOutputParser

@test.unit
class TestHHpredOutputParser(test.Case):

    def setUp(self):
        
        super(TestHHpredOutputParser, self).setUp()
        
        filename = self.config.getTestFile('d1ea0a1.hhr')
        content = open(filename).read()
        tmp = HHpredOutputParser(True)
        
        self.hitlist = tmp.parse_file(filename)
        self.hitlist2 = tmp.parse_string(content)

    def testParseFile(self):

        self.assertEqual(len(self.hitlist), 10)     
        self.assertEqual(self.hitlist.command, './hhsearch -i /ebio/abt1/kopec/projects/hhpred_on_cluster/split_fastas/d1ea0a1.hhm -d /ebio/abt1/kopec/projects/hhpred_on_cluster/astral_scop70_all_beta_only_v1.75_hhpred_database.hhm')
        self.assertEqual(self.hitlist.query_name, 'd1ea0a1 b.80.4.1 (A:1203-1472) Alpha subunit of glutamate synthase, C-terminal domain {Azospirillum brasilense [TaxId: 192]}')
        self.assertEqual(self.hitlist.neff, float('5.2'))
        self.assertEqual(self.hitlist.searched_hmms, 2936)   
        self.assertRaises(IndexError, self.hitlist.__getitem__, 20)
        
        segments = list(self.hitlist[9].alignment.segments)
        self.assertEqual(len(segments), 3)
        
        for s, l in zip(segments, [(9, 21), (22, 35), (37, 52)]):
            self.assertEqual(s.start, l[0])
            self.assertEqual(s.end, l[1])
            self.assertEqual(s.end - s.start, l[1] - l[0])
            
    def testParseString(self):

        self.assertEqual(len(self.hitlist2), 10)
        self.assertEqual(self.hitlist2.command, './hhsearch -i /ebio/abt1/kopec/projects/hhpred_on_cluster/split_fastas/d1ea0a1.hhm -d /ebio/abt1/kopec/projects/hhpred_on_cluster/astral_scop70_all_beta_only_v1.75_hhpred_database.hhm')
        self.assertEqual(self.hitlist2.query_name, 'd1ea0a1 b.80.4.1 (A:1203-1472) Alpha subunit of glutamate synthase, C-terminal domain {Azospirillum brasilense [TaxId: 192]}')
        self.assertEqual(self.hitlist2.neff, float('5.2'))
        self.assertEqual(self.hitlist2.searched_hmms, 2936)
        self.assertRaises(IndexError, self.hitlist2.__getitem__, 20)


if __name__ == '__main__':
    
    test.Console()
