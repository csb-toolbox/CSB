import csb.test as test

from csb.bio.io.whatif import WhatCheckParser


@test.functional
class TestWhatCheckParser(test.Case):


    def setUp(self):
        super(TestWhatCheckParser, self).setUp()
        self.file = self.config.getTestFile('pdbout.txt')
        self.parser = WhatCheckParser()


    def testParse(self):

        res = self.parser.parse(self.file)
        self.assertEqual(res['rama_z_score'], -4.617)
        self.assertEqual(res['bb_z_score'], -1.421)
        self.assertEqual(res['1st_packing_z_score'], -3.436)
        self.assertEqual(res['2nd_packing_z_score'], -4.424)
        self.assertEqual(res['rotamer_score'], -2.103)
        

if __name__ == '__main__':
    
    test.Console()

        
