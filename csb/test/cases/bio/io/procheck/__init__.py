import csb.test as test

from csb.bio.io.procheck import ProcheckParser

@test.functional
class TestProcheckParser(test.Case):

        
    def setUp(self):
        
        super(TestProcheckParser, self).setUp()
        self.file = self.config.getTestFile('2JZC.sum')
        self.parser =  ProcheckParser()

    def testParse(self):

        res = self.parser.parse(self.file)

        self.assertEqual(res['#residues'], 201)
        self.assertEqual(res['rama_core'], 69.5)
        self.assertEqual(res['rama_allow'], 22.6)
        self.assertEqual(res['rama_gener'], 5.6)
        self.assertEqual(res['rama_disall'], 2.3)

        self.assertEqual(res['g_dihedrals'], -0.1)
        self.assertEqual(res['g_bond'], 0.51)
        self.assertEqual(res['g_overall'], 0.14)

        self.assertEqual(res['badContacts'], 5581)
        

if __name__ == '__main__':
    
    test.Console()



        

