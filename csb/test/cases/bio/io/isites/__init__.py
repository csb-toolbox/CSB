import csb.test as test
from csb.bio.io.isites import ISitesParser


@test.functional
class TestISitesParser(test.Case):
    
    def setUp(self):
        
        super(TestISitesParser, self).setUp()
        
        self.file = self.config.getTestFile('ISL5.1.isl')
        self.parser = ISitesParser(self.file)
        
    def testParseAll(self):
        
        lib = self.parser.parseall()
        
        self.assertEqual(len(lib.clusters), 145)
        self.assertEqual(lib[11007], lib['7pcy'][0])
        
        c = lib[11007]
        self.assertEqual(c.motiflen, 11)
        self.assertEqual(c.profilelen, 15)
        self.assertEqual(c.overhang, 2)
        self.assertEqual(c.mda, 96.00)
        
        self.assertEqual(c.representative.accession, '7pcy')
        self.assertEqual(c.representative.angles.length, c.motiflen)
        self.assertEqual(c.representative.angles[1].phi, -149.42)
        self.assertEqual(c.representative.angles[1].psi, 151.65)
        self.assertEqual(c.representative.angles[1].omega, 176.37)
                
        

if __name__ == '__main__':
    
    test.Console()