import unittest
import csb.test as test

from csb.bio.io.wwpdb import StructureParser, RegularStructureParser,\
                             LegacyStructureParser, get
from csb.bio.sequence import SequenceAlphabets
from csb.bio.structure import ChemElements

@test.unit
class TestStrutureParser(test.Case):

    def setUp(self):
        
        super(TestStrutureParser, self).setUp()
        
        regular_file = self.config.getTestFile('1d3z.regular.pdb')
        legacy_file = self.config.getTestFile('1d3z.legacy.pdb')
        
        self.rp = StructureParser(regular_file)
        self.lp = StructureParser(legacy_file)

    def testFactory(self):
        
        self.assertTrue(isinstance(self.rp, RegularStructureParser))
        self.assertTrue(isinstance(self.lp, LegacyStructureParser))
        
    def testParseModels(self):
        
        ensemble = self.rp.parse_models()
        self.assertEquals(ensemble.models.length, 10)
        self.assertEquals(ensemble[0].model_id, 1)
        self.assertEquals(ensemble.models[1].model_id, 1)
        
    def testLegacyParseModels(self):
        
        ensemble = self.lp.parse_models()
        self.assertEquals(ensemble.models.length, 10)
        self.assertEquals(ensemble[0].model_id, 1)
        self.assertEquals(ensemble.models[1].model_id, 1)        

    def testParseStructure(self):

        structure = self.rp.parse(model=2)
        
        self.assertEquals(self.rp.parse_structure().model_id, 1)        

        self.assertEqual(structure.accession, '1d3z')
        self.assertEqual(structure.model_id, 2)
        
        # Chain level
        self.assertEqual(structure.chains.length, 1)
        self.assertEqual(len(structure.chains), 1)
        self.assertEqual(structure.chains['A'].sequence, 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG')        

        self.assertEqual(len(structure.chains['A']), 76)
        self.assertEqual(len(structure['A']), 76)

        # Residue level 
        self.assertEqual(len(structure['A'][1:10]), 9)
        self.assertEqual(structure['A'][0].type, SequenceAlphabets.Protein.MET)                     #@UndefinedVariable
        
    def testLegacyParseStructure(self):
        
        structure = self.lp.parse(model=1)
        
        self.assertEquals(self.lp.parse_structure().model_id, 1)        

        self.assertEqual(structure.accession, '1d3z')
        self.assertEqual(structure.model_id, 1)
        
        # Chain level        
        self.assertEqual(structure.chains.length, 1)
        self.assertEqual(len(structure.chains), 1)
        self.assertEqual(structure.first_chain.molecule_id, '1') 
        self.assertEqual(structure.chains['A'].sequence, 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG')       

        self.assertEqual(len(structure.chains['A']), 76)
        self.assertEqual(len(structure['A']), 76)

        # Residue level 
        self.assertEqual(len(structure['A'][1:10]), 9)
        self.assertEqual(structure['A'][0].type, SequenceAlphabets.Protein.MET)             #@UndefinedVariable
        
        # Atom level
        self.assertEqual(structure['A'][0].atoms['CA'].element, None)
        self.assertNotEqual(structure['A'][2].atoms['CA'].element, None)
        self.assertEqual(structure['A'][2].atoms['CA'].element, ChemElements.C)             #@UndefinedVariable

    def testGet(self):
        structure = get('1d3z')
        # Chain level
        self.assertEqual(structure.accession, '1d3z')
        # Chain level
        self.assertEqual(structure.chains.length, 1)
        self.assertEqual(len(structure.chains), 1)

        self.assertEqual(len(structure.chains['A']), 76)
        self.assertEqual(len(structure['A']), 76)

        # Residue level 
        self.assertEqual(len(structure['A'][1:10]), 9)
        self.assertEqual(structure['A'][0].type,SequenceAlphabets.Protein.MET)              #@UndefinedVariable


class PDBTestCase(test.Case):
    
    def runTest(self):
        
        try:
            StructureParser(self.entry).parse_structure()
        except:
            self.reRaise([self.entry])
        
@test.custom
def TestPDB():
    
    import glob
    
    suite = unittest.TestSuite()
    
    for entry in glob.glob('/media/DB/pdb/all/pdb1xx[a-f]*.ent'):
        
        case = PDBTestCase()
        case.entry = entry
        suite.addTest(case)
    
    return suite


if __name__ == '__main__':
    
    test.Console()
