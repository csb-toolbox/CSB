import unittest
from csb.bio.io.wwpdb import StructureParser, RegularStructureParser,\
                             LegacyStructureParser, get
from csb.bio.sequence import SequenceAlphabets
from csb.bio.structure import ChemElements

class TestStrutureParser(unittest.TestCase):

    def setUp(self):
        
        regular_file = 'data/test_pdb.pdb'
        legacy_file = 'data/test_legacy.pdb'
        
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
        self.assertEqual(structure['A'][0].type, SequenceAlphabets.Protein.MET)
        
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
        self.assertEqual(structure['A'][0].type, SequenceAlphabets.Protein.MET)     
        
        # Atom level
        self.assertEqual(structure['A'][0].atoms['CA'].element, None)
        self.assertNotEqual(structure['A'][2].atoms['CA'].element, None)
        self.assertEqual(structure['A'][2].atoms['CA'].element, ChemElements.C)      

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
        self.assertEqual(structure['A'][0].type,SequenceAlphabets.Protein.MET)


if __name__ == '__main__':
    unittest.main()
