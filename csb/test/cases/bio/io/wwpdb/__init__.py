import unittest
import csb.test as test

from csb.bio.io.wwpdb import StructureParser, RegularStructureParser,\
                             LegacyStructureParser, get, UnknownPDBResidueError
from csb.bio.sequence import SequenceAlphabets, SequenceTypes
from csb.bio.structure import ChemElements


@test.regression
class TestMappingRegressions(test.Case):
        
    def testHetMapping(self):
        """
        @see: [CSB 0000031]
        """
        pdbfile = self.config.getTestFile('1d3z.regular.pdb')
        structure = StructureParser(pdbfile).parse_structure(1)
                
        residue = structure.chains['A'].find(26)
        self.assertTrue(residue.has_structure)
        self.assertTrue(residue.atoms.length > 0)
        
        for an in residue.atoms:
            self.assertTrue(residue[an]._het)
            self.assertTrue(residue[an].vector.tolist())

    def testNonStandardResidueMapping(self):
        """
        @see: [CSB 0000052]
        """
        pdbfile = self.config.getTestFile('3p1u.pdb')
        chain = StructureParser(pdbfile).parse_structure().chains['A']
        
        for residue in chain.residues:
            
            if residue.rank < 15:
                self.assertFalse(residue.has_structure)
                
            elif residue.rank in (15, 16):
                self.assertTrue(residue.has_structure)
                
        self.assertEqual(chain.residues[2].sequence_number, None)
        self.assertEqual(chain.residues[15].sequence_number, 39)
        self.assertEqual(chain.residues[16].sequence_number, 40)

@test.unit
class TestStructureParser(test.Case):

    def setUp(self):
        
        super(TestStructureParser, self).setUp()
        
        regular_file = self.config.getTestFile('1d3z.regular.pdb')
        legacy_file = self.config.getTestFile('1d3z.legacy.pdb')
        
        self.rp = StructureParser(regular_file)
        self.lp = StructureParser(legacy_file)

    def testFactory(self):
        
        self.assertTrue(isinstance(self.rp, RegularStructureParser))
        self.assertTrue(isinstance(self.lp, LegacyStructureParser))
        
        
@test.unit
class TestLegacyStructureParser(test.Case):

    def setUp(self):
        
        super(TestLegacyStructureParser, self).setUp()
        
        self.pdb = self.config.getTestFile('1d3z.legacy.pdb')
        self.parser = LegacyStructureParser(self.pdb)
        
    def testParseModels(self):
        
        ensemble = self.parser.parse_models()
        self.assertEquals(ensemble.models.length, 10)
        self.assertEquals(ensemble[0].model_id, 1)
        self.assertEquals(ensemble.models[1].model_id, 1)        
        
    def testParseStructure(self):
        
        structure = self.parser.parse(model=1)
        
        self.assertEquals(self.parser.parse_structure().model_id, 1)        

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

        vector = [51.653, -89.304, 8.833]
        self.assertEqual(structure['A'][0]['CA'].vector.tolist(), vector)        

    def testParseResidue(self):
        
        self.assertEqual(self.parser.parse_residue('AGM'), SequenceAlphabets.Protein.ARG.name)                                  #@UndefinedVariable
        self.assertEqual(self.parser.parse_residue('AGM', as_type=SequenceTypes.Protein), SequenceAlphabets.Protein.ARG.name)   #@UndefinedVariable        
        self.assertRaises(UnknownPDBResidueError, self.parser.parse_residue, 'AGM', as_type=SequenceTypes.NucleicAcid)          #@UndefinedVariable                
    
    def testParseResidueSafe(self):
        
        self.assertEqual(self.parser.parse_residue_safe('AGM', as_type=None), SequenceAlphabets.Protein.ARG.name)                      #@UndefinedVariable
        self.assertEqual(self.parser.parse_residue_safe('AGM', as_type=SequenceTypes.Protein), SequenceAlphabets.Protein.ARG.name)     #@UndefinedVariable
        self.assertEqual(self.parser.parse_residue_safe('AGM', as_type=SequenceTypes.NucleicAcid), SequenceAlphabets.Nucleic.Any.name) #@UndefinedVariable                
        self.assertEqual(self.parser.parse_residue_safe('junk', as_type=SequenceTypes.Protein), SequenceAlphabets.Unknown.UNK.name)    #@UndefinedVariable 
    
    def testGuessSequenceType(self):
        
        self.assertEquals(self.parser.guess_sequence_type('AGM'), SequenceTypes.Protein)                                        #@UndefinedVariable
        self.assertEquals(self.parser.guess_sequence_type('DOC'), SequenceTypes.NucleicAcid)                                    #@UndefinedVariable  
        self.assertRaises(UnknownPDBResidueError, self.parser.guess_sequence_type, 'junk')
        
    def testFileName(self):
        self.assertEqual(self.parser.filename, self.pdb)
    
    def testModels(self):
        self.assertEqual(self.parser.models(), list(range(1, 11)))
        
        
@test.unit
class TestRegularStructureParser(test.Case):

    def setUp(self):
        
        super(TestRegularStructureParser, self).setUp()
        
        self.pdb = self.config.getTestFile('1d3z.regular.pdb')
        self.parser = RegularStructureParser(self.pdb)
        
    def testParseModels(self):
        
        ensemble = self.parser.parse_models()
        self.assertEquals(ensemble.models.length, 10)
        self.assertEquals(ensemble[0].model_id, 1)
        self.assertEquals(ensemble.models[1].model_id, 1)       

    def testParseStructure(self):

        structure = self.parser.parse(model=2)
        
        self.assertEquals(self.parser.parse_structure().model_id, 1)        

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
        self.assertEqual(structure['A'][0].type, SequenceAlphabets.Protein.MET)                                                 #@UndefinedVariable
        self.assertEqual(structure['A'][0].secondary_structure, None)
        
        # Atom
        vector = [52.647, -87.443, 9.674]
        self.assertEqual(structure['A'][0]['CA'].vector.tolist(), vector)
        self.assertEqual(structure['A'][0]['CA']._het, True) 

    def testParseResidue(self):
        
        self.assertEqual(self.parser.parse_residue('AGM'), SequenceAlphabets.Protein.ARG.name)                                  #@UndefinedVariable
        self.assertEqual(self.parser.parse_residue('AGM', as_type=SequenceTypes.Protein), SequenceAlphabets.Protein.ARG.name)   #@UndefinedVariable        
        self.assertRaises(UnknownPDBResidueError, self.parser.parse_residue, 'AGM', as_type=SequenceTypes.NucleicAcid)          #@UndefinedVariable                
    
    def testParseResidueSafe(self):
        
        self.assertEqual(self.parser.parse_residue_safe('AGM', as_type=None), SequenceAlphabets.Protein.ARG.name)                      #@UndefinedVariable
        self.assertEqual(self.parser.parse_residue_safe('AGM', as_type=SequenceTypes.Protein), SequenceAlphabets.Protein.ARG.name)     #@UndefinedVariable
        self.assertEqual(self.parser.parse_residue_safe('AGM', as_type=SequenceTypes.NucleicAcid), SequenceAlphabets.Nucleic.Any.name) #@UndefinedVariable                
        self.assertEqual(self.parser.parse_residue_safe('junk', as_type=SequenceTypes.Protein), SequenceAlphabets.Unknown.UNK.name)    #@UndefinedVariable 
    
    def testGuessSequenceType(self):
        
        self.assertEquals(self.parser.guess_sequence_type('AGM'), SequenceTypes.Protein)                                        #@UndefinedVariable
        self.assertEquals(self.parser.guess_sequence_type('DOC'), SequenceTypes.NucleicAcid)                                    #@UndefinedVariable  
        self.assertRaises(UnknownPDBResidueError, self.parser.guess_sequence_type, 'junk')
        
    def testFileName(self):
        self.assertEqual(self.parser.filename, self.pdb)
    
    def testModels(self):
        self.assertEqual(self.parser.models(), list(range(1, 11)))


@test.functional
class TestGet(test.Case):
    
    def runTest(self):
        
        structure = get('1d3z')
        self.assertEqual(structure.accession, '1d3z')
        
        # Chain level
        self.assertEqual(structure.chains.length, 1)
        self.assertEqual(len(structure.chains), 1)

        self.assertEqual(len(structure.chains['A']), 76)
        self.assertEqual(len(structure['A']), 76)

        # Residue level 
        self.assertEqual(len(structure['A'][1:10]), 9)
        self.assertEqual(structure['A'][0].type,SequenceAlphabets.Protein.MET)              #@UndefinedVariable        
       
       
@test.custom
def TestPDB():
    
    import glob

    mask = raw_input('PDB file(s) (enter file name or mask): ')
    suite = unittest.TestSuite()    
    
    class PDBTestCase(test.Case):
        
        def runTest(self):
        
            try:
                StructureParser(self.entry).parse_structure()
            except:
                self.reRaise([self.entry])    
        
    for entry in glob.glob(mask):
        
        case = PDBTestCase()
        case.entry = entry
        suite.addTest(case)
    
    return suite


if __name__ == '__main__':
    
    test.Console()
