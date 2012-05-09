import os
import unittest
import csb.test as test

from csb.bio.io.wwpdb import EntryID, StandardID, DegenerateID, SeqResID, InvalidEntryIDError
from csb.bio.io.wwpdb import StructureParser, RegularStructureParser, LegacyStructureParser, UnknownPDBResidueError
from csb.bio.io.wwpdb import get, find, FileSystemStructureProvider, RemoteStructureProvider, CustomStructureProvider, StructureNotFoundError
from csb.bio.sequence import SequenceAlphabets, SequenceTypes
from csb.bio.structure import ChemElements, SecStructures, Structure

@test.regression
class TestBiomoleculeRegressions(test.Case):

    def testCommaSplitting(self):
        """
        @see: [CSB 0000067]
        """
        pdbfile = self.config.getTestFile('3shm_ca.pdb')
        parser = LegacyStructureParser(pdbfile)

        s1 = parser.parse_biomolecule(1, True)

        self.assertEqual(len(s1.chains), 60)
        self.assertEqual(s1.first_chain.id, 'A')

@test.regression
class TestSecStructureRegressions(test.Case):
        
    def testSecStructureParsing(self):
        """
        @see: [CSB 0000045]
        """
        pdbfile = self.config.getTestFile('1d3z.regular.pdb')
        structure = StructureParser(pdbfile).parse_structure(1)
        
        self.assertTrue(structure.chains['A'].secondary_structure is not None)

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
        self.assertEqual(structure['A'][0].type, SequenceAlphabets.Protein.MET)             
        self.assertEqual(structure['A'][0]._pdb_name, 'MSE')
        self.assertEqual(structure['A'][1]._pdb_name, 'GLN')        
        
        # Atom level
        self.assertEqual(structure['A'][1].atoms['CA'].element, None)
        self.assertNotEqual(structure['A'][2].atoms['CA'].element, None)
        self.assertEqual(structure['A'][2].atoms['CA'].element, ChemElements.C)             

        vector = [51.653, -89.304, 8.833]
        self.assertEqual(structure['A'][0]['CA'].vector.tolist(), vector)        

    def testParseResidue(self):
        
        self.assertEqual(self.parser.parse_residue('AGM'), SequenceAlphabets.Protein.ARG.name)                                  #@UndefinedVariable
        self.assertEqual(self.parser.parse_residue('AGM', as_type=SequenceTypes.Protein), SequenceAlphabets.Protein.ARG.name)   #@UndefinedVariable        
        self.assertRaises(UnknownPDBResidueError, self.parser.parse_residue, 'AGM', as_type=SequenceTypes.NucleicAcid)                          
    
    def testParseResidueSafe(self):
        
        self.assertEqual(self.parser.parse_residue_safe('AGM', as_type=None), SequenceAlphabets.Protein.ARG.name)                      #@UndefinedVariable
        self.assertEqual(self.parser.parse_residue_safe('AGM', as_type=SequenceTypes.Protein), SequenceAlphabets.Protein.ARG.name)     #@UndefinedVariable
        self.assertEqual(self.parser.parse_residue_safe('AGM', as_type=SequenceTypes.NucleicAcid), SequenceAlphabets.Nucleic.Any.name) #@UndefinedVariable                
        self.assertEqual(self.parser.parse_residue_safe('junk', as_type=SequenceTypes.Protein), SequenceAlphabets.Unknown.UNK.name)    #@UndefinedVariable 
    
    def testGuessSequenceType(self):
        
        self.assertEquals(self.parser.guess_sequence_type('AGM'), SequenceTypes.Protein)                                        
        self.assertEquals(self.parser.guess_sequence_type('DOC'), SequenceTypes.NucleicAcid)                                      
        self.assertRaises(UnknownPDBResidueError, self.parser.guess_sequence_type, 'junk')
        
    def testFileName(self):
        self.assertEqual(self.parser.filename, self.pdb)
    
    def testModels(self):
        self.assertEqual(self.parser.models(), list(range(1, 11)))
        
    def testParseBiomolecule(self):

        pdbfile = self.config.getTestFile('3p1u.pdb')
        parser = LegacyStructureParser(pdbfile)

        s2 = parser.parse_biomolecule(2)

        self.assertEqual(len(s2.chains), 1)
        self.assertEqual(s2.first_chain.id, 'B1')
        self.assertRaises(KeyError, parser.parse_biomolecule, 3)

        
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
        
        ss = structure.chains['A'].secondary_structure
        self.assertEqual(ss.to_string(), '-EEEEE-----EEEEE-----HHHHHHHHHHHHHH-HHH-EEEEE--EE------HHHHH-----EEEEEE')
        self.assertEqual(len(ss.scan(1, 99, filter=SecStructures.Helix)), 3)
        self.assertEqual(ss[1].start, 2)
        self.assertEqual(ss[1].end, 6)
        
        self.assertEqual(len(structure.chains['A']), 76)
        self.assertEqual(len(structure['A']), 76)

        # Residue level 
        self.assertEqual(len(structure['A'][1:10]), 9)
        self.assertEqual(structure['A'][0].type, SequenceAlphabets.Protein.MET)                                                 
        self.assertEqual(structure['A'][0]._pdb_name, 'MSE')
        self.assertEqual(structure['A'][1]._pdb_name, 'GLN')
        
        # Atom
        vector = [52.647, -87.443, 9.674]
        self.assertEqual(structure['A'][0]['CA'].vector.tolist(), vector)
        self.assertEqual(structure['A'][0]['CA']._het, True) 

    def testParseResidue(self):
        
        self.assertEqual(self.parser.parse_residue('AGM'), SequenceAlphabets.Protein.ARG.name)                                  #@UndefinedVariable
        self.assertEqual(self.parser.parse_residue('AGM', as_type=SequenceTypes.Protein), SequenceAlphabets.Protein.ARG.name)   #@UndefinedVariable        
        self.assertRaises(UnknownPDBResidueError, self.parser.parse_residue, 'AGM', as_type=SequenceTypes.NucleicAcid)                          
    
    def testParseResidueSafe(self):
        
        self.assertEqual(self.parser.parse_residue_safe('AGM', as_type=None), SequenceAlphabets.Protein.ARG.name)                      #@UndefinedVariable
        self.assertEqual(self.parser.parse_residue_safe('AGM', as_type=SequenceTypes.Protein), SequenceAlphabets.Protein.ARG.name)     #@UndefinedVariable
        self.assertEqual(self.parser.parse_residue_safe('AGM', as_type=SequenceTypes.NucleicAcid), SequenceAlphabets.Nucleic.Any.name) #@UndefinedVariable                
        self.assertEqual(self.parser.parse_residue_safe('junk', as_type=SequenceTypes.Protein), SequenceAlphabets.Unknown.UNK.name)    #@UndefinedVariable 
    
    def testGuessSequenceType(self):
        
        self.assertEquals(self.parser.guess_sequence_type('AGM'), SequenceTypes.Protein)                                        
        self.assertEquals(self.parser.guess_sequence_type('DOC'), SequenceTypes.NucleicAcid)                                      
        self.assertRaises(UnknownPDBResidueError, self.parser.guess_sequence_type, 'junk')
        
    def testFileName(self):
        self.assertEqual(self.parser.filename, self.pdb)
    
    def testModels(self):
        self.assertEqual(self.parser.models(), list(range(1, 11)))


@test.unit
class TestFileSystemProvider(test.Case):
    
    def setUp(self):
        
        super(TestFileSystemProvider, self).setUp()
        
        self.path = self.config.data
        self.provider = FileSystemStructureProvider()
        self.provider.add(self.path)
    
    def testAdd(self):
        
        self.assertEqual(len(self.provider.paths), 1)
        
        self.provider.add('.')
        self.assertEqual(self.provider.paths[1], '.')
        self.assertEqual(len(self.provider.paths), 2)
        
        self.assertRaises(IOError, self.provider.add, 'non-exi$ting path')
    
    def testRemove(self):
        
        self.assertEqual(len(self.provider.paths), 1)
        self.provider.remove(self.path)
        self.assertEqual(len(self.provider.paths), 0)
        
        self.assertRaises(ValueError, self.provider.remove, 'non-exi$ting path')
    
    def testFind(self):
        
        f1 = self.provider.find('3p1u')
        f2 = self.config.getTestFile('3p1u.pdb')
        
        self.assertEqual(os.path.abspath(f1), os.path.abspath(f2))
        self.assertEqual(None, self.provider.find('$'))
        
    def testGet(self):
        
        s = self.provider.get('3p1u')
        self.assertEqual(s.accession, '3p1u')
        self.assertTrue(isinstance(s, Structure))
        self.assertRaises(StructureNotFoundError, self.provider.get, '$')


@test.unit
class TestCustomProvider(test.Case):
    
    def setUp(self):
        
        super(TestCustomProvider, self).setUp()
        
        self.path = self.config.getTestFile('3p1u.pdb')
        self.provider = CustomStructureProvider({'3p1u': self.path})
    
    def testAdd(self):
        
        self.assertEqual(len(self.provider.paths), 1)
        self.assertEqual(self.provider.paths[0], self.path)
                
        self.provider.add('test', self.config.getTestFile('d1nz0a_.pdb'))
        self.assertEqual(len(self.provider.paths), 2)
        
        self.assertRaises(IOError, self.provider.add, 'test', 'non-exi$ting path')
    
    def testRemove(self):
        
        self.assertEqual(len(self.provider.paths), 1)
        self.provider.remove('3p1u')
        self.assertEqual(len(self.provider.paths), 0)
        
        self.assertRaises(ValueError, self.provider.remove, '$')
    
    def testFind(self):
        
        f1 = self.provider.find('3p1u')
        f2 = self.config.getTestFile('3p1u.pdb')
        
        self.assertEqual(os.path.abspath(f1), os.path.abspath(f2))
        self.assertEqual(None, self.provider.find('$'))
        
    def testGet(self):
        
        s = self.provider.get('3p1u')
        self.assertEqual(s.accession, '3p1u')
        self.assertTrue(isinstance(s, Structure))
        self.assertRaises(StructureNotFoundError, self.provider.get, '$')
        
@test.unit
class TestRemoteProvider(test.Case):
    
    def setUp(self):
        
        super(TestRemoteProvider, self).setUp()
        self.provider = RemoteStructureProvider()
        
    def testGet(self):
        
        s = self.provider.get('3p1u')
        self.assertEqual(s.accession, '3p1u')
        self.assertTrue(isinstance(s, Structure))
        
        self.provider.prefix = 'http://NoSuchURL.test'
        self.assertRaises(StructureNotFoundError, self.provider.get, 'NoSuchFile')    
                

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
        self.assertEqual(structure['A'][0].type,SequenceAlphabets.Protein.MET)                      


@test.functional
class TestFind(test.Case):
    
    def runTest(self):
        
        f2 = find('3p1u', [self.config.data])
        f1 = self.config.getTestFile('3p1u.pdb')
        
        self.assertEqual(os.path.abspath(f1), os.path.abspath(f2))
        self.assertEqual(None, find('$', self.config.data))       


@test.unit
class TestEntryID(test.Case):

    def setUp(self):
        super(TestEntryID, self).setUp()
        
    def testFactory(self):
        self.assertTrue(isinstance(EntryID.create('abcdE'), StandardID))
        self.assertTrue(isinstance(EntryID.create('abcdddE'), DegenerateID))
        self.assertTrue(isinstance(EntryID.create('abcd_E'), SeqResID))
        
@test.unit
class TestStandardID(test.Case):

    def setUp(self):
        super(TestStandardID, self).setUp()
        
        self.id = StandardID('abCdE')
        self.accession = 'abcd'
        self.chain = 'E'
    
    def testAccession(self):
        self.assertEqual(self.id.accession, self.accession)

    def testChain(self):
        self.assertEqual(self.id.chain, self.chain)

    def testEntryID(self):
        self.assertEqual(self.id.entry_id, self.accession + self.chain)                       

    def testFormat(self):
        self.assertEqual(self.id.format(), self.accession + self.chain)
        
    def testConstructor(self):
        self.assertRaises(InvalidEntryIDError, StandardID, 'aE')

@test.unit
class TestDegenerateID(TestStandardID):

    def setUp(self):
        super(TestDegenerateID, self).setUp()
        
        self.id = DegenerateID('abCdddE')
        self.accession = 'abcddd'
        self.chain = 'E'  
        
@test.unit
class TestSeqResID(TestStandardID):

    def setUp(self):
        super(TestSeqResID, self).setUp()
        self.id = SeqResID('abCd_E')

    def testFormat(self):
        self.assertEqual(self.id.format(), 'abcd_E')
        
    def testConstructor(self):
        self.assertRaises(InvalidEntryIDError, SeqResID, 'abcdE')
                
                               
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
