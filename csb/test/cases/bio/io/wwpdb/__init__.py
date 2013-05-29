import os
import sys
import csb.test as test

from csb.bio.io.wwpdb import EntryID, StandardID, DegenerateID, SeqResID, InvalidEntryIDError, HeaderFormatError
from csb.bio.io.wwpdb import RobustResidueMapper, FastResidueMapper, CombinedResidueMapper, ResidueMappingError, SparseChainSequence
from csb.bio.io.wwpdb import StructureParser, RegularStructureParser, LegacyStructureParser, UnknownPDBResidueError
from csb.bio.io.wwpdb import get, find, FileSystemStructureProvider, RemoteStructureProvider, CustomStructureProvider, StructureNotFoundError
from csb.bio.sequence import SequenceAlphabets, ProteinAlphabet, SequenceTypes, RichSequence
from csb.bio.structure import ChemElements, SecStructures, Structure, Chain


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
        self.assertEqual(ensemble.models.length, 10)
        self.assertEqual(ensemble[0].model_id, 1)
        self.assertEqual(ensemble.models[1].model_id, 1)        
        
    def testParseStructure(self):
        
        structure = self.parser.parse(model=1)
        
        self.assertEqual(self.parser.parse_structure().model_id, 1)        

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
        self.assertEqual(structure['A'][0].label, 'MSE')
        self.assertEqual(structure['A'][1].label, 'GLN')
        self.assertTrue(structure['A'][0].is_modified)
        self.assertFalse(structure['A'][1].is_modified)   
        
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
        
        self.assertEqual(self.parser.guess_sequence_type('AGM'), SequenceTypes.Protein)                                        
        self.assertEqual(self.parser.guess_sequence_type('DOC'), SequenceTypes.NucleicAcid)                                      
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
        
    def testParseHetMolecules(self):
        
        with self.config.getTempStream() as tmp:
            
            tmp.write('HETATM    1  NA  BLM A   1     -14.575  27.241   3.310  1.00  0.00           N ')
            tmp.flush()
            
            parser = LegacyStructureParser(tmp.name)
            self.assertRaises(HeaderFormatError, parser.parse_structure)
            del parser

        
@test.unit
class TestRegularStructureParser(test.Case):

    def setUp(self):
        
        super(TestRegularStructureParser, self).setUp()
        
        self.pdb = self.config.getTestFile('1d3z.regular.pdb')
        self.mapping = self.config.getTestFile('mapping.pdb')
        self.parser = RegularStructureParser(self.pdb)
        
    def testMapper(self):
        
        p = RegularStructureParser(self.pdb, mapper=None)
        self.assertTrue(isinstance(p.mapper, CombinedResidueMapper))
        
        p.mapper = FastResidueMapper()
        self.assertTrue(isinstance(p.mapper, FastResidueMapper))
        
    def testCombinedMapping(self):

        # default mapper        
        c = self.parser.parse(self.mapping)['E']
        self.assertEqual(c.residues[14].type, ProteinAlphabet.GLU)
        self.assertEqual(c.residues[15].type, ProteinAlphabet.GLU)
        self.assertEqual(c.residues[16].type, ProteinAlphabet.THR)
        self.assertEqual(4, sum([1 for r in c if r.has_structure]))

        # explicit combined mapper        
        self.parser.mapper = CombinedResidueMapper()
        c = self.parser.parse(self.mapping)['E']
        self.assertEqual(4, sum([1 for r in c if r.has_structure]))

    def testFastMapping(self):
           
        self.parser.mapper = FastResidueMapper()
        self.assertRaises(ResidueMappingError, self.parser.parse, self.mapping)
        
        mapping2 = self.config.getTestFile('mapping2.pdb')
        
        c = self.parser.parse(mapping2)['E']
        self.assertEqual(2, sum([1 for r in c if r.has_structure]))
        
    def testRobustMapping(self):

        mapping3 = self.config.getTestFile('mapping3.pdb')
        
        self.parser.mapper = RobustResidueMapper()
        self.assertRaises(ResidueMappingError, self.parser.parse, mapping3)
        
        c = self.parser.parse(self.mapping)['E']                
        self.assertEqual(4, sum([1 for r in c if r.has_structure]))
                        
    def testParseModels(self):
        
        ensemble = self.parser.parse_models()
        self.assertEqual(ensemble.models.length, 10)
        self.assertEqual(ensemble[0].model_id, 1)
        self.assertEqual(ensemble.models[1].model_id, 1)
        
        self.assertRaises(ValueError, self.parser.parse_models, (999, 1000))
        
        pdb = self.config.getTestFile('3p1u.pdb')
        ensemble = RegularStructureParser(pdb).parse_models()
        self.assertEqual(ensemble.models.length, 1)    
        self.assertEqual(ensemble[0].model_id, 1)
        self.assertEqual(ensemble[0].resolution, 2.05)
    
    def testParseStructure(self):

        structure = self.parser.parse(model=2)
        
        self.assertEqual(self.parser.parse_structure().model_id, 1)
        self.assertEqual(structure.resolution, None)

        self.assertEqual(structure.accession, '1d3z')
        self.assertEqual(structure.model_id, 2)
        
        # Chain level
        self.assertEqual(structure.chains.length, 1)
        self.assertEqual(len(structure.chains), 1)
        self.assertEqual(structure.chains['A'].sequence, 'MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG')
        self.assertEqual(structure.chains['A'].sequence, ''.join([str(r.type) for r in structure.chains['A'] if r.has_structure]))        
        
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
        self.assertEqual(structure['A'][0].label, 'MSE')
        self.assertEqual(structure['A'][1].label, 'GLN')
        
        # Atom
        vector = [52.647, -87.443, 9.674]
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
        
        self.assertEqual(self.parser.guess_sequence_type('AGM'), SequenceTypes.Protein)                                        
        self.assertEqual(self.parser.guess_sequence_type('DOC'), SequenceTypes.NucleicAcid)                                      
        self.assertRaises(UnknownPDBResidueError, self.parser.guess_sequence_type, 'junk')
        
    def testFileName(self):
        self.assertEqual(self.parser.filename, self.pdb)
    
    def testModels(self):
        self.assertEqual(self.parser.models(), list(range(1, 11)))


@test.unit
class TestFastResidueMapper(test.Case):
    
    def setUp(self):
        super(TestFastResidueMapper, self).setUp()
        self.mapper = FastResidueMapper()
    
    def _build(self, string):
        
        id = str(hash(string))
        seq = RichSequence(id, "", string, SequenceTypes.Protein)
        
        return SparseChainSequence.create(Chain.from_sequence(seq))                

    def testMap(self):

        ref = self._build("ZABCD")        
        sparse = self._build("AC")
        
        self.assertRaises(ResidueMappingError, self.mapper.map, sparse, ref)
        
        sparse.residues[2].id = (22, None)
        result = self.mapper.map(sparse, ref)

        self.assertEqual(result.sequence, "-A-C-")
        
    def testModifiedResidueMapping(self):
        """
        Strictly speaking, this is a regression test. But it so essential that
        we should keep it here.
         
        @see: [csb: 19]
        """
        pdb = self.config.getTestFile('modified.pdb')
         
        structure = StructureParser(pdb, mapper=self.mapper).parse_structure()
        chain = structure.first_chain

        self.assertFalse(chain.residues[1].has_structure)
        self.assertEqual(chain.residues[1].label, "MET")
                        
        self.assertTrue(chain.residues[19].has_structure)
        self.assertEqual(chain.residues[19].label, "MSE")    
   
        
@test.unit
class TestRobustResidueMapper(TestFastResidueMapper):
    
    def setUp(self):
        super(TestRobustResidueMapper, self).setUp()
        self.mapper = RobustResidueMapper()              

    def testMap(self):

        ref = self._build("ABCD")        
        sparse = self._build("EF")
        
        self.assertRaises(ResidueMappingError, self.mapper.map, sparse, ref)
        
        ref = self._build("ZABCD")        
        sparse = self._build("AC")
        result = self.mapper.map(sparse, ref)

        self.assertEqual(result.sequence, "-A-C-")        
        
    def testModifiedResidueMapping(self):
        
        pdb = self.config.getTestFile('modified2.pdb')
         
        structure = StructureParser(pdb, mapper=self.mapper).parse_structure()
        chain = structure.first_chain

        self.assertTrue(chain.residues[1].has_structure)
        self.assertEqual(chain.residues[1].label, "MSE")
                        
        self.assertFalse(chain.residues[19].has_structure)
        self.assertEqual(chain.residues[19].label, "MET")           
        

@test.unit
class TestCombinedResidueMapper(TestFastResidueMapper):
    
    def setUp(self):
        super(TestCombinedResidueMapper, self).setUp()
        self.mapper = CombinedResidueMapper()              

    def testMap(self):

        ref = self._build("ZABCD")        
        sparse = self._build("AC")
        result = self.mapper.map(sparse, ref)

        self.assertEqual(result.sequence, "-A-C-") 
                
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
        
        self.provider.prefix = 'http://www.google.com/NotExisting'
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
    
    class PDBTestCase(test.Case):
        
        def runTest(self):
        
            try:
                StructureParser(self.entry).parse_structure()
            except:
                sys.stdout.write("\n{0}\n".format(self.entry))
                self.reRaise([self.entry])    

    var = 'PDBMASK'
    suite = test.unittest.TestSuite()
    
    if var in os.environ:
        for entry in glob.glob(os.environ[var]):
            case = PDBTestCase()
            case.entry = entry
            suite.addTest(case)
    
    return suite



if __name__ == '__main__':

    test.Console()
