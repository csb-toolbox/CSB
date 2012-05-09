import csb.test as test
import csb.pyutils

from csb.bio.io.wwpdb import FileSystemStructureProvider
from csb.bio.io.fasta import SequenceParser, PDBSequenceParser, SequenceAlignmentReader, StructureAlignmentFactory, OutputBuilder
from csb.bio.sequence import AbstractSequence, Sequence, SequenceTypes, ProteinAlphabet, DuplicateSequenceError, SequenceError
from csb.bio.sequence import AlignmentFormats, SequenceAlignment, A3MAlignment
         
        
@test.unit
class TestSequenceParser(test.Case):
    
    def setUp(self):
        
        super(TestSequenceParser, self).setUp()
        
        self.parser = SequenceParser(product=Sequence)
        self.file = self.config.getTestFile('test.fa')
        self.data = open(self.file).read() 
    
    def testParseFile(self):
        
        s = self.parser.parse_file(self.file)

        self.assertEquals(len(s), 3)

        self.assertTrue(isinstance(s[0], Sequence))
        self.assertEqual(s[0].id, 'gi|148654187')
        self.assertEqual(s[2].sequence, 'RLSRSHEFQRLRREGTRVRSGYLwCVMLQDPSLPGPAVAFAIGRPFGSAVRRNRLRRQLRSILSDRESAMGGGMFLIGVNNPHRDLPMPSFAQLTHDIDEILNK')        
        self.assertEqual(s[1].header, 'gi|142801636|gb|EDA68688.1 hypothetical protein GOS_1956086')
            
    def testParseString(self):

        s = self.parser.parse_string(self.data)
                                           
        self.assertEquals(len(s), 3)

        self.assertTrue(isinstance(s[0], Sequence))
        self.assertEqual(s[0].id, 'gi|148654187')
        self.assertEqual(s[2].sequence, 'RLSRSHEFQRLRREGTRVRSGYLwCVMLQDPSLPGPAVAFAIGRPFGSAVRRNRLRRQLRSILSDRESAMGGGMFLIGVNNPHRDLPMPSFAQLTHDIDEILNK')
        self.assertEqual(s[1].header, 'gi|142801636|gb|EDA68688.1 hypothetical protein GOS_1956086')        

    def testRead(self):
        
        for s in self.parser.read(self.file):
            self.assertTrue(s.id.startswith('gi'))
        

@test.unit
class TestPDBSequenceParser(test.Case):

    def setUp(self):
        
        super(TestPDBSequenceParser, self).setUp()
        
        self.protein = '>3bx6_A mol:protein length:192 ALPHA-1-ACID GLYCOPROTE\nAB-CD'
        self.nucleic = '>3bx6_A mol:na length:192 ALPHA-1-ACID GLYCOPROTE\nAB-CD'
        
        self.parser = PDBSequenceParser()
         
    def testReadJunk(self):
        
        self.assertRaises(ValueError, self.parser.read_sequence, '>junk')
        
    def testReadNucleicSequence(self):
        
        s = self.parser.read_sequence(self.nucleic)
        
        self.assertEqual(s.id, '3bx6A')
        self.assertEqual(s.header, '3bx6_A mol:na length:192 ALPHA-1-ACID GLYCOPROTE')
        self.assertEqual(s.sequence, 'AB-CD')
        self.assertEqual(s.type, SequenceTypes.NucleicAcid)
                
    def testReadProteinSequence(self):
                
        s = self.parser.read_sequence(self.protein)
        
        self.assertEqual(s.id, '3bx6A')
        self.assertEqual(s.header, '3bx6_A mol:protein length:192 ALPHA-1-ACID GLYCOPROTE')
        self.assertEqual(s.sequence, 'AB-CD')
        
    def testParseString(self):
        
        string = '{0}\n{1}'.format(self.nucleic, self.protein)
        seqs = self.parser.parse_string(string)
        
        self.assertEqual(len(seqs), 2)
        self.assertEqual(seqs[0].type, SequenceTypes.NucleicAcid)
        self.assertEqual(seqs[1].type, SequenceTypes.Protein)
        

@test.unit
class TestSequenceAlignmentReader(test.Case):

    def setUp(self):
        
        super(TestSequenceAlignmentReader, self).setUp()
        
        self.reader = SequenceAlignmentReader()
        self.a3m = self.config.getContent('d1nz0a_.a3m')
        self.fasta = self.config.getContent('d1nz0a_.mfasta')
        
    @property
    def size(self):
        return 9
    
    def testReadFASTA(self):
        
        ali = self.reader.read_fasta(self.fasta)
        
        self.assertEqual(ali.size, self.size)
        self.assertEqual(ali.length, 135)
        self.assertEqual(ali.rows[5].id, 'gi|139352214|gb|ECE59672.1|(37-150:150)')
        self.assertEqual(ali.rows[1].sequence, 'ERLRLRRDFLLIFKEG-KSLQNEYF-V---VLFRK--N------GMD---YSRLGIVV-KRK-FGKATRRNKLKRWVR---EIFRRNKGVI---PKGFDIVVIPRK--KLSEEFERVDFWTVREKLLNLLKRIEG')
        
        string = ali.format(AlignmentFormats.FASTA, headers=True).strip()
        self.assertEqual(string, self.fasta.strip())

    def testReadA3M(self):
        
        ali = self.reader.read_a3m(self.a3m)
        
        self.assertEqual(ali.size, self.size)
        self.assertEqual(ali.length, 135)
        self.assertEqual(ali.rows[5].id, 'gi|139352214|gb|ECE59672.1|(37-150:150)')
        self.assertEqual(ali.rows[3].sequence, '-RLTSSKDWKEVRTRG.RCSRSSFA.T...ICVLF..E......GES...E-KFGFAA.AKS.IGSVAKRNRAKRRLR...EAFRQTYKFG...SKPCLVIAIA--..--GPECLTMDFQELKSKL---------')
        
        string = ali.format(AlignmentFormats.A3M, headers=True).strip()
        self.assertEqual(string, self.a3m.strip())

        string = ali.format(AlignmentFormats.FASTA, headers=True).strip()
        self.assertEqual(string, self.fasta.strip())


@test.unit
class TestNonStrictAlignmentReader(TestSequenceAlignmentReader):

    def setUp(self):
        
        super(TestSequenceAlignmentReader, self).setUp()
        
        self.reader = SequenceAlignmentReader(strict=False)
        self.strict = SequenceAlignmentReader(strict=True)
        
        self.a3m = '{0}\n{0}'.format(self.config.getContent('d1nz0a_.a3m'))
        self.fasta = '{0}\n{0}'.format(self.config.getContent('d1nz0a_.mfasta'))

    @property
    def size(self):
        return super(TestNonStrictAlignmentReader, self).size * 2
    
    def testStrict(self):
        
        self.assertRaises(DuplicateSequenceError, self.strict.read_fasta, self.fasta)
        self.assertRaises(DuplicateSequenceError, self.strict.read_a3m, self.a3m) 


@test.unit
class TestStructureAlignmentFactory(test.Case):

    def setUp(self):
        
        super(TestStructureAlignmentFactory, self).setUp()
        
        self.provider = FileSystemStructureProvider(self.config.data)
        self.provider.add_template('{id}.regular.pdb')
        self.factory = StructureAlignmentFactory(self.provider)

        self.fasta = self.config.getContent('struct.ali.mfasta')
     
    def testMakeAlignment(self):
        
        ali = self.factory.make_alignment(self.fasta)
        
        self.assertEqual(ali.format().strip(), self.fasta.strip())
        self.assertEqual(ali.rows[1].columns[1].residue.type, ProteinAlphabet.GLU)
        self.assertEqual(ali.rows[2].columns[12].residue.atoms['CA'].vector[0], 58.911000000000001)
        
    def testMakeEntry(self):
        
        row = Sequence('1d3zA', '', 'DTI--ENVK', SequenceTypes.Protein)
        chain = self.provider.get('1d3z').chains['A']
        
        entry = self.factory.make_entry(row, chain)
        self.assertEqual(row.sequence, entry.sequence)
        self.assertEqual(row.residues[4].type, ProteinAlphabet.GAP)
        self.assertEqual(row.residues[5].type, ProteinAlphabet.GAP)
        
        row = Sequence('1d3zA', '', 'DBUG-ENVK', SequenceTypes.Protein)
        self.assertRaises(SequenceError, self.factory.make_entry, row, chain)
        

@test.unit
class TestOutputBuilder(test.Case):
    
    def testCreate(self):
        
        sequence = Sequence('s1', 's1', 'ABC', SequenceTypes.Protein) 
        
        with self.config.getTempStream() as tmp:
            for format in csb.pyutils.Enum.members(AlignmentFormats):
                
                builder = OutputBuilder.create(format, tmp, True)
                builder.add_sequence(sequence)
            
        
            self.assertRaises(ValueError, OutputBuilder.create, 'bug', tmp, True)
            
@test.unit
class TestFASTAOutputBuilder(test.Case):
    
    def setUp(self):
        super(TestFASTAOutputBuilder, self).setUp()
        
        self.s1 = Sequence('s1A', 's1A desc', 'A.-D', SequenceTypes.Protein)
        self.s2 = Sequence('s2A', 's2A desc', 'ABCD', SequenceTypes.Protein)
        self.ali = SequenceAlignment([self.s1, self.s2])
    
    @property    
    def format(self):
        return AlignmentFormats.FASTA
    
    @property
    def sequence_string(self):
        return '>s1A desc\nA--D'

    @property
    def alignment_string(self):
        return '>s1A desc\nA--D\n>s2A desc\nABCD'
            
    def testAddSequence(self):
        
        with self.config.getTempStream() as tmp:
            
            builder = OutputBuilder.create(self.format, tmp, True)
            builder.add_sequence(self.s1)
            
            tmp.flush()
            self.assertEqual(open(tmp.name).read().strip(), self.sequence_string)
                        
    def testAddAlignment(self):
        
        with self.config.getTempStream() as tmp:
            
            builder = OutputBuilder.create(self.format, tmp, True)
            builder.add_alignment(self.ali)
            
            tmp.flush()
            self.assertEqual(open(tmp.name).read().strip(), self.alignment_string)
            
@test.unit
class TestA3MOutputBuilder(TestFASTAOutputBuilder):
    
    def setUp(self):
        
        super(TestA3MOutputBuilder, self).setUp()
        self.ali = A3MAlignment([self.s1, self.s2])
        
    @property    
    def format(self):
        return AlignmentFormats.A3M
    
    @property
    def sequence_string(self):
        return '>s1A desc\nA.-D'

    @property
    def alignment_string(self):
        return '>s1A desc\nA-D\n>s2A desc\nAbCD'

@test.unit
class TestPIROutputBuilder(TestFASTAOutputBuilder):
        
    @property    
    def format(self):
        return AlignmentFormats.PIR
    
    @property
    def sequence_string(self):
        return '>P1;s1\nstructure:s1:.:A:.:A::::\nA--D*'

    @property
    def alignment_string(self):
        return r'''
>P1;s1
sequence:s1:.:A:.:A::::
A--D*
>P1;s2
structure:s2:.:A:.:A::::
ABCD*
        '''.strip()

@test.unit
class TestCrossAlignmentBuilding(test.Case):
    
    def setUp(self):
        super(TestCrossAlignmentBuilding, self).setUp()
        
        self.s1s = Sequence('s1A', 's1A desc', 'A-CD', SequenceTypes.Protein)
        self.s1a = Sequence('s1A', 's1A desc', 'A.CD', SequenceTypes.Protein)        
        self.s2  = Sequence('s2A', 's2A desc', 'AB-D', SequenceTypes.Protein)
        
        self.seq = SequenceAlignment([self.s1s, self.s2])
        self.a3m = A3MAlignment([self.s1a, self.s2])
        
    def testFastaToA3m(self):
        
        with self.config.getTempStream() as tmp:
            
            builder = OutputBuilder.create(AlignmentFormats.A3M, tmp, headers=False)
            builder.add_alignment(self.seq)
            
            tmp.flush()
            self.assertEqual(open(tmp.name).read().strip(), 'ACD\nAb-D')
            
    def testA3mToFasta(self):
        
        with self.config.getTempStream() as tmp:
            
            builder = OutputBuilder.create(AlignmentFormats.FASTA, tmp, headers=False)
            builder.add_alignment(self.a3m)
            
            tmp.flush()
            self.assertEqual(open(tmp.name).read().strip(), 'A-CD\nAB-D')
                        
            
if __name__ == '__main__':
    
    test.Console()
