import csb.bio.sequence as sequence
import csb.test as test


@test.unit
class TestSequence(test.Case):
    
    def setUp(self):
        
        super(TestSequence, self).setUp()
        
        self.data = '>id desc\nAB-CD'
        self.sequence = sequence.Sequence('id', '>id desc', 'AB-CD', sequence.SequenceTypes.Protein)     
        assert str(self.sequence) == self.data
        
    def testLength(self):
        self.assertEqual(self.sequence.length, 5)
        self.assertEqual(self.sequence.length, len(self.sequence))
        
    def testId(self):
        self.assertEqual(self.sequence.id, 'id')        
    
    def testHeader(self):
        self.assertEqual(self.sequence.header, '>id desc')
        
    def testSequence(self):
        self.assertEqual(self.sequence.sequence, 'AB-CD')
        
    def testType(self):        
        self.assertEqual(self.sequence.type, sequence.SequenceTypes.Protein)                            
    
    def testToFASTA(self):
        
        # in memory
        self.assertEqual(self.sequence.to_fasta(), self.data)
        self.assertEqual(str(self.sequence), self.data)
        
        # to stream
        temp = self.config.getTempStream()
        self.sequence.to_fasta(temp)
        temp.flush()
        self.assertEqual(open(temp.name).read(), self.data + '\n')
        
        # to file
        temp = self.config.getTempStream().name
        self.sequence.to_fasta(temp)
        self.assertEqual(open(temp).read(), self.data + '\n')
        
    def testExpand(self):
        self.assertEqual(self.sequence.expand(), [0, 1, 3, 4])   

    def testIsGap(self):
        self.assertEqual(self.sequence.is_gap(2), True)
        for i in [0, 1, 3, 4]:
            self.assertEqual(self.sequence.is_gap(i), False)

    def testSequenceIndex(self):
        for c, s in [(0, 0), (1, 1), (2, 2), (3, 2), (4, 3)]:
            self.assertEqual(self.sequence.sequence_index(c), s)        
    
    def testFromString(self):
        
        # normal FASTA
        s = sequence.Sequence.from_string('\n ' + self.data.replace('-', '\n\n -  \n') + ' \n\n')
        
        self.assertEqual(s.id, 'id')
        self.assertEqual(s.header, '>id desc')
        self.assertEqual(s.sequence, 'AB-CD')
        
        # plain text sequence
        s = sequence.Sequence.from_string('\n AB \n CD  \n')
        
        self.assertEqual(s.id, 'sequence')
        self.assertEqual(s.header, '>sequence')
        self.assertEqual(s.sequence, 'ABCD')
        

@test.unit
class TestPDBSequence(test.Case):
    
    def setUp(self):
        
        super(TestPDBSequence, self).setUp()
        
        self.protein = '>3bx6_A mol:protein length:192 ALPHA-1-ACID GLYCOPROTE\nAB-CD'
        self.nucleic = '>3bx6_A mol:na length:192 ALPHA-1-ACID GLYCOPROTE\nAB-CD'
        
        self.sequence = sequence.PDBSequence(
                    '3bx6A', 
                    '>3bx6_A mol:protein length:192 ALPHA-1-ACID GLYCOPROTE', 
                    'AB-CD', 
                    sequence.SequenceTypes.Protein)     
           
    def testAccession(self):
        self.assertEqual(self.sequence.accession, '3bx6')

    def testChain(self):
        self.assertEqual(self.sequence.chain, 'A')
        
    def testFromString(self):

        # protein
        s = sequence.PDBSequence.from_string('\n ' + self.protein + ' \n\n')
        
        self.assertEqual(s.id, '3bx6A')
        self.assertEqual(s.header, '>3bx6_A mol:protein length:192 ALPHA-1-ACID GLYCOPROTE')
        self.assertEqual(s.sequence, 'AB-CD')                       

        # nucleic
        s = sequence.PDBSequence.from_string('\n ' + self.nucleic + ' \n\n')
        
        self.assertEqual(s.id, '3bx6A')
        self.assertEqual(s.header, '>3bx6_A mol:na length:192 ALPHA-1-ACID GLYCOPROTE')
        self.assertEqual(s.sequence, 'AB-CD')
        self.assertEqual(s.type, sequence.SequenceTypes.NucleicAcid)        
        
        # junk
        self.assertRaises(ValueError, sequence.PDBSequence.from_string, '>junk')   
        
        
@test.unit
class TestA3MAlignment(test.Case):

    def setUp(self):
        
        super(TestA3MAlignment, self).setUp()
        
        self.file = self.config.getTestFile('d1nz0a_.a3m')
        self.a3m = open(self.file).read()
        self.ali = sequence.A3MAlignment.parse(self.a3m)
        
    def testConsensus(self):        
        self.assertEqual(self.ali.consensus.sequence, 'eRLkxxxdFxxvxxxgxxxxxxxxxlxxxxxxxxxxRxGxxvsKKvgxAVxRNriKRxlRexxrxxxxxlxxxxdivvixrxxxxxxxxxxxxxxlxxxlxxlxkkixg')

    def testMatchesCount(self):
        self.assertEqual(self.ali.matches_count, 109)

    def testColumnsCount(self):
        self.assertEqual(self.ali.cols_count, 135)
        
    def testRowsCount(self):
        self.assertEqual(self.ali.rows_count, 9)            

    def testRows(self):
        
        # row 1 (master)
        seqInfo, seq = self.ali.rows[1][0], self.ali.rows[1][1:]
        #  - sequence info: row[1][0]
        self.assertTrue(isinstance(seqInfo, sequence.Sequence))
        self.assertEqual(seqInfo.id, 'd1nz0a_')

        self.assertEqual(len(seq), 135)
        self.assertEqual(len(''.join(seq)), 109)        
        self.assertEqual(''.join(seq), 'ERLRLRRDFLLIFKEGKSLQNEYFVVLFRKNGMDYSRLGIVVKRKFGKATRRNKLKRWVREIFRRNKGVIPKGFDIVVIPRKKLSEEFERVDFWTVREKLLNLLKRIEG')
        self.assertEqual(seqInfo.sequence, ''.join(seq))
        
    def testToString(self):
 
        a3m = '{0}\n{1}'.format(self.ali.consensus.to_fasta(), self.ali.to_string())          
        self.assertEqual(self.a3m, a3m)
        
    def testToFASTA(self):
        
        fasta = self.ali.to_fasta(headers=False).splitlines()
        self.assertEqual(len(fasta), 9)
        
        for line in fasta:
            self.assertEqual(len(line), 135)
            
        ref = open(self.config.getTestFile('d1nz0a_.mfasta')).read()
        self.assertEqual(ref, self.ali.to_fasta())
    
    def testSubregion(self):
        
        sub = self.ali.subregion(2, 8)
        self.assertEqual(sub.rows[1][0].sequence, 'RLRLRRD')
        self.assertEqual(sub.rows_count, self.ali.rows_count)
        self.assertEqual(sub.matches_count, 8 - 2 + 1)                
        
        ref = open(self.config.getTestFile('d1nz0a_.mfasta')).read()
        self.assertEqual(ref, self.ali.subregion(1, 109).to_fasta())   
        
        self.assertRaises(IndexError, self.ali.subregion, -1, 2)
        self.assertRaises(IndexError, self.ali.subregion, 1, 110)                  
        
        
if __name__ == '__main__':

    test.Console()
