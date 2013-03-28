
import csb.bio.sequence as sequence
import csb.bio.io.fasta

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
        self.assertEqual(self.sequence.header, 'id desc')
        
    def testSequence(self):
        self.assertEqual(self.sequence.sequence, 'AB-CD')
        
    def testType(self):                
        
        def test(v):
            self.sequence.type = v
            
        self.assertRaises(TypeError, test, sequence.ProteinAlphabet.ALA)
        self.assertEqual(self.sequence.type, sequence.SequenceTypes.Protein)
        
    def testAlphabet(self):        
        self.assertEqual(self.sequence.alphabet, sequence.ProteinAlphabet)
                
    def testStrip(self):
        self.assertEqual(self.sequence.strip().sequence, 'ABCD')

    def testSubregion(self):
        
        s = self.sequence
        
        self.assertEqual(s.subregion(2, 3).sequence,  'B-')        
        self.assertRaises(sequence.SequencePositionError, s.subregion, -1, 2)
        self.assertRaises(sequence.SequencePositionError, s.subregion, 1, 6)
        
    def testExtract(self):
        
        s = self.sequence
        self.assertEqual(s.extract((2, 3, 5)).sequence, 'B-D')
        
        self.assertRaises(sequence.SequencePositionError, s.extract, [-1])
        self.assertRaises(sequence.SequencePositionError, s.extract, [6])
        
    def testIndexeres(self):
        
        s = self.sequence
        
        self.assertEqual(s.residues[2].type, sequence.ProteinAlphabet.ASX)
        self.assertEqual(s.residues[1].type, s[0].type)
        self.assertEqual(s[:].sequence, s.sequence)        
        self.assertEqual(s[2:].sequence, '-CD')
        self.assertEqual(s[2:3].sequence, '-')
        
        for rank in [-1, 0, 6]:
            self.assertRaises(sequence.SequencePositionError, lambda i: s.residues[i], rank)

        for index in [-1, 5]:
            self.assertRaises(IndexError, lambda i: s[i], index)
        
        self.assertRaises(IndexError, lambda: s[-1:])
            
    def testIterator(self):
        
        chars = [ str(r.type) for r in self.sequence ]
        seq = ''.join(chars)
        
        self.assertEqual(self.sequence.sequence.upper(), seq)
                    
    def testToString(self):
        
        self.assertEqual(str(self.sequence), self.data)
        self.assertEqual(str(self.sequence), self.data)
        
    def testDump(self):
        
        with self.config.getTempStream() as tmp:
            self.sequence.dump(tmp.name)
            tmp.flush()
            
            self.assertEqual(open(tmp.name).read().strip(), self.data)
            
        with self.config.getTempStream() as tmp:
            self.sequence.dump(tmp)
            tmp.flush()
            
            self.assertEqual(open(tmp.name).read().strip(), self.data)           

@test.unit
class TestRichSequence(TestSequence):
    
    def setUp(self):
        
        super(TestRichSequence, self).setUp()
        self.sequence = sequence.RichSequence.create(self.sequence)
                                              
@test.unit
class TestChainSequence(TestSequence):

    def setUp(self):
        
        super(TestChainSequence, self).setUp()
        
        from csb.bio.structure import ProteinResidue, Chain
        chain = Chain('A', name='desc', accession='accn')
        
        for rank, char in enumerate('AB-CD', start=1):
            chain.residues.append(ProteinResidue(rank, char))
         
        self.sequence = sequence.ChainSequence.create(chain)
        self.sequence.header = '>id desc'
        self.sequence.id = 'id'


@test.unit
class TestSequenceCollection(test.Case):

    def setUp(self):
        
        super(TestSequenceCollection, self).setUp()
        
        s1 = sequence.Sequence('id1', '>id1 desc', 'AB-CD', sequence.SequenceTypes.Protein)
        s2 = sequence.Sequence('id2', '>id2 desc', 'ABCDE', sequence.SequenceTypes.Protein)
        
        self.collection = sequence.SequenceCollection([s1, s2])
        self.data = '>id1 desc\nAB-CD\n>id2 desc\nABCDE' 
         
    def testToFASTA(self):
        
        with self.config.getTempStream() as tmp:
            self.collection.to_fasta(tmp.name)
            tmp.flush()
            
            self.assertEqual(open(tmp.name).read().strip(), self.data)
            
        with self.config.getTempStream() as tmp:
            self.collection.to_fasta(tmp)
            tmp.flush()
            
            self.assertEqual(open(tmp.name).read().strip(), self.data)        

    
@test.unit
class TestSequenceAlignment(test.Case):
    
    
    def _factory(self, sequences, strict=True):
        return sequence.SequenceAlignment(sequences, strict=strict)
    
    def setUp(self):
        
        super(TestSequenceAlignment, self).setUp()
        
        seq1 = sequence.Sequence('s1', 's1 desc1',          'AB-CD', sequence.SequenceTypes.Protein)
        seq2 = sequence.RichSequence('s2', 's2 desc2', list('ABX-D'), sequence.SequenceTypes.Protein)
            
        self.ali = self._factory([seq1, seq2])
        self.ali2 = self._factory([seq1, seq2, seq2, seq2], strict=False)

    def testAdd(self):
        
        # strict
        self.assertRaises(sequence.SequenceError, self.ali.add, sequence.Sequence('sn', 'sn','TOO-LONG-SEQ'))
        self.assertRaises(sequence.DuplicateSequenceError, self.ali.add, sequence.Sequence('s1', 's1','AB-CD'))
        
    def testLength(self):
        self.assertEqual(self.ali.length, 5)

    def testSize(self):
        self.assertEqual(self.ali.size, 2)
        self.assertEqual(self.ali2.size, 4)

    def testSubregion(self):
        
        sub = self.ali.subregion(2, 4)
        
        self.assertEqual(sub.length, 3)
        self.assertEqual(sub.size, 2)
        self.assertEqual(sub.rows[1].sequence, 'B-C')
        self.assertEqual(sub.rows[2].sequence, 'BX-')
        
        self.assertRaises(sequence.ColumnPositionError, self.ali.subregion, -1, 2)
        self.assertRaises(sequence.ColumnPositionError, self.ali.subregion, 1, 6)
        
        # should not raise DuplicateSequenceError
        self.ali2.subregion(1, 2)
        
    def testFormat(self):
        
        self.assertEqual(self.ali.format(headers=False).strip(), 'AB-CD\nABX-D')
        self.assertEqual(self.ali.format(headers=True).strip(), '>s1 desc1\nAB-CD\n>s2 desc2\nABX-D')
        
    def testRows(self):
        
        a = self.ali

        self.assertEqual(a.rows[2].id, 's2')
        self.assertRaises(sequence.SequenceNotFoundError, lambda i: a.rows[i], -1)
        self.assertRaises(sequence.SequenceNotFoundError, lambda i: a.rows[i], 3)
        
        # with duplicates
        self.assertEqual(self.ali2.rows['s2'].id, 's2')
        self.assertEqual(self.ali2.rows['s2:A1'].id, 's2')
        self.assertEqual(self.ali2.rows['s2:A2'].id, 's2')        

    def testColumns(self):
        
        a = self.ali
                
        self.assertEqual(a.columns[4][0].id, 's1')
        self.assertEqual(a.columns[4][0].column, 4)
        self.assertEqual(a.columns[4][0].residue.type, sequence.ProteinAlphabet.CYS)
        self.assertEqual(len(a.columns[4]), a.size)
        self.assertEqual(len(a.columns[4]), len(a.columns[3]))

    def testRowColumns(self):
        
        a = self.ali
        
        self.assertEqual(a.rows[1].columns[4].id, a.columns[4][0].id)
        self.assertEqual(a.rows[1].columns[4].column, a.columns[4][0].column)
        self.assertEqual(a.rows[1].columns[4].residue.type, a.columns[4][0].residue.type)       

    def testRowResidues(self):
        
        a = self.ali
        
        self.assertEqual(a.rows[1].residues[4].type, sequence.ProteinAlphabet.ASP)
        self.assertEqual(a.rows[1].residues[3].type, sequence.ProteinAlphabet.CYS)
        
        self.assertEqual(a.rows[2].residues[4].type, sequence.ProteinAlphabet.ASP)
        self.assertEqual(a.rows[2].residues[3].type, sequence.ProteinAlphabet.UNK)

    def testRowMap(self):
        
        a = self.ali
        
        self.assertEqual(a.rows[1].map_column(4), 3)
        self.assertEqual(a.rows[1].map_residue(3), 4)

    def testIndexer(self):
        
        a = self.ali

        self.assertEqual(a[0, 0].rows['s1'].id, 's1')
        self.assertEqual(a[1, 4].rows['s2'].columns[1].residue.type, sequence.ProteinAlphabet.ASP)
                
        self.assertEqual(a[0, 0].size, 1)
        self.assertEqual(a[0, 0].length, 1)

        self.assertEqual(a[:, 0].size, 2)
        self.assertEqual(a[:, 0].length, 1)

        self.assertEqual(a[0, :].size, 1)
        self.assertEqual(a[0, :].length, 5)

        self.assertEqual(a[0:2, 0:2].size, 2)
        self.assertEqual(a[0:2, 0:2].length, 2)

        self.assertEqual(a[(0, 1), (0, 1, 3)].size, 2)
        self.assertEqual(a[(0, 1), (0, 1, 3)].length, 3)

        self.assertRaises(IndexError, lambda: a[-1, :])
        self.assertRaises(IndexError, lambda: a[:, -1])
        
        self.assertRaises(TypeError, lambda: a['', :])
        self.assertRaises(TypeError, lambda: a[:, ''])

        self.assertRaises(ValueError, lambda: a[[], :])
        self.assertRaises(ValueError, lambda: a[:, []])

        self.assertRaises(IndexError, lambda: a[-1:, :])
        self.assertRaises(IndexError, lambda: a[:, -1:])        

    def testGapAt(self):
        
        self.assertFalse(self.ali.gap_at(1))
        self.assertTrue(self.ali.gap_at(3))
        

@test.unit
class TestA3MAlignmentSimple(TestSequenceAlignment):
    
    def _factory(self, sequences, strict=True):
        return sequence.A3MAlignment(sequences, strict=strict)
        
@test.unit
class TestA3MAlignment(test.Case):

    def setUp(self):
        
        super(TestA3MAlignment, self).setUp()
        
        self.file = self.config.getTestFile('d1nz0a_.a3m')
        self.a3m = open(self.file).read()
        self.ali = sequence.A3MAlignment.parse(self.a3m)

    def testMatches(self):
        self.assertEqual(self.ali.matches, 109)

    def testLength(self):
        self.assertEqual(self.ali.length, 135)
        
    def testSize(self):
        self.assertEqual(self.ali.size, 9)            

    def testRows(self):
        
        # row 1 (master)
        row = self.ali.rows[1]
        self.assertEqual(row.id, 'd1nz0a_')

        self.assertEqual(row.length, 135)
        self.assertEqual(row.strip().length, 109)        
        self.assertEqual(row.strip().sequence, 'ERLRLRRDFLLIFKEGKSLQNEYFVVLFRKNGMDYSRLGIVVKRKFGKATRRNKLKRWVREIFRRNKGVIPKGFDIVVIPRKKLSEEFERVDFWTVREKLLNLLKRIEG')
        
    def testToString(self):
         
        self.assertEqual(self.a3m.strip(), self.ali.format().strip())
        
    def testFormat(self):
        
        fasta = self.ali.format(sequence.AlignmentFormats.FASTA, headers=False).splitlines()
        self.assertEqual(len(fasta), 9)
        
        for line in fasta:
            self.assertEqual(len(line), 135)
            
        ref = open(self.config.getTestFile('d1nz0a_.mfasta')).read()
        self.assertEqual(ref.strip(), self.ali.format(sequence.AlignmentFormats.FASTA, headers=True).strip())
    
    def testHMMSubregion(self):
        
        sub = self.ali.hmm_subregion(2, 30)
        self.assertEqual(sub.rows['d1nz0a_'].strip().sequence, 'RLRLRRDFLLIFKEGKSLQNEYFVVLFRK')
        self.assertEqual(sub.size, self.ali.size)
        self.assertEqual(sub.matches, 30 - 2 + 1)                

        
        fasta = self.ali.hmm_subregion(1, 109).format(sequence.AlignmentFormats.FASTA, headers=True)       
        ref = open(self.config.getTestFile('d1nz0a_.mfasta')).read()
        self.assertEqual(ref, fasta.strip())   
        
        self.assertRaises(sequence.ColumnPositionError, self.ali.subregion, -1, 2)
        self.assertRaises(sequence.ColumnPositionError, self.ali.hmm_subregion, -1, 2)
        self.assertRaises(sequence.ColumnPositionError, self.ali.subregion, 1, 111110)
        self.assertRaises(sequence.ColumnPositionError, self.ali.hmm_subregion, 1, 110)
         
    def testInsertionAt(self):
        
        self.assertFalse(self.ali.insertion_at(1))
        self.assertTrue(self.ali.insertion_at(17))              
        
        
if __name__ == '__main__':

    test.Console()
