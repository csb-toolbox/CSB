import csb.test as test

from csb.bio.sequence import RichSequence, SequenceTypes
from csb.bio.sequence.alignment import IdentityMatrix, SimilarityMatrix
from csb.bio.sequence.alignment import GlobalAlignmentAlgorithm, LocalAlignmentAlgorithm, AlignmentResult


@test.unit
class TestIdentityMatrix(test.Case):
    
    def setUp(self):
        
        super(TestIdentityMatrix, self).setUp()
        self.matrix = IdentityMatrix(2, -3)
        
    def testScore(self):
        self.assertEqual(self.matrix.score("a", "a"), 2)
        self.assertEqual(self.matrix.score("a", "b"), -3)
        
@test.unit
class TestSimilarityMatrix(test.Case):
    
    def setUp(self):
        
        super(TestSimilarityMatrix, self).setUp()
        self.matrix = SimilarityMatrix(SimilarityMatrix.BLOSUM62)
        
    def testScore(self):
        self.assertEqual(self.matrix.score("A", "A"), 4)
        self.assertEqual(self.matrix.score("A", "R"), -1)  
        self.assertEqual(self.matrix.score("R", "A"), -1)
        
        
@test.unit
class TestGlobalAlignmentAlgorithm(test.Case):
    
    def setUp(self):
        
        super(TestGlobalAlignmentAlgorithm, self).setUp()
        
        self.seq1 = RichSequence('s1', '', 'CCABBBCBBCABAABCCEAAAAAAAAAAAAFAA', SequenceTypes.Protein)
        self.seq2 = RichSequence('s1', '', 'AZCBBABAABCCEF', SequenceTypes.Protein)
        self.algorithm = GlobalAlignmentAlgorithm(scoring=IdentityMatrix(1, -1), gap=0)  
        
    def testAlign(self):

        ali = self.algorithm.align(self.seq1, self.seq2)
        
        self.assertEqual(ali.query.sequence,   "CCA-BBBCBBCABAABCCEAAAAAAAAAAAAFAA")
        self.assertEqual(ali.subject.sequence, "--AZ---CBB-ABAABCCE------------F--")

        self.assertEqual(ali.query.residues[3], self.seq1.residues[3])
        self.assertTrue(ali.query.residues[3] is self.seq1.residues[3])
                
        self.assertEqual(ali.qstart, 1)
        self.assertEqual(ali.qend, 33)
        self.assertEqual(ali.start, 1)
        self.assertEqual(ali.end, 14)

        self.assertEqual(ali.length, 34)        
        self.assertEqual(ali.gaps, 21)
        self.assertEqual(ali.identicals, 13)
        self.assertEqual(ali.identity, 13 / 34.0 )
        self.assertEqual(ali.score, 13)
        
    def testEmptyAlignment(self):
        
        seq1 = RichSequence('s1', '', 'AAAA', SequenceTypes.Protein)
        seq2 = RichSequence('s2', '', 'BBBB', SequenceTypes.Protein)
        
        ali = self.algorithm.align(seq1, seq2)
        self.assertTrue(ali.is_empty)        

@test.unit
class TestLocalAlignmentAlgorithm(test.Case):
    
    def setUp(self):
        
        super(TestLocalAlignmentAlgorithm, self).setUp()
        
        self.seq1 = RichSequence('s1', '', 'CCABBBCBBCABAABCCEAAAAAAAAAAAAFAA', SequenceTypes.Protein)
        self.seq2 = RichSequence('s1', '', 'AZCBBABAACBCCEF', SequenceTypes.Protein)
        self.algorithm = LocalAlignmentAlgorithm(scoring=IdentityMatrix(1, -1), gap=-1)  
        
    def testAlign(self):

        ali = self.algorithm.align(self.seq1, self.seq2)
                
        self.assertEqual(ali.query.sequence,   "CBBCABAA-BCCE")
        self.assertEqual(ali.subject.sequence, "CBB-ABAACBCCE")         
        
        self.assertEqual(ali.qstart, 7)
        self.assertEqual(ali.qend, 18)
        self.assertEqual(ali.start, 3)
        self.assertEqual(ali.end, 14)

        self.assertEqual(ali.length, 13)        
        self.assertEqual(ali.gaps, 2)
        self.assertEqual(ali.identicals, 11)
        self.assertEqual(ali.identity, 11 / 13.0 )
        self.assertEqual(ali.score, 9)
        
    def testEmptyAlignment(self):
        
        seq1 = RichSequence('s1', '', 'AAAA', SequenceTypes.Protein)
        seq2 = RichSequence('s2', '', 'BBBB', SequenceTypes.Protein)
        
        ali = self.algorithm.align(seq1, seq2)
        self.assertTrue(ali.is_empty)
        

@test.unit
class TestAlignmentResult(test.Case):
    
    def setUp(self):
        
        super(TestAlignmentResult, self).setUp()        
        
        self.seq1 = RichSequence('s1', '', 'AB-D', SequenceTypes.Protein)
        self.seq2 = RichSequence('s2', '', 'A-CD', SequenceTypes.Protein)        
        self.ali = AlignmentResult(5.5, self.seq1, self.seq2, 10, 12, 20, 22)
        
        self.es = RichSequence('s1', '', '')
        self.empty = AlignmentResult(0, self.es, self.es, 0, 0, 0, 0)
        
    def testConstructor(self):
        
        self.assertRaises(ValueError, AlignmentResult, 1, self.es, self.es, 0, 0, 0, 0)
        self.assertRaises(ValueError, AlignmentResult, 0, self.es, self.es, 1, 0, 0, 0)
        self.assertRaises(ValueError, AlignmentResult, 0, self.es, self.es, 0, 1, 0, 0)
        self.assertRaises(ValueError, AlignmentResult, 0, self.es, self.es, 0, 0, 1, 0)
        self.assertRaises(ValueError, AlignmentResult, 0, self.es, self.es, 0, 0, 0, 1)
        
        self.assertRaises(ValueError, AlignmentResult, 1, self.seq1, self.seq2, 0, 0, 0, 0)
        
    def testStr(self):
        
        string = r"""
   10 AB-D 12   
   20 A-CD 22   """.strip("\r\n")
        self.assertEqual(string, str(self.ali))
        
    def testAlignment(self):
        
        ali = self.ali.alignment()
        self.assertEqual(ali.rows[1].sequence, self.seq1.sequence)
        self.assertEqual(ali.rows[2].sequence, self.seq2.sequence)
        
    def testQuery(self):
        self.assertEqual(self.ali.query.sequence, self.seq1.sequence)
        self.assertEqual(self.ali.query.residues[2], self.seq1.residues[2])
        self.assertTrue(self.ali.query.residues[2] is self.seq1.residues[2]) 
        
    def testSubject(self):
        self.assertEqual(self.ali.subject.sequence, self.seq2.sequence)
        self.assertEqual(self.ali.subject.residues[3], self.seq2.residues[3])
        self.assertTrue(self.ali.subject.residues[3] is self.seq2.residues[3])         
        
    def testQstart(self):
        self.assertEqual(self.ali.qstart, 10)

    def testQend(self):
        self.assertEqual(self.ali.qend, 12)
        
    def testStart(self):       
        self.assertEqual(self.ali.start, 20)
                
    def testEnd(self):
        self.assertEqual(self.ali.end, 22)
        
    def testLength(self):
        self.assertEqual(self.ali.length, 4)

    def testScore(self):
        self.assertEqual(self.ali.score, 5.5)
        
    def testGaps(self):
        self.assertEqual(self.ali.gaps, 2)
                                
    def testIdenticals(self):
        self.assertEqual(self.ali.identicals, 2)
        
    def testIdentity(self):
        self.assertEqual(self.ali.identity, 0.5) 
        
    def testIsEmpty(self):
        self.assertFalse(self.ali.is_empty)
        
        es = RichSequence('s1', '', '')
        empty = AlignmentResult(0, es, es, 0, 0, 0, 0)
        self.assertTrue(empty.is_empty)
                

if __name__ == '__main__':

    test.Console()
