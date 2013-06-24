"""
Collection of sequence alignment algorithms.

@note: The classes in this module have been optimized for performance.
       Think twice before switching a field to a generally nicer property
       access, because it turns out that these things often add significant
       constants to the running time of a dynamic programming algorithm.
"""

from csb.bio.sequence import AbstractSequence, SequenceAlignment, RichSequence, ResidueInfo
from abc import ABCMeta, abstractmethod


class ResidueNotFoundError(KeyError):
    pass


class AbstractScoringMatrix(object):
    """
    Defines a pairwise sequence alignment scoring function.
    """
    
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def score(self, x, y):
        """
        Return the pairwise score of residues C{x} and C{y}.
        C{x} and C{y} must be comparable (e.g. implement __eq__ and __hash__).
        
        @param x: first residue
        @type x: object
        @param y: second residue
        @type y: object
        
        @rtype: float
        
        @raise ResidueNotFoundError: if C{x} or C{y} cannot be handled by this
        scoring matrix          
        """
        pass
    
class IdentityMatrix(AbstractScoringMatrix):
    """
    Simple identity-based scoring matrix.
    
    @param match: score for a match
    @type match: float
    @param mismatch: score for a mismatch
    @type mismatch: float
    """

    def __init__(self, match=1, mismatch=-1):
        
        self._match = float(match)
        self._mismatch = float(mismatch)
        
    @property
    def match(self):
        """
        Score for a match
        @rtype: float
        """
        return self._match
    
    @property
    def mismatch(self):
        """
        Score for a mismatch
        @rtype: float
        """        
        return self._mismatch
    
    def score(self, x, y):
        
        if x == y:
            return self._match
        else:
            return self._mismatch      
        
class SimilarityMatrix(AbstractScoringMatrix):
    """
    Similarity scoring matrix.
    
    @param matrix: 
    """
    
    BLOSUM62 = { 
        'A': { 'A': 4.0, 'R':-1.0, 'N':-2.0, 'D':-2.0, 'C': 0.0, 'Q':-1.0, 'E':-1.0, 'G': 0.0, 'H':-2.0, 'I':-1.0, 'L':-1.0, 'K':-1.0, 'M':-1.0, 'F':-2.0, 'P':-1.0, 'S': 1.0, 'T': 0.0, 'W':-3.0, 'Y':-2.0, 'V': 0.0, 'B':-2.0, 'Z':-1.0, 'X': 0.0, '*':-4.0 },
        'R': { 'A':-1.0, 'R': 5.0, 'N': 0.0, 'D':-2.0, 'C':-3.0, 'Q': 1.0, 'E': 0.0, 'G':-2.0, 'H': 0.0, 'I':-3.0, 'L':-2.0, 'K': 2.0, 'M':-1.0, 'F':-3.0, 'P':-2.0, 'S':-1.0, 'T':-1.0, 'W':-3.0, 'Y':-2.0, 'V':-3.0, 'B':-1.0, 'Z': 0.0, 'X':-1.0, '*':-4.0 },
        'N': { 'A':-2.0, 'R': 0.0, 'N': 6.0, 'D': 1.0, 'C':-3.0, 'Q': 0.0, 'E': 0.0, 'G': 0.0, 'H': 1.0, 'I':-3.0, 'L':-3.0, 'K': 0.0, 'M':-2.0, 'F':-3.0, 'P':-2.0, 'S': 1.0, 'T': 0.0, 'W':-4.0, 'Y':-2.0, 'V':-3.0, 'B': 3.0, 'Z': 0.0, 'X':-1.0, '*':-4.0 },
        'D': { 'A':-2.0, 'R':-2.0, 'N': 1.0, 'D': 6.0, 'C':-3.0, 'Q': 0.0, 'E': 2.0, 'G':-1.0, 'H':-1.0, 'I':-3.0, 'L':-4.0, 'K':-1.0, 'M':-3.0, 'F':-3.0, 'P':-1.0, 'S': 0.0, 'T':-1.0, 'W':-4.0, 'Y':-3.0, 'V':-3.0, 'B': 4.0, 'Z': 1.0, 'X':-1.0, '*':-4.0 },
        'C': { 'A': 0.0, 'R':-3.0, 'N':-3.0, 'D':-3.0, 'C': 9.0, 'Q':-3.0, 'E':-4.0, 'G':-3.0, 'H':-3.0, 'I':-1.0, 'L':-1.0, 'K':-3.0, 'M':-1.0, 'F':-2.0, 'P':-3.0, 'S':-1.0, 'T':-1.0, 'W':-2.0, 'Y':-2.0, 'V':-1.0, 'B':-3.0, 'Z':-3.0, 'X':-2.0, '*':-4.0 },
        'Q': { 'A':-1.0, 'R': 1.0, 'N': 0.0, 'D': 0.0, 'C':-3.0, 'Q': 5.0, 'E': 2.0, 'G':-2.0, 'H': 0.0, 'I':-3.0, 'L':-2.0, 'K': 1.0, 'M': 0.0, 'F':-3.0, 'P':-1.0, 'S': 0.0, 'T':-1.0, 'W':-2.0, 'Y':-1.0, 'V':-2.0, 'B': 0.0, 'Z': 3.0, 'X':-1.0, '*':-4.0 },
        'E': { 'A':-1.0, 'R': 0.0, 'N': 0.0, 'D': 2.0, 'C':-4.0, 'Q': 2.0, 'E': 5.0, 'G':-2.0, 'H': 0.0, 'I':-3.0, 'L':-3.0, 'K': 1.0, 'M':-2.0, 'F':-3.0, 'P':-1.0, 'S': 0.0, 'T':-1.0, 'W':-3.0, 'Y':-2.0, 'V':-2.0, 'B': 1.0, 'Z': 4.0, 'X':-1.0, '*':-4.0 },
        'G': { 'A': 0.0, 'R':-2.0, 'N': 0.0, 'D':-1.0, 'C':-3.0, 'Q':-2.0, 'E':-2.0, 'G': 6.0, 'H':-2.0, 'I':-4.0, 'L':-4.0, 'K':-2.0, 'M':-3.0, 'F':-3.0, 'P':-2.0, 'S': 0.0, 'T':-2.0, 'W':-2.0, 'Y':-3.0, 'V':-3.0, 'B':-1.0, 'Z':-2.0, 'X':-1.0, '*':-4.0 },
        'H': { 'A':-2.0, 'R': 0.0, 'N': 1.0, 'D':-1.0, 'C':-3.0, 'Q': 0.0, 'E': 0.0, 'G':-2.0, 'H': 8.0, 'I':-3.0, 'L':-3.0, 'K':-1.0, 'M':-2.0, 'F':-1.0, 'P':-2.0, 'S':-1.0, 'T':-2.0, 'W':-2.0, 'Y': 2.0, 'V':-3.0, 'B': 0.0, 'Z': 0.0, 'X':-1.0, '*':-4.0 },
        'I': { 'A':-1.0, 'R':-3.0, 'N':-3.0, 'D':-3.0, 'C':-1.0, 'Q':-3.0, 'E':-3.0, 'G':-4.0, 'H':-3.0, 'I': 4.0, 'L': 2.0, 'K':-3.0, 'M': 1.0, 'F': 0.0, 'P':-3.0, 'S':-2.0, 'T':-1.0, 'W':-3.0, 'Y':-1.0, 'V': 3.0, 'B':-3.0, 'Z':-3.0, 'X':-1.0, '*':-4.0 },
        'L': { 'A':-1.0, 'R':-2.0, 'N':-3.0, 'D':-4.0, 'C':-1.0, 'Q':-2.0, 'E':-3.0, 'G':-4.0, 'H':-3.0, 'I': 2.0, 'L': 4.0, 'K':-2.0, 'M': 2.0, 'F': 0.0, 'P':-3.0, 'S':-2.0, 'T':-1.0, 'W':-2.0, 'Y':-1.0, 'V': 1.0, 'B':-4.0, 'Z':-3.0, 'X':-1.0, '*':-4.0 },
        'K': { 'A':-1.0, 'R': 2.0, 'N': 0.0, 'D':-1.0, 'C':-3.0, 'Q': 1.0, 'E': 1.0, 'G':-2.0, 'H':-1.0, 'I':-3.0, 'L':-2.0, 'K': 5.0, 'M':-1.0, 'F':-3.0, 'P':-1.0, 'S': 0.0, 'T':-1.0, 'W':-3.0, 'Y':-2.0, 'V':-2.0, 'B': 0.0, 'Z': 1.0, 'X':-1.0, '*':-4.0 },
        'M': { 'A':-1.0, 'R':-1.0, 'N':-2.0, 'D':-3.0, 'C':-1.0, 'Q': 0.0, 'E':-2.0, 'G':-3.0, 'H':-2.0, 'I': 1.0, 'L': 2.0, 'K':-1.0, 'M': 5.0, 'F': 0.0, 'P':-2.0, 'S':-1.0, 'T':-1.0, 'W':-1.0, 'Y':-1.0, 'V': 1.0, 'B':-3.0, 'Z':-1.0, 'X':-1.0, '*':-4.0 },
        'F': { 'A':-2.0, 'R':-3.0, 'N':-3.0, 'D':-3.0, 'C':-2.0, 'Q':-3.0, 'E':-3.0, 'G':-3.0, 'H':-1.0, 'I': 0.0, 'L': 0.0, 'K':-3.0, 'M': 0.0, 'F': 6.0, 'P':-4.0, 'S':-2.0, 'T':-2.0, 'W': 1.0, 'Y': 3.0, 'V':-1.0, 'B':-3.0, 'Z':-3.0, 'X':-1.0, '*':-4.0 },
        'P': { 'A':-1.0, 'R':-2.0, 'N':-2.0, 'D':-1.0, 'C':-3.0, 'Q':-1.0, 'E':-1.0, 'G':-2.0, 'H':-2.0, 'I':-3.0, 'L':-3.0, 'K':-1.0, 'M':-2.0, 'F':-4.0, 'P': 7.0, 'S':-1.0, 'T':-1.0, 'W':-4.0, 'Y':-3.0, 'V':-2.0, 'B':-2.0, 'Z':-1.0, 'X':-2.0, '*':-4.0 },
        'S': { 'A': 1.0, 'R':-1.0, 'N': 1.0, 'D': 0.0, 'C':-1.0, 'Q': 0.0, 'E': 0.0, 'G': 0.0, 'H':-1.0, 'I':-2.0, 'L':-2.0, 'K': 0.0, 'M':-1.0, 'F':-2.0, 'P':-1.0, 'S': 4.0, 'T': 1.0, 'W':-3.0, 'Y':-2.0, 'V':-2.0, 'B': 0.0, 'Z': 0.0, 'X': 0.0, '*':-4.0 },
        'T': { 'A': 0.0, 'R':-1.0, 'N': 0.0, 'D':-1.0, 'C':-1.0, 'Q':-1.0, 'E':-1.0, 'G':-2.0, 'H':-2.0, 'I':-1.0, 'L':-1.0, 'K':-1.0, 'M':-1.0, 'F':-2.0, 'P':-1.0, 'S': 1.0, 'T': 5.0, 'W':-2.0, 'Y':-2.0, 'V': 0.0, 'B':-1.0, 'Z':-1.0, 'X': 0.0, '*':-4.0 },
        'W': { 'A':-3.0, 'R':-3.0, 'N':-4.0, 'D':-4.0, 'C':-2.0, 'Q':-2.0, 'E':-3.0, 'G':-2.0, 'H':-2.0, 'I':-3.0, 'L':-2.0, 'K':-3.0, 'M':-1.0, 'F': 1.0, 'P':-4.0, 'S':-3.0, 'T':-2.0, 'W': 11.0,'Y': 2.0, 'V':-3.0, 'B':-4.0, 'Z':-3.0, 'X':-2.0, '*':-4.0 },
        'Y': { 'A':-2.0, 'R':-2.0, 'N':-2.0, 'D':-3.0, 'C':-2.0, 'Q':-1.0, 'E':-2.0, 'G':-3.0, 'H': 2.0, 'I':-1.0, 'L':-1.0, 'K':-2.0, 'M':-1.0, 'F': 3.0, 'P':-3.0, 'S':-2.0, 'T':-2.0, 'W': 2.0, 'Y': 7.0, 'V':-1.0, 'B':-3.0, 'Z':-2.0, 'X':-1.0, '*':-4.0 },
        'V': { 'A': 0.0, 'R':-3.0, 'N':-3.0, 'D':-3.0, 'C':-1.0, 'Q':-2.0, 'E':-2.0, 'G':-3.0, 'H':-3.0, 'I': 3.0, 'L': 1.0, 'K':-2.0, 'M': 1.0, 'F':-1.0, 'P':-2.0, 'S':-2.0, 'T': 0.0, 'W':-3.0, 'Y':-1.0, 'V': 4.0, 'B':-3.0, 'Z':-2.0, 'X':-1.0, '*':-4.0 },
        'B': { 'A':-2.0, 'R':-1.0, 'N': 3.0, 'D': 4.0, 'C':-3.0, 'Q': 0.0, 'E': 1.0, 'G':-1.0, 'H': 0.0, 'I':-3.0, 'L':-4.0, 'K': 0.0, 'M':-3.0, 'F':-3.0, 'P':-2.0, 'S': 0.0, 'T':-1.0, 'W':-4.0, 'Y':-3.0, 'V':-3.0, 'B': 4.0, 'Z': 1.0, 'X':-1.0, '*':-4.0 },
        'Z': { 'A':-1.0, 'R': 0.0, 'N': 0.0, 'D': 1.0, 'C':-3.0, 'Q': 3.0, 'E': 4.0, 'G':-2.0, 'H': 0.0, 'I':-3.0, 'L':-3.0, 'K': 1.0, 'M':-1.0, 'F':-3.0, 'P':-1.0, 'S': 0.0, 'T':-1.0, 'W':-3.0, 'Y':-2.0, 'V':-2.0, 'B': 1.0, 'Z': 4.0, 'X':-1.0, '*':-4.0 },
        'X': { 'A': 0.0, 'R':-1.0, 'N':-1.0, 'D':-1.0, 'C':-2.0, 'Q':-1.0, 'E':-1.0, 'G':-1.0, 'H':-1.0, 'I':-1.0, 'L':-1.0, 'K':-1.0, 'M':-1.0, 'F':-1.0, 'P':-2.0, 'S': 0.0, 'T': 0.0, 'W':-2.0, 'Y':-1.0, 'V':-1.0, 'B':-1.0, 'Z':-1.0, 'X':-1.0, '*':-4.0 },
        '*': { 'A':-4.0, 'R':-4.0, 'N':-4.0, 'D':-4.0, 'C':-4.0, 'Q':-4.0, 'E':-4.0, 'G':-4.0, 'H':-4.0, 'I':-4.0, 'L':-4.0, 'K':-4.0, 'M':-4.0, 'F':-4.0, 'P':-4.0, 'S':-4.0, 'T':-4.0, 'W':-4.0, 'Y':-4.0, 'V':-4.0, 'B':-4.0, 'Z':-4.0, 'X':-4.0, '*': 1.0 }
    }
    
    def __init__(self, matrix=BLOSUM62):
        self._matrix = matrix        
        
    def score(self, x, y):
        try:
            return self._matrix[x][y]
        except KeyError as ke:
            raise ResidueNotFoundError(ke.message)
      
    @staticmethod
    def parse(string):
        """
        Parse a standard scoring matrix file, where the first row and
        column are residue labels.
        
        @param string: scoring matrix string
        @type string: str
        
        @rtype: L{SimilarityMatrix}
        """
        
        residues = {}
        matrix = {}
        
        for line in string.splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            
            if not residues:
                residues = line.split()
                
            else:
                items = line.split()
                if len(items) != len(residues) + 1:
                    raise ValueError("{0} scoring columns expected".format(len(residues)))
                                    
                try:
                    aa, scores = items[0], map(float, items[1:])
                    matrix[aa] = dict((residues[i], s) for i, s in enumerate(scores))
                except (KeyError, ValueError):
                    raise ValueError("Corrupt scoring matrix")
        
        return SimilarityMatrix(matrix)
    
        
class AbstractAlignmentAlgorithm(object):
    """
    Base class for all sequence alignment algorithms.
    
    This class was designed with simple sequence alignment algorithms in mind.
    Implementors have full control over the behavior of the scoring function and
    the dynamic programming matrix, but implementing things that require
    additional matrices (such as affine gap penalties) might be trickier.
    
    @param scoring: scoring matrix 
    @type scoring: L{AbstractScoringMatrix}
    @param gap: simple gap penalty
    @type gap: float
    """
    
    __metaclass__ = ABCMeta
    
    def __init__(self, scoring=IdentityMatrix(), gap=0):
        
        if not isinstance(scoring, AbstractScoringMatrix):
            raise TypeError(scoring)
        
        self._gap = float(gap)
        self._scoring = scoring
        
    @property
    def gap(self):
        """
        Simple gap penalty
        @rtype: float 
        """
        return self._gap

    @property
    def scoring_matrix(self):
        """
        Scoring matrix 
        @rtype: L{AbstractScoringMatrix} 
        """        
        return self._scoring

    def align(self, query, subject):
        """
        Align two sequences and return the optimal alignment.
        
        @param query: first sequence
        @type query: L{AbstractSequence}
        @param subject: second sequence
        @type subject: L{AbstractSequence}
        
        @rtype: L{AlignmentResult}     
        """
        if query.length == 0 or subject.length == 0:
            raise ValueError("Can't align zero length sequence")        
        
        # working with string sequences results in a massive speed-up
        qseq = ["*"] + self._sequence(query)
        sseq = ["*"] + self._sequence(subject)
        
        # 1. create a dynamic programming matrix
        matrix = []        
        rows, cols = len(query), len(subject)
        
        for i in range(rows + 1):
            matrix.append([])
            for j in range(cols + 1):
                matrix[i].append(None)                    

        # 2. initialize the first row and column  
        self._initialize(matrix)
            
        # fill
        for i in range(1, rows + 1):
            for j in range(1, cols + 1):
                score = self._score(qseq[i], sseq[j])
                self._fill(matrix, i, j, score)
        
        # 3. compute alignment          
        return self._traceback(matrix, query, subject)    
    
    def _sequence(self, s):
        """
        Extract and return the string sequence of {s}.
        
        @param s: sequence object
        @type s: L{AbstractSequence}
        
        @rtype: list of str
        """
        return list(s.sequence)
    
    @abstractmethod
    def _initialize(self, matrix):
        """
        Initialize (typically the first row and column) of the dynamic
        programming matrix.
        
        @param matrix: list (2D)
        @type matrix: list
        """
        pass        
                                    
    def _score(self, residue1, residue2):
        """
        Retrieve the pairwise score of two residues using the current
        scoring matrix.
        
        @rtype: float
        """      
        return self._scoring.score(residue1, residue2)
    
    def _fill(self, matrix, i, j, score):
        """
        Compute and set the best score that leads to cell i,j in the dynamic
        programming matrix: right, down or diagonal.
        
        See also L{AbstractAlignmentAlgorithm._max}.
        
        @param score: pairwise score at matrix[i][j]
        @type score: float
        @return: the best score
        @rtype: float
        """
        
        match = matrix[i-1][j-1] + score
        insertion = matrix[i][j-1] + self._gap                
        deletion = matrix[i-1][j] + self._gap
        
        best = self._max(match, insertion, deletion)
        matrix[i][j] = best
        
        return best

    @abstractmethod    
    def _max(self, match, insertion, deletion):
        """
        Choose the best score among all given possibilities:
        scores for match, insertion and deletion. This will determine
        the direction taken at each step while building the dynamic programming
        matrix (diagonal, down or right). 

        This is an expected notable point of divergence for most sequence
        alignment algorithms.
        """
        pass

    def _traceback(self, m, seq1, seq2):
        """
        Trace back and return the optimal alignment.
        """

        query = []
        subject = []        

        # working with string sequences results in a massive speed-up
        qseq = ["*"] + self._sequence(seq1)
        sseq = ["*"] + self._sequence(seq2)

        i, j = self._terminus(m)
        qstart, start = i, j
        qend, end = i, j
        score = m[i][j]        
        
        while self._expandable(m, i, j):

            if i > 0 and j > 0 and m[i][j] == (m[i-1][j-1] + self._score(qseq[i], sseq[j])):
                query.append(seq1.residues[i])
                subject.append(seq2.residues[j])
                qstart, start = i, j
                i, j = i - 1, j - 1 
                
            elif i > 0 and  m[i][j] == (m[i-1][j] + self._gap):
                query.append(seq1.residues[i])
                subject.append(ResidueInfo(-1, seq2.alphabet.GAP))
                qstart = i
                i = i - 1
                
            elif j > 0 and  m[i][j] == (m[i][j-1] + self._gap):
                query.append(ResidueInfo(-1, seq1.alphabet.GAP))
                subject.append(seq2.residues[j])
                start = j
                j = j - 1
                
            else:
                assert False
                
        query.reverse()
        subject.reverse()
    
        aligned_query = RichSequence(seq1.id, seq1.header, query, seq1.type)
        aligned_subject = RichSequence(seq2.id, seq2.header, subject, seq2.type)
          
        return AlignmentResult(score, aligned_query, aligned_subject, qstart, qend, start, end)     

    @abstractmethod      
    def _terminus(self, matrix):
        """
        Find the coordinates of the optimal alignment's right endpoint in the
        dynamic programming matrix. This is the starting point of a traceback.
        
        @param matrix: the complete dynamic programming matrix
        @type matrix: 2D list
        
        @rtype: 2-tuple (i, j)
        """
        pass
    
    @abstractmethod
    def _expandable(self, i, j):
        """
        Return True if the traceback procedure must not terminate at
        position C{i,j} in the dynamic programming matrix.
        
        @rtype: bool
        """
        pass
    

class GlobalAlignmentAlgorithm(AbstractAlignmentAlgorithm):
    """
    Needleman-Wunsch global sequence alignment.
    """

    def __init__(self, scoring=IdentityMatrix(), gap=0):
        super(GlobalAlignmentAlgorithm, self).__init__(scoring, gap)
            
    def _initialize(self, matrix):
        
        for i in range(len(matrix)):
            matrix[i][0] = self._gap * i
        for j in range(len(matrix[0])):            
            matrix[0][j] = self._gap * j       

    def _max(self, match, insertion, deletion):
        return max(match, insertion, deletion)    
    
    def _terminus(self, matrix):
        
        i = len(matrix) - 1
        j = len(matrix[0]) - 1
        
        return (i, j)        

    def _expandable(self, matrix, i, j):
        return i > 0 or j > 0        
    
class LocalAlignmentAlgorithm(AbstractAlignmentAlgorithm):
    """
    Smith-Waterman local sequence alignment.
    """
        
    START = 0
    """
    Score for initiation of a new local alignment
    """

    def __init__(self, scoring=IdentityMatrix(), gap=-1):
        super(LocalAlignmentAlgorithm, self).__init__(scoring, gap)
            
    def _initialize(self, matrix):
        
        for i in range(len(matrix)):
            matrix[i][0] = LocalAlignmentAlgorithm.START
        for j in range(len(matrix[0])):            
            matrix[0][j] = LocalAlignmentAlgorithm.START      

    def _max(self, match, insertion, deletion):
        return max(match, insertion, deletion, LocalAlignmentAlgorithm.START)        

    def _terminus(self, matrix):
        
        maxi, maxj = 0, 0
        
        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                if matrix[i][j] > matrix[maxi][maxj]:
                    maxi, maxj = i, j
        
        return (maxi, maxj)
    
    def _expandable(self, matrix, i, j):
        return matrix[i][j] != LocalAlignmentAlgorithm.START  
    

class AlignmentResult(object):
    """
    Represents a pairwise sequence alignment result.
    
    @param score: raw alignment score
    @type score: float
    @param query: aligned query sequence (with gaps)
    @type query: L{AbstractSequence}
    @param subject: aligned subject sequence (with gaps)
    @type subject: L{AbstractSequence}     
    @param qstart: query start position
    @type qstart: int
    @param qend: query end position
    @type qend: int
    @param start: subject start position
    @type start: int
    @param end: subject end position
    @type end: int  
    """
    
    def __init__(self, score, query, subject, qstart, qend, start, end):
        
        if not isinstance(query, AbstractSequence):
            raise TypeError(query)
        if not isinstance(subject, AbstractSequence):
            raise TypeError(query)
        
        if not (len(query) == len(subject)):
            raise ValueError("Corrupt alignment")
                        
        self._score = float(score)
        self._query = query
        self._subject = subject
        self._qstart = int(qstart)
        self._qend = int(qend)
        self._start = int(start)
        self._end = int(end)
        self._identicals = 0
        self._gaps = 0        
        self._length = 0
        
        if query.length > 0 and subject.length > 0:

            if not 1 <= qstart <= qend:
                raise ValueError("Invalid query start/end positions")
            if not 1 <= start <= end:
                raise ValueError("Invalid subject start/end positions")
                        
            qgap = query.alphabet.GAP
            sgap = subject.alphabet.GAP
            
            for q, s in zip(query, subject):
                if q.type == qgap or s.type == sgap:
                    self._gaps += 1
                elif q.type == s.type:
                    self._identicals += 1
                    
            self._length = (self._gaps + (qend - qstart + 1) + (end - start + 1)) / 2
            
        else:
            if (score + qstart + qend + start + end) != 0:
                raise ValueError("Corrupt alignment")
            self._length = 0
                 
        
    def __str__(self):
        string = "{0.qstart:5} {0.query.sequence:} {0.qend:<5}\n"
        string += "{0.start:5} {0.subject.sequence:} {0.end:<5}"
        return string.format(self)        
    
    @property
    def is_empty(self):
        """
        Return True if this is an empty alignment (i.e. no matches)
        @rtype: bool
        """
        return self.length == 0 or self.gaps == self.length   
        
    @property
    def score(self):
        """
        Raw alignment score
        @rtype: float
        """
        return self._score
    
    @property
    def query(self):
        """
        Aligned query sequence (with gaps)
        @rtype: L{AbstractSequence}
        """
        return self._query
    
    @property
    def subject(self):
        """
        Aligned subject sequence (with gaps)
        @rtype: L{AbstractSequence}        
        """        
        return self._subject
    
    @property
    def qstart(self):
        """
        Query start position
        @rtype: int
        """        
        return self._qstart
    
    @property
    def qend(self):
        """
        Query end position
        @rtype: int
        """        
        return self._qend
    
    @property
    def start(self):
        """
        Subject start position
        @rtype: int
        """                
        return self._start
    
    @property
    def end(self):
        """
        Subject end position        
        @rtype: int
        """                
        return self._end
    
    @property
    def identicals(self):
        """
        Number of identical residues
        @rtype: int
        """                
        return self._identicals
    
    @property
    def identity(self):
        """
        Percentage of identical residues
        @rtype: int
        """                
        return float(self._identicals) / self._length
    
    @property
    def gaps(self):
        """
        Total number of gaps (query + subject)
        @rtype: int
        """               
        return self._gaps
    
    @property
    def length(self):
        """
        Alignment length (query + subject + gaps / 2)
        @rtype: int
        """               
        return self._length
    
    def alignment(self):
        """
        @return: as a sequence alignment object
        @rtype: L{SequenceAlignment}
        """
        return SequenceAlignment([self.query, self.subject])
                
    
