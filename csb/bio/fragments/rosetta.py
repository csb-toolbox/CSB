"""
Rosetta fragment libraries.

This module defines the L{RosettaFragmentMap} objects, which describes a
fragment library in Rosetta NNmake format. L{RosettaFragmentMap} has a 
static factory method for building a library from a fragment file:
    
    >>> RosettaFragmentMap.read('fragments.txt')
    <RosettaFragmentMap>
    
@note: Consider extracting L{RosettaFragmentMap.read} as a Rosetta
       fragment parser which naturally belongs to csb.bio.io.      
"""

from csb.bio.structure import TorsionAnglesCollection, TorsionAngles
from csb.core import AbstractContainer


class ResidueInfo(object):
    """
    Container struct for a single rosetta fragment residue.
    
    @param rank: residue position (in the source chain, 1-based)
    @type rank: int
    @param aa: amino acid
    @type aa: str
    @param ss: secondary structure
    @type ss: str
    @param torsion: torsion angles
    @type torsion: L{csb.bio.structure.TorsionAngles} 
    """
    
    def __init__(self, rank, aa, ss, torsion, calpha=[]):
        
        self.rank = rank
        self.aa = aa
        self.ss = ss
        self.torsion = torsion
        self.calpha = tuple(calpha)
        
    @property
    def phi(self):
        return self.torsion.phi or 0.

    @property
    def psi(self):
        return self.torsion.psi or 0.
    
    @property
    def omega(self):
        return self.torsion.omega or 0.
    
    def copy(self):
        """
        @return: a deep copy of the struct
        @rtype: L{ResidueInfo}
        """
        return ResidueInfo(self.rank, self.aa, self.ss, self.torsion.copy(), self.calpha)
    
class RosettaFragment(object):
    """
    Represents a single Rosetta fragment match.

    @param source_id: entry ID of the source PDB chain (in accnC format)
    @type source_id: str
    @param qstart: start position in target (rank)
    @type qstart: int
    @param qend: end position in target (rank)
    @type qend: int
    @param start: start position in C{source} (rank)
    @type start: int
    @param end: end position in C{source} (rank)
    @type end: int    
    @param score: score of the fragment
    @type score: float
    @param residues: fragment residue structs
    @type residues: iterable of L{ResidueInfo}    
    """
    
    def __init__(self, source_id, qstart, qend, start, end, score, residues):
        
        if not (qend - qstart + 1) == (end - start + 1) == len(residues):
            raise ValueError()
        if not len(source_id) == 5:
            raise ValueError(source_id)            
        
        self._source_id = str(source_id)
        self._qstart = int(qstart)
        self._qend = int(qend)
        self._start = int(start)
        self._end = int(end)
        self._score = float(score)
        self._residues = list(residues)
        
    def subregion(self, qstart, qend):
        """
        Extract a subregion from the fragment.
        
        @param qstart: start position in target
        @type qstart: int 
        @param qend: end position in target
        @type qend: int  
        
        @return: a new fragment (deep copy)
        @rtype: L{RosettaFragment}      
        """
        
        if not self.qstart <= qstart <= qend <= self.qend:
            raise ValueError('Invalid subregion')
        
        start = qstart - self.qstart + self.start
        end = qend - self.qend + self.end
        
        diff = qstart - self.qstart
        size = qend - qstart + 1
        assert 0 <= diff 
        
        residues = [ r.copy() for r in self.residues[diff : diff + size] ]
        assert len(residues) == size 
        
        return RosettaFragment(self.source_id, qstart, qend, start, end, self.score, residues)
        
    def __lt__(self, other):
        # lower score means a better fragment
        return self.score < other.score
    
    def __iter__(self):
        return iter(self._residues)
    
    def __len__(self):
        return len(self._residues)

    def __str__(self):

        out = []

        for residue in self.residues:
            line = ' {0.accession:4} {0.chain:1} {1.rank:>5} {1.aa:1} {1.ss:1} {1.phi:>8.3f} {1.psi:>8.3f} {1.omega:>8.3f} {0.score:>8.3f}'
            out.append(line.format(self, residue))
            
        return '\n'.join(out)
    
    @staticmethod
    def from_object(assignment):
        """
        Factory method: build a rosetta fragment from an assignment object.
        
        @param assignment: source assignment
        @type assignment: L{Assignment} 
        
        @rtype: L{RosettaFragment}
        """
        residues = []        
        a = assignment
        
        for rank, aa, torsion, calpha in zip(range(a.start, a.end + 1), a.sequence, a.torsion, a.backbone):
            residues.append(ResidueInfo(rank, aa, 'L', torsion, calpha))
            
        return RosettaFragment(a.source_id, a.qstart, a.qend, a.start, a.end, 1 - (a.probability or 0.0), residues)
    
    @property
    def length(self):
        return len(self)
    
    @property
    def source_id(self):
        return self._source_id

    @property
    def accession(self):
        return self.source_id[:4]
    
    @property
    def chain(self):
        return self.source_id[4:]
    
    @property
    def id(self):
        return '{0.source_id}:{0.start}-{0.end}'.format(self)

    @property
    def qstart(self):
        return self._qstart
    
    @property
    def qend(self):
        return self._qend
        
    @property
    def start(self):
        return self._start
    
    @property
    def end(self):
        return self._end

    @property
    def score(self):
        return self._score
    
    @property
    def residues(self):
        return tuple(self._residues)
    
    @property
    def torsion(self):
        return TorsionAnglesCollection([r.torsion for r in self._residues], start=0)    
        
class OutputBuilder(object):
    """
    Rosetta fragment file formatter.
    
    @param output: destination stream
    @type output: file
    """
    
    def __init__(self, output):
        self._out = output
        
    @property
    def output(self):
        return self._out
        
    def add_position(self, qstart, frags):
        """
        Write a new assignment origin.
        
        @param qstart: target position
        @type qstart: float
        @param frags: number of fragments, starting at that position
        @type frags: int 
        """
        self.output.write(' position: {0:>12} neighbors: {1:>12}\n\n'.format(qstart, len(frags)))
    
    def add_fragment(self, fragment):
        """
        Write a new fragment.
        @type fragment: L{RosettaFragment} 
        """
        for residue in fragment.residues:
            self.add_residue(fragment, residue)
            self.output.write('\n')
            
        self.output.write('\n')
        
    def add_residue(self, fragment, residue):
        """
        Write a new fragment residue.
        @type fragment: L{RosettaFragment}
        @type residue: L{ResidueInfo}
        """        
        line = ' {0.accession:4} {0.chain:1} {1.rank:>5} {1.aa:1} {1.ss:1} {1.phi:>8.3f} {1.psi:>8.3f} {1.omega:>8.3f} {0.score:>8.3f}'
        self.output.write(line.format(fragment, residue))
    
class ExtendedOutputBuilder(OutputBuilder):
    """
    Builds non-standard fragment files, which contain the CA coordinates of
    each residue at the end of each line.
    """
    
    def add_residue(self, fragment, residue):        
        
        super(ExtendedOutputBuilder, self).add_residue(fragment, residue)

        if residue.calpha:
            calpha = residue.calpha
        else:
            calpha = [0, 0, 0]
        
        self.output.write('        {0:>7.3f} {1:>7.3f} {2:>7.3f}'.format(*calpha))            
        
class RosettaFragmentMap(AbstractContainer):
    """
    Represents a Rosetta fragment library.
    
    @param fragments: library fragments
    @type fragments: iterable of L{RosettaFragment}
    @param length: target sequence's length. If not defined, the qend of the
                   last fragment will be used instead.
    @type length: int
    """
    
    def __init__(self, fragments, length=None):
        
        self._fragments = []

        self._unconf = set()
        self._sources = set()
        self._starts = set()
        self._ends = set()
        self._length = None
        
        for f in fragments:
            self.append(f)
            
        if length is not None:
            assert length >= self._maxend
            self._length = int(length)
        else:
            self._length = self._maxend
            
    @property
    def _maxend(self):
        return max(self._ends or [0])
            
    def append(self, fragment):
        """
        Append a new L{RosettaFragment}
        """
        
        if self._length and fragment.qend > self._length:
            raise ValueError('fragment out of range')
        
        self._fragments.append(fragment)
        self._sources.add(fragment.accession)
        self._starts.add(fragment.qstart)
        self._ends.add(fragment.qend)     
            
    def __len__(self):
        return len(self._fragments)
            
    @property
    def _children(self):
        return self._fragments
    
    @property
    def unconfident_positions(self):
        return tuple(sorted(self._unconf))
    
    @property
    def size(self):
        return len(self) 
        
    @property
    def sources(self):
        return tuple(self._sources)

    @property
    def start_positions(self):
        return tuple(sorted(self._starts))    

    def fromsource(self, accession):
        """
        @return: a tuple of all fragments, extracted from the specified C{source}.
        
        @param accession: source entry ID
        @type accession: str
        """
        return tuple(f for f in self._fragments if f.accession == accession)
    
    def starting_at(self, qrank):
        """
        @return: a tuple of all fragments, starting at the specified target position.
        
        @param qrank: fragment origin (in target, rank)
        @type qrank: int
        """        
        return tuple(f for f in self._fragments if f.qstart == qrank)
    
    def at(self, qrank):
        """
        @return: a tuple of all fragments, covering the specified position.
        
        @param qrank: position in target, rank
        @type qrank: int
        """            
        return tuple(f for f in self._fragments if f.qstart <= qrank <= f.qend)
    
    def mark_unconfident(self, rank):
        """
        Mark the specified position in the target as a low-confidence one.

        @param rank: position in target        
        @type rank: int 
        """
        if not 1 <= rank <= self._length:
            raise ValueError(rank)
        
        self._unconf.add(rank)
        
    def complement(self, fragment):
        """
        Append C{fragment} to the library, if the fragment is anchored
        around a low-confidence position.
        
        @type fragment: L{RosettaFragment} 
        """        
        if not self._unconf:
            raise ValueError('no unconfident regions to complement')
        
        f = fragment
        for rank in self._unconf:
            if f.qstart < rank < f.qend:
                if (rank - f.qstart + 1) > 0.4 * (f.qend - f.qstart + 1):
                    self.append(f)
                    break
    
    def sort(self, field='score', reverse=False):
        """
        Sort the fragments in the library.
        """
        
        self._fragments.sort(key=lambda i:getattr(i, field), reverse=reverse)
        
    def dump(self, file, builder=OutputBuilder):
        """
        Write the library to a Rosetta fragment file.
        
        @param file: destination file name
        @type file: str
        """
                       
        with open(file, 'w') as out:
            builder = builder(out)
            
            for qstart in self.start_positions:

                frags = self.starting_at(qstart)
                builder.add_position(qstart, frags)
                
                for fragment in frags:
                    builder.add_fragment(fragment)
    
    @staticmethod
    def read(file, top=None):
        """
        Read a standard fragment file.
        
        @param file: file name
        @type file: str
        @param top: if defined, read only C{top} fragments per start position
                    (default=all)
        @type top: int or None
        
        @return: parsed fragment library
        @rtype: L{RosettaFragmentMap}
        """
        # This is the format (rosetta_fragments/nnmake/makeoutput.f):
        # source chain rank residue ss phi psi omega    score dme dme_f best_nn_ss_type     dipolar+noe 'P' position 'F' fragment#
        
        def ang(a):
            a = float(a)
            if a < -180: return 360 + a
            elif a > 180: return -360 + a
            else: return a
        
        frags = []
    
        qstart, qend, start, end = 0, 0, 0, 0
        id = ''
        score = None
        count = 0
        residues = []
        
        in_frag = False
        
        for line in open(file):
            
            if line.startswith(' position:'):
                qstart = int(line.split()[1])
                count = 0
            
            elif not line.strip():
                if in_frag:
                    count += 1
                    if not top or count <= top:            
                        frags.append(RosettaFragment(id, qstart, qend, start, end, score, residues))
                in_frag = False                
                id = ''
                start = 0
                end = 0
                score = None
                residues = []
                
            else:
                fields = line.split()
                if not in_frag:
                    start = int(fields[2])
                    qend = qstart
                    in_frag = True
                else:
                    qend += 1
                                    
                end = int(fields[2])
                id = fields[0].lower() + fields[1]
                score = float(fields[8])
                rank = int(fields[2])
                aa, ss = fields[3:5]
                phi, psi, omega = map(ang, fields[5:8])
                residues.append(ResidueInfo(rank, aa, ss, TorsionAngles(phi, psi, omega)))
        
        return RosettaFragmentMap(frags)
