from csb.bio.structure import TorsionAnglesCollection, TorsionAngles
from csb.pyutils import AbstractContainer


class ResidueInfo(object):
    
    def __init__(self, rank, aa, ss, torsion):
        
        self.rank = rank
        self.aa = aa
        self.ss = ss
        self.torsion = torsion
        
    @property
    def phi(self):
        return self.torsion.phi or 0.

    @property
    def psi(self):
        return self.torsion.psi or 0.
    
    @property
    def omega(self):
        return self.torsion.omega or 0.       

class RosettaFragment(object):
    
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
        
    def __cmp__(self, other):
        # lower score means a better fragment
        return cmp(self.score, other.score)
    
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
        
        residues = []        
        a = assignment
        
        for rank, aa, torsion in zip(range(a.start, a.end + 1), a.sequence, a.torsion):
            residues.append(ResidueInfo(rank, aa, 'L', torsion))
            
        return RosettaFragment(a.source_id, a.qstart, a.qend, a.start, a.end, 1 - a.probability, residues)
    
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
        return '{0.source_id}.{0.start}-{0.end}'.format(self)

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
        return TorsionAnglesCollection(r.torsion for r in self._residues)    
        
class RosettaFragmentMap(AbstractContainer):
    
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
            assert length >= max(self._ends) 
            self._length = int(length)
        else:
            self._length = max(self._ends)
            
    def append(self, fragment):
        
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
    def size(self):
        return len(self) 
        
    @property
    def sources(self):
        return tuple(self._sources)

    @property
    def start_positions(self):
        return tuple(sorted(self._starts))    

    def fromsource(self, accession):
        return tuple(f for f in self._fragments if f.accession == accession)
    
    def starting_at(self, qrank):
        return tuple(f for f in self._fragments if f.qstart == qrank)
    
    def at(self, qrank):
        return tuple(f for f in self._fragments if f.qstart <= qrank <= f.qend)
    
    def mark_unconfident(self, rank):
        
        if not 1 <= rank <= self._length:
            raise ValueError(rank)
        
        self._unconf.add(rank)
        
    def complement(self, fragment):
        
        if not self._unconf:
            raise ValueError('no unconfident regions to complement')
        
        f = fragment
        for rank in self._unconf:
            if f.qstart < rank < f.qend:
                if (rank - f.qstart + 1) > 0.4 * (f.qend - f.qstart + 1):
                    self.append(f)
                    break
    
    def sort(self, field='score', reverse=False):
        
        self._fragments.sort(key=lambda i:getattr(i, field), reverse=reverse)
        
    def dump(self, file):
                
        with open(file, 'w') as out:

            for qstart in self.start_positions:
                
                frags = self.starting_at(qstart)
                out.write(' position: {0:>12} neighbors: {1:>12}\n\n'.format(qstart, len(frags)))
                
                for fragment in frags:
                    out.write(str(fragment))
                    out.write('\n\n')
    
    @staticmethod
    def read(file, top=None):
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
