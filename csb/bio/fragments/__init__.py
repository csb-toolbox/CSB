import numpy
import csb.pyutils
import csb.bio.utils


class FragmentMatch(object):
    
    def __init__(self, id, qstart, qend, probability, rmsd, tm_score, qlength):
        
        self._id = id
        self._qstart = qstart
        self._qend = qend
        self._probability = probability        
        self._rmsd = rmsd
        self._tm_score = tm_score
        self._qlength = qlength
        
    @property
    def id(self):
        return self._id

    @property
    def qstart(self):
        return self._qstart
    
    @property
    def qend(self):
        return self._qend
    
    @property
    def qlength(self):
        return self._qlength    
    
    @property
    def rmsd(self):
        return self._rmsd
    
    @property
    def tm_score(self):
        return self._tm_score
    
    @property
    def probability(self):
        return self._probability
            
    @property
    def length(self):
        return self.qend - self.qstart + 1
    
    @property
    def source_id(self):
        raise NotImplementedError()
        
    @property
    def start(self):
        raise NotImplementedError()

    @property
    def end(self):
        raise NotImplementedError()      
    
class PredictionContainer(object):
    
    def __init__(self, target, isites_prediction, hmm_prediction, combined_prediction):
        
        self.target = target
        
        self.isites = isites_prediction
        self.hmm = hmm_prediction
        self.combined = combined_prediction 
        
class Prediction(object):
    
    def __init__(self, alignment, coordinates):
        
        self.alignment = alignment
        self.coordinates = coordinates
        
class Target(csb.pyutils.AbstractNIContainer):
    
    def __init__(self, id, length, residues, overlap=None, segments=None):
    
        self._id = id
        self._accession = id[:4]
        self._chain_id = [4]
        self._length = length
        self._overlap = overlap
        
        self._assignments = csb.pyutils.ReadOnlyCollectionContainer(type=Assignment)
            
        resi = [TargetResidue(native) for native in residues]
        self._residues = csb.pyutils.ReadOnlyCollectionContainer(items=resi, 
                                            type=TargetResidue, start_index=1)
        
        if segments is not None:
            segments = dict([(s.start, s) for s in segments])
        self._segments = csb.pyutils.ReadOnlyDictionaryContainer(items=segments)
        
    @property
    def _children(self):
        return self._residues        
    
    @property
    def id(self):
        return self._id
    
    @property
    def accession(self):
        return self._accession
    
    @property
    def chain_id(self):
        return self._chain_id
    
    @property
    def max_overlap(self):
        return self._overlap
    
    @property
    def length(self):
        return self._length
    
    @property
    def sequence(self):
        return ''.join(str(r.native.type) for r in self)
    
    @property
    def matches(self):
        return self._assignments
    
    @property
    def residues(self):
        return self._residues

    @property
    def segments(self):
        return self._segments
        
    def assign(self, fragment):

        if not 1 <= fragment.qstart <= fragment.qend <= len(self._residues):
            raise ValueError("Fragment out of range")
                    
        for rank in range(fragment.qstart, fragment.qend + 1):
            ai = ResidueAssignmentInfo(fragment, rank)
            self._assignments._append_item(fragment)
            self._residues[rank].assign(ai)
            
        if fragment.segment is not None:
            try:
                self._segments[fragment.segment].assign(fragment)
            except KeyError:
                raise ValueError("Undefined segment starting at {0}".format(fragment.segment))
            
class TargetResidue(object):
    
    def __init__(self, native_residue):
        
        self._type = native_residue.type
        self._native = native_residue.clone()
        self._assignments = csb.pyutils.ReadOnlyCollectionContainer(type=ResidueAssignmentInfo)

    @property
    def type(self):
        return self._type
            
    @property
    def native(self):
        return self._native
    
    @property
    def assignments(self):
        return self._assignments
    
    def assign(self, assignment_info):
        self._assignments._append_item(assignment_info)
        
    
class TargetSegment(object):
    
    def __init__(self, start, end, count):
        
        self._start = start
        self._end = end
        self._count = count
        
        self._assignments = csb.pyutils.ReadOnlyCollectionContainer(type=Assignment)
    
    @property
    def count(self):
        return self._count
    
    @property
    def start(self):
        return self._start
    
    @property
    def end(self):
        return self._end

    @property
    def length(self):
        return (self._end - self._start + 1)
        
    @property
    def assignments(self):
        return self._assignments
    
    def assign(self, fragment):
        if fragment.segment != self.start:
            raise ValueError('Segment origin mismatch: {0} vs {1}'.format(fragment.segment, self.start))
        else:
            self._assignments._append_item(fragment)
            
    def pairwise_rmsd(self):
        
        rmsd  = []
        
        for q in self.assignments:
            for s in self.assignments:
                if q is not s:
                    rmsd.append(q.rmsd(s))
                else:
                    assert q.rmsd(s) < 0.01
        
        return [ r for r in rmsd if r is not None ]
            
class ResidueAssignmentInfo(object):
    
    def __init__(self, assignment, rank):
        
        if not assignment.qstart <= rank <= assignment.qend:
            raise ValueError('Rank {0} is not matched by this assignment')
        
        self._assignment = assignment
        self._rank = rank
        self._relrank = rank - assignment.qstart 
    
    @property
    def c_alpha(self):
        return self._assignment.backbone[self._relrank]
    
    @property
    def fragment(self):
        return self._assignment    

class Assignment(FragmentMatch):
    
    def __init__(self, source, start, end, id, qstart, qend, probability, rmsd, tm_score, 
                 score=None, neff=None, segment=None):

        sub = source.subregion(start, end, clone=True)
        calpha = [r.atoms['CA'].vector.copy() for r in sub.residues]
                    
        self._calpha = csb.pyutils.ReadOnlyCollectionContainer(items=calpha, type=numpy.ndarray)     
        
        self._source_id = source.entry_id
        self._start = start
        self._end = end

        self._score = score
        self._neff = neff
    
        self._segment_start = segment
        
        super(Assignment, self).__init__(id, qstart, qend, probability, rmsd, tm_score, None)

    @property
    def backbone(self):
        return self._calpha
        
    @property
    def source_id(self):
        return self._source_id
        
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
    def neff(self):
        return self._neff
        
    @property
    def segment(self):
        return self._segment_start    
    
    def transform(self, rotation, translation):
        
        for ca in self.backbone:
            newca = numpy.dot(ca, numpy.transpose(rotation)) + translation
            for i in range(3):
                ca[i] = newca[i]
                
    def backbone_at(self, qstart, qend):
        
        if not (self.qstart <= qstart <= qend <= self.qend):
            raise ValueError('Region {0}..{1} is out of range {2.qstart}..{2.qend}'.format(qstart, qend, self))
        
        relstart = qstart - self.qstart
        relend = qend - self.qstart + 1
        
        return self.backbone[relstart : relend]
    
    def rmsd(self, other, min_overlap=5):

        qranks = set(range(self.qstart, self.qend + 1))
        sranks = set(range(other.qstart, other.qend + 1))
        common = qranks.intersection(sranks)
        
        if len(common) >= min_overlap:
        
            qstart, qend = min(common), max(common)
            
            q = self.backbone_at(qstart, qend)
            s = other.backbone_at(qstart, qend)
            
            if len(q) > 0 and len(s) > 0:
                return csb.bio.utils.rmsd(q, s)
            
        return None
    
class FragmentTypes(object):
    
    ISites = 'IS'
    HMMFragments = 'HH'
    HHThread = 'TH'
    HHfrag = HHThread
    Rosetta = 'NN'        
    
class BenchmarkAdapter(object):
    
    class Connection(object):
    
        FACTORY = None
        DSN = None
        
        def __init__(self, factory=None, dsn=None):
    
            self.factory = factory or self.__class__.FACTORY
            self.cs = dsn or self.__class__.DSN
            self.connection = None
            self.cursor = None
    
        def __enter__(self):
    
            self.connection = self.factory(self.cs)
            try:
                self.cursor = self.connection.cursor()
            except:
                self.connection.close()
                raise 
            return self
    
        def __exit__(self, *args):
            try:
                if not self.cursor.closed:
                    self.cursor.close()
            finally:
                if not self.connection.closed:
                    self.connection.close()

    def __init__(self, pdb_paths, connection_string=None):
                
        self._pdb = pdb_paths
        self._connection = None
        
        from csb.bio.io.wwpdb import find, StructureParser
        self._parser = StructureParser
        self._find = find     
    
        try:    
            import psycopg2.extras
        except ImportError:
            raise RuntimeError('Please install the psycopg2 module first')
        
        if connection_string is None:
            connection_string = self.connection_string()
            
        BenchmarkAdapter.Connection.FACTORY = psycopg2.extras.DictConnection
        BenchmarkAdapter.Connection.DSN = connection_string
        
    @staticmethod
    def connection_string(database='FragmentBenchmarks', host='', username='', password=''):

        fields = ['dbname={0}'.format(database)]
        
        if host:
            fields.append('host={0}'.format(host))
        if username:
            fields.append('user={0}'.format(username))
            fields.append('password={0}'.format(password))

        return ' '.join(fields)
        
    def targets(self, benchmark_id):
            
        with BenchmarkAdapter.Connection() as db:
            
            db.cursor.callproc('reporting."GetTargets"', (benchmark_id,))
            return db.cursor.fetchall()            
            
    def target_details(self, target_id):
        
        with BenchmarkAdapter.Connection() as db:
            
            db.cursor.callproc('reporting."GetTargetDetails"', (target_id,))
            return db.cursor.fetchall()            
            
    def assignments(self, target_id, type):
        
        with BenchmarkAdapter.Connection() as db:
            
            db.cursor.callproc('reporting."GetAssignments"', (target_id, type))
            return db.cursor.fetchall()

    def scores(self, benchmark_id, type):
        
        with BenchmarkAdapter.Connection() as db:
            
            db.cursor.callproc('reporting."GetScores"', (benchmark_id, type))
            return db.cursor.fetchall()    
            
    def target_segments(self, target_id):
        
        with BenchmarkAdapter.Connection() as db:
            
            db.cursor.callproc('reporting."GetTargetSegments"', (target_id,))
            data = db.cursor.fetchall()
            
            return [ TargetSegment(row['Start'], row['End'], row['Count']) for row in data ]
            
    def structure(self, accession):

        pdbfile = self._find(accession, self._pdb)
        
        if not pdbfile:
            raise IOError('{0} not found here: {1}'.format(accession, self._pdb))
        
        return self._parser(pdbfile).parse_structure()
            
    def prediction(self, target_id, type):
        
        info = self.target_details(target_id)
        if not info:
            raise ValueError('No such Target ID in the database: {0}'.format(target_id))
        row = info[0]
        
        id = row["Accession"]
        length = float(row["Length"])
        overlap = float(row["MaxOverlap"]) / (length or 1.)
        
        native = self.structure(id[:4]).chains[id[4]]
        segments = self.target_segments(target_id)
        target = Target(id, length, native.residues, overlap, segments)        
        
        source = None
        
        for row in self.assignments(target_id, type):
            
            src_accession = row['Source'][:4]
            src_chain = row['Source'][4]
            
            if source is None or source.accession != src_accession:
                source = self.structure(src_accession)
                
            frag_chain = source.chains[src_chain]
            
            fragment = Assignment(source=frag_chain, 
                                  start=row['SourceStart'], 
                                  end=row['SourceEnd'], 
                                  id=row['FragmentName'], 
                                  qstart=row['Start'], 
                                  qend=row['End'], 
                                  probability=row['Probability'], 
                                  score=row['Score'],
                                  neff=row['Neff'],
                                  rmsd=row['RMSD'], 
                                  tm_score=row['TMScore'],
                                  segment=row['SegmentStart'])
            
            target.assign(fragment)
        
        return target
            
    
  
if __name__ == '__main__':
     
    fb = BenchmarkAdapter('/media/DB/pdb/all')
    target = fb.prediction(361, FragmentTypes.HHfrag)
    
    print target.segments.length
    
    for seg_start in target.segments:
        print '    ', seg_start, target.segments[seg_start].count
        rmsd = target.segments[seg_start].pairwise_rmsd()
    
            