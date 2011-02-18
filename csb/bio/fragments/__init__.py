import os
import cPickle
import numpy
import csb.pyutils
import csb.bio.utils

class FragmentTypes(object):
    
    ISites = 'IS'
    HMMFragments = 'HH'
    HHThread = 'TH'
    HHfrag = HHThread
    Rosetta = 'NN'   
    
class Metrics(object):
    
    RMSD = 'rmsd_to'
    MDA = 'mda_to' 
    
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
        self._chain_id = id[4]
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
        
    def verybest(self):
        
        best = None
        
        for ai in self.assignments:
            a = ai.fragment
            if a.length < 5:
                continue
            if best is None or a.rmsd < best.rmsd:
                best = a
            elif a.rmsd == best.rmsd and a.length > best.length:
                best = a
        
        return best
                
    def filter(self, method=Metrics.RMSD):
        
        try:
            assignments = [ ai.fragment for ai in self.assignments ]
            cluster = FragmentCluster(assignments, threshold=1.5,
                                      connectedness=0.5, method=method)
          
            centroid = cluster.shrink(minitems=0)
            return centroid
        
        except ClusterExhaustedError:
            return None
    
    def longest(self):
    
        best = None        
        
        for q in self.assignments:
            if best is None or (q.fragment.length > best.length):
                best = q.fragment
                
        return best     
    
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
                
    def verybest(self):
        
        best = None
        
        for a in self.assignments:
            if a.length < 5:
                continue
            if best is None or a.rmsd < best.rmsd:
                best = a
            elif a.rmsd == best.rmsd and a.length > best.length:
                best = a
        
        return best
                
    def best(self, method=Metrics.RMSD):
        
        try:
            cluster = FragmentCluster(self.assignments, threshold=1.5,
                                      connectedness=0.5, method=method)          
            centroid = cluster.shrink(minitems=1)
            return centroid
                  
        except ClusterExhaustedError:
            return None
        finally:
            del cluster
    
    def longest(self):
    
        best = None        
        
        for q in self.assignments:
            if best is None or (q.length > best.length):
                best = q
                
        return best
            
    def pairwise_rmsd(self, min_overlap=5):
        
        rmsds  = []
        
        for q in self.assignments:
            for s in self.assignments:
                if q is not s:
                    r = q.rmsd_to(s, min_overlap)
                    if r is not None:
                        rmsds.append(r)
                else:
                    assert q.rmsd_to(s, 1) < 0.01
        
        return rmsds
    
    def pairwise_mda(self, min_overlap=5):
        
        mdas  = []
        
        for q in self.assignments:
            for s in self.assignments:
                if q is not s:
                    m = q.mda_to(s, min_overlap)
                    if m is not None:
                        mdas.append(m)
        return mdas

    def pairwise_sa_rmsd(self, profiles='.', min_overlap=5):
        
        from csb.bio.hmm import RELATIVE_SA
        from csb.bio.io.hhpred import ScoreUnits, HHProfileParser

        def convert_sa(sa):
            return numpy.array([ RELATIVE_SA[i] for i in sa ])
            
        sources = {}        
        scores  = []
        
        for q in self.assignments:
            for s in self.assignments:
                
                if s.source_id not in sources:
                    hmm = HHProfileParser(os.path.join(profiles, s.source_id + '.hhm')).parse()
                    sources[s.source_id] = hmm.dssp_solvent
                    
                if q is not s:
                    
                    common = q.overlap(s)
                    if len(common) >= min_overlap:
                        
                        qsa = q.solvent_at(sources[q.source_id], min(common), max(common))
                        ssa = s.solvent_at(sources[s.source_id], min(common), max(common))
                        
                        if '-' in qsa + ssa:
                            continue
                        
                        qsa = convert_sa(qsa)
                        ssa = convert_sa(ssa)
                        assert len(qsa) == len(ssa)
                        sa_rmsd = numpy.sqrt(numpy.sum((qsa - ssa) ** 2) / float(len(qsa)))
                        
                        scores.append(sa_rmsd)        
        return scores           
    
    def pairwise_scores(self, profiles='.', min_overlap=5):
        
        from csb.bio.hmm import BACKGROUND
        back = numpy.sqrt(numpy.array(BACKGROUND))

        sources = {}        
        scores  = []
        
        for q in self.assignments:
            for s in self.assignments:
                
                if s.source_id not in sources:
                    # hmm = HHProfileParser(os.path.join(hmm_path, s.source_id + '.hhm')).parse(ScoreUnits.Probability)
                    sources[s.source_id] = cPickle.load(open(os.path.join(profiles, s.source_id + '.pkl')))
                    
                if q is not s:
                    
                    common = q.overlap(s)
                    if len(common) >= min_overlap:
                        
                        qprof = q.profile_at(sources[q.source_id], min(common), max(common))
                        sprof = s.profile_at(sources[s.source_id], min(common), max(common))
                        
                        #score = qhmm.emission_similarity(shmm)
                        assert len(qprof) == len(sprof)
                        dots = [ numpy.dot(qprof[i] / back, sprof[i] / back) for i in range(len(qprof)) ]
                        score = numpy.log(numpy.prod(dots))
                        if score is not None:
                            scores.append(score)        
        return scores       
    
    def _entropy(self, data, binsize):
        
        binsize = float(binsize)
        bins = numpy.ceil(numpy.array(data) / binsize)

        hist = dict.fromkeys(bins, 0)
        for bin in bins:
            hist[bin] += (1.0 / len(bins))
        
        freq = numpy.array(hist.values())
        return -numpy.sum(freq * numpy.log(freq))
    
    def rmsd_entropy(self, binsize=0.1):
        
        rmsds = self.pairwise_rmsd()
        return self._entropy(rmsds, binsize)
    
    def score_entropy(self, profiles='.', binsize=1):
        
        scores = self.pairwise_scores(profiles)
        return self._entropy(scores, binsize)    
    
    def rmsd_consistency(self, threshold=1.5):

        rmsds = self.pairwise_rmsd()
        
        if len(rmsds) < 1:
            return None
        
        return sum([1 for i in rmsds if i <= threshold]) / float(len(rmsds))

    def sa_rmsd_consistency(self, threshold=0.4, profiles='.'):

        sa_rmsds = self.pairwise_sa_rmsd(profiles=profiles)
        
        if len(sa_rmsds) < 1:
            return None
        
        return sum([1 for i in sa_rmsds if i <= threshold]) / float(len(sa_rmsds))
        
    def true_positives(self, threshold=1.5):
        
        if self.assignments.length < 1:
            return None
        
        return sum([1 for i in self.assignments if i.rmsd <= threshold]) / float(self.assignments.length)
    
    def confidence(self):
        
        cons = self.rmsd_consistency()
        
        if cons is None:
            return 0
        else:
            return numpy.log10(self.count) * cons
            
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

        assert source.has_torsion
        sub = source.subregion(start, end, clone=True)
        calpha = [r.atoms['CA'].vector.copy() for r in sub.residues]
        torsion = [r.torsion.copy() for r in sub.residues]
                    
        self._calpha = csb.pyutils.ReadOnlyCollectionContainer(items=calpha, type=numpy.ndarray)
        self._torsion = torsion     
        
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
    def torsion(self):
        return self._torsion
            
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

    def _check_range(self, qstart, qend):

        if not (self.qstart <= qstart <= qend <= self.qend):
            raise ValueError('Region {0}..{1} is out of range {2.qstart}..{2.qend}'.format(qstart, qend, self))
                                
    def backbone_at(self, qstart, qend):
        
        self._check_range(qstart, qend)
        
        relstart = qstart - self.qstart
        relend = qend - self.qstart + 1
        
        return self.backbone[relstart : relend]
    
    def torsion_at(self, qstart, qend):
        
        self._check_range(qstart, qend)
        
        relstart = qstart - self.qstart
        relend = qend - self.qstart + 1
        
        return self.torsion[relstart : relend]    
    
    def solvent_at(self, sa_string, qstart, qend):
        
        self._check_range(qstart, qend)
        
        relstart = qstart - self.qstart
        relend = qend - self.qstart + 1
        
        return sa_string[relstart : relend]    
    
    def profile_at(self, source, qstart, qend):
        
        self._check_range(qstart, qend)
        
        start = qstart - self.qstart + self.start
        end = qend - self.qstart + self.start
        
        if hasattr(source, 'subregion'):
            return source.subregion(start, end)
        else:
            return source[start - 1 : end]
    
    def overlap(self, other):

        qranks = set(range(self.qstart, self.qend + 1))
        sranks = set(range(other.qstart, other.qend + 1))
        
        return qranks.intersection(sranks)
            
    def rmsd_to(self, other, min_overlap=5):

        common = self.overlap(other)
        
        if len(common) >= min_overlap:
        
            qstart, qend = min(common), max(common)
            
            q = self.backbone_at(qstart, qend)
            s = other.backbone_at(qstart, qend)
            
            if len(q) > 0 and len(s) > 0:
                return csb.bio.utils.rmsd(q, s)
            
        return None
    
    def mda_to(self, other, min_overlap=5):

        common = self.overlap(other)
        
        if len(common) >= min_overlap:
        
            qstart, qend = min(common), max(common)
            
            q = self.torsion_at(qstart, qend)
            s = other.torsion_at(qstart, qend)
            
            if len(q) > 0 and len(s) > 0:
                
                maxphi = max( numpy.abs(i.phi-j.phi) for i, j in zip(q, s)[1:] )   # phi: 2 .. L
                maxpsi = max( numpy.abs(i.psi-j.psi) for i, j in zip(q, s)[:-1] )  # psi: 1 .. L-1
                
                return max(maxphi, maxpsi)
            
        return None  

class ClusterExhaustedError(ValueError):
    pass    
class ClusterEmptyError(ClusterExhaustedError):
    pass      
    
class FragmentCluster(object):

    def __init__(self, items, threshold=1.5, connectedness=0.5, method=Metrics.RMSD):

        items = set(i for i in items if i.length > 5)
        
        self._matrix = {}        
        self.threshold = float(threshold)
        self.connectedness = float(connectedness)
           
        for i in items:
            
            self._matrix[i] = {}
            conn = 0.0
            
            for j in items:
                distance = getattr(i, method)(j)
                if distance is not None:
                    conn += 1
                    self._matrix[i][j] = distance
            
            if conn / len(items) < self.connectedness:
                del self._matrix[i]
                
                
        self._items = set(self._matrix.keys())
                     
        if len(self._items) < 1:
            raise ClusterEmptyError()           
               
    @property
    def count(self):
        return len(self._items)
    
    @property    
    def items(self):
        return tuple(self._items)
    
    def _distances(self, skip=None):

        d = []
        
        for i in self._matrix:
            if skip is i:
                continue

            for j in self._matrix[i]:
                if skip is not j:
                    d.append(self._matrix[i][j])
                    
        return d

    def mean(self, skip=None):
        
        d = self._distances(skip=skip) 
            
        if len(d) > 0:
            return numpy.mean(d)
        else:
            raise ClusterExhaustedError()

    def centroid(self):

        cen = None
        avg = None

        for i in self._matrix:
            
            curravg = numpy.mean(self._matrix[i].values())
            
            if avg is None or curravg < avg:
                avg = curravg
                cen = i
            elif curravg == avg:
                if i.length > cen.length:
                    cen = i
    
        d = self._distances()
        mean = numpy.mean(d)
        cons = sum(1.0 for i in d if i <= self.threshold) / len(d)
        
        return CentroidInfo(cen, mean, cons, self.count)

    def reject(self, item):
        
        if self.count == 1:
            raise ClusterExhaustedError()
        
        for i in self._matrix:
            if item in self._matrix[i]:
                del self._matrix[i][item]
            
        del self._matrix[item]
        self._items.remove(item)

    def shrinkone(self):
        
        mean = self.mean()
        if mean <= self.threshold or self.count == 1:
            return False                # already shrunk enough

        m = {}
        
        for i in self._matrix:
            newmean = self.mean(skip=i)
            m[newmean] = i

        newmean = min(m)
        assert newmean <= mean

        if newmean < mean:            
            junk = m[newmean]
            self.reject(junk)
            return True                 # successful shrink
        else:
            return False                # converged
        
    def shrink(self, minitems=2):

        if self.count > minitems:
            
            while self.shrinkone():
                if self.count <= minitems:
                    raise ClusterExhaustedError()
        else:
            raise ClusterExhaustedError()
            
        return self.centroid()
    
class CentroidInfo(object):
    
    def __init__(self, centroid, mean, consistency, count):
        
        self.centroid = centroid
        self.mean = mean
        self.consistency = consistency
        self.count = count
    
    @property
    def confidence(self):
        if self.count <= 0 or self.count is None or self.consistency is None:
            return 0
        else:
            return numpy.log10(self.count) * self.consistency
            
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
            
    def structure(self, accession, chain=None):

        pdbfile = self._find(accession, self._pdb)
        
        if not pdbfile and chain:
            pdbfile = self._find(accession + chain, self._pdb)
                    
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
        
        native = self.structure(id[:4], id[4]).chains[id[4]]
        segments = self.target_segments(target_id)
        target = Target(id, length, native.residues, overlap, segments)        
        
        source = None
        
        for row in self.assignments(target_id, type):
            
            src_accession = row['Source'][:4]
            src_chain = row['Source'][4]
            
            if source is None or source.accession != src_accession:
                source = self.structure(src_accession, src_chain)
                
            frag_chain = source.chains[src_chain]
            if not frag_chain.has_torsion:
                frag_chain.compute_torsion()
            
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

    