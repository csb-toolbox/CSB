import os
import cPickle
import numpy

import csb.pyutils
import csb.bio.utils
import csb.bio.structure
import csb.bio.sequence


class FragmentTypes(object):
    
    ISites = 'IS'
    HMMFragments = 'HH'
    HHThread = 'TH'
    HHfrag = HHThread
    Rosetta = 'NN'   
    
class Metrics(object):
    
    RMSD = 'rmsd_to'
    NORMALIZED_RMSD = 'nrmsd_to'
    MDA = 'mda_to'
    
RANDOM_RMSD = {  5: 1.8749005857255376,  6: 2.4314283686276261,  7: 2.9021135267789608,  8: 3.2477716200172715,  9: 3.5469606556031708, 10: 3.8295465524456329, 
                11: 4.1343107114131783, 12: 4.3761697929053014, 13: 4.6707299668248394, 14: 4.9379016881069733, 15: 5.1809028645084911, 16: 5.4146957142595662, 
                17: 5.7135948448156988, 18: 5.9597935432566782, 19: 6.1337340535741962, 20: 6.3962825155503271, 21: 6.6107937773415166, 22: 6.8099096274123401, 
                23: 7.0435583846849639, 24: 7.2160956482560970, 25: 7.4547896324594962, 26: 7.6431870072434211, 27: 7.8727812194173836, 28: 8.0727393298443637, 
                29: 8.2551450998965326, 30: 8.4413583511786587, 31: 8.5958719774122052, 32: 8.7730435506242408, 33: 8.9970648837941649, 34: 9.1566521405105163, 
                35: 9.2828620878454728, 36: 9.4525824357923405, 37: 9.6322126445253300, 38: 9.7851684750961176, 39: 9.9891454649821476, 40: 10.124373939352028, 
                41: 10.284348528344765, 42: 10.390457305096271, 43: 10.565792044674239, 44: 10.676532740033737, 45: 10.789537132283652, 46: 11.004475543757550, 
                47: 11.064541647783571, 48: 11.231219875286985, 49: 11.319222637391441, 50: 11.485478165340824, 51: 11.607522494435521, 52: 11.700268836069840, 
                53: 11.831245255954073, 54: 11.918975893263905 }     
    
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
        
class TorsionAnglesPredictor(object):
    """
    Fragment-based phi/psi angles predictor.
    
    @param target: target protein, containing fragment assignments
    @type target: L{Target}
    @param threshold: RMSD distance threshold for L{FragmentCluster}-based filtering
    @type threshold: float
    @param extend: pick alternative, longer cluster reps, if possible
    @type extend: bool
    @param init: populate all L{FragmentCluster}s on instantiation. If False, this step
                 will be performed on demand (the first time C{predictor.compute()} is invoked)
                 
    @note: if C{init} is False, the first call to C{predictor.compute()} might take a long
           time. Subsequent calls will be very fast.
    """
    
    def __init__(self, target, threshold=1.5, extend=False, init=False):
        
        if not isinstance(target, Target):
            raise TypeError(target)
        if target.matches.length == 0:
            raise ValueError('This target has no fragment assignments')
        
        self._target = target
        self._threshold = float(threshold)
        self._extend = bool(extend)
        
        self._initialized = False
        self._reps = {}
            
        if init:
            self.init()
        
    @property
    def target(self):
        return self._target

    @property
    def threshold(self):
        return self._threshold
    
    @property
    def extend(self):
        return self._extend
    
    def init(self):
        """
        Compute and cache all L{FragmentCluster}s.
        """

        self._reps = {}
                
        for residue in self.target.residues:
            rep = residue.filter(threshold=self.threshold, extend=self.extend)
            
            if rep is not None:
                self._reps[residue.native.rank] = rep
                
        self._initialized = True      
       
    def _residue(self, rank):
        
        for r in self._target.residues:
            if r.native.rank == rank:
                return r
        
        raise ValueError('Rank {0} is out of range'.format(rank))
    
    def compute_single(self, rank):
        """
        Extract torsion angles from the L{ClusterRep} at residue C{#rank}.
        
        @param rank: target residue rank
        @type rank: int
        
        @rtype: L{TorsionPredictionInfo} 
        """
        
        residue = self._residue(rank)
        rep = residue.filter(threshold=self.threshold, extend=self.extend)
        
        if rep is None:
            return None
        
        else:
            fragment = rep.centroid
            torsion = fragment.torsion_at(rank, rank)[0]
            
            return TorsionPredictionInfo(rank, rep.confidence, torsion, primary=True)    
            
    def compute(self, rank):
        """
        Extract torsion angles from all L{ClusterRep}s, covering residue C{#rank}.
        
        @param rank: target residue rank
        @type rank: int
        
        @return: a tuple of L{TorsionPredictionInfo}, sorted by confidence  
        @rtype: tuple  
        """        
        
        if not self._initialized:
            self.init()
        
        residue = self._residue(rank)
        prediction = []
        
        for rep in self._reps.values():
            
            if rep.centroid.qstart <= residue.native.rank <= rep.centroid.qend:
                
                fragment = rep.centroid
                torsion = fragment.torsion_at(rank, rank)[0]
                info = TorsionPredictionInfo(rank, rep.confidence, torsion)

                if rep is self._reps.get(rank, None):
                    info.primary = True
                    
                prediction.append(info)
        
        prediction.sort(reverse=True)
        return tuple(prediction)

class TorsionPredictionInfo(object):
    """
    Struct container for a single torsion angle prediction.
    
    @param rank: target residue rank
    @type rank: int
    @param confidence: confidence of prediction
    @type confidence: float
    @param torsion: assigned phi/psi/omega angles
    @type torsion: L{TorsionAngles}
    @param primary: if True, designates that the assigned angles are extracted
                    from the L{ClusterRep} at residue C{#rank}; otherwise: the
                    angles are coming from another, overlapping L{ClusterRep}
    
    """
    
    def __init__(self, rank, confidence, torsion, primary=False):
        
        self.rank = rank
        self.confidence = confidence
        self.torsion = torsion
        self.primary = primary
            
    def as_tuple(self):
        """
        @return: convert this prediction to a tuple: (confidence, phi, psi, omega)
        @rtype: tuple
        """        
        return tuple([self.confidence, self.torsion.phi, self.torsion.psi, self.torsion.omega])
    
    def __str__(self):
        return '<TorsionPredictionInfo: {0.confidence:6.3f} at #{0.rank}>'.format(self)
    
    def __cmp__(self, other):
        return cmp(self.confidence, other.confidence)


class AssignmentFactory(object):
    
    def target(self, *a, **k):
        return Target(*a, **k)
    
    def residue(self, *a, **k):
        return TargetResidue(*a, **k)
    
    def assignment(self, *a, **k):
        return Assignment(*a, **k)
            
class Target(csb.pyutils.AbstractNIContainer):
    
    def __init__(self, id, length, residues, overlap=None, segments=None, factory=AssignmentFactory()):
    
        self._id = id
        self._accession = id[:-1]
        self._chain_id = id[-1]
        self._length = length
        self._overlap = overlap
        self._factory = factory
        
        self._assignments = csb.pyutils.ReadOnlyCollectionContainer(type=Assignment)
        self._errors = csb.pyutils.CollectionContainer()
            
        resi = [factory.residue(native) for native in residues]
        self._residues = csb.pyutils.ReadOnlyCollectionContainer(items=resi,
                                            type=TargetResidue, start_index=1)
        
        if segments is not None:
            segments = dict([(s.start, s) for s in segments])
        self._segments = csb.pyutils.ReadOnlyDictionaryContainer(items=segments)
        
    @staticmethod
    def from_sequence(id, sequence):
        
        if isinstance(sequence, csb.bio.sequence.Sequence):
            sequence = sequence.sequence
        
        residues = []
        
        for rn, aa in enumerate(sequence, start=1):
            residue = csb.bio.structure.ProteinResidue(rank=rn, type=aa)
            residues.append(residue)
            
        return Target(id, len(residues), residues)
    
    @staticmethod
    def from_profile(hmm):
        
        residues = [ r.clone() for r in hmm.residues ]
        return Target(hmm.id, hmm.layers.length, residues)
    
    @staticmethod
    def deserialize(pickle):
        
        with open(pickle) as stream:
            return cPickle.load(stream)
    
    @property
    def _children(self):
        return self._residues       
    
    @property
    def errors(self):
        return self._errors 
    
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
        
        self._assignments._append_item(fragment)
                    
        for rank in range(fragment.qstart, fragment.qend + 1):
            ai = ResidueAssignmentInfo(fragment, rank)
            self._residues[rank].assign(ai) 
            
        if fragment.segment is not None:
            try:
                self._segments[fragment.segment].assign(fragment)
            except KeyError:
                raise ValueError("Undefined segment starting at {0}".format(fragment.segment))
            
    def assignall(self, fragments):
        
        for frag in fragments:
            self.assign(frag)         
    
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
            if a.length < FragmentCluster.MIN_LENGTH:
                continue
            if best is None or a.rmsd < best.rmsd:
                best = a
            elif a.rmsd == best.rmsd and a.length > best.length:
                best = a
        
        return best
                
    def filter(self, method=Metrics.RMSD, threshold=1.5, extend=False):
        
        try:
            nodes = []
            for ai in self.assignments:
                if ai.fragment.probability > 0.7 and ai.fragment.length >= FragmentCluster.MIN_LENGTH:
                    node = ClusterNode(ai.fragment, distance=method, fixed=extend)
                else:
                    node = ClusterNode(ai.fragment, distance=method, fixed=False)
                nodes.append(node)
            cluster = FragmentCluster(nodes, threshold=threshold)
          
            center = cluster.shrink(minitems=0)
            if center.has_alternative:
                center.exchange()
            return center
        
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
            if a.length < FragmentCluster.MIN_LENGTH:
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
        
        rmsds = []
        
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
        
        mdas = []
        
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
        scores = []
        
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
        scores = []
        
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
        return - numpy.sum(freq * numpy.log(freq))
    
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
                 score=None, neff=None, segment=None, internal_id=None):

        assert source.has_torsion
        sub = source.subregion(start, end, clone=True)
        try:
            calpha = [r.atoms['CA'].vector.copy() for r in sub.residues]
        except csb.pyutils.ItemNotFoundError:
            raise csb.bio.structure.Broken3DStructureError()
        torsion = [r.torsion.copy() for r in sub.residues]
                    
        self._calpha = csb.pyutils.ReadOnlyCollectionContainer(items=calpha, type=numpy.ndarray)
        self._torsion = torsion
        self._sequence = sub.sequence
        
        self._source_id = source.accession[:4] + source.id 
        self._start = start
        self._end = end

        self._score = score
        self._neff = neff
    
        self._segment_start = segment
        self.internal_id = internal_id
        
        super(Assignment, self).__init__(id, qstart, qend, probability, rmsd, tm_score, None)

    @property
    def backbone(self):
        return self._calpha
    
    @property
    def sequence(self):
        return self._sequence

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
        
    def anchored_around(self, rank):
        
        if self.qstart < rank < self.qend:
            if (rank - self.qstart + 1) > 0.4 * (self.qend - self.qstart + 1):
                return True
        
        return False
                                
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
        
    def chain_at(self, source, qstart, qend):
        
        self._check_range(qstart, qend)
        
        start = qstart - self.qstart + self.start
        end = qend - self.qstart + self.start
        
        return source.subregion(start, end)
    
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
    
    def nrmsd_to(self, other, min_overlap=5):
        
        common = self.overlap(other)
        
        if len(common) >= min_overlap:
        
            qstart, qend = min(common), max(common)
            
            q = self.backbone_at(qstart, qend)
            s = other.backbone_at(qstart, qend)
            
            if len(q) > 0 and len(s) > 0:
                return csb.bio.utils.rmsd(q, s) / RANDOM_RMSD[ len(common) ]
            
        return None
    
    def mda_to(self, other, min_overlap=5):

        common = self.overlap(other)
        
        if len(common) >= min_overlap:
        
            qstart, qend = min(common), max(common)
            
            q = self.torsion_at(qstart, qend)
            s = other.torsion_at(qstart, qend)
            
            if len(q) > 0 and len(s) > 0:
                
                maxphi = max(numpy.abs(i.phi - j.phi) for i, j in zip(q, s)[1:])   # phi: 2 .. L
                maxpsi = max(numpy.abs(i.psi - j.psi) for i, j in zip(q, s)[:-1])  # psi: 1 .. L-1
                
                return max(maxphi, maxpsi)
            
        return None  
    
    def to_rosetta(self, source, qstart=None, qend=None, weight=None):
        
        from cStringIO import StringIO
        stream = StringIO()
        
        if weight is None:
            weight = self.probability
        if not qstart:
            qstart = self.qstart
        if not qend:
            qend = self.qend            
        
        source.compute_torsion()
        chain = self.chain_at(source, qstart, qend)
        
        for i, r in enumerate(chain.residues):
                
            acc = self.source_id[:4]
            ch = self.source_id[4].upper()

            start = qstart - self.qstart + self.start + i
            aa = r.type
            ss = 'L'
            phi, psi, omega = 0, 0, 0
            if r.torsion.phi:
                phi = r.torsion.phi
            if r.torsion.psi:
                psi = r.torsion.psi
            if r.torsion.omega:
                omega = r.torsion.omega
            
            stream.write(' {0:4} {1:1} {2:>5} {3!s:1} {4!s:1} {5:>8.3f} {6:>8.3f} {7:>8.3f} {8:>8.3f}\n'.format(acc, ch, start, aa, ss, phi, psi, omega, weight))            
        
        return stream.getvalue()    

class ClusterExhaustedError(ValueError):
    pass    
class ClusterEmptyError(ClusterExhaustedError):
    pass
class ClusterDivergingError(RuntimeError):
    pass
    
class FragmentCluster(object):
    
    MIN_LENGTH = 6

    def __init__(self, items, threshold=1.5, connectedness=0.5):

        items = set(i for i in items if i.fragment.length >= FragmentCluster.MIN_LENGTH)
        
        self._matrix = {}        
        self.threshold = float(threshold)
        self.connectedness = float(connectedness)
           
        for i in items:
            
            self._matrix[i] = {}
            conn = 0.0
            
            for j in items:
                distance = i.distance(j)
                if distance is not None:
                    conn += 1
                    self._matrix[i][j] = distance
            
            if conn / len(items) < self.connectedness:
                # reject i as a first class node
                del self._matrix[i]
                
        self._items = set(self._matrix.keys())
                     
        if len(self._items) < 1:
            raise ClusterEmptyError()
        
        self._initcount = self.count 
               
    @property
    def count(self):
        return len(self._items)

    @property    
    def items(self):
        return tuple(self._items)
        
    @property    
    def fragments(self):
        return [i.fragment for i in self._items]
    
    def _distances(self, skip=None):

        d = []
        
        for i in self._matrix:
            if skip is i:
                continue

            for j in self._matrix[i]:
                if skip is not j:
                    d.append(self._matrix[i][j])
                    
        return d
    
    def _distance(self, i, j):
        
        if j in self._matrix[i]:
            return self._matrix[i][j]
        else:
            return None

    def mean(self, skip=None):
        
        d = self._distances(skip=skip) 
            
        if len(d) > 0:
            return numpy.mean(d)
        else:
            raise ClusterExhaustedError()

    def centroid(self):

        alt = None
        cen = None
        avg = None

        for i in self._matrix:
            
            curravg = numpy.mean(self._matrix[i].values())
            
            if avg is None or curravg < avg:
                avg = curravg
                cen = i
            elif curravg == avg:
                if i.fragment.length > cen.fragment.length:
                    cen = i
    
        d = self._distances()
        mean = numpy.mean(d)
        cons = sum(1.0 for i in d if i <= self.threshold) / len(d)
        
        for i in self._matrix:
            if i is not cen and i.fixed and i.fragment.length > cen.fragment.length:
                distance = self._distance(i, cen)
                if distance is not None and distance < 0.5 * self.threshold:
                    if alt is None or alt.fragment.length < i.fragment.length:
                        alt = i

        return ClusterRep(cen, mean, cons, len(self._matrix[cen]), alternative=alt,
                            rejections=(self._initcount - self.count))

    def reject(self, item):
        
        if self.count == 1:
            raise ClusterExhaustedError()
        
        assert not item.fixed
        
        for i in self._matrix:
            if item in self._matrix[i]:
                del self._matrix[i][item]
            
        del self._matrix[item]
        self._items.remove(item)

    def shrinkone(self):
        
        mean = self.mean()
        if mean <= self.threshold or self.count == 1:
            return False                  # already shrunk enough

        m = {}
        
        for i in self._matrix:
            if not i.fixed:
                newmean = self.mean(skip=i)
                m[newmean] = i

        if len(m) == 0:                   # only fixed items remaining
            raise ClusterExhaustedError()

        newmean = min(m)

        if newmean > mean:
            raise ClusterDivergingError() # can't converge, usually when fixed items are too far away from the average         
        elif newmean < mean:            
            junk = m[newmean]
            self.reject(junk)
            return True                   # successful shrink
        else:
            return False                  # converged
        
    def shrink(self, minitems=2):

        if self.count > minitems:

            while self.shrinkone():
                if self.count <= minitems:
                    raise ClusterExhaustedError()
        else:
            raise ClusterExhaustedError()
            
        return self.centroid()
    
class ClusterNode(object):
    
    def __init__(self, fragment, distance=Metrics.RMSD, fixed=False):
        
        if fixed and fragment.length < FragmentCluster.MIN_LENGTH:
            raise ValueError("Can't fix a short fragment")
        
        self.fragment = fragment
        self.fixed = bool(fixed)
                
        self._distance = getattr(self.fragment, distance)
        
    def distance(self, other):
        return self._distance(other.fragment)
    
class ClusterRep(object):
    
    def __init__(self, centroid, mean, consistency, count, rejections=0, alternative=None):
        
        if isinstance(centroid, ClusterNode):
            centroid = centroid.fragment
        if isinstance(alternative, ClusterNode):
            alternative = alternative.fragment
                    
        self._centroid = centroid
        self._alternative = alternative
        self._mean = mean
        self._consistency = consistency
        self._count = count
        self._rejections = rejections
    
    @property
    def confidence(self):
        if self.count <= 0 or self.count is None or self.consistency is None:
            return 0
        else:
            return numpy.log10(self.count) * self.consistency
        
    @property
    def centroid(self):
        return self._centroid
    
    @property
    def alternative(self):
        return self._alternative
    
    @property    
    def has_alternative(self):
        return self._alternative is not None        
    
    @property
    def mean(self):
        return self._mean
    
    @property
    def consistency(self):
        return self._consistency
    
    @property
    def count(self):
        return self._count
    
    @property
    def rejections(self):
        return self._rejections
    
    def exchange(self):
        
        if self._alternative is not None:

            centroid = self._centroid
            self._centroid = self._alternative
            self._alternative = centroid

    def to_rosetta(self, source):
        return self.centroid.to_rosetta(source, weight=self.confidence)
            
class AdaptedAssignment(object):
    
    @staticmethod
    def with_overhangs(center, start, end, overhang=1):
        
        if center.centroid.qstart <= (start - overhang):
            start -= overhang
        elif center.centroid.qstart < start:
            start = center.centroid.qstart
            
        if center.centroid.qend >= (end + overhang):
            end += overhang
        elif center.centroid.qend > end:
            end = center.centroid.end
                        
        return AdaptedAssignment(center, start, end)
    
    def __init__(self, center, qstart, qend):
        
        if qstart < center.centroid.qstart:
            raise ValueError(qstart)
        if qend > center.centroid.qend:
            raise ValueError(qend)
                
        self._qstart = qstart
        self._qend = qend
        self._center = center
        
    @property
    def fragment(self):
        return self._center.centroid

    @property
    def center(self):
        return self._center
        
    @property
    def confidence(self):
        return self._center.confidence
    
    @property
    def qstart(self):
        return self._qstart
    
    @property
    def qend(self):
        return self._qend
    
    @property
    def backbone(self):
        return self.fragment.backbone_at(self.qstart, self.qend)
    
    def chain(self, source):
        return self.fragment.chain_at(source, self.qstart, self.qend) 
    
    def to_rosetta(self, source):
        return self.fragment.to_rosetta(source, self.qstart, self.qend, self.confidence)
        
class SmoothFragmentMap(csb.pyutils.AbstractContainer):
    
    def __init__(self, length, centroids):
        
        if not length > 0:
            raise ValueError(length)
        
        self._length = int(length)  
        self._slots = set(range(1, self._length + 1))
        self._map = {}
    
        centers = list(centroids)
        centers.sort(key=lambda i: i.confidence, reverse=True)
        
        for c in centers:
            self.assign(c)
        
    @property
    def _children(self):
        return self._map
    
    def assign(self, center):
        
        for r in range(center.centroid.qstart, center.centroid.qend + 1):
            if r in self._slots:
                self._map[r] = center
                self._slots.remove(r)
    
    def patches(self):            
            
        center = None
        start = None
        end = None
        
        for r in range(1, self._length + 1):
            
            if center is None:
                if r in self._map:
                    center = self._map[r]
                    start = end = r
                else:
                    center = None
                    start = end = None
            else:
                if r in self._map:
                    if self._map[r] is center:
                        end = r
                    else:
                        yield AdaptedAssignment(center, start, end)
                        center = self._map[r]
                        start = end = r
                else:
                    yield AdaptedAssignment(center, start, end)
                    center = None
                    start = end = None
    

class ResidueEventInfo(object):
    
    def __init__(self, rank, confidence=None, count=None, confident=True, rep=None):
        
        self.rank = rank
        self.confidence = confidence
        self.confident = confident
        self.count = count
        self.rep = rep
    
class RosettaFragsetFactory(object):
    
    def __init__(self):
        import csb.bio.fragments.rosetta as rosetta
        self.rosetta = rosetta
    
    def make_fragset(self, target):
        
        frag_factory = self.rosetta.RosettaFragment
        fragments = map(frag_factory.from_object, target.matches)
        #fragments = [ frag_factory.from_object(f) for f in target.matches if f.length >= 6 ]
        fragments.sort()
                
        return self.rosetta.RosettaFragmentMap(fragments, target.length) 
    
    def make_chopped(self, fragments, window):
        
        frags = []
        
        for f in fragments:
            for qs in range(f.qstart, f.qend - window + 1):
                frags.append(f.subregion(qs, qs + window - 1))
        
        return self.rosetta.RosettaFragmentMap(frags)
    
    def make_combined(self, target, filling, threshold=0.5, callback=None):
        
        fragmap = self.make_fragset(target)
        covered = set()
        
        for r in target.residues:
            
            if r.assignments.length == 0:
                if callback:
                    callback(ResidueEventInfo(r.native.rank, None, 0, False))
                continue
            
            cluster = r.filter()
            if cluster is None:
                if callback:
                    callback(ResidueEventInfo(r.native.rank, 0, 0, False))                
                continue

            if cluster.confidence >= threshold:
                covered.add(r.native.rank)
            elif callback:
                callback(ResidueEventInfo(r.native.rank, cluster.confidence, cluster.count, False))
                
        for r in target.residues:
            if r.native.rank not in covered:               # true for gaps and low-conf residues
                fragmap.mark_unconfident(r.native.rank)
                
        for frag in filling:
            fragmap.complement(frag)
            
        return fragmap
    
    def make_filtered(self, target, extend=False, callback=None):
        
        fragments = []
        
        for r in target.residues:
            if r.assignments.length == 0:
                continue    
            
            cluster = r.filter(extend=extend)
            if cluster is None:
                continue
            
            if extend and cluster.has_alternative:
                best = cluster.alternative
            else:
                best = cluster.centroid
                
            fragment = self.rosetta.RosettaFragment.from_object(best)
            fragments.append(fragment)
            if callback:
                callback(ResidueEventInfo(r.native.rank, cluster.confidence, cluster.count, rep=cluster.centroid))
        
        fragments.sort()
        return self.rosetta.RosettaFragmentMap(fragments, target.length)
    
    def mix(self, *fragsets):
        
        fragments = []
        length = 0
        
        for fragset in fragsets:
            if fragset._length > length:
                length = fragset._length  
            
            for fragment in fragset:
                fragments.append(fragment)
                
        return self.rosetta.RosettaFragmentMap(fragments, length)        
                    
            
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

    def __init__(self, pdb_paths, connection_string=None, factory=AssignmentFactory()):
                
        self._pdb = pdb_paths
        self._connection = None
        
        from csb.bio.io.wwpdb import find, StructureParser
        self._parser = StructureParser
        self._find = find
        self._factory = factory
    
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

    def centroids(self, benchmark_id):
        
        with BenchmarkAdapter.Connection() as db:
            
            db.cursor.callproc('reporting."GetCentroids"', (benchmark_id,))
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
        target = self._factory.target(id, length, native.residues, overlap, segments)        
        
        source = None
        
        for row in self.assignments(target_id, type):
            
            src_accession = row['Source'][:4]
            src_chain = row['Source'][4]
            
            if source is None or source.accession != src_accession:
                try:
                    source = self.structure(src_accession, src_chain)
                except IOError as ioe:
                    target.errors.append(ioe)
                    continue
            
            if src_chain == '_':
                frag_chain = source.first_chain
            else:
                frag_chain = source.chains[src_chain]
            if not frag_chain.has_torsion:
                frag_chain.compute_torsion()
            
            fragment = self._factory.assignment(
                                        source=frag_chain,
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
                                        segment=row['SegmentStart'],
                                        internal_id=row['InternalID'])
            
            target.assign(fragment)
        
        return target


    
if __name__ == '__main__':
    
    class A(object):
        
        def __init__(self, r):
            self.r = r
            self.length = 8
        
        def __repr__(self):
            return 'A' + str(self.r)
        
        def rmsd_to(self, other):
            return abs(self.r - other.r)
        
    c = FragmentCluster([ ClusterNode(A(0.4)), ClusterNode(A(1.3)), ClusterNode(A(1.2)), ClusterNode(A(1.2)), ClusterNode(A(0.1), fixed=0), ClusterNode(A(0.09), fixed=1) ], threshold=0.1)
    print c.mean()
    ci = c.shrink(minitems=0)
    print '==============='
    print ci.mean, ci.rejections, ci.centroid
    print c.fragments
    