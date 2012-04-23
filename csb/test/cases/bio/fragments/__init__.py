import random
import csb.test as test
import csb.test.cases.bio.hmm

from csb.bio.fragments import Target, Assignment, TorsionAnglesPredictor
from csb.bio.fragments import FragmentCluster, ClusterNode, ClusterDivergingError, ClusterExhaustedError
from csb.bio.fragments import RosettaFragsetFactory
from csb.bio.fragments.rosetta import RosettaFragment, RosettaFragmentMap

from csb.bio.structure import Atom
from csb.bio.sequence import ProteinAlphabet


class SampleTarget(object):
    
    _instance = None
    
    def __init__(self):
        raise NotImplementedError('you are doing it wrong')
    
    @staticmethod
    def get(fraglength=10):
        
        if SampleTarget._instance is None:
            SampleTarget._instance = SampleTarget.build_target(fraglength)
        
        return SampleTarget._instance

    @staticmethod
    def build_target(fraglength=10):
        
        config = test.Config()
        chain = config.getPickle('1nz9.model1.pickle').chains['A']
        chain.compute_torsion()
        
        target = Target('testA', chain.length, chain.residues)
        
        for start in range(1, target.length + 1, 5):
            
            end = min(start + fraglength - 1,  target.length)
            if end - start + 1 < 2:
                continue
            
            source = chain.clone()
            for atom in source.components(klass=Atom):
                atom.vector[0] += 0.005 * start * atom.vector[0]
                
            rmsd = chain.subregion(start, end).rmsd(source.subregion(start, end))
            
            a = Assignment(source, start, end, id, start, end, probability=1.0, rmsd=rmsd)
            target.assign(a)
        
        return target


@test.unit
class TestTarget(test.Case):
    
    FRAG_LENGTH = 10

    def setUp(self):
        
        super(TestTarget, self).setUp()
        self.target = SampleTarget.get(self.FRAG_LENGTH)
        
    def testFromSequence(self):
        
        t = Target.from_sequence('testA', 'ABC')
        
        self.assertEqual(t.length, 3)
        self.assertEqual(t.sequence, 'ABC')
        self.assertEqual(t.accession, 'test')
        self.assertEqual(t.residues[1].native.type, ProteinAlphabet.ALA)

    def testFromProfile(self):
        
        hmm = csb.test.cases.bio.hmm.build_hmm()
        t = Target.from_profile(hmm)

        self.assertEqual(t.length, 2)
        self.assertEqual(t.sequence, 'AC')
        self.assertEqual(t.id, 'id')
        self.assertEqual(t.residues[1].native.type, ProteinAlphabet.ALA)
        
    def testAccession(self):
        self.assertEqual(self.target.accession, 'test')

    def testEntryID(self):
        self.assertEqual(self.target.id, 'testA')        

    def testLength(self):
        self.assertEqual(self.target.length, 58)
        
    def testSequence(self):
        self.assertEqual(self.target.sequence, 'AQVAFREGDQVRVVSGPFADFTGTVTEINPERGKVKVMVTIFGRETPVELDFSQVVKA')
        
    def testResidues(self):
        self.assertEqual(self.target.residues.length, 58)     
        self.assertEqual(self.target.residues[1].native.type, ProteinAlphabet.ALA)
        self.assertRaises(IndexError, lambda:self.target.residues[0])
    
    def testMatches(self):
        self.assertEqual(self.target.matches.length, 12)     
        self.assertEqual(self.target.matches[0].qstart, 1)
        self.assertEqual(self.target.matches[-1].qend, 58)
        
    def testAssign(self):
        
        a = self.target.matches[0]
        
        for r in range(a.qstart, a.qend + 1):
            
            frags = [rai.fragment for rai in self.target.residues[r].assignments]
            self.assertTrue(a in frags)
            
    def testFilter(self):
        
        t = self.target.filter(threshold=1.5, extend=False)
        self.assertEqual(t.matches.length, self.target.length)
        

@test.unit
class TestTargetResidue(test.Case):

    def setUp(self):
        
        super(TestTargetResidue, self).setUp()
        
        self.target = SampleTarget.get()
        self.residue = self.target.residues[6]
        
    def testType(self):
        self.assertEqual(self.residue.type, ProteinAlphabet.ARG)
        
    def testNative(self):
        self.assertEqual(self.residue.native.type, ProteinAlphabet.ARG)
    
    def testAssignments(self):
        
        self.assertEqual(self.residue.assignments.length, 2)
        self.assertEqual(self.residue.assignments[0].fragment.qstart, 1)
        self.assertEqual(self.residue.assignments[0].fragment.qend, 10)   
        self.assertEqual(self.residue.assignments[1].fragment.qstart, 6)
        self.assertEqual(self.residue.assignments[1].fragment.qend, 15)           
        
    def testAssign(self):
        pass # implicitly covered
        
    def testVeryBest(self):
        b = self.residue.verybest()
        self.assertEqual(b, self.residue.assignments[0].fragment)
        
    def testLongest(self):
        b = self.residue.longest()
        self.assertEqual(b, self.residue.assignments[0].fragment)        
    
    def testPrecision(self):
        self.assertEqual(self.residue.precision(), 100.0)
        self.assertEqual(self.residue.precision(threshold=0.05), 50.0)
        
    def testFilter(self):
        
        rep = self.residue.filter()
        self.assertTrue(rep.centroid in [i.fragment for i in self.residue.assignments])

        
@test.unit
class TestTorsionAnglesPredictor(test.Case):

    def setUp(self):
        
        super(TestTorsionAnglesPredictor, self).setUp()
        
        self.target = SampleTarget.get()
        self.predictor = TorsionAnglesPredictor(self.target)
        
    def testTarget(self):
        self.assertEqual(self.predictor.target, self.target)

    def testThreshold(self):
        self.assertEqual(self.predictor.threshold, 1.5)
        
    def testComputeSingle(self):
        
        native = self.target.residues[3].native
        pred = self.predictor.compute_single(3)
        
        self.assertEqual(pred.torsion.phi, native.torsion.phi)
        self.assertEqual(pred.torsion.psi, native.torsion.psi)
        self.assertEqual(pred.as_tuple(), (0.0, -143.77110677043459, 148.81948663098206, -179.98643191250056))
                
        self.assertEqual(pred.rank, 3)
        self.assertEqual(pred.confidence, self.target.residues[3].filter().confidence)
        self.assertEqual(pred.confidence, 0)
        self.assertTrue(pred.primary)
    
    def testCompute(self):
        
        pred = self.predictor.compute(3)
        self.assertTrue(len(pred) > 1)

        for pi in pred:
            self.assertEqual(pi.as_tuple()[1:], (-143.77110677043459, 148.81948663098206, -179.98643191250056))            
            self.assertTrue(pi.confidence <= pred[0].confidence)
            
            
@test.unit
class TestAssignment(test.Case):

    def setUp(self):
        
        super(TestAssignment, self).setUp()
        
        self.target = SampleTarget.get()
        self.a = self.target.matches[0]
        
    def testBackbone(self):
        
        self.assertEqual(self.a.backbone.length, 10)
        
        x = self.target.residues[1].native.atoms['CA'].vector[0]
        self.assertEqual(self.a.backbone[0][0], x + 0.005 * x)
        self.assertEqual(list(self.a.backbone[0]), [-0.96882, -17.865, -6.46])
        
    def testSequence(self):
        self.assertEqual(self.a.sequence, self.target.sequence[self.a.qstart - 1 : self.a.qend])
        
    def testTorsion(self):
        self.assertEqual(self.a.torsion[0].psi, self.target.residues[1].native.torsion.psi)
        
    def testSourceID(self):
        self.assertEqual(self.a.source_id, '1nz9A')
    
    def testStart(self):
        self.assertEqual(self.a.start, 1)

    def testEnd(self):
        self.assertEqual(self.a.end, 10)
            
    def testQStart(self):
        self.assertEqual(self.a.qstart, 1)
        self.assertEqual(self.target.residues[7].assignments[1].fragment.qstart, 6)

    def testQEnd(self):
        self.assertEqual(self.a.qend, 10)
        self.assertEqual(self.target.residues[7].assignments[1].fragment.qend, 15)
        
    def testRMSD(self):    
        self.assertAlmostEqual(self.a.rmsd, 0.0128, places=4)
    
    def testAnchoredAround(self):
        self.assertTrue(self.a.anchored_around(5))
        self.assertFalse(self.a.anchored_around(1))
        
    def testBackboneAt(self):
        self.assertEqual(list(self.a.backbone_at(1, 1)[0]), [-0.96882, -17.865, -6.46])
        
    def testTorsionAt(self):
        self.assertEqual(self.a.torsion_at(1, 1)[0].psi, self.target.residues[1].native.torsion.psi)
        
    def testOverlap(self):
        overlap = self.a.overlap(self.target.residues[7].assignments[1].fragment)
        self.assertEqual(overlap, set([6, 7, 8, 9, 10]))    
        
    def testRMSDTo(self):
        rmsd = self.a.rmsd_to(self.target.residues[7].assignments[1].fragment)
        self.assertAlmostEqual(rmsd, 0.038, places=4)
        

@test.unit
class TestFragmentCluster(test.Case):
    
    class A(object):
        
        def __init__(self, r):
            self.r = r
            self.length = 8
        
        def rmsd_to(self, other):
            return abs(self.r - other.r)    

    def setUp(self):
        
        super(TestFragmentCluster, self).setUp()
        self.target = SampleTarget.get()
        
        A = TestFragmentCluster.A
        self.c = FragmentCluster([ ClusterNode(A(0.4)), ClusterNode(A(1.3)), ClusterNode(A(1.2)), ClusterNode(A(1.2)), ClusterNode(A(0.1)), ClusterNode(A(0.09), fixed=False) ], threshold=0.1)
        self.fc = FragmentCluster([ ClusterNode(A(0.4)), ClusterNode(A(1.3)), ClusterNode(A(1.2)), ClusterNode(A(1.2)), ClusterNode(A(0.1)), ClusterNode(A(0.09), fixed=True) ], threshold=0.1)
        
    def testCount(self):
        self.assertEqual(self.c.count, 6)

    def testShrink(self):
        
        # covers also ClusterRep and ClusterNode; also the methods: mean, reject, shrinkone
        
        self.assertRaises(ClusterDivergingError, self.fc.shrink)
        self.assertRaises(ClusterExhaustedError, self.c.shrink, minitems=100)
        
        rep = self.c.shrink(minitems=1)
        
        self.assertEqual(rep.centroid.r, 1.2)
        self.assertEqual(rep.count, 3)
        self.assertEqual(rep.rejections, 3)
        self.assertAlmostEqual(rep.mean, 0.0444, places=4)
        self.assertAlmostEqual(self.c.mean(), 0.0444, places=4)
        self.assertAlmostEqual(rep.confidence, 0.2651, places=4)
        self.assertEqual(rep.centroid, self.c.centroid().centroid)
        
        for item in self.c.items:
            self.assertGreater(item.fragment.r, 1.1)
        for item in self.c.fragments:
            self.assertGreater(item.r, 1.1)            
    
    def testMean(self):
        self.assertAlmostEqual(self.c.mean(), 0.5639, places=4)
        
    def testShrinkOne(self):
        
        self.assertTrue(self.c.shrinkone())
        self.assertEqual(self.c.count, 5)
        
        self.c.shrink()
        self.assertFalse(self.c.shrinkone())

@test.unit
class TestRosettaFragsetfactory(test.Case):

    def setUp(self):
        
        super(TestRosettaFragsetfactory, self).setUp()
        
        self.target = SampleTarget.get()
        self.factory = RosettaFragsetFactory()
        
    def testMakeFragset(self):
        
        fragset = self.factory.make_fragset(self.target)
        
        self.assertEqual(len(fragset), self.target.matches.length)
        self.assertEqual(len(fragset.at(3)), 1)
        self.assertEqual(len(fragset.at(7)), 2)
        
    def testMakeChopped(self):
        
        fragset = self.factory.make_fragset(self.target)
        chopped = self.factory.make_chopped(fragset, window=5)
        
        for fragment in chopped:
            self.assertEqual(fragment.length, 5)

    def testMakeCombined(self):
        
        assert self.target.residues[1].assignments.length == 1
                
        fragment = RosettaFragment.from_object(self.target.matches[0])
        filling = RosettaFragmentMap([fragment], self.target.length)
        combined = self.factory.make_combined(self.target, filling)
        
        self.assertEqual(len(combined.starting_at(1)), 2)
        self.assertEqual(len(combined.starting_at(6)), 1)
        
    def testMakeFiltered(self):
        
        fragset = self.factory.make_filtered(self.target)
        self.assertGreater(len(fragset.at(1)), 0)

    def testMix(self):
        
        fragset = self.factory.make_fragset(self.target)
        mixed = self.factory.mix(fragset, fragset, fragset)
        
        self.assertEqual(len(mixed), len(fragset) * 3)


@test.unit
class TestRosettaFragment(test.Case):                   

    def setUp(self):
        
        super(TestRosettaFragment, self).setUp()
        
        self.target = SampleTarget.get()
        self.fragment = RosettaFragment.from_object(self.target.matches[0])
        
    def testSubregion(self):
        
        sub = self.fragment.subregion(2, 4)
        self.assertEqual(sub.length, 3)
        self.assertEqual(sub.qstart, 2)
        self.assertEqual(sub.qend, 4)
        self.assertEqual(sub.residues[0].rank, self.fragment.residues[1].rank)
        self.assertEqual(sub.residues[-1].rank, self.fragment.residues[3].rank)
        
    def testLength(self):
        self.assertEqual(self.fragment.length, 10)
        
    def testFomObject(self):
        fragment = RosettaFragment.from_object(self.target.matches[0])
        self.assertEqual(fragment.length, 10)
        self.assertEqual(fragment.qstart, 1)
        self.assertEqual(fragment.qend, 10)
        
    def testSourceID(self):
        self.assertEqual(self.fragment.source_id, '1nz9A')

    def testAccession(self):
        self.assertEqual(self.fragment.accession, '1nz9')
        
    def testQStart(self):
        self.assertEqual(self.fragment.qstart, 1)
        self.assertEqual(self.target.residues[6].assignments[1].fragment.qstart, 6)
                
    def testStart(self):
        self.assertEqual(self.fragment.start, 1)

    def testQEnd(self):
        self.assertEqual(self.fragment.qend, 10)
        self.assertEqual(self.target.residues[6].assignments[1].fragment.qend, 15)
                
    def testEnd(self):
        self.assertEqual(self.fragment.end, 10)   
    
    def testScore(self):
        self.assertEqual(self.fragment.score, 0.0)
        
    def testTorsion(self):
        self.assertAlmostEqual(self.fragment.torsion[0].psi, 112.12957, places=4)
        
    def testResidues(self):        
        self.assertEqual(len(self.fragment.residues), 10)
        
        first = self.fragment.residues[0]
        self.assertEqual(first.aa, self.target.residues[1].type)
        self.assertEqual(first.rank, 1)
        self.assertAlmostEqual(first.torsion.psi, 112.12957, places=4)
        
                
@test.unit
class TestRosettaFragmentMap(test.Case):                   

    def setUp(self):
        
        super(TestRosettaFragmentMap, self).setUp()
        
        self.target = SampleTarget.get()
        self.map = RosettaFragsetFactory().make_fragset(self.target)

    def testSources(self):
        self.assertEqual(self.map.sources, tuple(['1nz9']))
        
    def testSize(self):
        self.assertEqual(self.map.size, 12)
        
    def testMarkUnconfident(self):
     
        self.map.mark_unconfident(5)
        self.map.mark_unconfident(8)
        
        self.assertEqual(self.map.unconfident_positions, tuple([5, 8]))
        
    def testFromSource(self):
        self.assertEqual(self.map.fromsource('1nz9'), tuple(self.map))
        
    def testStartingAt(self):
        self.assertEqual(len(self.map.starting_at(2)), 0)        
        self.assertEqual(len(self.map.starting_at(6)), 1)    

    def testAt(self):
        self.assertEqual(len(self.map.at(2)), 1)        
        self.assertEqual(len(self.map.at(6)), 2)        

    def testComplement(self):
        
        self.map.mark_unconfident(7)
        for fragment in self.map.at(6):
            self.map.complement(fragment)
        
        self.assertEqual(self.map.size, 13)
    
    def testSort(self):
        
        fragments = [ RosettaFragment.from_object(f) for f in self.target.matches ]
        random.shuffle(fragments)
        
        fragmap = RosettaFragmentMap(fragments)
        fragmap.sort(field='length')
        
        for i, fragment in enumerate(fragmap):
            if i >= 1:
                self.assertTrue(fragment.length >= fragmap[i-1].length)
                
    def testIndexer(self):
        
        self.assertEqual(self.map[0].qstart, 1)
        self.assertEqual(self.map[0].qend, 10)
        
    def testRead(self):
        
        file = self.config.getTestFile('1nz9A.frags')
        fragmap = RosettaFragmentMap.read(file)
        
        self.assertEqual(fragmap.size, self.map.size)
        self.assertEqual(fragmap[1].residues[0].torsion.phi, -150.513)
        self.assertEqual(fragmap[1].residues[1].torsion.psi, 147.798)
        self.assertEqual(fragmap[1].residues[1].rank, 7)
        self.assertEqual(fragmap[1].residues[1].aa, 'E')

    def testDump(self):
        
        with self.config.getTempStream() as tmp:
            
            self.map.dump(tmp.name)
            tmp.flush()
            
            ref = self.config.getContent('1nz9A.frags').rstrip()
            new = open(tmp.name).read().rstrip() 
            
            self.assertEqual(ref, new) 
    
                                 
if __name__ == '__main__':
    
    test.Console()
    