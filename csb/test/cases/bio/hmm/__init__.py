"""
@note: HHpredHitList and HHpredHit are pretty much data containers rather than
       API objects. They are already covered quite well by csb.test.cases.bio.io.hhpred.
"""

import csb.test as test

from csb.pyutils import Enum

from csb.bio.hmm import State, Transition, ProfileHMM, HMMLayer, ProfileLength, StateFactory, ProfileHMMSegment,\
    StateExistsError, UnobservableStateError, EmissionExistsError
from csb.bio.hmm import HHpredHitList, HHpredHit, ScoreUnits, States, EVDParameters, BACKGROUND
from csb.bio.sequence import ProteinAlphabet
from csb.bio.structure import ProteinResidue



def build_hmm():
    
        hmm = ProfileHMM(units=ScoreUnits.Probability)
        
        factory = StateFactory()
        background = { ProteinAlphabet.ALA: 0.02457563, ProteinAlphabet.CYS: 0.00325358, ProteinAlphabet.GLU: 0.01718016 }
        emission = dict( (aa, 1.0 / i) for i, aa in enumerate(background, start=1) )
        
        # States
        #  Start / End
        hmm.start = State(States.Start)
        hmm.start_insertion = factory.create_insertion(background)
        hmm.end = State(States.End)
        #  L1
        match1 = factory.create_match(emission, background)
        insertion1 = factory.create_insertion(background)
        deletion1 = factory.create_deletion()
        states = {States.Match: match1, States.Insertion: insertion1, States.Deletion: deletion1}
        layer1 = HMMLayer(1, ProteinResidue(1, ProteinAlphabet.ALA), states)
        hmm.layers.append(layer1)
        #  L2
        match2 = factory.create_match(emission, background)
        insertion2 = factory.create_insertion(background)
        deletion2 = factory.create_deletion()
        states = {States.Match: match2, States.Insertion: insertion2, States.Deletion: deletion2}
        layer2 = HMMLayer(2, ProteinResidue(2, ProteinAlphabet.CYS), states)
        hmm.layers.append(layer2)
        
        # Transitions
        #  start
        hmm.start.transitions.append(Transition(hmm.start, match1, 0.8))
        hmm.start.transitions.append(Transition(hmm.start, deletion1, 0.1))
        hmm.start.transitions.append(Transition(hmm.start, hmm.start_insertion, 0.1))
        hmm.start_insertion.transitions.append(Transition(hmm.start_insertion, hmm.start_insertion, 0.5))
        hmm.start_insertion.transitions.append(Transition(hmm.start_insertion, match1, 0.5))
        #  L1
        match1.transitions.append(Transition(match1, match2, 0.8))
        match1.transitions.append(Transition(match1, insertion1, 0.1))
        match1.transitions.append(Transition(match1, deletion2, 0.1))
        insertion1.transitions.append(Transition(insertion1, insertion1, 0.5))
        insertion1.transitions.append(Transition(insertion1, match2, 0.5))
        deletion1.transitions.append(Transition(deletion1, deletion2, 0.5))
        deletion1.transitions.append(Transition(deletion1, match2, 0.5))        
        #  L2
        match2.transitions.append(Transition(match2, hmm.end, 0.9))
        match2.transitions.append(Transition(match2, insertion2, 0.1))
        insertion2.transitions.append(Transition(insertion2, insertion2, 0.5)) 
        insertion2.transitions.append(Transition(insertion2, hmm.end, 0.5))
        deletion2.transitions.append(Transition(deletion2, hmm.end, 1.0))
        
        hmm.effective_matches = 10
        hmm.version = 1.5
        hmm.name = 'name'
        hmm.id = 'id'
        hmm.family = 'fam'
        hmm.length = ProfileLength(2, 2)

        return hmm 
        
        
@test.unit
class TestProfileHMM(test.Case):

    def setUp(self):
        
        super(TestProfileHMM, self).setUp()
        
        self.hmm = build_hmm()
        self.hmm2 = build_hmm()
        self.string = self.config.getContent('test.hhm')
        
    def testName(self):
        self.assertEqual(self.hmm.name, 'name')
    
    def testID(self):
        self.assertEqual(self.hmm.id, 'id')
        
    def testLayers(self):
        
        self.assertEqual(self.hmm.layers.length, 2)
        self.assertTrue(isinstance(self.hmm.layers[1], HMMLayer))
        self.assertEqual(set(self.hmm.layers[1]), set([States.Match, States.Insertion, States.Deletion]))
        
        for i in [0, 3]:
            self.assertRaises(IndexError, lambda:self.hmm.layers[i])
            
    def testStart(self):
        self.assertEqual(self.hmm.start.type, States.Start)

    def testStartInsertion(self):
        self.assertEqual(self.hmm.start_insertion.type, States.Insertion)
        
    def testEnd(self):
        self.assertEqual(self.hmm.end.type, States.End)               

    def testAllLayers(self):
        self.assertEqual(len(self.hmm.all_layers), 3)
        
    def testScoreUnits(self):
        self.assertEqual(self.hmm.score_units, ScoreUnits.Probability)
        
        self.hmm2.convert_scores(units=ScoreUnits.LogScales)
        self.assertEqual(self.hmm2.score_units, ScoreUnits.LogScales)
                
    def testConvertScores(self):
        
        self.hmm2.convert_scores(units=ScoreUnits.LogScales)
        self.assertEqual(self.hmm2.score_units, ScoreUnits.LogScales)
        self.assertEqual(self.hmm2.layers[1][States.Match].emission[ProteinAlphabet.ALA], 0)
        
        self.hmm2.convert_scores(units=ScoreUnits.Probability)
        self.assertEqual(self.hmm2.score_units, ScoreUnits.Probability)
        self.assertEqual(self.hmm2.layers[1][States.Match].emission[ProteinAlphabet.ALA], 1)
        
    def testSegment(self):
        
        segment = self.hmm.segment(1, 2)
        self.assertEqual(segment.layers.length, 2)
        self.assertEqual(segment.layers[2][States.Match].transitions[States.End].probability, 1)
        self.assertRaises(IndexError, self.hmm.segment, 1, 3)

    def testSubregion(self):
        
        segment = self.hmm.subregion(1, 2)
        self.assertEqual(segment.layers.length, 2)
        self.assertEqual(segment.layers[1], self.hmm.layers[1])
        self.assertRaises(IndexError, self.hmm.subregion, 1, 3)        
    
    def testStructure(self):
        s = self.hmm.structure('$', 'abcd')
        
        self.assertEqual(s.chains.length, 1)
        self.assertEqual(s.chains['$'].entry_id, 'abcd$')
        self.assertEqual(s.chains['$'].sequence, 'AC')
        
    def testChain(self):
        c = self.hmm.chain('$')
        
        self.assertEqual(c.id, '$')
        self.assertEqual(c.sequence, 'AC')
        
    def testEmissionSimilarity(self):
        
        score = self.hmm.emission_similarity(self.hmm)
        self.assertAlmostEqual(score, 9.6405, places=3)                 
    
    def testToHMM(self):
        
        self.hmm2.convert_scores(units=ScoreUnits.LogScales)
        self.assertEqual(self.string, self.hmm2.to_hmm().strip())
    

@test.unit
class TestProfileHMMSegment(test.Case):

    def setUp(self):
        
        super(TestProfileHMMSegment, self).setUp()
        
        self.hmm = build_hmm()
        self.segment = ProfileHMMSegment(self.hmm, 1, 2)
        self.string = self.config.getContent('test2.hhm')
    
    def testConstructor(self):
        
        self.segment.convert_scores(ScoreUnits.Probability)
        
        self.assertEqual(self.hmm.segment(1, 1).layers.length, 1)
        self.assertEqual(self.segment.layers.length, 2)
        self.assertEqual(self.segment.start.transitions[States.Match].probability, 1)
        self.assertEqual(self.segment.layers[1][States.Match].emission[ProteinAlphabet.CYS], 0.5)
        self.assertRaises(IndexError, self.hmm.segment, 1, 3)
                    
    def testToHMM(self):
        
        self.segment.convert_scores(ScoreUnits.LogScales) 
        self.assertEqual(self.string, self.segment.to_hmm().strip())


@test.unit
class TestHMMLayer(test.Case):

    def setUp(self):
        
        super(TestHMMLayer, self).setUp()        
        self.hmm = build_hmm()
        self.layer = self.hmm.layers[1] 
        
        self.layer.effective_matches = 5
        self.layer.effective_insertions = 4
        self.layer.effective_deletions = 3
        
    def testRank(self):
        self.assertEqual(self.hmm.layers[1].rank, 1)
        self.assertEqual(self.hmm.layers[2].rank, 2)
        
    def testNeff(self):
        self.assertEqual(self.layer.effective_matches, 5)
        self.assertEqual(self.layer.effective_insertions, 4)
        self.assertEqual(self.layer.effective_deletions, 3)
        
    def testResidue(self):
        def test(type):
            self.layer.residue = ProteinResidue(1, type=type)
            
        self.assertEqual(self.layer.residue.type, ProteinAlphabet.ALA)
        self.assertRaises(ValueError, test, ProteinAlphabet.GAP)
        test(ProteinAlphabet.ALA)
        self.assertEqual(self.layer.residue.type, ProteinAlphabet.ALA)
        
    def testAppend(self):
        
        s = State(States.End)
        self.layer.append(s)
                
        self.assertRaises(StateExistsError, self.layer.append, s)
        self.assertEqual(self.layer[States.End], s)

    def testIndexer(self):
        self.assertEqual(self.layer[States.Match].type, States.Match)
        

@test.unit
class TestState(test.Case):

    def setUp(self):
        
        super(TestState, self).setUp()
                
        self.hmm = build_hmm()
        self.layer = self.hmm.layers[1]
        self.m = self.hmm.layers[1][States.Match]
        self.d = self.hmm.layers[1][States.Deletion]
        self.i = self.hmm.layers[1][States.Insertion]
        self.s = self.hmm.start
        
    def testType(self):
        self.assertEqual(self.m.type, States.Match)
        
    def testTransitions(self):
        # also covers TransitionTable and Transition
        self.assertEqual(self.m.transitions.length, 3)
        self.assertEqual(self.m.transitions[States.Match].predecessor, self.m)
        self.assertEqual(self.m.transitions[States.Match].successor, self.hmm.layers[2][States.Match])
        self.assertEqual(self.m.transitions[States.Match].probability, 0.8)
        self.assertEqual(self.m.transitions[States.Insertion].type.source_state, States.Match)
        self.assertEqual(self.m.transitions[States.Insertion].type.target_state, States.Insertion)
        
    def testEmission(self):
        # also covers EmissionTable
        self.assertEqual(self.m.emission.length, 3)
        self.assertEqual(list(self.m.emission), [ProteinAlphabet.ALA, ProteinAlphabet.CYS, ProteinAlphabet.GLU])
        self.assertEqual(self.m.emission[ProteinAlphabet.CYS], 0.5)
        self.assertRaises(lambda: self.m.emission[ProteinAlphabet.GAP])
        self.assertRaises(UnobservableStateError, lambda: self.d.emission)
        
        self.m.emission.append(ProteinAlphabet.ILE, 0.5)
        self.assertEqual(self.m.emission[ProteinAlphabet.ILE], 0.5)
        self.m.emission.append('L', 0)
        self.assertRaises(EmissionExistsError, lambda:self.m.emission.append('I', 0))
        self.assertRaises(KeyError, lambda:self.m.emission.append('BUG', 0))
     
    def testSilent(self):   
        self.assertTrue(self.d.silent)
        self.assertFalse(self.m.silent)
        
    def testBackground(self):
        # also covers EmissionTable
        self.assertEqual(self.m.background.length, 3)
        self.assertEqual(list(self.m.background), [ProteinAlphabet.ALA, ProteinAlphabet.CYS, ProteinAlphabet.GLU])
        self.assertEqual(self.m.background[ProteinAlphabet.CYS], 0.00325358)
        self.assertRaises(lambda: self.m.background[ProteinAlphabet.GAP])
        
        self.m.background.append(ProteinAlphabet.ILE, 0.5)
        self.assertEqual(self.m.background[ProteinAlphabet.ILE], 0.5)
        self.m.background.append('L', 0)
        self.assertRaises(EmissionExistsError, lambda:self.m.background.append('I', 0))
        self.assertRaises(KeyError, lambda:self.m.background.append('BUG', 0))



if __name__ == '__main__':
    
    test.Console()
        
        
        