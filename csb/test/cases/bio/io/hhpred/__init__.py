import csb.test as test

from csb.bio.hmm import ScoreUnits, States
from csb.bio.io.hhpred import HHProfileParser, HHOutputParser
from csb.bio.sequence import SequenceAlphabets
from csb.bio.structure import SecStructures


@test.regression
class TestHHOutputRegressions(test.Case):

    def setUp(self):
        
        super(TestHHOutputRegressions, self).setUp()
        
        filename = self.config.getTestFile('d1ea0a1.hhr')
        self.hitlist = HHOutputParser().parse_file(filename)
        
    def testLargeHitlist(self):
        """
        @see: [CSB 0000056]
        """
        hit = self.hitlist[-1]
        
        self.assertEqual(hit.rank, 10000)
        self.assertEqual(hit.id, 'd1v8qa_')

@test.regression
class TestHHProfileRegressions(test.Case):
        
    def testParseRedundantA3M(self):
        """
        @see: [CSB 0000068]
        """
        
        def strip(s):
            return s.replace('\t', '').replace(' ', '')
        
        filename = self.config.getTestFile('d1b24a2.hhm')
        content = self.config.getContent('d1b24a2.hhm')
        
        hmm = HHProfileParser(filename).parse()
        
        with self.config.getTempStream() as tmp:
            
            hmm.to_hmm(tmp.name)
            tmp.flush()
            generated = open(tmp.name).read()

            self.assertEqual(strip(content), strip(generated))

        
@test.unit
class TestHHProfileParser(test.Case):
    
    def setUp(self):
        
        super(TestHHProfileParser, self).setUp()
        
        self.hhm = self.config.getTestFile('d1nz0a_.hhm')
        self.pdb = self.config.getTestFile('d1nz0a_.pdb')        
        #self.profile = HHProfileParser(self.hhm, self.pdb).parse()

    def testParseLogScales(self):
        
        hmm = HHProfileParser(self.hhm, self.pdb).parse()
        self.assertEqual(hmm.score_units, ScoreUnits.LogScales) 
        
        self.assertEqual(hmm.name, 'd1nz0a_ d.14.1.2 (A:) RNase P protein {Thermotoga maritima [TaxId: 2336]}')
        self.assertEqual(hmm.id, 'd1nz0a_')
        self.assertEqual(hmm.family, 'd.14.1.2')
        self.assertEqual(hmm.length.matches, 109)
        self.assertEqual(hmm.length.layers, 109) 
               
        self.assertEqual(hmm.alignment.size, 9)
        self.assertEqual(hmm.consensus.sequence, 'eRLkxxxdFxxvxxxgxxxxxxxxxlxxxxxxxxxxRxGxxvsKKvgxAVxRNriKRxlRexxrxxxxxlxxxxdivvixrxxxxxxxxxxxxxxlxxxlxxlxkkixg')        
        self.assertEqual(hmm.alignment.rows[1][0].id, 'd1nz0a_')
        self.assertEqual(hmm.alignment.rows[1].strip().sequence, 'ERLRLRRDFLLIFKEGKSLQNEYFVVLFRKNGMDYSRLGIVVKRKFGKATRRNKLKRWVREIFRRNKGVIPKGFDIVVIPRKKLSEEFERVDFWTVREKLLNLLKRIEG')
                                        
        self.assertEqual(hmm.dssp[1].type, SecStructures.Coil)
        self.assertEqual(hmm.psipred[1].type, SecStructures.Coil)        
        self.assertEqual(hmm.dssp.to_string(), 'CCCCHHHHHHHHHHHSEEEECSSEEEEEEECSSSSCEEEECCCGGGCSHHHHHHHHHHHHHHHHHHTTTSCSSEEEEEEECHHHHHHGGGSCHHHHHHHHHHHHTTCCC')
        self.assertEqual(hmm.psipred.to_string(), 'CCCCCHHHHHHHHHCCCEEECCCEEEEEECCCCCCCEEEEEEECCCHHHHHHHHHHHHHHHHHHHHHHHCCCCCEEEEEECCCCCCCCCCCCHHHHHHHHHHHHHHHCC')        
        
        self.assertEqual(hmm.effective_matches, 7.4)
        self.assertEqual(hmm.pseudocounts, False)
        
        self.assertEqual(hmm.start.silent, True)
        self.assertEqual(hmm.start.type, States.Start)
        self.assertEqual(set(hmm.start.transitions), set([States.Match]))
        
        self.assertEqual(hmm.end.silent, True)
        self.assertEqual(hmm.end.type, States.End)
        self.assertEqual(set(hmm.end.transitions), set([]))      

        self.assertEqual(hmm.has_structure, True)
        self.assertEqual(hmm.chain().sequence, 'ERLRLRRDFLLIFKEGKSLQNEYFVVLFRKNGMDYSRLGIVVKRKFGKATRRNKLKRWVREIFRRNKGVIPKGFDIVVIPRKKLSEEFERVDFWTVREKLLNLLKRIEG')        
                
        layer = hmm.layers[2]
        self.assertEqual(layer.rank, 2)
        self.assertAlmostEqual(layer.effective_matches, 6.9939, places=3)
        self.assertAlmostEqual(layer.effective_insertions, 0)
        self.assertAlmostEqual(layer.effective_deletions, 0)
        self.assertEqual(set(layer), set([States.Match, States.Insertion]))
        
        self.assertEqual(layer.residue.type, SequenceAlphabets.Protein.ARG)             
        self.assertEqual(layer.residue.has_structure, True)
        self.assertEqual(layer.residue.secondary_structure.type, SecStructures.Coil)                            
        self.assertEqual(layer.residue.secondary_structure.end, 4)
                
        match = layer[States.Match]
        self.assertEqual(match.type, States.Match)        
        self.assertEqual(match.rank, 2)
        self.assertEqual(match.background[SequenceAlphabets.Protein.ALA], 3706.0)                                               
        self.assertEqual(match.emission[SequenceAlphabets.Protein.LYS], 3316.0)         
        self.assertEqual(set(match.transitions), set([States.Match]))
        self.assertEqual(match.transitions[States.Match].predecessor, match)
        self.assertEqual(match.transitions[States.Match].successor.rank, 3)                
        self.assertEqual(match.transitions[States.Match].successor.type, States.Match)   
        self.assertEqual(repr(match.transitions[States.Match].type), 'M->M')
        self.assertEqual(match.transitions[States.Match].probability, 0.0)        
                                                
    def testParseProbability(self):
        
        hmm = HHProfileParser(self.hhm, self.pdb).parse(ScoreUnits.Probability)
        self.assertEqual(hmm.score_units, ScoreUnits.Probability)

        layer = hmm.layers[2]        
        self.assertAlmostEqual(layer.effective_matches, 6.9939, places=3)
        self.assertAlmostEqual(layer.effective_insertions, 0.0)
        self.assertAlmostEqual(layer.effective_deletions, 0.0)
                
        match = layer[States.Match]
        self.assertAlmostEqual(match.background[SequenceAlphabets.Protein.ALA], 0.07662, places=4)                                               
        self.assertAlmostEqual(match.emission[SequenceAlphabets.Protein.LYS], 0.10041, places=4)                 
        self.assertAlmostEqual(match.transitions[States.Match].probability, 1.0)           
        
@test.unit
class TestHHOutputParser(test.Case):

    def setUp(self):
        
        super(TestHHOutputParser, self).setUp()
        
        filename = self.config.getTestFile('d1ea0a1.hhr')
        content = open(filename).read()
        tmp = HHOutputParser(True)
        
        self.hitlist = tmp.parse_file(filename)
        self.hitlist2 = tmp.parse_string(content)

    def testParseFile(self):

        # header  
        self.assertEqual(self.hitlist.command, './hhsearch -i /ebio/abt1/kopec/projects/hhpred_on_cluster/split_fastas/d1ea0a1.hhm -d /ebio/abt1/kopec/projects/hhpred_on_cluster/astral_scop70_all_beta_only_v1.75_hhpred_database.hhm')
        self.assertEqual(self.hitlist.query_name, 'd1ea0a1 b.80.4.1 (A:1203-1472) Alpha subunit of glutamate synthase, C-terminal domain {Azospirillum brasilense [TaxId: 192]}')
        self.assertEqual(self.hitlist.neff, float('5.2'))
        self.assertEqual(self.hitlist.searched_hmms, 2936)   
        
        # alignments
        segments = list(self.hitlist[9].alignment.segments)
        self.assertEqual(len(segments), 3)
        
        for s, l in zip(segments, [(9, 21), (22, 35), (37, 52)]):
            self.assertEqual(s.start, l[0])
            self.assertEqual(s.end, l[1])
            self.assertEqual(s.end - s.start, l[1] - l[0])
            
        # hits
        self.assertEqual(len(self.hitlist), 11)
        self.assertRaises(IndexError, self.hitlist.__getitem__, 20)                   
        self.assertEqual(self.hitlist[0].rank, 1)

        hit = self.hitlist[4]
        self.assertEqual(hit.rank, 5)                
        self.assertEqual(hit.id, 'd1g3pa1')
        self.assertEqual(hit.start, 34)
        self.assertEqual(hit.end, 49)
        self.assertEqual(hit.qstart, 161)
        self.assertEqual(hit.qend, 176)
        self.assertEqual(hit.qlength, 270)
        self.assertEqual(hit.probability * 100, 35.1)
        self.assertEqual(hit.length, 16)       
        self.assertEqual(hit.slength, 65)
        self.assertEqual(hit.evalue, 1.2)
        self.assertEqual(hit.pvalue, 0.00041)
        self.assertEqual(hit.score, 27.0)
        self.assertEqual(hit.identity, 31.0)
        self.assertEqual(hit.similarity, 0.659)
        self.assertEqual(hit.prob_sum, 13.1)                                                                                           
            
    def testParseString(self):

        # header  
        self.assertEqual(self.hitlist.command, './hhsearch -i /ebio/abt1/kopec/projects/hhpred_on_cluster/split_fastas/d1ea0a1.hhm -d /ebio/abt1/kopec/projects/hhpred_on_cluster/astral_scop70_all_beta_only_v1.75_hhpred_database.hhm')
        self.assertEqual(self.hitlist.query_name, 'd1ea0a1 b.80.4.1 (A:1203-1472) Alpha subunit of glutamate synthase, C-terminal domain {Azospirillum brasilense [TaxId: 192]}')
        self.assertEqual(self.hitlist.neff, float('5.2'))
        self.assertEqual(self.hitlist.searched_hmms, 2936)   
        
        # alignments
        segments = list(self.hitlist[9].alignment.segments)
        self.assertEqual(len(segments), 3)
        
        for s, l in zip(segments, [(9, 21), (22, 35), (37, 52)]):
            self.assertEqual(s.start, l[0])
            self.assertEqual(s.end, l[1])
            self.assertEqual(s.end - s.start, l[1] - l[0])
            self.assertEqual(self.hitlist[9].score, s.score)
            self.assertEqual(self.hitlist[9].ss_score, s.ss_score)
            self.assertEqual(self.hitlist[9].probability, s.probability)
            
        # hits
        self.assertEqual(len(self.hitlist), 11)
        self.assertRaises(IndexError, self.hitlist.__getitem__, 20)                   
        self.assertEqual(self.hitlist[0].rank, 1)

        hit = self.hitlist[4]
        self.assertEqual(hit.rank, 5)                
        self.assertEqual(hit.id, 'd1g3pa1')
        self.assertEqual(hit.start, 34)
        self.assertEqual(hit.end, 49)
        self.assertEqual(hit.qstart, 161)
        self.assertEqual(hit.qend, 176)
        self.assertEqual(hit.qlength, 270)
        self.assertEqual(hit.probability * 100, 35.1)
        self.assertEqual(hit.length, 16)       
        self.assertEqual(hit.slength, 65)
        self.assertEqual(hit.evalue, 1.2)
        self.assertEqual(hit.pvalue, 0.00041)
        self.assertEqual(hit.score, 27.0)
        self.assertEqual(hit.identity, 31.0)
        self.assertEqual(hit.similarity, 0.659)
        self.assertEqual(hit.prob_sum, 13.1)  


if __name__ == '__main__':
    
    test.Console()
