import csb.bio.nmr as nmr
import csb.bio.sequence as sequence
import csb.test as test


@test.unit
class TestRandomCoil(test.Case):
    
    def setUp(self):
                
        super(TestRandomCoil, self).setUp()
        
        self.rc = nmr.RandomCoil.get()
        self.chain = self.config.getPickle('1nz9.model1.pickle').first_chain
        self.residue = self.chain.residues[2]
        
    def testFactory(self):
        self.assertTrue(nmr.RandomCoil.get() is nmr.RandomCoil.get())
        
    def testSimpleSecondaryShift(self):
        
        raw = 200.0
        
        for r in ['A', 'ALA', sequence.SequenceAlphabets.Protein.ALA]:          #@UndefinedVariable
            
            self.assertEqual(
                    self.rc.simple_secondary_shift(r, 'N', raw),
                    raw - 125)            
        
        self.assertRaises(
                nmr.InvalidResidueError,
                lambda:self.rc.simple_secondary_shift('$', 'N', 0))
        self.assertRaises(
                nmr.EntityNotSupportedError,
                lambda:self.rc.simple_secondary_shift('A', '$', 0))        
        
    def testSecondaryShift(self):
        
        raw = 200.0 
        
        for r in [self.residue, self.residue.rank, self.residue.id]:
            
            self.assertEqual(
                    self.rc.secondary_shift(self.chain, r, 'H', raw),
                    raw - (8.44 + 0.07 - 0.05 - 0.01))
        
if __name__ == '__main__':

    test.Console()
