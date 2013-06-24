import csb.test as test

from csb.bio.nmr import RandomCoil, AtomConnectivity, ContactMap, Filters, Label
from csb.bio.nmr import ChemShiftInfo, ChemShiftScoringModel, NOESpectrum, NOEPeak
from csb.bio.nmr import InvalidResidueError, EntityNotSupportedError
from csb.bio.sequence import ProteinAlphabet
from csb.bio.structure import ChemElements, Atom


def get_chain():
    s = test.Config().getPickle('1nz9.model1.pickle')
    return s.first_chain.subregion(1, 5, clone=True)     


@test.unit
class TestRandomCoil(test.Case):
    
    def setUp(self):
        
        super(TestRandomCoil, self).setUp()
        
        self.rc = RandomCoil.get()
        self.chain = get_chain()
        self.residue = self.chain.residues[2]
        
    def testFactory(self):
        self.assertTrue(RandomCoil.get() is RandomCoil.get())
        
    def testSimpleSecondaryShift(self):
        
        raw = 200.0
        
        for r in ['A', 'ALA', ProteinAlphabet.ALA]:
            self.assertEqual(
                    self.rc.simple_secondary_shift(r, 'N', raw),
                    raw - 125)            
        
        self.assertRaises(
                InvalidResidueError,
                lambda:self.rc.simple_secondary_shift('$', 'N', 0))
        self.assertRaises(
                EntityNotSupportedError,
                lambda:self.rc.simple_secondary_shift('A', '$', 0))        
        
    def testSecondaryShift(self):
        
        raw = 200.0 
        
        for r in [self.residue, self.residue.rank, self.residue.id]:
            self.assertEqual(
                    self.rc.secondary_shift(self.chain, r, 'H', raw),
                    raw - (8.44 + 0.07 - 0.05 - 0.01))
            

@test.unit
class TestAtomConnectivity(test.Case):

    def setUp(self):
        super(TestAtomConnectivity, self).setUp()        
        self.ac = AtomConnectivity.get()
    
    def testFactory(self):
        self.assertTrue(AtomConnectivity.get() is AtomConnectivity.get())
        
    def testConnected(self):
        
        self.assertTrue(self.ac.connected(ProteinAlphabet.CYS, "SG", "HG"))
        self.assertTrue(self.ac.connected(ProteinAlphabet.CYS, "SG", "CB"))
        self.assertFalse(self.ac.connected(ProteinAlphabet.CYS, "SG", "H"))
        
        self.assertTrue(self.ac.connected(ProteinAlphabet.CYS, "CA", "C"))
        self.assertTrue(self.ac.connected(ProteinAlphabet.CYS, "CA", "HA"))
        self.assertTrue(self.ac.connected(ProteinAlphabet.CYS, "CA", "CB"))
        self.assertFalse(self.ac.connected(ProteinAlphabet.CYS, "CA", "HB"))
    
    def testConnectedAtoms(self):

        partners = self.ac.connected_atoms(ProteinAlphabet.CYS, "CA") 
        self.assertTrue(partners, ("N", "C", "CB", "HA"))
        
        partners = self.ac.connected_atoms(ProteinAlphabet.CYS, "SG")
        self.assertTrue(partners, ("CB"))
        
    def testContains(self):
        
        self.assertTrue(self.ac.contains(ProteinAlphabet.CYS, "SG"))
        self.assertTrue(self.ac.contains(ProteinAlphabet.ASP, "CG"))
        self.assertFalse(self.ac.contains(ProteinAlphabet.ALA, "CG"))
        
    def testGetAtoms(self):
        
        atoms = frozenset(['C', 'HD2', 'CB', 'CA', 'CG', 'O', 'N', 
                           'H', 'OD1', 'HA', 'OD2', 'HB3', 'HB2'])
        self.assertEqual(self.ac.get_atoms(ProteinAlphabet.ASP), atoms)


@test.unit
class TestContactMap(test.Case):

    CHAIN = get_chain()
    MAP = ContactMap(CHAIN, cutoff=1.75, filter=Filters.HYDROGENS)
    MAP.build()

    def setUp(self):
        super(TestContactMap, self).setUp()
        self.chain = TestContactMap.CHAIN
        self.cm = TestContactMap.MAP
                
    def testBuild(self): 
        self.assertEqual(len(self.cm.atoms), 38)
        
    def testContains(self):
        self.assertFalse(self.chain[0]['C'] in self.cm)
        self.assertTrue(self.chain[0]['HA'] in self.cm)
        
    def testIterator(self):
        self.assertEqual(list(self.cm), list(self.cm.contacts))
        
    def testCutoff(self):
        self.assertEqual(self.cm.cutoff, 1.75)
        
    def testChain(self):
        self.assertTrue(self.cm.chain is self.chain)
        
    def testAtoms(self):
        self.assertTrue(self.cm.atoms[0] is self.chain[0]['H1'])
        self.assertTrue(self.cm.atoms[-1] is self.chain[-1]['HZ'])
        
    def testContacts(self):
        
        c = self.chain        
        contacts = set([
                        (c[0]['H2'], c[0]['H3']),
                        (c[0]['H2'], c[0]['H1']),
                        (c[0]['H1'], c[0]['H3']),
                        (c[1]['HE22'], c[1]['HE21']) ])

        self.assertEqual(len(contacts), len(set(self.cm.contacts)))
        for a1, a2 in contacts:
            self.assertTrue((a1, a2) in contacts or (a2, a1) in contacts)
            
    def testConnect(self):
        # covers also cm.connected()
        
        c = self.chain
        cm = ContactMap(TestContactMap.CHAIN, cutoff=1.75, filter=Filters.CALPHAS)
        
        self.assertFalse(cm.connected(c[0]['CA'], c[1]['CA']))
        cm.connect(c[0]['CA'], c[1]['CA'])
        self.assertTrue(cm.connected(c[0]['CA'], c[1]['CA']))
        self.assertTrue(cm.connected(c[1]['CA'], c[0]['CA']))
        
    def testConnected(self):
        c = self.chain
        self.assertTrue(self.cm.connected(c[0]['H3'], c[0]['H2']))
        self.assertTrue(self.cm.connected(c[0]['H2'], c[0]['H3']))
            
    def testAtomContacts(self):
        atoms = set([self.chain[0]['H2'], self.chain[0]['H3']])
        self.assertEqual(self.cm.atom_contacts(self.chain[0]['H1']), atoms)
        
    def testResidueContacts(self):
        residues = set([self.chain[0]])
        self.assertEqual(self.cm.residue_contacts(self.chain[0]), residues)
        
    def testPosition(self):
        self.assertEqual(self.cm.position(1, 'H1'), 1.0)
        self.assertEqual(self.cm.position(2, 'HA'), 2.125)
        
    def testCompare(self):
        
        ci = ContactMap.compare(self.cm, self.cm, 0)
        self.assertEqual(ci.precision, 1)
        self.assertEqual(ci.coverage, 1)
        
        ci = ContactMap.compare(self.cm, self.cm, 1)
        self.assertEqual(ci.precision, 0)
        self.assertEqual(ci.coverage, 0)        


@test.unit
class TestLabel(test.Case):
    
    def testBuild(self):
        self.assertEqual(Label.build(ProteinAlphabet.ALA, 2, 'CA'), "A#2:CA")

    def testFromShift(self):
        shift = ChemShiftInfo(2, ProteinAlphabet.ALA, 'CA', ChemElements.C, 12)
        self.assertEqual(Label.from_shift(shift), "A#2:CA")
        
    def testFromAtom(self):
        atom = get_chain()[1]['CA']
        self.assertEqual(Label.from_atom(atom), "Q#2:CA")
        
    def testGetAtom(self):
        chain = get_chain()
        self.assertEqual(Label.get_atom(chain, "Q#2:CA"), chain[1]['CA'])
        
    def testMatch(self):        
        atom = get_chain()[1]['CA']
        shift = ChemShiftInfo(2, ProteinAlphabet.GLN, 'CA', ChemElements.C, 12)
        self.assertTrue(Label.match(shift, atom))
        
    def testParse(self):
        self.assertEqual(Label.parse("Q#2:CA"), ("Q", 2, 'CA'))
        
    def testFromString(self):
        
        label = Label.from_string("Q#2:CA")
        
        self.assertEqual(label.residue, ProteinAlphabet.GLN)
        self.assertEqual(label.rank, 2)
        self.assertEqual(label.atom_name, 'CA')


@test.unit
class TestChemShiftInfo(test.Case):
    
    def testConstructor(self):
        
        shift = ChemShiftInfo(2, 'ALA', 'CA', 'C', 12)
        
        self.assertEqual(shift.element, ChemElements.C)
        self.assertEqual(shift.residue, ProteinAlphabet.ALA)
        self.assertEqual(shift.position, 2)
        self.assertEqual(shift.shift, 12)
        
    def testLabel(self):

        shift = ChemShiftInfo(2, 'ALA', 'CA', 'C', 12)        
        self.assertEqual(shift.label, Label.from_shift(shift))


@test.unit
class TestChemShiftScoringModel(test.Case):
    
    def setUp(self):
        super(TestChemShiftScoringModel, self).setUp()
        self.model = ChemShiftScoringModel()
    
    def testPositive(self):
        self.assertAlmostEqual(self.model.positive('CA', 1), 0.191, places=3)
        self.assertAlmostEqual(self.model.positive('CA', [1, 1])[1], 0.191, places=3)

    def testNegative(self):
        self.assertAlmostEqual(self.model.negative('CA', 1), 0.127, places=3)
        self.assertAlmostEqual(self.model.negative('CA', [1, 1])[1], 0.127, places=3)    
    
    def testScore(self):
        self.assertAlmostEqual(self.model.score('CA', 1), 0.588, places=3)
        self.assertAlmostEqual(self.model.score('CA', [1, 1])[1], 0.588, places=3)
        

@test.unit
class TestNOESpectrum(test.Case):
    
    def setUp(self):
        super(TestNOESpectrum, self).setUp()
        
        self.elements = (ChemElements.H, ChemElements.C, ChemElements.H)
        self.spectrum = NOESpectrum(self.elements)
        self.spectrum.connect(1, 2)
        self.spectrum.add(11, [1, 2, 3])
        self.spectrum.add(12, [11, 22, 33])
        
    def testGetitem(self):
        self.assertEqual(self.spectrum[0].intensity, 11)
        self.assertRaises(IndexError, lambda i: self.spectrum[i], 3) 
                          
    def testLen(self):
        self.assertEqual(len(self.spectrum), 2)

    def testMinIntensity(self):        
        self.assertEqual(self.spectrum.min_intensity, 11)
        
    def testMaxIntensity(self):        
        self.assertEqual(self.spectrum.max_intensity, 12)

    def testElement(self):        
        self.assertEqual(self.spectrum.element(0), self.elements[0])
        self.assertEqual(self.spectrum.element(1), self.elements[1])
        self.assertRaises(IndexError, self.spectrum.element, 3) 

    def testDimensions(self):        
        self.assertEqual(self.spectrum.dimensions, self.elements)
        
    def testProtonDimensions(self):        
        self.assertEqual(self.spectrum.proton_dimensions, (0, 2))
        
    def testNumDimensions(self):        
        self.assertEqual(self.spectrum.num_dimensions, 3)
        
    def testNumProtonDimensions(self):        
        self.assertEqual(self.spectrum.num_proton_dimensions, 2)

    def testHasElement(self):
        self.assertFalse(self.spectrum.has_element(ChemElements.Ca))
        self.assertTrue(self.spectrum.has_element(ChemElements.C))

    def testHasConnectedDimensions(self):
        self.assertFalse(self.spectrum.has_connected_dimensions(0))
        self.assertTrue(self.spectrum.has_connected_dimensions(1))        
        self.assertTrue(self.spectrum.has_connected_dimensions(2))        

    def testConnectedDimensions(self):
        self.assertEqual(self.spectrum.connected_dimensions(0), ())
        self.assertEqual(self.spectrum.connected_dimensions(1), (2,))        
        self.assertEqual(self.spectrum.connected_dimensions(2), (1,))    
        
    def testConnect(self):
        self.assertRaises(ValueError, self.spectrum.connect, 0, 2) # H-H
        self.assertRaises(IndexError, self.spectrum.connect, 0, 3)            
    
    def testIterator(self):
        
        peaks = list(self.spectrum)
        
        self.assertEqual(peaks[0].intensity, 11)
        self.assertEqual(peaks[0].get(0), 1)
        self.assertEqual(peaks[0].get(1), 2)
        self.assertEqual(peaks[0].get(2), 3)        
        

@test.unit
class TestNOEPeak(test.Case):
    
    def setUp(self):
        super(TestNOEPeak, self).setUp()
        
        self.elements = (ChemElements.H, ChemElements.C, ChemElements.H)
        self.spectrum = NOESpectrum(self.elements)
        self.spectrum.connect(1, 2)
        self.spectrum.add(11, [1, 2, 3])
        self.spectrum.add(12, [11, 22, 33])
        
        self.peaks = list(self.spectrum)
        
    def testIntensity(self):
        self.assertEqual(self.peaks[0].intensity, 11)
        self.assertEqual(self.peaks[1].intensity, 12)
        
    def testNumDimensions(self):        
        self.assertEqual(self.peaks[0].num_dimensions, 3)
        
    def testHasConnectedDimensions(self):
        self.assertFalse(self.peaks[0].has_connected_dimensions(0))
        self.assertTrue(self.peaks[0].has_connected_dimensions(1))        
        self.assertTrue(self.peaks[0].has_connected_dimensions(2))        

    def testConnectedDimensions(self):
        self.assertEqual(self.peaks[0].connected_dimensions(0), ())
        self.assertEqual(self.peaks[0].connected_dimensions(1), (2,))        
        self.assertEqual(self.peaks[0].connected_dimensions(2), (1,))              

    def testElement(self):        
        self.assertEqual(self.peaks[0].element(0), self.elements[0])
        self.assertEqual(self.peaks[0].element(1), self.elements[1])          

    def testHasElement(self):        
        self.assertTrue(self.peaks[0].has_element(self.elements[0]))
        self.assertTrue(self.peaks[0].has_element(self.elements[1]))         
        
    def testGet(self):
        self.assertEqual(self.peaks[0].get(0), 1)
        self.assertEqual(self.peaks[0][0], 1)
        self.assertEqual(self.peaks[1].get(1), 22)         
        self.assertEqual(self.peaks[1][1], 22)
        
        self.assertRaises(IndexError, lambda i: self.peaks[0][i], 4)
        self.assertRaises(IndexError, self.peaks[0].get, -1)
            

if __name__ == '__main__':

    test.Console()
