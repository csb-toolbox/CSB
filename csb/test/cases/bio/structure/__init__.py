import numpy

import csb.test as test
import csb.bio.structure as structure
import csb.core

import csb.bio.sequence as sequence
    
    
@test.regression
class ChainRegressions(test.Case):

    def setUp(self):
        
        super(ChainRegressions, self).setUp()

        self.structure = self.config.getPickle('1nz9.model1.pickle')        
        self.chain = self.structure['A']
            
    def testClone(self):
        """
        @see: [CSB 0000029]
        """
        self.assertEqual(self.chain.residues[1].clone()._container, None)
        
    def testApplyTransformation(self):
        """
        @see: [CSB 0000030]
        """
        r = numpy.eye(3)
        t = numpy.array([1, 1, 1])
        
        residue = self.chain.residues[1]
        
        atom1 = structure.Atom(99999, 'Cx', structure.ChemElements.C, numpy.array([2, 2, 2]), alternate=True)    
        atom1.occupancy = 0.2
        atom2 = structure.Atom(99999, 'Cx', structure.ChemElements.C, numpy.array([2, 2, 2]), alternate=True)    
        atom2.occupancy = 0.8
                        
        residue.atoms.append(atom1)
        residue.atoms.append(atom2)
        
        self.chain.transform(r, t)
        
        disatom = residue.atoms['Cx']
        self.assertTrue(isinstance(disatom, structure.DisorderedAtom))
        self.assertEqual(disatom.vector.tolist(), [3, 3, 3])
        
        # now the regression itself
        for altatom in disatom:
            self.assertEqual(altatom.vector.tolist(), [3, 3, 3])
    
    @test.skip("may fail on slower machines")
    def testListCoordinates(self):
        """
        Performance tweaks in AbstractEntity.get_coordinates()
        """
        
        def test():
            for dummy in range(1000):
                self.structure.get_coordinates(what=['CA'])
                
        self.assertFasterThan(0.50, test)
            
    def testAppendAtom(self):
        """
        @see: [CSB 0000122]
        """
        r = structure.ProteinResidue(1, sequence.ProteinAlphabet.ALA)
        
        a1 = structure.Atom(1, 'CA', structure.ChemElements.C, [1, 1, 1], alternate='A')
        a2 = structure.Atom(1, 'CA', structure.ChemElements.C, [1, 1, 1], alternate='B')
        
        r.atoms.append(a1)
        r.atoms.append(a2)
        
        self.assertEqual(a1.residue, r)
        self.assertEqual(a2.residue, r)
        
        
@test.unit
class TestAbstractEntity(test.Case): 
    
    def setUp(self):
        
        super(TestAbstractEntity, self).setUp()
        
        self.e = self.config.getPickle('1nz9.full.pickle')
        self.s = self.e.models[1]
        self.c = self.s.chains['A']
        self.r = self.c.residues[1]
        self.a = self.r.atoms['CA']
        
    def testItems(self):
        
        self.assertEqual(list(self.e.items), list(self.e.models))
        self.assertEqual(list(self.s.items), list(self.s.chains[c] for c in self.s.chains))
        self.assertEqual(list(self.c.items), list(self.c.residues))
        self.assertEqual(list(self.r.items), list(self.r.atoms[a] for a in self.r.atoms))
        self.assertEqual(list(self.a.items), list([]))    

    def testComponents(self):
        
        # also covers the composite iterators

        items = list(structure.CompositeEntityIterator.create(self.e, leaf=structure.Chain))
        self.assertEqual(items, [self.e[0], self.e[0]['A'], self.e[1], self.e[1]['A']])
        
        self.assertEqual(list(self.e.components())[:4], [self.s, self.c, self.r, list(self.r.items)[0]])        
        self.assertEqual(list(self.e.components(structure.Structure)), list(self.e.models))
        self.assertEqual(list(self.e.components(structure.Chain)), [self.e[0]['A'], self.e[1]['A']])

        self.assertEqual(list(self.r.components()), [self.r[a] for a in self.r.atoms])
        self.assertEqual(list(self.r.components(structure.Atom)), [self.r[a] for a in self.r.atoms])        
        self.assertEqual(list(self.r.components(structure.Residue)), [])                

        self.assertEqual(list(self.a.components()), [])
        self.assertEqual(list(self.a.components(structure.Atom)), [])
      
    def testApplyTransformation(self):
        
        R = numpy.eye(3)
        t = numpy.array([1, 2, 3])
        
        original = self.s['A'].get_coordinates(('CA',))
        assert len(original) > 0

        self.e.transform(R, t)
        translated = self.s['A'].get_coordinates(('CA',))
        
        self.e.transform(R, -t)
        restored = self.s['A'].get_coordinates(('CA',))                
        
        self.assertTrue(len(original) == len(translated) == len(restored))
        
        for i in range(len(original)):
            for j in range(3):
                self.assertEqual((original[i] + t)[j], translated[i][j])                
                self.assertAlmostEqual(original[i][j], restored[i][j], places=13)  
                
    def testListCoordinates(self):
        
        # all atoms
        original, listed = [], []
        
        for s in self.e.items:
            for c in s.items:
                for r in c.items:
                    for a in r.items:
                        original.append(list(a.vector))
        
        for i in self.e.get_coordinates(what=None):
            listed.append(list(i))
            
        self.assertEqual(original, listed)
        self.assertEqual([list(i) for i in self.r.get_coordinates()], [list(a.vector) for a in self.r.items])
        self.assertEqual(list(self.a.get_coordinates()[0]), list(self.a.vector))
        
        # specific atoms
        original, listed = [], []
        
        for s in self.e.items:
            for c in s.items:
                for r in c.items:
                    original.append(list(r.atoms['CA'].vector))
        
        for i in self.e.get_coordinates(what=['CA']):
            listed.append(list(i))
            
        self.assertEqual(original, listed)
        self.assertEqual([list(i) for i in self.r.get_coordinates(['CA'])], [list(self.r.atoms['CA'].vector)])
        self.assertEqual(list(self.a.get_coordinates(['CA'])[0]), list(self.a.vector))
                        
        # other stuff
        self.assertRaises(structure.Broken3DStructureError, lambda: self.c.get_coordinates(what=['BUG']))
        self.assertRaises(structure.Broken3DStructureError, lambda: self.r.get_coordinates(what=['BUG']))        
        self.assertRaises(structure.Broken3DStructureError, lambda: self.a.get_coordinates(what=['BUG']))
                    
@test.unit
class TestEnsemble(test.Case):
    
    def setUp(self):
        
        super(TestEnsemble, self).setUp()
        self.ensemble = self.config.getPickle('1nz9.full.pickle')

    def testItems(self):
        
        models = list(self.ensemble.items)
        self.assertEqual(len(models), self.ensemble.models.length)
        
    def testComponents(self):
        
        models = list(self.ensemble.components(klass=structure.Structure))
        self.assertEqual(len(models), self.ensemble.models.length)
                    
    def testGetitem(self):
        
        self.assertEqual(self.ensemble.models.length, 2)
        self.assertEqual(self.ensemble.models[1].model_id, 1)
        self.assertEqual(self.ensemble[1].model_id, 2)

    def testFirstModel(self):
        
        self.assertEqual(self.ensemble.models[1], self.ensemble.first_model)
        self.assertEqual(self.ensemble[0], self.ensemble.first_model)  
        
    def testAppend(self):
        
        testModel = structure.Structure('test')
        testModel.model_id = 3
        rank = self.ensemble.models.append(testModel)
        
        self.assertEqual(rank, 3)
        self.assertEqual(self.ensemble.models[3], testModel)        
        self.assertEqual(self.ensemble.models[3].model_id, 3)        
        
        self.assertRaises(Exception, self.ensemble.models.append, 1)
        
    def testToPDB(self):

        pdbraw = self.ensemble.to_pdb()
        pdb = pdbraw.splitlines()
        
        model = 0
        endmdl = 0
        
        for line in pdb:
            if line.startswith('MODEL '):
                model += 1
            if line.startswith('ENDMDL '):
                endmdl += 1
                
        self.assertEqual(model, self.ensemble.models.length)
        self.assertEqual(endmdl, self.ensemble.models.length)
        
        self.assertTrue(pdb[0].startswith, 'HEADER')
        self.assertTrue(pdb[0].endswith('1NZ9              '))
        self.assertEqual(pdb[-3].rstrip(), 'TER')
        self.assertEqual(pdb[-2].rstrip(), 'ENDMDL')
        self.assertEqual(pdb[-1].rstrip(), 'END')                

@test.unit
class TestStructure(test.Case):
    
    def setUp(self):
        
        super(TestStructure, self).setUp()
        self.structure = self.config.getPickle('1nz9.model1.pickle')
        self.structure2 = self.config.getPickle('1nz9.model1.pickle')        
        
    def testGetitem(self):
        
        self.assertEqual(self.structure.chains.length, 1)
        self.assertTrue('A' in self.structure.chains)
        self.assertTrue('A' in self.structure)
        self.assertEqual(self.structure.chains['A'], self.structure['A'])        
        self.assertEqual(self.structure.chains['A'].id, 'A')        

    def testFirstChain(self):
        
        self.assertEqual(self.structure.chains['A'], self.structure.first_chain)
        self.assertEqual(self.structure['A'], self.structure.first_chain)
        
    def testItems(self):
        
        self.assertEqual(self.structure['A'], tuple(self.structure.items)[0])          
        
    def testAppend(self):
        
        testChain = structure.Chain('$', accession='____')
        self.structure2.chains.append(testChain)
        
        self.assertTrue('$' in self.structure2.chains)
        self.assertEqual(self.structure2.chains['$'], testChain)
        self.assertEqual(self.structure2.chains['$'].accession, self.structure2.accession)        
        
        self.assertRaises(structure.DuplicateChainIDError, self.structure2.chains.append, testChain)
        self.assertRaises(structure.InvalidOperation, self.structure2.chains.append, self.structure.first_chain)                
    
    def testRemove(self):
        
        self.assertTrue('A' in self.structure2.chains)
        self.structure2.chains.remove('A')
        self.assertFalse('A' in self.structure2.chains)
        
        for c in ['.', 'A']:
            self.assertRaises(csb.core.ItemNotFoundError, self.structure2.chains.remove, c)       
            
    def testAccession(self):
        
        self.assertEqual(self.structure.accession, '1nz9')
        
        self.structure.accession = '  TeST  '
        self.assertEqual(self.structure.accession, 'test')
        
        for chain in self.structure.items:
            self.assertEqual(chain.accession, self.structure.accession)

        self.structure.accession = '1nz9'            
                    
    def testToFASTA(self):
        
        fasta = self.structure.to_fasta().splitlines()
        self.assertTrue(fasta[0].startswith('>1nz9_A'))        
        self.assertEqual(fasta[1], 'AQVAFREGDQVRVVSGPFADFTGTVTEINPERGKVKVMVTIFGRETPVELDFSQVVKA')

    def testToPDB(self):
        
        pdbraw = self.structure.to_pdb()
        pdb = pdbraw.splitlines()
        
        self.assertTrue(pdb[0].startswith, 'HEADER')
        self.assertTrue(pdb[0].endswith('1NZ9              '))
        self.assertEqual(pdb[-2].rstrip(), 'TER')  
        self.assertEqual(pdb[-1].rstrip(), 'END')
        
        for line in pdb:
            self.assertTrue(not line.startswith('MODEL '))
            self.assertTrue(not line.startswith('ENDMDL '))
                
        has_seqres = False
        has_atom = False 
        # check first SEQRES line
        for l in pdb:
            if l.startswith('SEQRES') and not has_seqres:
                has_seqres = True                
                self.assertEqual(l, 'SEQRES   1 A   58  ALA GLN VAL ALA PHE ARG GLU GLY ASP GLN VAL ARG VAL          ')
        
            if l.startswith('SEQRES   3 A   58'):          
                self.assertEqual(l, 'SEQRES   3 A   58  GLU ILE ASN PRO GLU ARG GLY LYS VAL LYS VAL MSE VAL          ')
                break                   
        
        # check first ATOM line
        for l in pdb:
            if l.startswith('ATOM') and not has_atom:
                has_atom = True                
                self.assertEqual(l, 'ATOM      1  N   ALA A 127      -0.977 -18.702  -7.689  1.00  0.00           N  ')
            
            if l.startswith('ATOM    572'):
                self.assertEqual(l, 'ATOM    572 SE   MSE A 164       8.107   5.802   4.753  1.00  0.00          Se  ')
                break
                
        self.assertTrue(has_seqres)
        self.assertTrue(has_atom)
        
        # test also the content written to a file stream 
        with self.config.getTempStream() as tmp:
            
            self.structure.to_pdb(tmp)
            tmp.flush()
                        
            self.assertEqual(pdbraw, open(tmp.name).read())
            
        # and to a file 
        with self.config.getTempStream() as tmp:   
            self.structure.to_pdb(tmp.name)
            self.assertEqual(pdbraw, open(tmp.name).read())            
            
@test.unit
class TestChain(test.Case):

    def setUp(self):
        
        super(TestChain, self).setUp()

        self.structure = self.config.getPickle('1nz9.model1.pickle')        
        self.chain = self.structure['A']
        
    def testFromSequence(self):
        
        seq = sequence.Sequence('testA', 'test seq', 'ABCDE', type=sequence.SequenceTypes.Protein)
        chain = structure.Chain.from_sequence(seq, 'Z')
        
        self.assertEqual(chain.sequence, 'ABCDE')
        self.assertEqual(chain.id, 'Z')
        self.assertEqual(chain.residues[2].sequence_number, 2)
        self.assertEqual(chain.residues[2].rank, 2)
        self.assertEqual(chain.residues[2].type, sequence.ProteinAlphabet.ASX)
    
    def testId(self):
        
        self.assertEqual(self.chain.id, 'A')
        self.chain.id = '$'
        self.assertEqual(self.chain.id, '$')    
        self.assertEqual(tuple(self.structure.items)[0].id, '$')
        self.assertTrue('$' in self.structure)        
        self.chain.id = 'A'    
    
    def testAccession(self):
        
        def assign(chain, value):
            chain.accession = value
            
        self.assertRaises(structure.InvalidOperation, assign, self.chain, 'xxxx')
        chain = self.chain.clone()
        chain.accession = 'xxxx'
        self.assertEqual(chain.accession, 'xxxx')

    def testType(self):
        self.assertEqual(self.chain.type, structure.SequenceTypes.Protein)              

    def testEntryID(self):
        self.assertEqual(self.chain.entry_id, '1nz9A')
        
    def testSequence(self):
        self.assertEqual(self.chain.sequence, 'AQVAFREGDQVRVVSGPFADFTGTVTEINPERGKVKVMVTIFGRETPVELDFSQVVKA')
        
    def testHeader(self):
        self.assertTrue(self.chain.header.startswith('1nz9_A'))
        
    def testLength(self):
        self.assertEqual(self.chain.length, 58)
        self.assertEqual(self.chain.length, len(self.chain))
        
    def testAlphabet(self):
        
        self.assertEqual(self.chain.alphabet, structure.SequenceAlphabets.Protein)

        chain = self.chain.clone()
        chain.type = structure.SequenceTypes.DNA                                                    
        self.assertEqual(chain.alphabet, structure.SequenceAlphabets.Nucleic)
        
    def testSecStructure(self):
        
        self.assertEqual(self.chain.secondary_structure.length, 6)
        self.assertEqual(self.chain.secondary_structure[1].type, structure.SecStructures.Strand)    
        self.assertEqual(self.chain.secondary_structure[5].type, structure.SecStructures.Helix)      
        self.assertEqual(self.chain.secondary_structure.to_string(), '---------EEEE--------EEEEEEEE----EEEEEEE----EEEEEE-HHHEEE')        
    
    def testClone(self):
        
        chain = self.chain
        clone = chain.clone()
        
        self.assertEqual(chain.id, clone.id)
        self.assertEqual(chain.length, clone.length)
        self.assertNotEqual(id(chain), id(clone))               

        # residues
        self.assertEqual(chain[0].type, clone[0].type)
        self.assertNotEqual(id(chain[0]), id(clone[0]))
        # atoms
        self.assertEqual(set(chain[0].atoms), set(clone[0].atoms))        
        self.assertNotEqual(id(chain[0]['CA']), id(clone[0]['CA']))        
        
        # detached?
        self.assertEqual(chain._structure, self.structure)
        self.assertEqual(clone._structure, None) 
        
    def testSubregion(self):
        
        sub = self.chain.subregion(1, self.chain.length)
        self.assertEqual(list(self.chain.residues), list(sub.residues))
        
        self.assertRaises(IndexError, self.chain.subregion, -1, 5)
        self.assertRaises(IndexError, self.chain.subregion, 1, self.chain.length + 1)
        
    def testFind(self):
        
        self.assertEqual(self.chain.find('127'), self.chain[0])
        
        self.chain[0].id = 127, 'X'
        self.assertEqual(self.chain.find('127', 'X'), self.chain[0])
        
        self.chain[0].id = 127, None
        self.assertRaises(csb.core.ItemNotFoundError, self.chain.find, 127, 'X')
        self.assertRaises(csb.core.ItemNotFoundError, self.chain.find, 999999)
        
    def testComputeTorsion(self):
        
        self.assertRaises(structure.InvalidOperation, lambda:self.chain.torsion)
        self.chain.compute_torsion()

        self.assertEqual(self.chain.residues[1].torsion.phi, None)     
        self.assertNotEqual(self.chain.residues[1].torsion.psi, None)        
        
        # http://www.fz-juelich.de/nic/cbb/wang/dssp/dssp/1nz9.dssp
        self.assertAlmostEqual(self.chain.residues[2].torsion.phi, -134.597, places=3)      
        self.assertAlmostEqual(self.chain.residues[2].torsion.psi, 155.738, places=3)        
        self.assertAlmostEqual(self.chain.residues[2].torsion.omega, 179.916, places=3)        
    
    def testListCoordinates(self):
        
        coords = self.chain.get_coordinates(('CA',))
        self.assertEqual(coords[0].tolist(), self.chain[0]['CA'].vector.tolist())
        self.assertEqual(len(coords), self.chain.length)
                
        self.assertRaises(structure.Broken3DStructureError, self.chain.get_coordinates, ('H',))
        
    def testSuperimpose(self):
        
        clone = self.chain.clone()
        si = self.chain.superimpose(clone)
        
        self.assertAlmostEqual(si.rmsd, 0.0, places = 6)
        for coord in si.translation:
            self.assertAlmostEqual(coord, 0.0)
            
        clone.transform(si.rotation, si.translation)
        self.assertAlmostEqual(self.chain.rmsd(clone), 0, 5)
        
        # not testable for now, but ensure at least no exception is raised
        si = self.chain.superimpose(clone, what=('C', 'CA',))
        si = self.chain.superimpose(clone, how=structure.AlignmentTypes.Local)          
        
        # make sure the subject is not moving (as it is the case in .align)
        #  - grab a different structure   
        diff = self.chain.clone()
        diff[1]['CA'].vector += 2.0
        x = diff[1]['CA'].vector[0]
        #  - superimpose them
        si = self.chain.superimpose(diff)
        #  - assert diff is untouched
        assert si.rmsd > 0, 'test with different structures' 
        self.assertEqual(diff[1]['CA'].vector[0], x)
        
        self.assertRaises(ValueError, clone.superimpose, clone.subregion(1, 1))        
        
    def testAlign(self):
        
        clone = self.chain.clone()
        si = self.chain.align(clone)
        
        self.assertAlmostEqual(si.rmsd, 0.0, places = 6)
        for coord in si.translation:
            self.assertAlmostEqual(coord, 0.0)
            
        self.assertAlmostEqual(self.chain.rmsd(clone), 0, 5)
        
        # mutate clone a bit to achieve RMSD != 0
        clone[0]['CA'].vector[0] += 2.0
        si = self.chain.align(clone)
        self.assertNotEqual(si.rmsd, 0.0)   
        
        # ensure clone has been transformed; self.chain holds the original coords of clone
        self.assertNotEqual(clone[3]['CA'].vector[0], self.chain[3]['CA'].vector[0])
        
        self.assertRaises(ValueError, clone.align, clone.subregion(1, 1))        
                
        # not testable for now, but ensure at least no exception is raised
        si = self.chain.align(clone, what=('C', 'CA',))
        si = self.chain.align(clone, how=structure.AlignmentTypes.Local)            
        
    def testRMSD(self):
        
        chain = self.chain
        self.assertAlmostEqual(chain.rmsd(chain), 0.0, places=6)  
        self.assertAlmostEqual(chain.rmsd(chain, what=['N', 'CA', 'C', 'O']), 0.0, places=6)
        
        self.assertRaises(ValueError, chain.rmsd, chain.subregion(1, 1))
        
    def testTMSuperimpose(self):
        
        clone = self.chain.clone()
        si = self.chain.tm_superimpose(clone)
        
        self.assertEqual(si.tm_score, 1.0)
        for coord in si.translation:
            self.assertAlmostEqual(coord, 0.0)
            
        clone.transform(si.rotation, si.translation)
        self.assertEqual(self.chain.tm_score(clone), 1.0)
        
        # not testable for now, but ensure at least no exception is raised
        si = self.chain.tm_superimpose(clone, what=('C', 'CA',))
        si = self.chain.tm_superimpose(clone, how=structure.AlignmentTypes.Local)       
        
        # make sure the subject is not moving (as it is the case in .align)
        #  - grab a different structure   
        diff = self.chain.clone()
        diff[1]['CA'].vector += 2.0
        x = diff[1]['CA'].vector[0]
        #  - superimpose them
        si = self.chain.tm_superimpose(diff)
        #  - assert diff is untouched
        assert si.tm_score < 1, 'test with different structures' 
        self.assertEqual(diff[1]['CA'].vector[0], x)
        
        self.assertRaises(ValueError, clone.tm_superimpose, clone.subregion(1, 1)) 
        
    def testTMScore(self):
        
        chain = self.chain
        self.assertAlmostEqual(chain.tm_score(chain), 1.0, places=6)  
        self.assertAlmostEqual(chain.tm_score(chain, what=['N', 'CA', 'C', 'O']), 1.0, places=6)
        
        self.assertRaises(ValueError, chain.tm_score, chain.subregion(1, 1))        
              
    def testGetitem(self):
        
        self.assertEqual(self.chain.residues[1], self.chain[0])
        
        self.assertRaises(csb.core.CollectionIndexError, lambda:self.chain.residues[0])
        self.assertRaises(csb.core.CollectionIndexError, lambda:self.chain.residues[59])
        self.assertRaises(IndexError, lambda:self.chain[58])

        self.assertEqual(self.chain[0].type, structure.SequenceAlphabets.Protein.ALA)       
        self.assertEqual(self.chain[0].rank, 1)                                   
        self.assertEqual(self.chain[0].sequence_number, 127)        
        self.assertEqual(self.chain[0].id, '127')
        
    def testAppend(self):
        
        chain = self.chain.clone()
        chain.compute_torsion()
        assert chain.torsion is not None
        residue = structure.ProteinResidue(9999, structure.SequenceAlphabets.Protein.ALA, 9999, 'A')        
        
        rank = chain.residues.append(residue)        
        self.assertTrue(chain.residues._contains('9999A'))
        self.assertRaises(structure.InvalidOperation, lambda:chain.torsion)        
        
        self.assertEqual(chain[-1], residue)
        self.assertEqual(chain.residues[59], residue)        
        self.assertEqual(rank, 59)        
        
        self.assertRaises(structure.DuplicateResidueIDError, chain.residues.append, residue)     
        # assert NOT raises DuplicateResidueIDError if id is None:
        residue = structure.ProteinResidue(99999, structure.SequenceAlphabets.Protein.ALA)              
        chain.residues.append(residue)
        chain.residues.append(residue)           
        
    def testContains(self):
        self.assertTrue(self.chain.residues._contains('127'))

@test.unit  
class TestResidue(test.Case):              
    
    def setUp(self):
        
        super(TestResidue, self).setUp()

        self.structure = self.config.getPickle('1nz9.model1.pickle')        
        self.chain = self.structure['A']
        self.chain.compute_torsion()
        self.residue = self.chain[0]
        
    def testGetitem(self):
        
        self.assertTrue('CA' in self.residue)
        self.assertEqual(set(self.residue), set(['C', 'H2', 'CB', 'CA', 'H1', 'O', 'N', 'H3', 'HA', 'HB1', 'HB3', 'HB2']))
        
        self.assertEqual(self.residue['CA'], self.residue.atoms['CA'])
        self.assertEqual(self.residue['CA'], tuple(self.residue.items)[1])        
        
        self.assertEqual(self.residue['CA'].element, structure.ChemElements.C)              
        self.assertEqual(self.residue['CA'].serial_number, 2)    
        
    def testItems(self):
        self.assertEqual(tuple(self.residue.items)[0], self.residue['N'])
        self.assertTrue(tuple(self.residue.items)[0] < tuple(self.residue.items)[1])  
    
    def testHasStructure(self):
        
        self.assertTrue(self.residue.has_structure)

        residue = structure.ProteinResidue(111, 'ALA')
        self.assertFalse(residue.has_structure)
        
    def testId(self):
        
        def assign(sn, ins=None):
            self.residue.id = sn, ins
        
        # should raise if the residue is part of a chain and SN=128 is already there (true)
        self.assertRaises(structure.DuplicateResidueIDError, assign, 128)
        self.assertRaises(structure.DuplicateResidueIDError, assign, '128 ')
        # should not raise (same residue)
        self.residue.id = 127, None 
        
        # sequence_number must be defined when an insertion_code is specified
        self.assertRaises(structure.InvalidOperation, assign, None, 'X')
                
        self.residue.id = 111127, 'X'
        self.assertEqual(self.chain.find(111127, 'X'), self.residue)           
        self.residue.id = 127, None
        
    def testSequenceNumber(self):
        self.assertEqual(self.residue.sequence_number, 127)
                
    def testInsertionCode(self):
        self.assertEqual(self.residue.insertion_code, None)    
        
    def testType(self):
        self.assertEqual(self.residue.type, structure.SequenceAlphabets.Protein.ALA)        

    def testTorsion(self):
        self.assertEqual(self.residue.torsion.phi, None)
        self.assertNotEqual(self.residue.torsion.psi, None)        

    def testFactory(self):
        
        factory = structure.Residue.create

        nucleotide = factory(structure.SequenceTypes.NucleicAcid, 1, structure.SequenceAlphabets.Nucleic.Adenine)   
        aminoacid = factory(structure.SequenceTypes.Protein, 1, structure.SequenceAlphabets.Protein.ALA)                  
        unknown = factory(structure.SequenceTypes.Unknown, 1, structure.SequenceAlphabets.Unknown.UNK)              
        
        self.assertTrue(isinstance(nucleotide, structure.NucleicResidue))
        self.assertTrue(isinstance(aminoacid, structure.ProteinResidue))
        self.assertTrue(isinstance(unknown, structure.UnknownResidue))
        
    def testAppend(self):
        
        # should raise if appending an atom with atom._residue (already part of another residue)        
        self.assertRaises(structure.InvalidOperation, self.residue.atoms.append, self.chain[2]['CA'])
        # should not raise InvalidOperation (same residue)
        self.assertRaises(structure.DuplicateAtomIDError, self.residue.atoms.append, self.residue['CA'])
                
        atom = structure.Atom(999999, 'Cx', structure.ChemElements.C, numpy.array([1, 1, 1]))                   
        self.residue.atoms.append(atom)
        self.assertTrue(atom.name in self.residue)
        self.assertEqual(atom.residue, self.residue)
        
        # test alternate handling
        atom2 = structure.Atom(999999, 'Cx', structure.ChemElements.C, numpy.array([2, 2, 2]), alternate=True)  
        self.residue.atoms.append(atom2)
        self.assertTrue(isinstance(self.residue['Cx'], structure.DisorderedAtom))
        self.assertEqual(self.residue['Cx'].length, 2)

        atom3 = structure.Atom(999999, 'Cx', structure.ChemElements.C, numpy.array([3, 3, 3]))             
        self.residue.atoms.append(atom3)
        self.assertEqual(self.residue['Cx'].length, 3)  

@test.unit
class TestAtom(test.Case):

    def setUp(self):
        
        super(TestAtom, self).setUp()

        self.structure = self.config.getPickle('1nz9.model1.pickle')        
        self.chain = self.structure['A']
        self.residue = self.chain[0]
        self.atom = self.residue['CA']
        
    def testName(self):
        self.assertEqual(self.atom.name, 'CA')
    
    def testSerialNumber(self):
        """
        @bug: check for overlapping SNs in a Chain
        """
        
        def assign(sn):
            self.atom.serial_number = sn
            
        self.assertEqual(self.atom.serial_number, 2)
        self.assertRaises(TypeError, assign, -2)            
        
    def testElement(self):
        self.assertEqual(self.atom.element, structure.ChemElements.C)                   
        
    def testResidue(self):
        """
        @todo: more checks are needed here (e.g. append an atom to different residues)
        """    
        self.assertEqual(self.atom.residue, self.residue) 
        
    def testSorting(self):
        self.assertTrue(self.residue['N'] < self.residue['CA'])
        
    def testVector(self):
        
        def assign(v):
            self.atom.vector = v
            
        self.assertEqual(self.atom.vector.tolist(), [-0.96399999999999997, -17.864999999999998, -6.46])
        self.assertRaises(ValueError, assign, numpy.array([1]))  

@test.unit
class TestDisorderedAtom(TestAtom):

    def setUp(self):
        
        super(TestDisorderedAtom, self).setUp()

        self.structure = self.config.getPickle('1nz9.model1.pickle')        
        self.chain = self.structure['A']
        self.residue = self.chain[0]
        self.atom = self.residue['CA']
        
        self.alt = structure.Atom(999, 'CA', structure.ChemElements.C,
                                  self.atom.vector.copy(), alternate='B')
        
        self.atom.occupancy = 0.9
        self.alt.occupancy = 0.1
        self.residue.atoms.append(self.alt)
        
        self.disatom = self.residue['CA']
        
    def testAppend(self):

        self.assertEqual(self.residue.atoms['CA'].serial_number, 2)
        
        alt = structure.Atom(1000, 'CA', structure.ChemElements.C,
                                  self.atom.vector.copy(), alternate='C')
        
        alt.occupancy = 1.0
        self.residue.atoms.append(alt)
        
        self.assertEqual(self.residue.atoms['CA'].serial_number, 1000)
        self.assertTrue(isinstance(self.residue.atoms['CA'], structure.DisorderedAtom))
        
    def testIteration(self):
        
        atoms = list(self.residue.atoms['CA'])
        self.assertEqual(len(atoms), 2)
    
    def testLength(self):
        self.assertEqual(self.disatom.length, 2)
        
    def testFind(self):
        self.assertTrue(self.disatom.find('B') is self.alt)
        self.assertRaises(structure.EntityNotFoundError, self.disatom.find, 'X')

            
@test.unit
class TestSecondaryStructure(test.Case):
    """
    @todo: implement secondary structure tests
    """
    
    def setUp(self):
        
        super(TestSecondaryStructure, self).setUp()

        self.structure = self.config.getPickle('1nz9.model1.pickle')        
        self.chain = self.structure['A']
        self.residue = self.chain[0]    
    
@test.unit
class TestTorsionAngles(test.Case):
        
    def setUp(self):
        
        super(TestTorsionAngles, self).setUp()

        self.structure = self.config.getPickle('1nz9.model1.pickle')        
        self.chain = self.structure['A']
        self.chain.compute_torsion()
        
        self.angles1 = structure.TorsionAnglesCollection()
        self.angles2 = structure.TorsionAnglesCollection()
        self.angles3 = structure.TorsionAnglesCollection()
        self.angles4 = structure.TorsionAnglesCollection()
        
        for dummy in range(10):
            self.angles1.append(structure.TorsionAngles(20, 30, 180))
            self.angles2.append(structure.TorsionAngles(20, 30, 180))
            self.angles3.append(structure.TorsionAngles(20, 30, 180))
            self.angles4.append(structure.TorsionAngles(20, 30, 180))
            
        self.angles2[5].phi = 40
        self.angles2[7].psi = 10
        self.angles3[4].phi = None
        self.angles4[1].phi = None

    def testRMSD(self):
        
        self.assertEqual(self.chain.torsion.rmsd(self.chain.torsion), 0.0)
        self.assertNotEqual(self.angles1.rmsd(self.angles2), 0.0)
        self.assertEqual(self.angles1.rmsd(self.angles4), 0.0)
        self.assertAlmostEqual(self.angles1.rmsd(self.angles2), 0.005442, places=5)
        
        self.assertRaises(structure.Broken3DStructureError, self.angles1.rmsd, self.angles3)
        self.assertRaises(ValueError, self.angles1.rmsd, structure.TorsionAnglesCollection())

    def testConversion(self):
        
        ta = structure.TorsionAngles(0, 90, 180)
        
        ta.to_radians()
        self.assertEqual(ta.phi, 0)
        self.assertAlmostEqual(ta.psi, 1.57, places=2)
        self.assertAlmostEqual(ta.omega, 3.14, places=2)

        ta.to_degrees()
        self.assertEqual(ta.phi, 0)
        self.assertEqual(ta.psi, 90)
        self.assertEqual(ta.omega, 180)
        
    def testBool(self):
        self.assertTrue(structure.TorsionAngles(20, None, None))
        self.assertFalse(structure.TorsionAngles(None, None, None))



if __name__ == '__main__':
    
    test.Console()  
