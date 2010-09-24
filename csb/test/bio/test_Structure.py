import unittest



class TestStructure(unittest.TestCase):

    def setUp(self):
        from csb.bio.io.wwpdb import StructureParser
        
        regular_file = 'io/data/test_pdb.pdb'
        self.rp = StructureParser(regular_file)
        self.structure = self.rp.parse()


    def test_apply_transformation(self):
        from numpy import eye, zeros
        
        R = eye(3)
        t = zeros((3,))
        position_before = zeros((3,))
        position_after = zeros((3,))
        
        position_before[:]= self.structure['A'][1]['CA'].vector[:]         
        self.structure.apply_transformation(R,t)
        position_after[:]= self.structure['A'][1]['CA'].vector[:]
        self.assertEqual(position_after[0], position_before[0])
        self.assertEqual(position_after[1], position_before[1])
        self.assertEqual(position_after[2], position_before[2])

    def test_Chain_compute_torsion(self):
        
        c = self.structure['A']
        c.compute_torsion()
        self.assertTrue(c.torsion)
        self.assertTrue(len(c.torsion) > 1) 
        
if __name__ == '__main__':
    unittest.main()

        
    

        
