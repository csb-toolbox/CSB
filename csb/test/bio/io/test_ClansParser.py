import unittest
from csb.bio.io import ClansParser


class TestClansParser(unittest.TestCase):

    def setUp(self):
        filename = 'data/test_ClansParser.clans'
        cp = ClansParser()
        self.clans_instance = cp.parse_file(filename)

    def testParseFile(self):
        from numpy import array
        self.assertEqual(len(self.clans_instance), 41)
        self.assertRaises(IndexError, self.clans_instance.__getitem__, 41)

        correct_rotmtx = array([[0.75614862, 0.65439992, 0.],
                                [-0.65439992, 0.75614862, 0.],
                                [0., 0., 1.]])
        self.assertEqual(self.clans_instance.rotmtx.shape, (3, 3))
        self.assertTrue(
            (self.clans_instance.rotmtx - correct_rotmtx < 1e-6).all())

if __name__ == '__main__':
    unittest.main()
