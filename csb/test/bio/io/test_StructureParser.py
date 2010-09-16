import unittest
from csb.bio.io.wwpdb import StructureParser,get
from csb.pyutils import CollectionIndexError
from csb.bio.sequence import SequenceAlphabets

class TestStrutureParser(unittest.TestCase):

    def setUp(self):
        filename = 'data/test_pdb.pdb'
        self.cp = StructureParser(filename)


    def testParseFile(self):
        from numpy import array
        structure = self.cp.parse()

        self.assertEqual(structure.accession, '1d3z')
        # Chain level
        self.assertEqual(structure.chains.length, 1)
        self.assertEqual(len(structure.chains), 1)

        self.assertEqual(len(structure.chains['A']),76)
        self.assertEqual(len(structure['A']),76)

        # Residue level 
        self.assertEqual(len(structure['A'][1:10]),9)
        self.assertEqual(structure['A'][1].type,SequenceAlphabets.Protein.MET)

    def test_get(self):
        structure = get('1d3z')
        # Chain level
        self.assertEqual(structure.accession, '1d3z')
        # Chain level
        self.assertEqual(structure.chains.length, 1)
        self.assertEqual(len(structure.chains), 1)

        self.assertEqual(len(structure.chains['A']),76)
        self.assertEqual(len(structure['A']),76)

        # Residue level 
        self.assertEqual(len(structure['A'][1:10]),9)
        self.assertEqual(structure['A'][1].type,SequenceAlphabets.Protein.MET)


if __name__ == '__main__':
    unittest.main()
