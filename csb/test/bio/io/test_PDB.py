import os
import sys
import glob
import unittest

from csb.bio.io import StructureParser


class PDBEntryTestCase(unittest.TestCase):
    
    def runTest(self):
        parser = StructureParser(self.file)
        try:
            parser.parse_structure()
        except:
            ex = sys.exc_info()
            ex[1].args = list(ex[1].args) + [self.file]
            raise ex[0], ex[1], ex[2]


def buildSuite(path_mask):
    """
    @param path_mask: search path and file mask for the PDB files
    @type path_mask: str
    """
    
    suite = unittest.TestSuite()
    
    for file in glob.glob(path_mask):
        
        case = PDBEntryTestCase()
        case.file = file
        suite.addTest(case)
    
    return suite


if __name__ == '__main__':
    
    usage = r"""
 python {0} path_mask
     
 Test StructureParser over the entire PDB. Options:
 
     path_mask - search path and file mask for the PDB files
 
    """.format(os.path.basename(__file__))
    
    if len(sys.argv) != 2:
        print usage
        sys.exit(1)        
        
    path = sys.argv[1]    

    if not os.path.isdir(os.path.dirname(path)):
        print '\n ERROR: directory not found: {0}\n'.format(os.path.dirname(path))
        print usage
        sys.exit(1)

    print 'Loading tests...',     
    sys.stdout.flush()
       
    suite = buildSuite(path)

    print '{0} tests loaded.\n'.format(suite.countTestCases())
    
    unittest.TextTestRunner().run(suite)
