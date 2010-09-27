import os
import sys
import glob
import unittest

from csb.bio.io import AsyncStructureParser


class PDBEntryTestCase(unittest.TestCase):
    
    def runTest(self):
        parser = AsyncStructureParser()
        try:
            parser.parse_structure(self.file, 10)
        except:
            ex = sys.exc_info()
            ex[1].args = list(ex[1].args) + [self.file]
            print '\n{0}    {1}'.format(ex[0].__name__, ex[1].args)
            raise ex[0], ex[1], ex[2]


def buildSuite(path_mask, top=None):
    """
    @param path_mask: search path and file mask for the PDB files
    @type path_mask: str
    """
    
    suite = unittest.TestSuite()
    
    for i, file in enumerate(glob.glob(path_mask), start=1):
        
        if top and i > top:
            break
        
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
    
    if len(sys.argv[1:]) not in (1, 2):
        print usage
        sys.exit(1)        
        
    path = sys.argv[1]
    top = None
    if len(sys.argv[1:]) == 2:
        top = int(sys.argv[2])

    if not os.path.isdir(os.path.dirname(path)):
        print '\n ERROR: directory not found: {0}\n'.format(os.path.dirname(path))
        print usage
        sys.exit(1)

    print 'Loading tests...',     
    sys.stdout.flush()
       
    suite = buildSuite(path, top)

    print '{0} tests loaded.\n'.format(suite.countTestCases())
    
    unittest.TextTestRunner().run(suite)
