import csb.test as test

from csb.bio.io.cs import ChemShiftReader, ChemShift3Reader, ChemShiftFormatError
from csb.bio.structure import ChemElements
from csb.bio.sequence import ProteinAlphabet



@test.unit
class TestChemShiftReader(test.Case):
    
    def setUp(self):
        
        super(TestChemShiftReader, self).setUp()
        
        self.parser = self.klass()
        self.file2 = self.config.getTestFile('2l01.v2.str')
        self.file3 = self.config.getTestFile('2l01.v3.str')
        
    @property
    def file(self):
        return self.file2
    
    @property
    def klass(self):
        return ChemShiftReader
    
    def testCreate(self):
        
        klass = self.klass
        
        self.assertTrue(isinstance(klass.create(version=2), ChemShiftReader))
        self.assertTrue(isinstance(klass.create(version=3), ChemShift3Reader))
        
        self.assertRaises(ValueError, klass.create, version=1)
        
    def testGuess(self):
        
        klass = self.klass
                
        self.assertTrue(isinstance(klass.guess(self.file2), ChemShiftReader))
        self.assertTrue(isinstance(klass.guess(self.file3), ChemShift3Reader))
        
        dummy = self.config.getTestFile("2JZC.sum")
        self.assertRaises(ChemShiftFormatError, klass.guess, dummy)
        
    def testReadShifts(self):
        
        content = open(self.file).read()
        cs = self.parser.read_shifts(content)
        
        self.assertEqual(len(cs), 11)
        
        self.assertEqual(cs[0].name, "HA")
        self.assertEqual(cs[0].element, ChemElements.H)
        self.assertEqual(cs[0].shift, 3.977)
        
        self.assertEqual(cs[1].name, "HB2")
        self.assertEqual(cs[1].shift, 2.092)
         
        self.assertEqual(cs[7].element, ChemElements.C)
        self.assertEqual(cs[7].residue, ProteinAlphabet.MET)
        
        self.assertEqual(cs[10].residue, ProteinAlphabet.LYS)
        self.assertEqual(cs[10].shift, 4.423)
        
    def testReadFile(self):
        
        cs = self.parser.read_file(self.file)   
        self.assertEqual(len(cs), 11)
        
@test.unit
class TestChemShift3Reader(TestChemShiftReader):
    
    @property
    def file(self):
        return self.file3

    @property
    def klass(self):
        return ChemShift3Reader        
    
    
if __name__ == '__main__':
    
    test.Console()
    