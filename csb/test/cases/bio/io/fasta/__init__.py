import csb.test as test

from csb.bio.io import SequenceParser
from csb.bio.sequence import Sequence


@test.functional
class TestSequenceParser(test.Case):
    
    def runTest(self):
        
        parser = SequenceParser()
        file = self.config.getTestFile('test.fa')
        
        s = parser.parse_file(file)
        self.assertEquals(len(s), 1)        
        
        
@test.unit
class TestSequenceParserMethods(test.Case):
    
    def setUp(self):
        
        super(TestSequenceParserMethods, self).setUp()
        
        self.parser = SequenceParser()
        self.file = self.config.getTestFile('test.fa')        
        self.data = open(self.file).read()
    
    def testParseFile(self):
        
        s = self.parser.parse_file(self.file)
        self.assertEquals(len(s), 1)
        self.assertTrue(isinstance(s[0], Sequence))
    
    def testParseString(self):
        
        s = self.parser.parse_string(self.data)
        self.assertEquals(len(s), 1)
        self.assertTrue(isinstance(s[0], Sequence))
        
        
if __name__ == '__main__':
    
    test.Console()
