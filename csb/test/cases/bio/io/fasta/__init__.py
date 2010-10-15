import csb.test as test

from csb.bio.io import SequenceParser
from csb.bio.sequence import Sequence
        
        
@test.unit
class TestSequenceParser(test.Case):
    
    def setUp(self):
        
        super(TestSequenceParser, self).setUp()
        
        self.parser = SequenceParser()
        self.file = self.config.getTestFile('test.fa')        
        self.data = open(self.file).read()
    
    def testParseFile(self):
        
        s = self.parser.parse_file(self.file)

        self.assertEquals(len(s), 3)

        self.assertTrue(isinstance(s[0], Sequence))
        self.assertTrue(s[0].id, 'gi|148654187')
        self.assertTrue(s[0].sequence, 'RLSRSHEFQRLRREGTRVRSGYLwCVMLQDPSLPGPAVAFAIGRPFGSAVRRNRLRRQLRSILSDRESAMGGGMFLIGVNNPHRDLPMPSFAQLTHDIDEILNK')        
        self.assertTrue(s[0].header, '>gi|142801636|gb|EDA68688.1 hypothetical protein GOS_1956086')
            
    def testParseString(self):

        s = self.parser.parse_string(self.data)
                                           
        self.assertEquals(len(s), 3)

        self.assertTrue(isinstance(s[0], Sequence))
        self.assertTrue(s[0].id, 'gi|148654187')
        self.assertTrue(s[0].sequence, 'RLSRSHEFQRLRREGTRVRSGYLwCVMLQDPSLPGPAVAFAIGRPFGSAVRRNRLRRQLRSILSDRESAMGGGMFLIGVNNPHRDLPMPSFAQLTHDIDEILNK')
        self.assertTrue(s[0].header, '>gi|142801636|gb|EDA68688.1 hypothetical protein GOS_1956086')        

        
if __name__ == '__main__':
    
    test.Console()
