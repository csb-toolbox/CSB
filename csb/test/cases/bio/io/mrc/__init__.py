import csb.test as test

from csb.io import MemoryStream
from csb.bio.io.mrc import DensityMapReader, DensityMapWriter, DensityMapFormatError, HeaderInfo, ByteOrder

        
@test.unit
class TestDensityMapReader(test.Case):
    
    def setUp(self):
        
        super(TestDensityMapReader, self).setUp()
        
        self.file = self.config.getTestFile('1C3W_10.mrc')
        self.reader = DensityMapReader(self.file)
        self.rawheader = None
        
        with open(self.file, 'rb') as stream:
            self.rawheader = self.reader._rawheader(stream) 
        
    def testReadRawHeader(self):
        self.assertEqual(len(self.rawheader), DensityMapReader.HEADER_SIZE)
        
    def testReadHeader(self):
        
        density = self.reader.read_header()

        self.assertEqual(density.data, None)
        self.assertEqual(density.header, self.rawheader)
        self.assertEqual(density.origin, [-36.0, -36.0, -36.0])
        self.assertEqual(density.shape, (72, 72, 72))
        self.assertEqual(density.spacing, (1.0, 1.0, 1.0))
        
    def testRead(self):
        
        density = self.reader.read()

        self.assertNotEqual(density.data, None)
        self.assertEqual(density.header, self.rawheader)
        self.assertEqual(density.origin, [-36.0, -36.0, -36.0])
        self.assertEqual(density.shape, (72, 72, 72))
        self.assertEqual(density.spacing, (1.0, 1.0, 1.0))
        
        
@test.unit
class TestDensityMapWriter(test.Case):
    
    def setUp(self):
        
        super(TestDensityMapWriter, self).setUp()
        
        self.file = self.config.getTestFile('1C3W_10.mrc')
        self.writer = DensityMapWriter()
        self.reader = DensityMapReader(self.file)
        self.density = self.reader.read()        
        
    def testWriteDensity(self):
        
        with self.config.getTempStream(mode='b') as temp:
            with open(self.file, 'rb') as source:
                self.writer.write(temp, self.density)
                temp.flush()
                if temp.content() != source.read(): 
                    self.fail('binary strings differ')
                    
    def testReconstructHeader(self):
        
        raw = self.density.header
        self.density.header = None
        
        new = self.writer.reconstruct_header(self.density)
        
        original = self.reader._inspect(raw, ByteOrder.NATIVE)
        generated = self.reader._inspect(new, ByteOrder.NATIVE)
        
        for o, g in zip(original, generated):
            self.assertAlmostEqual(o, g, places=4)

        
                        
if __name__ == '__main__':
    
    test.Console()
