import csb.test as test
import csb.io

from csb.bio.io.noe import SparkyPeakListReader, XeasyPeakListReader, XeasyFileBuilder
from csb.bio.structure import ChemElements
from csb.bio.sequence import ProteinAlphabet


@test.unit
class TestSparkyPeakListReader(test.Case):
    
    def setUp(self):
        
        super(TestSparkyPeakListReader, self).setUp()
        
        self.elements = (ChemElements.H, ChemElements.C, ChemElements.H)
        self.parser = SparkyPeakListReader(self.elements, [(1, 2)])
        self.file = self.config.getTestFile('Sparky.peaks')
    
    def testRead(self):
        
        content = open(self.file).read()
        spectrum = self.parser.read(content)
        
        self.assertEqual(len(spectrum), 3)
        
        self.assertEqual(spectrum.min_intensity, 147454)
        self.assertEqual(spectrum.max_intensity, 204746)
        
        self.assertEqual(spectrum.element(0), self.elements[0])
        self.assertEqual(spectrum.element(1), self.elements[1])
        
        self.assertEqual(spectrum.dimensions, self.elements)
        self.assertEqual(spectrum.proton_dimensions, (0, 2))
        self.assertEqual(spectrum.num_dimensions, 3)
        self.assertEqual(spectrum.num_proton_dimensions, 2)
        
        self.assertFalse(spectrum.has_element(ChemElements.Ca))
        self.assertTrue(spectrum.has_element(ChemElements.C))

        self.assertFalse(spectrum.has_connected_dimensions(0))
        self.assertEqual(spectrum.connected_dimensions(0), ())
        self.assertTrue(spectrum.has_connected_dimensions(1))        
        self.assertEqual(spectrum.connected_dimensions(1), (2,))
        self.assertTrue(spectrum.has_connected_dimensions(2))        
        self.assertEqual(spectrum.connected_dimensions(2), (1,))
        
        peaks = list(spectrum)
        self.assertEqual(peaks[0].intensity, 157921)
        self.assertEqual(peaks[0].get(0), 3.418)
        self.assertEqual(peaks[0].get(1), 114.437)
        self.assertEqual(peaks[0].get(2), 7.440)
    
    def testReadFile(self):

        spectrum = self.parser.read_file(self.file)
        self.assertEqual(len(spectrum), 3)
    
    def testReadAll(self):
        
        spectrum = self.parser.read_all([self.file, self.file])
        self.assertEqual(len(spectrum), 6)        


@test.unit
class TestXeasyPeakListReader(test.Case):
    
    def setUp(self):
        
        super(TestXeasyPeakListReader, self).setUp()

        self.elements = (ChemElements.H, ChemElements.C, ChemElements.H)        
        self.parser = XeasyPeakListReader()
        self.file = self.config.getTestFile('Xeasy1.peaks')
    
    def testRead(self):

        content = open(self.file).read()
        spectrum = self.parser.read(content)
        
        self.assertEqual(len(spectrum), 3)
        
        self.assertEqual(spectrum.min_intensity, 1.291120e05)
        self.assertEqual(spectrum.max_intensity, 4.243830e05)
        
        self.assertEqual(spectrum.element(0), self.elements[0])
        self.assertEqual(spectrum.element(1), self.elements[1])
        
        self.assertEqual(spectrum.dimensions, self.elements)
        self.assertEqual(spectrum.proton_dimensions, (0, 2))
        self.assertEqual(spectrum.num_dimensions, 3)
        self.assertEqual(spectrum.num_proton_dimensions, 2)
        
        self.assertFalse(spectrum.has_element(ChemElements.Ca))
        self.assertTrue(spectrum.has_element(ChemElements.C))

        self.assertFalse(spectrum.has_connected_dimensions(0))
        self.assertEqual(spectrum.connected_dimensions(0), ())
        self.assertTrue(spectrum.has_connected_dimensions(1))        
        self.assertEqual(spectrum.connected_dimensions(1), (2,))
        self.assertTrue(spectrum.has_connected_dimensions(2))        
        self.assertEqual(spectrum.connected_dimensions(2), (1,))        

        peaks = list(spectrum)
        self.assertEqual(peaks[0].intensity, 1.565890e05)
        self.assertEqual(peaks[0].get(0), 7.050)
        self.assertEqual(peaks[0].get(1), 10.374)
        self.assertEqual(peaks[0].get(2), 0.889)        

@test.unit
class TestXeasyPeakListReader2(TestXeasyPeakListReader):
    
    def setUp(self):
        
        super(TestXeasyPeakListReader2, self).setUp()

        self.elements = (ChemElements.H, ChemElements.C, ChemElements.H)        
        self.parser = XeasyPeakListReader()
        self.file = self.config.getTestFile('Xeasy2.peaks')


@test.unit
class TestXeasyFileBuilder(test.Case):
    
    def setUp(self):
        super(TestXeasyFileBuilder, self).setUp()
        
    def testBuild(self):
        
        content = self.config.getContent("Xeasy1.peaks")
        spectrum = XeasyPeakListReader().read(content)
        
        stream = csb.io.MemoryStream()
        
        builder = XeasyFileBuilder(stream)
        builder.add_header(spectrum)
        builder.add_peaks(spectrum)
        
        self.assertEqual(stream.getvalue().strip(), content.strip())
        
        


if __name__ == '__main__':
    
    test.Console()
    
            