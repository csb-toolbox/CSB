"""
Simple XEASY and Sparky peak list parsers.
"""

from abc import ABCMeta, abstractmethod
from csb.bio.nmr import NOESpectrum


class PeakListFormatError(ValueError):
    pass

class BasePeakListReader(object):
    
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def read(self, table):
        """
        Parse a peak list table.
        
        @param table: input peak list table
        @type table: str
        @rtype: L{NOESpectrum}
        """
        pass
    
    def read_file(self, filename):
        """
        Parse a peak list file.
        
        @param filename: input file name
        @type filename: str
        @rtype: L{NOESpectrum}
        """
        with open(filename) as input:
            return self.read(input.read())
        
    def read_all(self, filenames):
        """
        Parse a list of peak list files and merge the resulting spectra.
        All spectra must have identical dimensions.
        
        @param filenames: input file names
        @type filenames: iterable of str
        
        @return: joint spectrum
        @rtype: L{NOESpectrum}
        """
        spectra = [self.read_file(f) for f in filenames]
        return NOESpectrum.join(*spectra)
        
class SparkyPeakListReader(BasePeakListReader):
    """
    Sparky NOE peak list parser.
    
    @param elements: list of element names for each dimension
    @type elements: list of (str or L{EnumItem})
    @param connected: list of covalently connected dimension indices in the
                      format: [(i1,i2),...]
    @type connected: list of (int,int) tuples
    """
    
    def __init__(self, elements, connected):
        
        self._elements = list(elements)
        self._connected = [(d1, d2) for d1, d2 in connected]
        
        if len(self._elements) < 1:
            raise ValueError("Can't parse a 0-dimensional peak list")
    
    def read(self, table):
        """
        Parse a Sparky peak list table.
        
        @param table: input peak list
        @type table: str
        @rtype: L{NOESpectrum}
        """
        offset = 0                
        spectrum = NOESpectrum(self._elements)
        
        for d1, d2 in self._connected:
            spectrum.connect(d1, d2)
        
        for l in table.splitlines():
            if not l.strip() or ('w1' in l and 'w2' in l):
                if l.lstrip().lower().startswith('assignment'):
                    offset = 1
                continue
            
            line = l.split()[offset:]
            try:
                float(line[-1])             # last item may or may not be a comment
            except ValueError:
                if len(line) > 0:
                    line.pop()
            
            items = list(map(float, line))
            intensity = items[-1]
            dimensions = items[:-1]
            
            if len(dimensions) != len(self._elements):
                raise PeakListFormatError("Expected {0} dimensional spectrum, got {1}".format(
                                                                    len(self._elements), len(dimensions)))
            
            spectrum.add(intensity, dimensions)
        
        return spectrum       
                    
class XeasyPeakListReader(BasePeakListReader):
    """
    XEASY NOE peak list parser.
    """
        
    def __init__(self):
        pass
            
    def read(self, table):
        """
        Parse an XEASY peak list table.
        
        @param table: input peak list
        @type table: str
        @rtype: L{NOESpectrum}
        """      
        lines = table.splitlines()  
        spectrum = self._read_header(lines)
        
        for l in lines:
            if not l.strip() or l.startswith('#'):
                continue
            
            parts = l.split()[1:]
            peak = parts[:spectrum.num_dimensions]
            height = parts[spectrum.num_dimensions + 2] 
            
            intensity = float(height)
            dimensions = map(float, peak)
            
            spectrum.add(intensity, dimensions)
        
        return spectrum           
    
    
    def _read_header(self, lines):
        
        num = 0
        dim = {}
        el = {}
        el2 = {}
        connectivity = None
        
        for l in lines:
            if l.startswith('#'):
                if l[1:].lstrip().lower().startswith('number of dimensions'):
                    num = int(l.split()[-1])
                
                if l.startswith('#INAME'):
                    parts = l.split()[1:]
                    if len(parts) != 2:
                        raise PeakListFormatError("Invalid Xeasy header")                    
                    
                    index = int(parts[0]) - 1
                    if index < 0:
                        raise PeakListFormatError("Invalid Xeasy header: dimension index < 1")
                    
                    element = ''.join(i for i in parts[1] if i.isalpha())
                    el[parts[1]] = index
                    el2[element] = index
                     
                    dim[index] = element
                    
                if l.startswith('#CYANAFORMAT'):
                    connectivity = l.split()[1]
      
        if len(dim) != num or num == 0:
            raise PeakListFormatError("Invalid Xeasy header")
        
        elements = tuple(dim[i] for i in sorted(dim))
        spectrum = NOESpectrum(elements)
        
        if connectivity:
            # cyanaformat - explicitly defines connected dimensions:
            # upper case dimensions are connected, e.g. "#CYANAFORMAT hCH" => 2-3 
            if connectivity.upper() != ''.join(elements).upper():
                raise ValueError("Invalid XEASY/CYANA header") 
            for i1 in range(len(connectivity)):
                for i2 in range(len(connectivity)):
                    e1, e2 = connectivity[i1], connectivity[i2] 
                    if i1 != i2 and e1.isupper() and e2.isupper():
                        spectrum.connect(i1, i2)                        
        else:
            # dimension labels starting with a number are connected, e.g. "1A B2 3C" => 1-3
            if len(el) != num:
                raise PeakListFormatError("Invalid XEASY header")                
            for e1 in el:
                for e2 in el:
                    if e1 != e2:
                        element1 = dim[el[e1]]
                        element2 = dim[el[e2]]
                        
                        num1 = e1.replace(element1, '')
                        num2 = e2.replace(element2, '')
                        
                        if e1.startswith(num1) and e2.startswith(num2):
                            spectrum.connect(el[e1], el[e2])
                        
        return spectrum


class XeasyFileBuilder(object):
    """
    XEASY output format builder.
    
    @param stream: destination stream, were the output is written
    @type stream: file
    """
    
    def __init__(self, stream):
        self._out = stream
        
    def add_spectrum(self, spectrum):
        
        self.add_header(spectrum)
        self.add_peaks(spectrum)
        
    def add_header(self, spectrum):
        """
        Write the XEASY header.
        
        @param spectrum: NOE spectrum
        @type spectrum: L{NOESpectrum} 
        """
        
        self._out.write(
            '# Number of dimensions {0}\n'.format(spectrum.num_dimensions))
        
        conn = ''
        
        for en, e in enumerate(spectrum.dimensions, start=1):
            element = repr(e).upper()
            self._out.write('#INAME {0} {1}{0}\n'.format(en, element))
            
            if spectrum.has_connected_dimensions(en - 1):
                conn += element.upper()
            else:
                conn += element.lower()
                
        self._out.write(
            '#CYANAFORMAT {0}\n'.format(conn))
    
    def add_peaks(self, spectrum):
        """
        Write all peaks from C{spectrum}.
        
        @param spectrum: NOE spectrum
        @type spectrum: L{NOESpectrum} 
        """
        
        for pn, peak in enumerate(spectrum, start=1):
            self._out.write("{0:5} ".format(pn))
            
            for dim in range(spectrum.num_dimensions):
                data = "{0:7.3f} ".format(peak.get(dim))
                self._out.write(data)
                            
            self._out.write("2 U ")
            self._out.write("{0:18e} ".format(peak.intensity))
            self._out.write("0.00e+00 m   0    0    0    0 0\n")
