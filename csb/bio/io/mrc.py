"""
Cryo-EM density map I/O

@warning: dragons ahead, this module is experimental
"""

import numpy
import struct


class DensityMapFormatError(ValueError):
    pass


class ByteOrder(object):

    NATIVE = '='    
    LITTLE = '<'
    BIG = '>'


class DensityInfo(object):

    def __init__(self, data, spacing, origin, shape=None, header=None, axes=None):
        
        self.data = data
        self.spacing = spacing
        self.origin = origin
        self.header = header
        self.shape = shape
        self.axes = axes
        
        if shape is None and data is not None:
            self.shape = self.data.shape
        
class HeaderInfo(object):
    
    def __init__(self, fields):
        
        fields = tuple(fields)
        if not len(fields) == 25:
            raise ValueError(fields)
        
        self._fields = fields
        
    def __getitem__(self, i):
        return self._fields[i]
    
    def __iter__(self):
        return iter(self._fields)

    @property
    def nc(self):
        return self._fields[0]
    @property
    def nr(self):
        return self._fields[1]
    @property
    def ns(self):
        return self._fields[2]
    @property
    def mode(self):
        return self._fields[3]
    @property
    def ncstart(self):
        return self._fields[4]
    @property
    def nrstart(self):
        return self._fields[5]
    @property
    def nsstart(self):
        return self._fields[6]
    @property
    def nx(self):
        return self._fields[7]
    @property
    def ny(self):
        return self._fields[8]
    @property
    def nz(self):
        return self._fields[9]
    @property
    def x(self):
        return self._fields[10]
    @property
    def y(self):
        return self._fields[11]
    @property
    def z(self):
        return self._fields[12]
    @property
    def alpha(self):
        return self._fields[13]
    @property
    def beta(self):
        return self._fields[14]
    @property
    def gamma(self):
        return self._fields[15]
    @property
    def mapc(self):
        return self._fields[16]
    @property
    def mapr(self):
        return self._fields[17]
    @property
    def maps(self):
        return self._fields[18]
    @property
    def amin(self):
        return self._fields[19]
    @property
    def amax(self):
        return self._fields[20]
    @property
    def amean(self):
        return self._fields[21]
    @property
    def ispg(self):
        return self._fields[22]
    @property
    def nsymbt(self):
        return self._fields[23]
    @property
    def lskflg(self):
        return self._fields[24]
        
    
class DensityMapReader(object):
    """
    Binary MRC density map reader.
    
    @param filename: input MRC file name
    @type filename: str
    """
    
    HEADER_SIZE = 1024
    
    def __init__(self, filename):
        self._filename = filename
        
    @property
    def filename(self):
        """
        Input MRC file name
        """
        return self._filename
    
    def _rawheader(self, stream):
        """
        Read and return the raw binary header.
        """
        raw = stream.read(DensityMapReader.HEADER_SIZE)
        return bytes(raw)

    def _inspect(self, rawheader, order):
        """
        Parse a raw binary header.
        """
        format = '{0}10l6f3l3f3l'.format(order)
        fields = struct.unpack(format, rawheader[:4 * 25])
        
        return HeaderInfo(fields)
    
    def _spacing(self, header):

        if header.nx != 0 and header.ny != 0 and header.nz != 0:
            return (header.x / header.nx, header.y / header.ny, header.z / header.nz)
        else:
            return (0, 0, 0)
    
    def _origin(self, header, spacing=None):
        
        if spacing is None:
            spacing = self._spacing(header)
            
        origin = header.ncstart, header.nrstart, header.nsstart
        return [origin[i] * spacing[i] for i in range(3)]
    
    def _shape(self, header):
        return (header.ns, header.nr, header.nc)
    
    def inspect_header(self, order=ByteOrder.NATIVE):
        """
        Parse the raw binary header of the density map.
        
        @param order: byte order (defaults to L{ByteOrder.NATIVE})
        @type order: str
        
        @return: header information
        @rtype: L{HeaderInfo}
        """
        with open(self.filename, 'rb') as stream:        
            raw = self._rawheader(stream)
            return self._inspect(raw, order)

    def read_header(self):
        """
        Read the header of the density map only.
        
        @return: density info without any actual data (density.data is None)
        @rtype: L{DensityInfo}
        """
        with open(self.filename, 'rb') as stream:
    
            raw = self._rawheader(stream)

            header = self._inspect(raw, ByteOrder.NATIVE)
            spacing = self._spacing(header)
            origin = self._origin(header, spacing)
            shape = self._shape(header)

            return DensityInfo(None, spacing, origin, shape=shape, header=raw)
        
    def read(self):
        """
        Read the entire density map.
        
        @return: complete density info
        @rtype: L{DensityInfo}
        """
        with open(self.filename, 'rb') as stream:
    
            raw = self._rawheader(stream)
            header = self._inspect(raw, ByteOrder.NATIVE)
                
            if header.mode == 2 or header.mode == 1:
                byte_order = ByteOrder.NATIVE
        
            elif header.mode == 33554432:
                header = self._inspect(raw, ByteOrder.BIG)
                byte_order = ByteOrder.BIG
                        
                if header.mode == 33554432:
                    header = self._inspect(raw, ByteOrder.LITTLE)
                    byte_order = ByteOrder.LITTLE
        
            else:
                raise DensityMapFormatError("Not a mode 2 CCP4 map file")
        
            stream.read(header.nsymbt)  # symmetry_data
            
            count = header.ns * header.nr * header.nc
            map_data = stream.read(4 * count)
            
            if byte_order == ByteOrder.NATIVE:
                array = numpy.fromstring(map_data, numpy.float32, count)
            else:
                array = numpy.zeros((count,), numpy.float32)
                index = 0
                while len(map_data) >= 4 * 10000:
                    values = struct.unpack(byte_order + '10000f', map_data[:4 * 10000])
                    array[index:index + 10000] = numpy.array(values, numpy.float32)
                    index += 10000
                    map_data = map_data[4 * 10000:]
                values = struct.unpack(byte_order + '%df' % (len(map_data) / 4), map_data)
                array[index:] = numpy.array(values, numpy.float32)
        
            del map_data
        
            array.shape = self._shape(header)
            data = array.T
        
            spacing = self._spacing(header)
            origin = self._origin(header, spacing)
            
            return DensityInfo(data, spacing, origin, header=raw)


class DensityMapWriter(object):
    """
    Binary MRC density map writer.
    """
        
    def reconstruct_header(self, density):
        """
        Attempt to reconstruct the header, given L{DensityInfo}'s
        data shape, spacing and origin.
        
        @param density: density info
        @type density: L{DensityInfo}
        
        @return: reconstructed binary header
        @rtype: bytes
        """
        N = list(density.data.shape)
        MODE = 2
    
        if isinstance(density.spacing, float):
            spacing = 3 * [density.spacing]
        else:
            spacing = density.spacing
            
        if density.origin is None:
            origin = 3 * [0.]
        else:
            origin = density.origin

        if density.axes is None:
            MAP = list(range(1, 4))
        else:
            MAP = list(density.axes)            
            
        start = [int(round(origin[i] / spacing[i], 0)) for i in range(3)]
        M = [density.data.shape[i] for i in range(3)]
        cella = [density.data.shape[i] * spacing[i] for i in range(3)]
        cellb = 3 * [90.]

        stats = [density.data.min(), density.data.max(), density.data.mean()]
    
        ISPG = 0 
        NSYMBT = 0
        LSKFLG = 0
    
        JUNK = [0] * 25
        ORIGIN = [0., 0., 0.]
        MACHST = 0
        
        args = N + [MODE] + start + M + cella + cellb + \
               MAP + stats + [ISPG, NSYMBT, LSKFLG] + JUNK + \
               ORIGIN + [0, MACHST, 0., 0] + [b' ' * 796]

        return struct.pack('=10l6f3l3f3l25l3f2l1f1l796s', *args)
    
    def write(self, stream, density):
        """
        Write C{density} to a binary C{stream}.
        
        @param stream: destination binary stream
        @type stream: stream        
        @param density: input density info
        @type density: L{DensityInfo}
        """
        if density.header is not None:
            stream.write(density.header)
        else:
            stream.write(self.reconstruct_header(density))
    
        data = density.data.T.flatten()
        x = struct.pack('=%df' % len(data), *data.tolist())
        stream.write(x)

    def write_file(self, filename, density):
        """
        Write C{density} to a binary file.
        
        @param filename: destination file name
        @type filename: str        
        @param density: input density info
        @type density: L{DensityInfo}             
        """
        with open(filename, 'wb') as stream:
            self.write(stream, density)
