"""
Build SVG diagrams from various csb objects.
"""

import math

from csb.bio.structure import SecondaryStructure, SecStructures


class SSCartoonBuilder(object):
    """
    Creates a 2D vector diagram from a L{SecondaryStructure} object.
    
    @param ss: source secondary structure (either a SS string or a SS object)
    @type ss: str or L{SecondaryStructure}
    @param width: output width of the diagram in pixels
    @type width: int
    @param height: output height of the diagram in pixels
    @type height: int
    
    @param thickness: stroke-width (2px by default)
    @param helix: SVG color for helicies (red by default)
    @param strand: SVG color for strands (blue by default)
    @param coil: SVG color for coils (orange by default)
    @param gap: SVG color for gaps (grey by default)
    @param cap: stroke-linecap (round by default)
    """
    
    def __init__(self, ss, width, height, thickness='2px', 
                 helix='#C24641', strand='#6698FF', coil='#FF8C00', gap='#E0E0E0', 
                 cap='round'):
         
        if ss:
            if isinstance(ss, basestring):
                self._ss = SecondaryStructure(ss)
            else:
                self._ss = ss.clone()
            self._ss.to_three_state()
            self._residues = sum(e.length for e in self._ss)
            if self._residues == 0:
                raise ValueError('Zero-length secondary structure')
        else:
            raise ValueError('Invalid secondary structure')
        
        self.thickness = thickness
        self.helixcolor = helix
        self.strandcolor = strand
        self.coilcolor = coil
        self.gapcolor = gap
        self.cap = cap
        
        self._realwidth = float(width)
        self._width = self._realwidth - 2 # this is to compensate for antialiasing and rounded caps
        self._height = float(height)
        self._x = 0
        self._y = 0        
        
        self._svg = ''
        
    def build(self):
        """
        Build a SVG image using the current size and color settings.
        
        @return: SVG diagram
        @rtype: str (SVG document) 
        """
        
        self._x = 0
        self._y = 0
        self._svg = [r'''<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" 
    width="{0._realwidth}" height="{0._height}">
    
    <g transform="translate(0, {1})">'''.format(self, self._height / 2.0)]
        
        for e in self._ss:
            
            if e.type == SecStructures.Helix:
                cartoon = self._helix(e.length)
                color = self.helixcolor
                
            elif e.type == SecStructures.Strand:
                cartoon = self._strand(e.length)
                color = self.strandcolor
                
            elif e.type == SecStructures.Coil:
                cartoon = self._coil(e.length)
                color = self.coilcolor
                
            elif e.type == SecStructures.Gap:
                cartoon = self._gap(e.length)
                color = self.gapcolor
                
            else:
                assert False, "Unhandled SS Type: {0!r}".format(e.type)
        
            path = r'''        <path fill="none" stroke="{0}" stroke-width="{1.thickness}" stroke-linecap="{1.cap}"
            d="{2}" />'''.format(color, self, cartoon)
        
            self._svg.append(path)

        self._svg.append('    </g>')        
        self._svg.append('</svg>')        
        return '\n'.join(self._svg)
    
    def _format(self, path):
        
        formatted = []
        
        for i in path:
            
            if i == -0:
                i = 0
            
            if isinstance(i, float):
                i = round(i, ndigits=7)
                if i == -0:
                    i = 0
                formatted.append('{0:.7f}'.format(i))
            else:
                formatted.append(str(i))
        
        return ' '.join(formatted)
    
    def _helix(self, length, arc_width=3.0):
        
        if length < 1:
            return ''
        
        helix_width = float(length) * self._width / self._residues
        helix_end = self._x + helix_width
        path = ['M', self._x, self._y, 'Q']
        
        arcs = int(helix_width / arc_width)
        for i in range(1, arcs + 1):

            # quadratic bezier control points: sine curve's min, max and inflection points (0, 1, 0, -1, 0, 1 ...)
            # one arc is the curve from 0 to pi/2                                    
            if i < arcs:
                # inner arc
                self._x += arc_width            
                self._y = math.sin(math.pi * i / 2) * (self._height / 2.0)
                path.append(self._x)
                path.append(self._y)                 
            else:
                # last arc; stretch it to make the helix pixel-precise, ending also at y=0
                # also the number of arcs/controlpoints must be even, otherwise the path is broken
                
                # remaining pixels on x
                remainder = helix_end - self._x 
                
                if i % 2 == 0:
                    # even number of arcs, just extend the last arc with the remainder
                    self._x += remainder
                    self._y = 0                
                    path.append(self._x)
                    path.append(self._y)
                else:
                    # odd number of arcs
                    
                    #  1) keep this arc at the expected y, but stretch it half of the x remainder
                    self._x += remainder / 2.0     
                    self._y = math.sin(math.pi * i / 2) * (self._height / 2.0)
                    path.append(self._x)
                    path.append(self._y)
                    
                    #  2) append a final arc, ending at [helix_end, 0]
                    self._x += remainder / 2.0
                    self._y = 0
                    path.append(self._x)
                    path.append(self._y)  
                    
        return self._format(path)
                
    def _strand(self, length, arrow_width=3.0):
        
        offset = 1.0
        strand_width = float(length) * self._width / self._residues
        path = ['M', self._x, self._y, 'H']

        self._x += strand_width         
        path.append(self._x)
        
        if offset < arrow_width < strand_width:
            arrow_start = self._x - offset - arrow_width
            path.extend(['M', self._x - offset, self._y])
            path.extend(['L', arrow_start, self._y + self._height / 9])
            path.extend(['L', arrow_start, self._y - self._height / 9])
            path.extend(['L', self._x - offset, self._y])
                
        return self._format(path)
        
    def _coil(self, length):
        
        coil_width = float(length) * self._width / self._residues
        path = ['M', self._x, self._y, 'Q']
    
        # first control point    
        self._x += coil_width / 2.0 
        self._y = self._height / -2.0
        path.append(self._x)
        path.append(self._y)
        
        # second
        self._x += coil_width / 2.0 
        self._y = 0
        path.append(self._x)
        path.append(self._y)
        
        return self._format(path)               
    
    def _gap(self, length):
        
        return self._strand(length, arrow_width=0)
