"""
Simple NMR STAR chemical shift readers.
"""

from csb.bio.nmr import ChemShiftInfo


class ChemShiftFormatError(ValueError):
    pass


class ChemShiftReader(object):
    """
    Simple NMR STAR v2 chemical shift reader. 
    
    @note: This is not a full-fledged, semantic NMR STAR parser. It handles
           only the chemical shift table. 
    """
    
    FRAME = 'save_assigned_chemical_shifts'
    
    RANK = '_Residue_seq_code'
    RESIDUE = '_Residue_label'
    ATOM = '_Atom_name'
    ELEMENT = '_Atom_type'
    SHIFT = '_Chem_shift_value'
    
    @staticmethod
    def create(frame=FRAME, version=2):
        """
        Parser factory: create a new parser, given a saveframe name 
        and format verison.
        
        @param frame: name of the saveframe to read
        @type frame: str
        @param version: NMR STAR format version
        @type version: int
        
        @return: an instance of any L{ChemShiftReader} class
        @rtype: L{ChemShiftReader}
        """
        
        if version == 3:
            return ChemShift3Reader(frame=frame)
        elif version == 2:
            return ChemShiftReader(frame=frame)
        else:
            raise ValueError('Unknown NMR-STAR version')
        
    @staticmethod
    def guess(file, frame=FRAME):
        """
        Parser factory: try to guess the correct NMR STAR version from a given
        file and create an appropriate parser.
        
        @param file: NMR STAR path and file name
        @type file: str 
        @param frame: name of the saveframe to read
        @type frame: str
        
        @return: an instance of any L{ChemShiftReader} class
        @rtype: L{ChemShiftReader}
        
        @raise ChemShiftFormatError: on failure to determine the NMR STAR version
        """
                
        with open(file) as cs:
            content = cs.read()
            
            if not content.strip():
                return ChemShiftReader.create()
            elif ChemShift3Reader.SHIFT3 in content:
                return ChemShiftReader.create(frame, version=3)
            elif ChemShiftReader.SHIFT in content:
                return ChemShiftReader.create(frame, version=2)
            else:
                raise ChemShiftFormatError("Can't guess NMR-STAR version")
            
    def __init__(self, frame=FRAME):
        self._frame = frame
    
    def read_file(self, filename):
        """
        Parse the specified file.
        
        @param filename: file path and name
        @type filename: str
        
        @rtype: tuple of L{ChemShiftInfo}     
        """
        with open(filename) as input:
            return self.read_shifts(input.read())
    
    def read_shifts(self, star_table):
        """
        Parse a given NMR STAR chemical shift table.
        
        @param star_table: NMR STAR chemical shift table
        @type star_table: str
        
        @rtype: tuple of L{ChemShiftInfo} 
        @raise ChemShiftFormatError: on parse error    
        """
                
        shifts = []
        
        init = False
        in_shifts = False
        fields = []
        lines = iter(star_table.splitlines())
        
        if self._frame in star_table:
            self._scroll(lines, self._frame)

        
        for l in lines:
            ls = l.strip()
            
            if not in_shifts:

                if ls == 'loop_':
                    assert in_shifts is False and not fields and init is False
                    init = True
                    continue

                elif init and ls.startswith('_'):
                    assert in_shifts is False
                    fields.append(l.strip())
                    continue
                
                elif init and not ls:
                    if len(fields) < 1:
                        raise ChemShiftFormatError("No fields found in the CS table")             
                    in_shifts = True
                    continue
                    
            else:
                
                if ls == 'stop_':
                    break
                
                elif ls.startswith('#'):
                    continue
                
                elif ls:
                    values = l.split()
                    if len(values) < len(fields):
                        raise ChemShiftFormatError("Insufficient number of values: {0}".format(l))
                    data = dict(zip(fields, values))
                                        
                    shifts.append(self._create_shift(data))
                    
        return tuple(shifts)
    
    def _scroll(self, iterator, field):
        
        for line in iterator:
            if line.lstrip().startswith(field):
                break
            
    def _create_shift(self, data):
        
        try:
            position = int(data[ChemShiftReader.RANK])
            residue = data[ChemShiftReader.RESIDUE]
            name = data[ChemShiftReader.ATOM]
            element = data[ChemShiftReader.ELEMENT]
            shift = float(data[ChemShiftReader.SHIFT])
            
        except KeyError as ke:
            raise ChemShiftFormatError("Required field {0} not found".format(str(ke)))
        except ValueError as ve:
            raise ChemShiftFormatError("Can't parse value: {0}".format(str(ve)))
        
        return ChemShiftInfo(position, residue, name, element, shift)


class ChemShift3Reader(ChemShiftReader):
    """
    Simple NMR STAR v3 chemical shift reader. 
    
    @note: This is not a full-fledged, semantic NMR STAR parser. It handles
           only the chemical shift table. 
    """    
    
    RANK3 = '_Atom_chem_shift.Seq_ID'
    RESIDUE3 = '_Atom_chem_shift.Comp_ID'
    ATOM3 = '_Atom_chem_shift.Atom_ID'
    ELEMENT3 = '_Atom_chem_shift.Atom_type'
    SHIFT3 = '_Atom_chem_shift.Val'
    
    def _create_shift(self, data):

        try:        
            position = data[ChemShift3Reader.RANK3]
            residue = data[ChemShift3Reader.RESIDUE3]
            name = data[ChemShift3Reader.ATOM3]
            element = data[ChemShift3Reader.ELEMENT3]
            shift = data[ChemShift3Reader.SHIFT3]
            
        except KeyError as ke:
            raise ChemShiftFormatError("Required field {0} not found".format(str(ke)))
        except ValueError as ve:
            raise ChemShiftFormatError("Can't parse value: {0}".format(str(ve)))
                
        return ChemShiftInfo(position, residue, name, element, shift)
