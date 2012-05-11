import csb.pyutils


class ShiftInfo(object):

    def __init__(self, residue_id, amino_acid, nucleus,
                 shift, element, secondary_structure):

        self.residue_id = residue_id
        self.nucleus = nucleus
        self.element = element
        self.amino_acid = amino_acid
        self.shift = shift
        self.secondary_structure = secondary_structure

    def __str__(self):
        return '{0.amino_acid} {0.nucleus} {0.shift}'.format(self)

    __repr__ = __str__

        
class ChemicalShiftContainer(csb.pyutils.DictionaryContainer):

    def __init__(self, bmrb_id='', pdb_id='', sequence='',
                 chain='', exptype=''):
        
        self.bmrb_id = bmrb_id
        self.pdb_id = pdb_id
        self.sequence = sequence
        self.chain = chain
        self.exptype = exptype
        
        super(ChemicalShiftContainer, self).__init__()

class VascoStructureParser(object):
    """
    Simple Vasco Parser
    """

    def __init__(self):
        self._stream  = None

    def parse(self, file_name, ignore_outliers=True):
        """
        @param file_name: source  file to parse
        @type file_name: str
        @return: a L{ChemicalShiftContainer} of L{ShiftInfo} objects
        @rtype: dict
        """
        self._stream  = open(file_name)
        shifts = self._parse_header()

        self._parse_shifts(shifts, ignore_outliers=ignore_outliers)
        self._stream.close()

        return shifts
        
    def _parse_header(self):

        bmrb_id = ''
        pdb_id = ''
        sequence = ''
        chain = ''
        exptype = ''
        self._stream.seek(0)
        
        while True:
            try:
                line = next(self._stream)
            except StopIteration :
                break

            if line.startswith('#'):
                if line[2:].startswith('BMRB ORIGIN'):
                    bmrb_id = line[20:].strip()
                elif line[2:].startswith('PDB ORIGIN'):
                    pdb_id = line[20:].strip()
                elif line[2:].startswith('SEQUENCE PDB'):
                    sequence = line[20:].strip()
                    chain = line[17]
                elif line[2:].startswith('PDB EXPTYPE'):
                    exptype = line[20:].strip()
            else:
                break

         
        return ChemicalShiftContainer(bmrb_id, pdb_id, chain,
                                      sequence, exptype )
            
    
    def _parse_shifts(self, data, ignore_outliers=True):

        while True:
            try:
                line = next(self._stream)
            except StopIteration:
                break

            if ignore_outliers and "Shift outlier" in line:
                continue
      
            chain_id = line[7]
            res_code = line[9:14].strip()
            res_label = line[16:19].strip()
            res_ss = line[21]
            nucleus_name = line[23:28].strip()
            nucleus_element = line[41]
            shift = float(line[43:52])

            info = ShiftInfo(res_code, res_label,
                             nucleus_name, shift,
                             nucleus_element, res_ss)

            if not chain_id in data:
                data.append(chain_id, csb.pyutils.OrderedDict())

            if not res_code in data[chain_id]:
                data[chain_id][res_code] = {}

            
            data[chain_id][res_code][nucleus_name] = info

    