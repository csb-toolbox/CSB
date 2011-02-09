"""
DSSP Parser
"""

import csb.pyutils

from csb.bio.structure import SecStructures, UnknownSecStructureError


class DSSPParseError(ValueError):
    pass


class ResidueAssignmentInfo(object):
    
    def __init__(self, residue_id, accession, chain, secondary_structure, phi, psi):
        
        self.residue_id = residue_id
        self.accession = accession
        self.chain = chain
        self.secondary_structure = secondary_structure
        self.phi = phi
        self.psi = psi
        

class DSSPParser(object):
    """
    Simple DSSP Secondary Structure Parser.
    """
    
    def parse(self, dssp_file):
        """
        @param dssp_file: source DSSP file to parse
        @type dssp_file: str
        @return: a dictionary of L{ResidueAssignmentInfo} objects
        @rtype: dict
        """
        
        data = {}
        start = False
        offset = 0                  # assume old DSSP format
        
        for line in open(dssp_file):
            
            if not start:
                
                if line.startswith('HEADER'):
                    accession = line[62:66].strip().lower()
                     
                elif line.startswith('  #  RESIDUE'):
                    if len(line) >= 140:
                        offset = 4  # the new DSSP format           
                    start = True
            else:
                if line[13] == '!':
                    continue                
                
                residue_id = line[6:11].strip()
                chain = line[11]         
                try:
                    ss = line[16].strip()
                    if ss == '':
                        ss = SecStructures.Gap
                    else:
                        ss = csb.pyutils.Enum.parse(SecStructures, ss)  
                except csb.pyutils.EnumValueError as e:
                    raise UnknownSecStructureError(str(e)) 
                phi = float(line[104 + offset : 109 + offset])
                psi = float(line[110 + offset : 115 + offset])
                
                
                if chain not in data:
                    data[chain] = {}
                
                data[chain][residue_id] = ResidueAssignmentInfo(residue_id, accession, chain, ss, phi, psi)
                
        return data

class StrideParser(object):
    """
    Simple STRIDE Secondary Structure Parser.
    """
    
    def parse(self, stride_file):
        """
        @param stride_file: source STRIDE file to parse
        @type stride_file: str
        @return: a dictionary of L{ResidueAssignmentInfo} objects
        @rtype: dict
        """
        
        data = {}
        
        for line in open(stride_file):
            if line.startswith('ASG '):
                
                fields = line.split()
                
                residue_id = fields[3]
                chain = fields[2]
                accession = fields[-1].lower()
                try:
                    ss = csb.pyutils.Enum.parse(SecStructures, fields[5])  
                except csb.pyutils.EnumValueError as e:
                    raise UnknownSecStructureError(str(e)) 
                phi = float(fields[7])
                psi = float(fields[8])
                
                if chain not in data:
                    data[chain] = {}
                
                data[chain][residue_id] = ResidueAssignmentInfo(residue_id, accession, chain, ss, phi, psi)
                
        return data

    
def get(accession, prefix='http://www.pdb.org/pdb/files/'):
    """
    Download and parse a DSSP entry.

    @param accession: accession number of the entry
    @type accession: str
    @param prefix: download URL prefix
    @type prefix: str

    @return: see L{DSSPParser.parse}
    @rtype: dict
    """
    import urllib2
    from tempfile import NamedTemporaryFile

    dssp = NamedTemporaryFile()

    browser = urllib2.urlopen(prefix + accession.lower() + '.dssp')
    dssp.write(browser.read())
    dssp.flush()

    return DSSPParser().parse(dssp.name)

