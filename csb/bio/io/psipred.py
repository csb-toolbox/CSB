"""
PSIPRED Parser
"""

import csb.core

from csb.bio.structure import SecondaryStructure, SecStructures, UnknownSecStructureError


class PSIPredParseError(ValueError):
    pass


class PSIPredResidueInfo(object):
    
    def __init__(self, rank, residue, sec_structure, helix, strand, coil):
        
        self.rank = rank
        self.residue = residue
        self.sec_structure = sec_structure
        self.helix = helix
        self.coil = coil
        self.strand = strand    


class PSIPredParser(object):
    """
    Simple PSI-PRED Secondary Structure Parser.
    """
    
    def parse(self, psipred_file):
        """
        @param psipred_file: source PSI-PRED *.horiz file to parse
        @type psipred_file: str
        @rtype: L{SecondaryStructure}
        """
        
        ss = []
        conf = []
        
        for line in open(psipred_file):
            
            if line.startswith('Conf:'):
                conf.extend(line[6:].strip())
                
            elif line.startswith('Pred:'):
                ss.append(line[6:].strip())
        
        ss = ''.join(ss)
        conf = ''.join(conf)
        
        if len(ss) != len(conf):
            raise PSIPredParseError('Invalid PSI-PRED output file')
        
        if ss:
            return SecondaryStructure(ss, conf)
        else:
            return SecondaryStructure(None)

    def parse_scores(self, scores_file):
        """
        @param scores_file: source PSI-PRED *.ss2 file to parse
        @type scores_file: str
        @rtype: list of L{PSIPredResidueInfo}
        """
        residues = [] 
        
        for line in open(scores_file):
            
            if line.startswith('#') or not line.strip():
                continue
            else:
                line = line.split()         

                rank = int(line[0])
                residue = line[1]
                                
                try:
                    ss = csb.core.Enum.parse(SecStructures, line[2])  
                except csb.core.EnumValueError as e:
                    raise UnknownSecStructureError(str(e))
                
                coil, helix, strand = map(float, line[3:6])
                
                residues.append(PSIPredResidueInfo(rank, residue, ss, helix, strand, coil))
        
        return tuple(residues)
                