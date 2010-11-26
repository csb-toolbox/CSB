"""
PSIPRED Parser
"""

from csb.bio.structure import SecondaryStructure


class PSIPredParseError(ValueError):
    pass


class PSIPredParser(object):
    """
    Simple PSI-PRED Secondary Structure Parser.
    """
    
    def parse(self, psipred_file):
        """
        @param psipred_file: source PSI-PRED file to parse
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
