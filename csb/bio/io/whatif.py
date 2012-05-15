"""
Simple WhatIf/WhatCheck Summary parser
"""

import re
import os
import shutil

from csb.pyutils import Shell
from csb.io import TempFolder

class WhatCheckParser(object):
    """
    Simple WhatIf/WhatCheck Summary parser
    """

    def __init__(self, binary='DO_WHATCHECK.COM'):
        self.binary = binary
    
    def parse_summary(self, fn):
        """
        @param fn: whatif pdbout.txt file to parse
        @type fn: str

        @return: A dict containing some of the WhatCheck results
        @rtype: a dict
        """
        f_handler = open(os.path.expanduser(fn))
        text = f_handler.read()

        info = dict()
        re_ramachandran = re.compile('Ramachandran\s*Z-score\s*:\s*([0-9.Ee-]+)')
        re_1st = re.compile('1st\s*generation\s*packing\s*quality\s*:\s*([0-9.Ee-]+)')
        re_2nd = re.compile('2nd\s*generation\s*packing\s*quality\s*:\s*([0-9.Ee-]+)')
        re_backbone = re.compile('Backbone\s*conformation\s*Z-score\s*:\s*([0-9.Ee-]+)')
        re_rotamer = re.compile('chi-1\S*chi-2\s*rotamer\s*normality\s*:\s*([0-9.Ee-]+)')
        

        info['rama_z_score'] = float(re_ramachandran.search(text).groups(0)[0])
        info['bb_z_score'] = float(re_backbone.search(text).groups(0)[0])
        info['1st_packing_z_score'] = float(re_1st.search(text).groups(0)[0])
        info['2nd_packing_z_score'] = float(re_2nd.search(text).groups(0)[0])
        info['rotamer_score'] = float(re_rotamer.search(text).groups(0)[0])

        f_handler.close()
        return info

    parse = parse_summary


    def run(self, pdb_file):
        """
        Runs WhatCheck for the given pdbfile and parses the output.
        Will fail if the WhatCheck binary is not in the path.
        
        @param pdb_file: file to parse
        @return: dict of parsed values
        """
        wd = os.getcwd()
        base = os.path.basename(pdb_file)

        with TempFolder() as tmp:
            shutil.copy(os.path.expanduser(pdb_file), tmp.name)
            os.chdir(tmp.name)
            Shell.run('{0} {1}'.format(self.binary,
                                       os.path.join(tmp.name, base)))
            out = self.parse_summary(os.path.join(tmp.name, 'pdbout.txt'))
            os.chdir(wd)

        return out
                    
        

        
