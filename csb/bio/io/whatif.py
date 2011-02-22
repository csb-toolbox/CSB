"""
Simple WhatIf/WhatCheck Summary parser
"""

import re
import os
import tempfile
import subprocess

def grep(lines, what):
    for i in range(len(lines)):
        if lines[i].find(what) <> -1:
            return i


class WhatCheckSummaryInfo(dict):
    pass


class WhatCheckParser(object):
    """
    Simple WhatIf/WhatCheck Summary parser
    """

    def __init__(self, binary = '/home/mechelke/install/whatcheck/DO_WHATCHECK.COM'):
        self.binary = binary
    
    def parse_summary(self, fn):
        """
        @param fn: source whatif file to parse
        @type fn: str

        @return: a dictinary
        """
        f_handler = open(os.path.expanduser(fn))
        text = f_handler.read()

        info = WhatCheckSummaryInfo()
        reRamachandran = re.compile('Ramachandran\s*Z-score\s*:\s*([0-9.Ee-]+)')
        re1st = re.compile('1st\s*generation\s*packing\s*quality\s*:\s*([0-9.Ee-]+)')
        re2nd = re.compile('2nd\s*generation\s*packing\s*quality\s*:\s*([0-9.Ee-]+)')
        reBB = re.compile('Backbone\s*conformation\s*Z-score\s*:\s*([0-9.Ee-]+)')
        reRotamer = re.compile('chi-1\S*chi-2\s*rotamer\s*normality\s*:\s*([0-9.Ee-]+)')
        

        info['rama_z_score'] = float(reRamachandran.search(text).groups(0)[0])
        info['bb_z_score'] = float(reBB.search(text).groups(0)[0])
        info['1st_packing_z_score'] = float(re1st.search(text).groups(0)[0])
        info['2nd_packing_z_score'] = float(re2nd.search(text).groups(0)[0])
        info['rotamer_score'] = float(reRotamer.search(text).groups(0)[0])

        f_handler.close()
        return info

    parse = parse_summary

    def run(self, pdb_file):
        wd = os.getcwd()
        tmp = tempfile.mkdtemp()
        base = os.path.basename(pdb_file)
        
        p = subprocess.call('cp %s %s/.' %(os.path.expanduser(pdb_file),
                                           tmp), shell = True)
        os.chdir(tmp)

        p = subprocess.Popen('%s %s' %(self.binary, os.path.join(tmp, base)),
                             stdout=open("/dev/null", "w"),
                             stderr=open("/dev/null", "w"),
                             shell = True)
        
        sts = os.waitpid(p.pid, 0)[1]
        
        out = self.parse_summary(os.path.join(tmp,'pdbout.txt'))
        
        os.chdir(wd)
        p = subprocess.call('rm -rf %s' %tmp, shell = True)
        return out
                    
        

        


if __name__ == "__main__":

    filename =  '~/data/whatif/pdbout.txt'
    wp = WhatCheckParser()
    inf = wp.parse_summary(filename)

    pdbfile =  '~/data/whatif/2JZC.pdb'
#    pdbfile =  '~/2kd7.pdb'
    inf2 = wp.run(pdbfile)




    
