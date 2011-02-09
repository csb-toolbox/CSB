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

    def __init__(self, checks = None):
        self.binary = '/home/mechelke/install/whatcheck/DO_WHATCHECK.COM'
        if checks is None:
            self.checks  = ['NAMCHK', 'HNDCHK', 'WGTCHK',
                            'MISCHK', 'MO2CHK', 'BNDCHK',
                            'ANGCHK', 'PLNCHK', 'PL2CHK',
                            'PUCCHK', 'PC2CHK', 'CHICHK',
                            'RAMCHK', 'C12CHK', 'INOCHK',
                            'BMPCHK', 'NQACHK', 'FLPCHK',
                            'ROTCHK', 'BBCCHK', 'HNQCHK',
                            'BH2CHK', 'BA2CHK', 'SCOLST',
                            'BVALST', 'ACCLST', 'PDBLST',
                            'AXACHK', 'H2OCHK', 'H2OCHK',
                            'H2OHBO']
        else:
            self.checks = checks

    
    def parse_summary(self, fn):
        """
        @param whatif_file: source whatif file to parse
        @type whatif_file: str

        @return a dictinary
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

    def run(self, pdb_file):
        wd = os.getcwd()
        tmp = tempfile.mkdtemp()
        base = os.path.basename(pdb_file)
        
        p = subprocess.call('cp %s %s/.' % (os.path.expanduser(pdbfile), tmp), shell = True)
        os.chdir(tmp)

        p = subprocess.Popen('%s %s' %(self.binary, os.path.join(tmp,base)),
                             stdout=open("/dev/null", "w"),
                             stderr=open("/dev/null", "w"),
                             shell = True)
        
        sts = os.waitpid(p.pid, 0)[1]
        
        out = self.parse_summary(os.path.join(tmp,'pdbout.txt'))
        
        os.chdir(wd)
        p = subprocess.call('rm -rf %s' %tmp, shell = True)
        return out


    def parse_checkdb(self, fn):
        """
        @param whatif_file: source whatif file to parse
        @type whatif_file: str

        """

        x = open(os.path.expanduser(fn)).read().split('CheckID')[1:]

        d = {}
        
        
        for xx in x:
            dd = {}

            check, ids, data = self.parse_check([s.strip() for s in xx.split('\n')])

            for id, datum in zip(ids,data):

                model, chain, residue = self._parse_id(id)
                if not model in dd:
                    dd[model] = {}

                if not chain in dd[model]:
                    dd[model][chain] = {}

                try:
                    dd[model][chain][residue] = float(data)
                    
                except ValueError:
                    dd[model][chain][residue] = data
                    

    def parse_check(self, entry):

        id = entry[0][-6:]

        values = {}

        if not id in CHECKIDS:
            print id
            return id, values

        i = grep(entry, 'Type')

        t = entry[i][-5:].strip().lower()

        if t == 'text':
            t = str
        else:
            t = eval(t)

        entry = entry[i:]

        i = grep(entry, 'Name')

        if i is None:
            return id, values

        entry = entry[i:]

        while len(entry) > 1:

            name = entry[0].split(':')[1].strip()
            value = t(entry[1].split(':')[1].strip())
            values[name] = value

            i = grep(entry[1:], 'Name')

            if i is None: break

            entry = entry[i+1:]

        return id, values

        
                    
    def parse_id(self, id):
        a,b,c = id.split(';')[:3]

        return int(a), b.strip(), int(c)
        
     

        

        


if __name__ == "__main__":

    filename =  '~/data/whatif/pdbout.txt'
    wp = WhatCheckParser()
    inf = wp.parse_summary(filename)

    pdbfile =  '~/data/whatif/2JZC.pdb'
    pdbfile =  '~/2kd7.pdb'
    inf2 = wp.run(pdbfile)




    
