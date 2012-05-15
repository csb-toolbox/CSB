"""
Procheck parser
"""
import os
import re
import shutil

from csb.pyutils import Shell
from csb.io import TempFolder

class ProcheckParser():
    """
    Simple Prochceck Summary parser
    """
    def __init__(self):
        self.binary = 'procheck.scr'
        self.acc = 2.0
        
    def parse(self, fn):
        """
        @param fn: source  file to parse
        @type fn: str

        @return: dicttionary of parsed quality indicatiors
        """
        info = dict()
        
        f_handler = open(os.path.expanduser(fn))
        text = f_handler.read()
        
        input_file_name = re.compile('>>>-----.*?\n.*?\n\s*\|\s*(\S+)\s+')
        residues = re.compile('(\d+)\s*residues\s\|')
        ramachandran_plot = re.compile('Ramachandran\splot:\s*(\d+\.\d+)' + 
                                      '%\s*core\s*(\d+\.\d+)%\s*allow\s*(\d+\.\d+)' + 
                                      '%\s*gener\s*(\d+\.\d+)%\s*disall')
        labelled_all = re.compile('Ramachandrans:\s*(\d+)\s*.*?out\sof\s*(\d+)')
        labelled_chi = re.compile('Chi1-chi2\splots:\s*(\d+)\s*.*?out\sof\s*(\d+)')
        bad_contacts = re.compile('Bad\scontacts:\s*(\d+)')
        g_factors = re.compile('G-factors\s*Dihedrals:\s*([0-9-+.]+)' + 
                              '\s*Covalent:\s*([0-9-+.]+)\s*Overall:\s*([0-9-+.]+)')

        info['input_file'] = input_file_name.search(text).groups()[0]
        info['#residues'] = int(residues.search(text).groups()[0])
        info['rama_core'], info['rama_allow'], info['rama_gener'], info['rama_disall'] = \
                           [float(g) for g in ramachandran_plot.search(text).groups()]
        info['g_dihedrals'], info['g_bond'], info['g_overall'] = \
                             [float(g) for g in g_factors.search(text).groups()]
        info['badContacts'] = int(bad_contacts.search(text).groups()[0])
        info['labelledAll'] = float(labelled_all.search(text).groups()[0]) / \
                              float(labelled_all.search(text).groups()[1])
        info['labelledChi'] = float(labelled_chi.search(text).groups()[0]) / \
                              float(labelled_chi.search(text).groups()[0])

        f_handler.close()
        
        return info


    def run(self, pdb_file):
        """
        Runs procheck for the given pdbfile and parses the output.
        Will fail if the procheck binary is not in the path.
        
        @param pdb_file: file to parse
        @return: dict of parsed values
        """
        wd = os.getcwd()
        base = os.path.basename(pdb_file)

        with TempFolder() as tmp:
            shutil.copy(os.path.expanduser(pdb_file), tmp.name)
            os.chdir(tmp.name)
            Shell.run('{0} {1} {2}'.format(self.binary,
                                           os.path.join(tmp.name, base),
                                           self.acc))
            summary = '.'.join([os.path.splitext(base)[0], 'sum'])
            out = self.parse(os.path.join(tmp.name, summary))
            os.chdir(wd)

        return out
