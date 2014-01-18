"""
Build an HMM from a FASTA sequence. This program is a proxy to hhblits/addss.pl
and hhmake from the HHpred package.

@note: assuming you have the full HHpred package installed and configured.
"""


import os
import abc

import csb.apps
import csb.core
import csb.io

from csb.bio.io.wwpdb import StructureParser
from csb.bio.io.hhpred import HHProfileParser
from csb.bio.io.fasta import FASTAOutputBuilder
from csb.bio.sequence import ChainSequence



class ExitCodes(csb.apps.ExitCodes):
    
    IO_ERROR = 2
    INVALID_DATA = 3
    EXT_TOOL_FAILURE = 4

class AppRunner(csb.apps.AppRunner):
    
    @property
    def target(self):
        return BuildProfileApp
    
    def command_line(self):
        
        cmd = csb.apps.ArgHandler(self.program, __doc__)
        
        cmd.add_scalar_option('query-id', 'q', str, 'ID of the query, in PDB-like format (accessionCHAIN).'
                              'Used for naming the output files. Also, if the input is a PDB file with '
                              'multiple chains, CHAIN is used to pull the required chain from the file.',
                              required=True)
        cmd.add_scalar_option('tk-root', 't', str, 'path to the ToolkitRoot folder in your HHsuite setup', default='/ebio/abt1_toolkit/share/wye')
        cmd.add_scalar_option('database', 'd', str, 'custom HHblits database; if not defined, toolkit\'s unirpto20 will be used', default=None)        
        cmd.add_scalar_option('tk-config', 'c', str, 'path to a folder containing custom HHsuite configs (e.g. HHPaths.pm)', default='.')
        cmd.add_scalar_option('cpu', None, int, 'maximum degree of parallelism', default=1)        

        cmd.add_boolean_option('no-ss', None, 'do not include secondary structure', default=False)        
        cmd.add_boolean_option('no-pseudo', None, 'do not add emission and transition pseudocounts', default=False)
        cmd.add_boolean_option('no-calibration', None, 'do not calibrate the profile', default=False)

        cmd.add_positional_argument('query', str, 'input sequence (FASTA or PDB file)')        
                
        return cmd
    

class BuildProfileApp(csb.apps.Application):
    
    def main(self):            
        
        if os.path.isfile(self.args.query_id + '.hhm'):
            BuildProfileApp.exit('# Profile "{0}" already exists, skipping'.format(self.args.query_id),
                                 ExitCodes.CLEAN)  
        
        try:
            self.log('# Building profile HMM for {0}...'.format(self.args.query))
            pb = ProfileBuilder.create(self.args.query, self.args.query_id, self.args.database, self.args.tk_root, self.args.tk_config,
                                       pseudo=not self.args.no_pseudo, ss=not self.args.no_ss, cpu=self.args.cpu)
            
            pb.build_alignment()
            pb.make_hmm()
            
            if not self.args.no_calibration:                
                pb.calibrate_hmm()

        except BuildArgError as ae:
            BuildProfileApp.exit(str(ae), ExitCodes.INVALID_DATA)
                        
        except BuildIOError as ioe:
            BuildProfileApp.exit(str(ioe), ExitCodes.IO_ERROR)
            
        except csb.io.InvalidCommandError as ose:
            msg = '{0!s}: {0.cmd}'.format(ose)
            BuildProfileApp.exit(msg, ExitCodes.IO_ERROR)            

        except NoOutputError as noe:
            msg = 'Expected file {0} not produced by: {1.cmd}.\nSTDERR: {1.stderr}\nSTDOUT: {1.stdout}'.format(noe.expected, noe.context)
            BuildProfileApp.exit(msg, ExitCodes.EXT_TOOL_FAILURE)  
            
        except csb.io.ProcessError as pe:
            msg = 'Bad exit code #{0.code} from: {0.cmd}.\nSTDERR: {0.stderr}\nSTDOUT: {0.stdout}'.format(pe.context)
            BuildProfileApp.exit(msg, ExitCodes.EXT_TOOL_FAILURE)              

        self.log('  successfully created profile "{0}"'.format(self.args.query_id))            


class BuildError(Exception):
    pass

class BuildIOError(BuildError):
    pass

class BuildArgError(BuildError):
    pass

class NoOutputError(BuildError):

    def __init__(self, expected, context, *args):
        
        self.expected = expected
        self.context = context
        super(NoOutputError, self).__init__(*args)


class ProfileBuilder(object):
    
    __metaclass__ = abc.ABCMeta
    
    EMISSION_PSEUDO = '-pcm 4 -pca 2.5 -pcb 0.5 -pcc 1.0'
    TRANSITION_PSEUDO = '-gapb 1.0 -gapd 0.15 -gape 1.0 -gapf 0.6 -gapg 0.6 -gapi 0.6'
    
    @staticmethod
    def create(query, target_id, database, tk_root, tk_config, pseudo=True, ss=True, cpu=1):
        
        if database is None:
            database = os.path.join(tk_root, "databases", "hhblits", "uniprot20")        

        if not os.path.isfile(query):
            raise BuildIOError('File not found: ' + query)
                
        for line in open(query):
            
            if not line.strip():
                continue
            
            if line.startswith('>'):
                return FASTAProfileBuilder(query, target_id, database, tk_root, tk_config, pseudo, ss, cpu)
            elif line.startswith('HEADER') or line.startswith('ATOM'): 
                return PDBProfileBuilder(query, target_id, database, tk_root, tk_config, pseudo, ss, cpu)
            else:
                raise BuildArgError('Unknown input file format')
                
    def __init__(self, query, target_id, database, tk_root, tk_config, pseudo=True, ss=True, cpu=1):
        
        self.tk_root = tk_root
        self.tk_config = tk_config
        self.hhlib = os.path.join(tk_root, "bioprogs", "hhsuite")
        
        if 'TK_ROOT' not in os.environ or not os.environ['TK_ROOT']:
            os.putenv('TK_ROOT', self.tk_root)
        if 'HHLIB' not in os.environ or not os.environ['HHLIB']:
            os.putenv('HHLIB', self.hhlib)
        os.environ["PATH"] += os.pathsep + os.path.join(self.hhlib, "bin")
                                
        self.query = query
        self.accession = target_id[:-1]
        self.chain = target_id[-1]
        self.database = database
        self.pseudo = bool(pseudo)
        self.ss = bool(ss)
        self.cpu = cpu
        
        self._input = None
        self._a3m = None
        self._hhm = None
        
        self.configure_input()
            
    def run(self):

        self.build_alignment()
        self.make_hmm()
        self.calibrate_hmm()        
    
    @property
    def target_id(self):
        return self.accession + self.chain
    
    @abc.abstractmethod
    def configure_input(self):
        pass
    
    def build_alignment(self):
        assert self._input is not None
        
        program = os.path.join(self.tk_root, 'bioprogs', 'hhsuite', 'bin', 'hhblits')
            
        ali = self.target_id + '.a3m'                    
        cmd = '{0} -cpu {1} -i {2} -d {3} -nodiff -oa3m {4}'.format(
                                program, self.cpu, self._input, self.database, ali)
        bali = csb.io.Shell.run(cmd)
        
        if bali.code != 0:
            raise csb.io.ProcessError(bali)
        if not os.path.isfile(ali):
            raise NoOutputError(ali, bali)
        
        if self.ss:
            program2 = os.path.join(self.tk_root, 'bioprogs', 'hhsuite', 'scripts', 'addss.pl')
            
            with csb.io.TempFile() as patch:
                for l in open(program2):
                    if l.lstrip().startswith("use HHPaths"):
                        patch.write('use lib "{0}";\n'.format(self.tk_config))
                    patch.write(l);
                patch.flush()
                
                cmd2 = "perl {0} {1}".format(patch.name, ali)
                addss = csb.io.Shell.run(cmd2)            
                if addss.code != 0:
                    raise csb.io.ProcessError(addss)       
        
        self._ali = ali
        return ali
        
    def make_hmm(self):
        assert self._ali is not None
        
        program = os.path.join(self.tk_root, 'bioprogs', 'hhpred', 'hhmake')
        hhm = self.target_id + '.hhm'
        cmd = '{0} -i {1} -o {2}'.format(program, self._ali, hhm)
        
        if self.pseudo:
            cmd = '{0} {1} {2}'.format(cmd, ProfileBuilder.EMISSION_PSEUDO, ProfileBuilder.TRANSITION_PSEUDO)
                    
        nnmake = csb.io.Shell.run(cmd)
        if nnmake.code != 0:
            raise csb.io.ProcessError(nnmake)
        if not os.path.isfile(hhm):
            raise NoOutputError(hhm, nnmake)
        
        self._hhm = hhm
        return hhm
    
    def calibrate_hmm(self):
        assert self._hhm is not None
        
        program = os.path.join(self.tk_root, 'bioprogs', 'hhpred', 'hhsearch')
        caldb = os.path.join(self.tk_root, 'databases', 'hhpred', 'cal.hhm')
          
        cmd = '{0} -i {1}.hhm -d {2} -cal -cpu {3}'.format(program, self.target_id, caldb, self.cpu)
        csb.io.Shell.runstrict(cmd)


class FASTAProfileBuilder(ProfileBuilder):
    
    def configure_input(self):
        
        if not os.path.isfile(self.query):
            raise BuildIOError('File not found: ' + self.query)

        fasta = self.target_id + '.fa'
        
        with csb.io.EntryWriter(fasta) as f:
            with open(self.query) as q:
                f.write(q.read())
                    
        self._input = fasta
        return fasta


class PDBProfileBuilder(ProfileBuilder):

    def configure_input(self):
        
        try:
            s = StructureParser(self.query).parse()            
            chain = s.chains[self.chain]
        except csb.core.ItemNotFoundError:
            raise BuildArgError('Chain {0.chain} not found in {0.query}'.format(self))
        except IOError as ioe:
            raise BuildIOError(str(ioe))
        
        fasta = self.target_id + '.fa'
        
        with open(fasta, 'w') as f:
            sequence = ChainSequence.create(chain)
            FASTAOutputBuilder(f).add_sequence(sequence)
            
        self._input = fasta
        return fasta
        
    def make_hmm(self):
        
        super(PDBProfileBuilder, self).make_hmm()
        self.format_structure()                    
        
    def format_structure(self):
        assert self._hhm is not None
        
        pdb = self.target_id + '.pdb'
                
        parser = HHProfileParser(self._hhm)
        parser.format_structure(self.query, self.chain, pdb)
        
        self._pdb = pdb
        return pdb             


def main():
    AppRunner().run()
    
    
if __name__ == '__main__':
    main()