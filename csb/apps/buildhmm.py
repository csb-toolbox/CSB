"""
Build an HMM from a FASTA sequence. This program is a proxy to buildali.pl
and hhmake from the HHpred package.

@note: assuming you have the full HHpred package installed and configured.
"""


import os
import abc

import csb.apps
import csb.io
import csb.bio.io
import csb.core


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
        cmd.add_scalar_option('tk_root', 't', str, 'path to the ToolkitRoot folder in your HHpred setup', default='/ebio/abt1_toolkit/share/wye')
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
            pb = ProfileBuilder.create(self.args.query, self.args.query_id, self.args.tk_root,
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
    def create(query, target_id, tk_root, pseudo=True, ss=True, cpu=1):

        if not os.path.isfile(query):
            raise BuildIOError('File not found: ' + query)
                
        for line in open(query):
            
            if not line.strip():
                continue
            
            if line.startswith('>'):
                return FASTAProfileBuilder(query, target_id, tk_root, pseudo, ss, cpu)
            elif line.startswith('HEADER') or line.startswith('ATOM'): 
                return PDBProfileBuilder(query, target_id, tk_root, pseudo, ss, cpu)
            else:
                raise BuildArgError('Unknown input file format')
                
    def __init__(self, query, target_id, tk_root, pseudo=True, ss=True, cpu=1):
        
        self.tk_root = tk_root
        if 'TK_ROOT' not in os.environ or not os.environ['TK_ROOT']:
            os.putenv('TK_ROOT', tk_root)
                    
        self.query = query
        self.accession = target_id[:-1]
        self.chain = target_id[-1]
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
        
        program = os.path.join(self.tk_root, 'bioprogs', 'hhpred', 'buildali.pl')
        
        if not self.ss:
            noss = '-noss'
        else:
            noss = ''        
        cmd = 'perl {0} {1} -cpu {2} {3}'.format(program, noss, self.cpu, self._input)
        bali = csb.io.Shell.run(cmd)
        
        ali = self.target_id + '.a3m'
        if bali.code != 0:
            raise csb.io.ProcessError(bali)
        if not os.path.isfile(ali):
            raise NoOutputError(ali, bali)        
        
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
            s = csb.bio.io.StructureParser(self.query).parse()            
            chain = s.chains[self.chain]
        except csb.core.ItemNotFoundError:
            raise BuildArgError('Chain {0.chain} not found in {0.query}'.format(self))
        except IOError as ioe:
            raise BuildIOError(str(ioe))
        
        fasta = self.target_id + '.fa'
        
        with csb.io.EntryWriter(fasta) as f:
            f.writeline(chain.header)
            f.writeline(chain.sequence)
            
        self._input = fasta
        return fasta
        
    def make_hmm(self):
        
        super(PDBProfileBuilder, self).make_hmm()
        self.format_structure()                    
        
    def format_structure(self):
        assert self._hhm is not None
        
        pdb = self.target_id + '.pdb'
                
        parser = csb.bio.io.HHProfileParser(self._hhm)
        parser.format_structure(self.query, self.chain, pdb)
        
        self._pdb = pdb
        return pdb             



if __name__ == '__main__':
    
    AppRunner().run()    
