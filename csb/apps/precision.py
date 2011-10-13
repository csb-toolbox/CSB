"""
Measure the precision and coverage of a fragment library stored in Rosetta
NNmake format.
"""

import os
import multiprocessing

import csb.apps

import csb.bio.io.wwpdb  as wwpdb
import csb.bio.structure as structure
import csb.bio.fragments
import csb.bio.fragments.rosetta as rosetta
import csb.pyutils

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot


class ExitCodes(csb.apps.ExitCodes):
    
    IO_ERROR = 2
    INVALID_DATA = 3

class AppRunner(csb.apps.AppRunner):
    
    @property
    def target(self):
        return PrecisionApp
     
    def command_line(self):
        
        cpu = multiprocessing.cpu_count()
        cmd = csb.apps.ArgHandler(self.program, __doc__)

        cmd.add_scalar_option('pdb', 'p', str, 'the PDB database (a directory containing all PDB files)', required=True)
        cmd.add_scalar_option('native', 'n', str, 'native structure of the target (PDB file)', required=True)
        cmd.add_scalar_option('chain', 'c', str, 'chain identifier (if not specified, the first chain)', default=None)

        cmd.add_scalar_option('top', 't', int, 'read top N fragments per position', default=25)
        cmd.add_scalar_option('cpu', 'C', int, 'maximum degree of parallelism', default=cpu)
        cmd.add_scalar_option('rmsd', 'r', float, 'RMSD cutoff for precision and coverage', default=1.5)         
        cmd.add_scalar_option('output', 'o', str, 'output directory', default='.')
        
        cmd.add_positional_argument('library', str, 'Fragment library file in Rosetta NNmake format')
                        
        return cmd

class PrecisionApp(csb.apps.Application):
    
    def main(self):
        
        if not os.path.isdir(self.args.output):
            PrecisionApp.exit('Output directory does not exist', code=ExitCodes.INVALID_DATA, usage=True)
        
        for file in [self.args.native, self.args.library]:
            if not os.path.isfile(file):
                PrecisionApp.exit('File not found: {0}'.format(file), code=ExitCodes.INVALID_DATA, usage=True)            
        
        self.log('\nLoading {0.library} (top {0.top} per position)... '.format(self.args), ending='')
        
        try:
            library = rosetta.RosettaFragmentMap.read(self.args.library, top=self.args.top)
            self.log('{0} fragments'.format(len(library)))
            
            try:
                native = wwpdb.RegularStructureParser(self.args.native).parse_structure()
                if self.args.chain:
                    native = native.chains[self.args.chain]
                else:
                    native = native.first_chain
            except (structure.Broken3DStructureError, wwpdb.PDBParseError) as se:
                raise ArgumentError("Can't parse native structure: " + str(se))
            
            except (csb.pyutils.ItemNotFoundError) as ce:
                raise ArgumentError("Chain {0!s} not found in native structure".format(ce))            

            self.log('Superimposing fragments...\n')
            lsi = LibrarySuperimposer(native, library, self.args.pdb, self.args.output, self.args.rmsd)
            matches = lsi.run(self.args.cpu)
            
            info = lsi.plot(matches)
            self.log(' RMSD Cutoff    {0:>6.2f} A'.format(self.args.rmsd))            
            self.log(' AVG Precision  {0.precision:>6.2f} %'.format(info))
            self.log(' Coverage       {0.coverage:>6.2f} %'.format(info))
            self.log(' RW  Precision  {0.figure}'.format(lsi))


            
        except ArgumentIOError as ioe:
            PrecisionApp.exit(str(ioe), code=ExitCodes.INVALID_DATA)
        
        except ArgumentError as ae:
            PrecisionApp.exit(str(ae), code=ExitCodes.INVALID_DATA)

        self.log('\nDone.')      

class ArgumentError(ValueError):
    pass

class ArgumentIOError(ArgumentError):
    pass         

class GlobalInfo(object):
    
    def __init__(self, precision, coverage):
        self.precision = precision
        self.coverage = coverage
        
class LibrarySuperimposer(object):
    
    def __init__(self, native, library, pdb, output, cutoff=1.5):
    
        if not isinstance(native, structure.Chain):
            raise TypeError(native)
        elif native.length < 1:
            raise ArgumentError('The native chain has no residues')
        
        if not isinstance(library, rosetta.RosettaFragmentMap):
            raise TypeError(library)

        for i in [pdb, output]:
            if not os.path.isdir(i):
                raise ArgumentIOError("{0} is not a valid directory".format(i))

        self._pdb = pdb
        self._native = native
        self._library = library
        self._output = os.path.abspath(output)
        self._tab = os.path.join(self._output, native.entry_id + '.fragments.tab')
        self._figure = os.path.join(self._output, native.entry_id + '.precision.png')
        self._out = open(self._tab, 'w')
        self._cutoff = float(cutoff)

    def __del__(self):
        try:
            if self._out and not self._out.closed:
                self._out.close()
        except:
            pass
        
    @property
    def table(self):
        return self._tab
    
    @property
    def figure(self):
        return self._figure
            
    def run(self, cpu=1):
        
        tasks = []
        matches = []
        pool = multiprocessing.Pool(cpu)
        
        for source in self._library.sources:
            fragments = self._library.fromsource(source)
            task = pool.apply_async(rmsd, [self._native, source, fragments, self._pdb])
            tasks.append(task)
        
        for task in tasks:
            for match in task.get():
                if isinstance(match, basestring):       # error
                    self._out.write(' [E] ')                    
                    self._out.write(match.rstrip())
                    self._out.write('\n')
                else:
                    line = '{0.id}\t{0.qstart}\t{0.qend}\t{0.length}\t{0.rmsd}\n'.format(match)
                    self._out.write(line)
                    matches.append(match)
        
        return matches
                    
    def plot(self, matches):
        
        residues = range(1, self._native.length + 1)
        precision = []
        background = []
        covered = set()
        
        all = {}
        positive = {}
        
        for rank in residues:
            all[rank] = 0
            positive[rank] = 0
        
        for match in matches:
            for rank in range(match.qstart, match.qend + 1):                
                all[rank] += 1
                
                if match.rmsd <= self._cutoff:    
                    positive[rank] += 1
                    covered.add(rank)
        
        assert sorted(all.keys()) == sorted(positive.keys()) == residues
        
        for rank in residues:
            if all[rank] == 0:
                precision.append(0)
                background.append(0)
            else:
                precision.append((positive[rank] * 100.0 / all[rank]))
                background.append(100)
                
        coverage = len(covered) * 100.0 / len(residues)
        avg_precision = sum(precision) / len(precision)                 
        
        
        pyplot.ion()
        pyplot.bar(residues, background, color='#FFB0B0', linewidth=None, edgecolor='#FFB0B0')
        pyplot.bar(residues, precision, color='#50A6DA', linewidth=None, edgecolor='#50A6DA')
        
        pyplot.title(self._native.entry_id)
        pyplot.xlabel('Residue')
        pyplot.xlim(1, len(residues))
        pyplot.ylabel('Precision, %')
        pyplot.ylim(0, 100)
        pyplot.axes().xaxis.set_minor_locator(matplotlib.ticker.IndexLocator(1, 0))
        pyplot.axes().xaxis.set_major_locator(matplotlib.ticker.IndexLocator(5, 0))        
        
        try:
            pyplot.gcf().set_size_inches(15, 5.5)
            pyplot.savefig(self._figure)
        except IOError as io:
            raise ArgumentIOError("Can't save figure: " + str(io))
        
        return GlobalInfo(avg_precision, coverage)
            
def rmsd(target, source, fragments, pdb):

    matches = []
    
    try:
        src_file = wwpdb.find(source, [pdb])
        if src_file is None:
            raise IOError("Can't find structure {0} in {1}".format(source, pdb))

        src_structure = wwpdb.RegularStructureParser(src_file).parse_structure()

    except (IOError, structure.Broken3DStructureError, wwpdb.PDBParseError) as ex:
        matches.append('Error parsing {0:5}: {1!s}'.format(source, ex))
        return matches
            
    for fragment in fragments:            
            try:                
                if fragment.chain not in ('_', '', None):
                    src_chain = src_structure.chains[fragment.chain]
                else:
                    src_chain = src_structure.first_chain

                try: 
                    query = target.subregion(fragment.qstart, fragment.qend)
                    subject = src_chain.subregion(fragment.start, fragment.end, clone=True)
                    
                except IndexError:
                    matches.append('Fragment {1.source_id:>5} {0.start:>4}-{0.end:>4} is out of range'.format(fragment))                 
                    continue                    
            
                rmsd = query.rmsd(subject)
                match = csb.bio.fragments.FragmentMatch(fragment.id, qstart=fragment.qstart, qend=fragment.qend,
                                                        probability=None, rmsd=rmsd, tm_score=None, qlength=target.length)              
                matches.append(match)

            except (structure.Broken3DStructureError, IOError) as ex:
                matches.append("Can't superimpose fragment {0}: {1!s}".format(fragment.id, ex))          
                continue
                
    return matches         
     


if __name__ == '__main__':
    
    AppRunner().run()
