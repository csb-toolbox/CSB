"""
CSfrag: build a dynamic library of analogous fragments, given a list 
of assigned chemical shifts.
"""

import os
import numpy
import multiprocessing

import csb.io
import csb.apps

from csb.bio.io.wwpdb import FileSystemStructureProvider, StructureNotFoundError, PDBParseError
from csb.bio.nmr import RandomCoil, ChemShiftScoringModel
from csb.bio.structure import Chain, Broken3DStructureError
from csb.bio.fragments import ChemShiftTarget, ChemShiftAssignment, RosettaFragsetFactory
from csb.bio.io.cs import ChemShiftReader, ChemShiftFormatError
from csb.bio.io.fasta import SequenceParser, SequenceFormatError


class ExitCodes(csb.apps.ExitCodes):
    
    IO_ERROR = 2
    INVALID_DATA = 3
    NO_OUTPUT = 5    

class AppRunner(csb.apps.AppRunner):
    
    @property
    def target(self):
        return CSfragApp
     
    def command_line(self):
        
        cmd = csb.apps.ArgHandler(self.program, __doc__)
        cpu = multiprocessing.cpu_count()
        
        cmd.add_scalar_option('database', 'd', str, 'PDBS25 database directory (containing PDBS25cs.scs)', required=True)
        cmd.add_scalar_option('shifts', 's', str, 'assigned chemical shifts table (NMR STAR file fragment)', required=True)    

        cmd.add_scalar_option('window', 'w', int, 'sliding window size', default=8)
        cmd.add_scalar_option('top', 't', int, 'maximum number per starting position', default=25)                
        cmd.add_scalar_option('cpu', 'c', int, 'maximum degree of parallelism', default=cpu)

        cmd.add_scalar_option('verbosity', 'v', int, 'verbosity level', default=1)        
        cmd.add_scalar_option('output', 'o', str, 'output directory', default='.')        
        cmd.add_boolean_option('filtered-map', 'f', 'make an additional filtered fragment map of centroids', default=False)
                
        cmd.add_positional_argument('QUERY', str, 'query sequence  (FASTA file)')
                        
        return cmd
    
class CSfragApp(csb.apps.Application):
    
    def main(self):
        if not os.path.isdir(self.args.output):
            CSfragApp.exit('Output directory does not exist', code=ExitCodes.INVALID_DATA, usage=True)

        try:
            csf = CSfrag(self.args.QUERY, self.args.shifts, self.args.database, self.args.window, logger=self)
            output = os.path.join(self.args.output, csf.query.accession)                
                    
            frags = csf.extract_fragments(self.args.top, self.args.cpu)
            
            if len(frags) == 0:
                CSfragApp.exit('No fragments found!', code=ExitCodes.NO_OUTPUT)                
            
            fragmap = csf.build_fragment_map()
            fragmap.dump(output + '.csfrags.08')
            
            if self.args.filtered_map:
                fragmap = csf.build_filtered_map()
                fragmap.dump(output + '.filtered.08')
                
            self.log('\nDONE.')

        except ArgumentIOError as ae:
            CSfragApp.exit(str(ae), code=ExitCodes.IO_ERROR)
                    
        except ArgumentError as ae:
            CSfragApp.exit(str(ae), code=ExitCodes.INVALID_DATA, usage=True)

        except ChemShiftFormatError as ce:
            msg = "Can't parse input chemical shifts: " + str(ce)
            CSfragApp.exit(msg, code=ExitCodes.INVALID_DATA)
            

    def log(self, message, ending='\n', level=1):
        
        if level <= self.args.verbosity:
            super(CSfragApp, self).log(message, ending)
            
    
class SecondaryShiftConverter(object):
    """
    Helper, which reads assigned shifts from NMR STAR files and calculates
    corrected secondary shifts.
    """
    
    def convert(self, file, chain):
        """
        Compute secondary shofts.
        
        @param file: NMR STAR path and file name
        @type file: str
        @param chain: the protein chain, containing the chemical shifts
                      (L{Chain.from_sequence} may be useful)
        @type chain: L{Chain}
        
        @return: dictionary of the form: [rank: [nucleus: sec shift]]
        @rtype: dict
        """
        rc = RandomCoil.get()
        cs = {}
        
        for ni in ChemShiftReader().guess(file).read_file(file):
            
            if ni.name in ChemShiftScoringModel.NUCLEI:    
                ni.shift = rc.secondary_shift(chain, ni.position, ni.name, ni.shift)
                    
                cs.setdefault(ni.position, {})
                cs[ni.position][ni.name] = ni.shift
                
        return cs
    
class SecondaryShiftReader(object):
    """
    Reads secondary shifts from files in CSfrag format.
    """
    
    DB = 'pdbs25cs.scs'
    
    def read_shifts(self, string):
        """
        Read secondary shifts.
        @param string: complete secondary shift block
        @type string: str
        
        @return: dictionary of the form: [rank: [nucleus: sec shift]]
        @rtype: dict
        """
        
        shifts = {}
        
        for l in string.splitlines():
            
            if l.startswith('#') or not l.strip():
                continue
            
            l = l.split('\t')
            rank = int(l[0])
            
            for n, cs in zip(ChemShiftScoringModel.NUCLEI, l[1:]):
                if cs != '':
                    shifts.setdefault(rank, {})[n] = float(cs)
                    
        return shifts
        
    def load_database(self, path, file=DB):
        """
        Read the entire PDBS25CS database.
        
        @return: dictionary of the form: [entry ID: [rank: [nucleus: sec shift]]]
        @rtype: dict        
        """
        
        db = {}
        file = os.path.join(path, file)
        
        with open(file) as stream:
            er = csb.io.EntryReader(stream, '#', None)
            
            for e in er.entries():
                entry = e[10:15]
                db[entry] = self.read_shifts(e)
                
        return db         

class ScoringHelper(object):
    
    def __init__(self, window):
        
        self._window = window
        self._model = ChemShiftScoringModel()
        
    @property
    def window(self):
        return self._window

    def score(self, qcs, scs, qstart, qend, sstart, send):
        
        window = self._window
        
        if window is None:
            window = min(qend - qstart + 1, send - sstart + 1)
            
        off_start, off_end = self.offsets(qstart, qend, window=window)
        qs = qstart + off_start
        qe = qend - off_end
        ss = sstart + off_start
        se = send - off_end
        
        assert qe - qs + 1 == se - ss + 1 == window
                
        score = 0

        for nucleus in ChemShiftScoringModel.NUCLEI:
            query = []
            subject = []

            for qr, sr in zip(range(qs, qe + 1), range(ss, se + 1)):
                try:
                    qshift = qcs[qr][nucleus]
                    sshift = scs[sr][nucleus]
                    
                    if qshift is not None and sshift is not None:
                        query.append(qshift)
                        subject.append(sshift)

                except KeyError:
                    continue
        
            if query and subject:
                deltas = numpy.array(query) - numpy.array(subject)
                score += self._model.score(nucleus, deltas).sum()

        return score
        
    def offsets(self, start, end, window=6):
        
        if end - start + 1 <= window:
            return 0, 0
                
        d1 = ((end - start + 1) - window) / 2
        ns = start + d1
        ne = ns + window - 1
        d2 = end - ne
        
        return d1, d2
    

class ArgumentError(ValueError):
    pass

class ArgumentIOError(ArgumentError):
    pass 

class InvalidOperationError(ValueError):
    pass

    
class CSfrag(object):
    """
    @param query: query FASTA sequence path and file name
    @type query: str
    @param cstable: file, containing the table of assigned experimental chemical shifts 
    @type cstable: str
    @param database: path to the PDBS25 directory
    @type database: str
    @param logger: logging client (needs to have a C{log} method)
    @type logger: L{Application}
    """
    
    def __init__(self, query, cstable, database, window=8, logger=None):
        
        self._query = None
        self._qcs = None  
        self._matches = None
        self._helper = ScoringHelper(window)
        self._database = None
        self._window = None
        self._app = logger
        self._pdb = None
        
        try:
            fasta = SequenceParser().parse_file(query)
            if len(fasta) != 1:
                raise ArgumentError("The input FASTA file should contain one sequence")
            elif fasta[0].length < 1:
                raise ArgumentError("Zero-length query sequence")                
            
            self._query = Chain.from_sequence(fasta[0], 'A')
            self._query.accession = fasta[0].id
            self._qcs = SecondaryShiftConverter().convert(cstable, self._query)
            
            if len(self._qcs) == 0:
                raise ArgumentError("No chemical shifts read; check your input")    
        
        except IOError as io:
            raise ArgumentIOError(str(io))
        
        except SequenceFormatError as se:
            raise ArgumentError("Can't parse FASTA file: {0}".format(str(se)))
           
                
        self.database = database
        self.window = window

    @property
    def query(self):
        return self._query
            
    @property
    def database(self):
        return self._database
    @database.setter
    def database(self, value):
        database = value
        pdbs25cs = os.path.join(value, SecondaryShiftReader.DB)
        if not os.path.isfile(pdbs25cs):
            raise ArgumentError('PDBS25CS not found here: ' + pdbs25cs)
        self._database = database
        self._pdb = FileSystemStructureProvider(database)
        
    @property
    def window(self):
        return self._window
    @window.setter
    def window(self, value):
        value = int(value)
        if value < 1:
            raise ValueError("Invalid sliding window: {0}".format(value)) 
        self._window = value        
        
    def log(self, *a, **ka):
        if self._app:
            self._app.log(*a, **ka)        

    def extract_fragments(self, top=25, cpu=2):
        """
        Extract fragments with matching chemical shifts using a sliding window.
        
        @param top: L{MatchTable} capacity per starting position
        @type top: int
        @param cpu: degree of parallelism
        @type cpu: int
                
        @rtype: tuple of L{ChemShiftAssignment}s
        """        
        self.log("# Reading chemical shifts...", level=1)
        db = SecondaryShiftReader().load_database(self.database)
        matches = MatchTable(self.query.length, capacity=top)
                        
        slices = []
        fragments = []
            
        for qs in range(1, self.query.length + 1):
            qe = qs + self.window - 1
            if qe > self.query.length:
                break
            
            slices.append((qs, qe))
            
        self.log("\n# Processing target {0}...".format(self.query.accession), level=1)
        pool = multiprocessing.Pool(cpu)
                
        try:
            for subject in db:
                tasks = []
                
                for qs, qe in slices:
                    task = pool.apply_async(_task, [self._helper, subject, qs, qe, self._qcs, db[subject]])
                    tasks.append(task)
                    
                for task in tasks:
                    for result in task.get():
                        if result.score > ChemShiftAssignment.BIT_SCORE_THRESHOLD * self.window:
                            matches.add(result)
                            
        except KeyboardInterrupt:
            pass
        finally:
            pool.terminate()
        
        for rank in matches:
            msg = '{0:3} {1:3} ({2:2} aa)    {3:3} fragments'
            self.log(msg.format(rank, rank + self.window - 1, self.window, len(matches[rank])), 
                     level=1)                      

        
        self.log("\n# Extracting fragments...")
        
        for group in matches.by_source:
            try:
                source_id = group[0].entry_id
                source = self._pdb.get(source_id).first_chain
                source.compute_torsion()
                         
                for match in group:
                    try:
                        row = ' {0.entry_id:5}   L{0.qs:3} {0.qe:3}  {1}aa   S:{0.score:5.1f}'
                        self.log(row.format(match, self.window), ending='', level=2)
                        
                        fragment = ChemShiftAssignment(source=source, start=match.ss, end=match.se, 
                                                       qstart=match.qs, qend=match.qe, 
                                                       window=self.window, rmsd=None, score=match.score)                
                        fragments.append(fragment)
                        self.log('', level=2)                
                                                
                    except Broken3DStructureError:
                        self.log('    corrupt', level=2)
                        continue
            except PDBParseError:
                continue
            except StructureNotFoundError:
                self.log("  Warning: Template {0} is missing!".format(source_id))

        self._matches = fragments
        return tuple(fragments)

    def build_fragment_map(self):
        """
        Build a full Rosetta fragset.
        @rtype: L{RosettaFragmentMap}
        """

        if self._matches is None:
            self.extract_fragments()   

        self.log('\n# Building fragment map...')
        
        target = ChemShiftTarget(self.query.accession, self.query.length, self.query.residues)
        target.assignall(self._matches)
             
        factory = RosettaFragsetFactory()
        return factory.make_fragset(target)
    
    def build_filtered_map(self):
        """
        Build a filtered fragset of centroids.
        @rtype: L{RosettaFragmentMap}
        """
        
        if self._matches is None:
            self.extract_fragments()  
            
        self.log('\n# Building filtered map...')
                
        target = ChemShiftTarget(self.query.accession, self.query.length, self.query.residues)
        target.assignall(self._matches)

        factory = RosettaFragsetFactory()
        return factory.make_filtered(target, extend=False)        
        
class MatchInfo(object):
    
    def __init__(self, entry_id, qs, qe, ss, se, score):
        
        self.entry_id = entry_id
        self.qs = qs
        self.qe = qe
        self.ss = ss
        self.se = se
        self.score = score

    def __str__(self):
        return '{0.qs:4} {0.qe:4} {0.ss:4} {0.se:4} {0.score:10.3f}'.format(self)
    
    def __cmp__(self, other):
        return cmp(self.score, other.score)        
    
class MatchTable(object):
    
    def __init__(self, length, capacity=25):
        
        if capacity < 1:
            capacity = 1
            
        self._capacity = capacity
        self._length = length
        self._t = {}
        
        for i in range(1, length + 1):
            self._t[i] = []
            
    def add(self, m):
        
        matches = self._t[m.qs]
        
        if len(matches) < self._capacity:
            
            matches.append(m)
            matches.sort()
                        
        elif m.score > matches[-1].score:
            
            matches.pop()
            matches.append(m)
            matches.sort()
    
    def __getitem__(self, rank):
        return tuple(self._t[rank])
    
    def __iter__(self):
        return iter(self._t)
    
    @property
    def by_source(self):
        
        matches = {}
        
        for rank in self:
            for m in self[rank]:
                if m.entry_id not in matches:
                    matches[m.entry_id] = []
                
                matches[m.entry_id].append(m)
            
        for entry_id in matches:
            yield tuple(matches[entry_id])
                
def _task(helper, subject, qs, qe, qcs, scs):
    
    try:
        results = []
        slength = max(scs or [0])
        
        for ss in range(1, slength + 1, 3):
            se = ss + helper.window - 1
            
            if se > slength:
                break
            
            score = helper.score(qcs, scs, qs, qe, ss, se)
            
            if score is not None:
                info = MatchInfo(subject, qs, qe, ss, se, score)
                results.append(info)
        
        return results
    
    except KeyboardInterrupt:
        return []




if __name__ == '__main__':
    
    args = "cs.py -v 1 -f -t 12 -d /home/ivan/Desktop/cstest/db -s /home/ivan/Desktop/cstest/t.str -o /home/ivan/Desktop/cstest /home/ivan/Desktop/cstest/t.fa".split() 
    AppRunner().run()
    