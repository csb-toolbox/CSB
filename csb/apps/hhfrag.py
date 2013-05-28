"""
HHfrag: build a dynamic variable-length fragment library for protein structure
prediction with Rosetta AbInitio.
"""

import os
import multiprocessing

import csb.apps
import csb.apps.hhsearch as hhsearch

import csb.bio.io.hhpred
import csb.bio.fragments
import csb.bio.fragments.rosetta as rosetta
import csb.bio.structure

import csb.io.tsv
import csb.core


class ExitCodes(csb.apps.ExitCodes):
    
    IO_ERROR = 2
    INVALID_DATA = 3
    HHSEARCH_FAILURE = 4
    NO_OUTPUT = 5


class AppRunner(csb.apps.AppRunner):
    
    @property
    def target(self):
        return HHfragApp
     
    def command_line(self):
        
        cmd = csb.apps.ArgHandler(self.program, __doc__)
        cpu = multiprocessing.cpu_count()
        
        cmd.add_scalar_option('hhsearch', 'H', str, 'path to the HHsearch executable', default='hhsearch')
        cmd.add_scalar_option('database', 'd', str, 'database directory (containing PDBS25.hhm)', required=True)        
        
        cmd.add_scalar_option('min', 'm', int, 'minimum query segment length', default=6)
        cmd.add_scalar_option('max', 'M', int, 'maximum query segment length', default=21)
        cmd.add_scalar_option('step', 's', int, 'query segmentation step', default=3)
        cmd.add_scalar_option('cpu', 'c', int, 'maximum degree of parallelism', default=cpu)
        
        cmd.add_scalar_option('gap-filling', 'g', str, 'path to a fragment file (e.g. CSfrag or Rosetta NNmake), which will be used '
                              'to complement low-confidence regions (when specified, a hybrid fragment library be produced)')
        cmd.add_scalar_option('filtered-filling', 'F', str, 'path to a filtered fragment file (e.g. filtered CSfrag-ments), which will '
                              'be mixed with the  HHfrag-set and then filtered, resulting in a double-filtered library')
        cmd.add_boolean_option('filtered-map', 'f', 'make an additional filtered fragment map of centroids and predict torsion angles', default=False)        
        cmd.add_boolean_option('c-alpha', None, 'include also C-alpha vectors in the output', default=False)
        cmd.add_scalar_option('confidence-threshold', 't', float, 'confidence threshold for gap filling', default=0.7)

        cmd.add_scalar_option('verbosity', 'v', int, 'verbosity level', default=2)        
        cmd.add_scalar_option('output', 'o', str, 'output directory', default='.')
                        
        cmd.add_positional_argument('QUERY', str, 'query profile HMM (e.g. created with csb.apps.buildhmm)')
                        
        return cmd
    
    
class HHfragApp(csb.apps.Application):
    
    def main(self):
        if not os.path.isdir(self.args.output):
            HHfragApp.exit('Output directory does not exist', code=ExitCodes.INVALID_DATA, usage=True)
            
        if self.args.c_alpha:
            builder = rosetta.ExtendedOutputBuilder
        else:
            builder = rosetta.OutputBuilder
            
        try:
            hhf = HHfrag(self.args.QUERY, self.args.hhsearch, self.args.database, logger=self)
            output = os.path.join(self.args.output, hhf.query.id)                
                    
            hhf.slice_query(self.args.min, self.args.max, self.args.step, self.args.cpu)
            frags = hhf.extract_fragments()
            
            if len(frags) == 0:
                HHfragApp.exit('No fragments found!', code=ExitCodes.NO_OUTPUT)                
            
            fragmap = hhf.build_fragment_map()
            fragmap.dump(output + '.hhfrags.09', builder)
            
            if self.args.filtered_map:
                fragmap, events = hhf.build_filtered_map()
                fragmap.dump(output + '.filtered.09', builder)
                tsv = PredictionBuilder.create(events).product
                tsv.dump(output + '.centroids.tsv')

            if self.args.filtered_filling:
                fragmap, events = hhf.build_hybrid_filtered_map(self.args.filtered_filling)
                fragmap.dump(output + '.hybrid.filtered.09', builder)
                tsv = PredictionBuilder.create(events).product
                tsv.dump(output + '.hybrid.centroids.tsv')
                                                
            if self.args.gap_filling:
                fragmap = hhf.build_combined_map(self.args.gap_filling, self.args.confidence_threshold)
                fragmap.dump(output + '.complemented.09', builder)
                
                                
            self.log('\nDONE.')

        except ArgumentIOError as ae:
            HHfragApp.exit(str(ae), code=ExitCodes.IO_ERROR)
                    
        except ArgumentError as ae:
            HHfragApp.exit(str(ae), code=ExitCodes.INVALID_DATA)

        except csb.io.InvalidCommandError as ose:
            msg = '{0!s}: {0.program}'.format(ose)
            HHfragApp.exit(msg, ExitCodes.IO_ERROR)

        except csb.bio.io.hhpred.HHProfileFormatError as hfe:
            msg = 'Corrupt HMM: {0!s}'.format(hfe)
            HHfragApp.exit(msg, code=ExitCodes.INVALID_DATA)           
                              
        except csb.io.ProcessError as pe:
            message = 'Bad exit code from HHsearch: #{0.code}.\nSTDERR: {0.stderr}\nSTDOUT: {0.stdout}'.format(pe.context)
            HHfragApp.exit(message, ExitCodes.HHSEARCH_FAILURE)
            
              
        
    def log(self, message, ending='\n', level=1):
        
        if level <= self.args.verbosity:
            super(HHfragApp, self).log(message, ending)

class ArgumentError(ValueError):
    pass

class ArgumentIOError(ArgumentError):
    pass 

class InvalidOperationError(ValueError):
    pass


class HHfrag(object):
    """
    The HHfrag dynamic fragment detection protocol.
    
    @param query: query HMM path and file name
    @type query: str
    @param binary: the HHsearch binary
    @type binary: str
    @param database: path to the PDBS25 directory
    @type database: str
    @param logger: logging client (needs to have a C{log} method)
    @type logger: L{Application}
    """
    
    PDBS = 'pdbs25.hhm'
    
    def __init__(self, query, binary, database, logger=None):
        
        try:
            self._query = csb.bio.io.HHProfileParser(query).parse()
        except IOError as io:
            raise ArgumentIOError(str(io))
        self._hsqs = None
        self._matches = None
        
        self._app = logger
        self._database = None
        self._pdbs25 = None
        self._aligner = None
        
        self.database = database
        self.aligner = hhsearch.HHsearch(binary, self.pdbs25, cpu=2)
        
        if self.query.layers.length < 1:
            raise ArgumentError("Zero-length sequence profile")

    @property
    def query(self):
        return self._query
    
    @property
    def pdbs25(self):
        return self._pdbs25
             
    @property
    def database(self):
        return self._database
    @database.setter
    def database(self, value):
        database = value
        pdbs25 = os.path.join(value, HHfrag.PDBS)
        if not os.path.isfile(pdbs25):
            raise ArgumentError('PDBS25 not found here: ' + pdbs25)
        self._database = database
        self._pdbs25 = pdbs25
    
    @property
    def aligner(self):
        return self._aligner
    @aligner.setter
    def aligner(self, value):
        if hasattr(value, 'run') and hasattr(value, 'runmany'):
            self._aligner = value
        else:
            raise TypeError(value)
        
    def log(self, *a, **ka):
        
        if self._app:
            self._app.log(*a, **ka)
        
    def slice_query(self, min=6, max=21, step=3, cpu=None):
        """
        Run the query slicer and collect the optimal query segments.
        
        @param min: min segment length
        @type min: int
        @param max: max segment length
        @type max: int
        @param step: slicing step
        @type step: int
        @param cpu: degree of parallelism
        @type cpu: int
                
        @rtype: tuple of L{SliceContext}
        """

        if not 0 < min <= max:
            raise ArgumentError('min and max must be positive numbers, with max >= min')
        if not 0 < step:
            raise ArgumentError('step must be positive number')        

        self.log('\n# Processing profile HMM "{0}"...'.format(self.query.id))
        self.log('', level=2)                
        qp = self.query
        hsqs = []
              
        if not cpu:
            cpu = max - min + 1
        
        for start in range(1, qp.layers.length - min + 1 + 1, step):
            
            self.log('{0:3}.      '.format(start), ending='', level=1)
            probes = []
            
            for end in range(start + min - 1, start + max):
                if end > qp.layers.length:
                    break
                context = SliceContext(qp.segment(start, end), start, end)
                probes.append(context)
                
            probes = self.aligner.runmany(probes, workers=cpu)
            probes.sort()
            
            if len(probes) > 0:
                rep = probes[-1]
                hsqs.append(rep)
                self.log('{0.start:3} {0.end:3} ({0.length:2} aa)    {0.recurrence:3} hits'.format(rep), level=1)
            else:
                self.log('                     no hits', level=1)
        
        self._hsqs = hsqs
        return tuple(hsqs)
    
    def extract_fragments(self):
        """
        Extract all matching fragment instances, given the list of optimal 
        query slices, generated during the first stage.
        
        @rtype: tuple of L{Assignment}s
        """
                        
        if self._hsqs is None:
            self.slice_query()

        self.log('\n# Extracting fragments...')        
        fragments = []
        
        for si in self._hsqs:
            self.log('\nSEGMENT:  {0.start:3} {0.end:3}   ({0.recurrence})'.format(si), level=2)
            
            for hit in si.hits:
                cn = 0
                for chunk in hit.alignment.segments:
                    chunk.qstart = chunk.qstart + si.start - 1
                    chunk.qend = chunk.qend + si.start - 1
                    cn += 1
                    self.log(' {0.id:5}   L{0.qstart:3} {0.qend:3}  {0.length:2}aa   P:{0.probability:5.3f}'.format(chunk), ending='', level=2)
                                    
                    sourcefile = os.path.join(self.database, hit.id + '.pdb')
                    if not os.path.isfile(sourcefile):
                        self.log('    missing', level=2)
                        continue
                    source = csb.bio.io.StructureParser(sourcefile).parse().first_chain
                    assert hit.id[-1] == source.id
                    
                    source.compute_torsion()
                    try:
                        fragment = csb.bio.fragments.Assignment(source, chunk.start, chunk.end,  
                                                                chunk.qstart, chunk.qend, chunk.probability, 
                                                                rmsd=None, tm_score=None)
                        fragments.append(fragment)
                        if cn > 1:
                            self.log('    (chunk #{0})'.format(cn), level=2)
                        else:
                            self.log('', level=2)
                                          
                    except csb.bio.structure.Broken3DStructureError:
                        self.log('    corrupt', level=2)
                        continue                    
        
        self._matches = fragments
        return tuple(fragments)
    
    def _plot_lengths(self):
        
        self.log('\n  {0} ungapped assignments'.format(len(self._matches)))
        self.log('', level=2)
                
        histogram = {}
        for f in self._matches:
            histogram[f.length] = histogram.get(f.length, 0) + 1
            
        for length in sorted(histogram):
            
            percent = histogram[length] * 100.0 / len(self._matches)
            bar = '{0:3} |{1} {2:5.1f}%'.format(length, 'o' * int(percent), percent)
            self.log(bar, level=2)
    
    def build_fragment_map(self):
        """
        Build a full Rosetta fragset.
        @rtype: L{RosettaFragmentMap}
        """

        if self._matches is None:
            self.extract_fragments()        

        self.log('\n# Building dynamic fragment map...')        
        self._plot_lengths()
        
        target = csb.bio.fragments.Target.from_profile(self.query)
        target.assignall(self._matches)
            
        factory = csb.bio.fragments.RosettaFragsetFactory()
        return factory.make_fragset(target)

    def _filter_event_handler(self, ri):    
               
        if ri.gap is True or ri.confident is False:
            self.log('{0.rank:3}.     {0.confidence:5.3f}         {0.count:3}'.format(ri), level=2)            
        
        else:
            phi = PredictionBuilder.format_angle(ri.torsion.phi)
            psi = PredictionBuilder.format_angle(ri.torsion.psi)
            omega = PredictionBuilder.format_angle(ri.torsion.omega)
        
            pred = "{0.source_id:5}  {0.start:3} {0.end:3}    {1} {2} {3}".format(ri.rep, phi, psi, omega)
            self.log('{0.rank:3}.     {0.confidence:5.3f}         {0.count:3}    {1}'.format(ri, pred), level=2)
            
    def build_filtered_map(self):
        """
        Build a filtered fragset of centroids.
        @return: filtered fragset and a list of residue-wise predictions
                 (centroid and torsion angles) 
        @rtype: L{RosettaFragmentMap}, list of L{ResidueEventInfo}
        """
        
        if self._matches is None:
            self.extract_fragments()  
            
        self.log('\n# Building filtered map...')
        self.log('\n    Confidence  Recurrence    Representative       Phi    Psi  Omega', level=2)
        
        events = []
        def logger(ri):
            events.append(ri)
            self._filter_event_handler(ri)
                
        target = csb.bio.fragments.Target.from_profile(self.query)
        target.assignall(self._matches)

        factory = csb.bio.fragments.RosettaFragsetFactory()
        fragset = factory.make_filtered(target, extend=True, callback=logger)
        
        return fragset, events

    def _merge_event_handler(self, ri):
        
        marked = ""
        
        if ri.gap is True or ri.confident is False:
            marked = "*"
    
        self.log('{0.rank:3}.     {0.confidence:5.3f}         {0.count:3}     {1:>3}'.format(ri, marked), level=2)        
                    
    def build_combined_map(self, fragfile, threshold=0.7, top=25):
        """
        Build a hybrid map, where low-confidence regions are complemented
        with the specified filling.
        
        @param threshold: confidence threshold
        @type threshold: float
        @param fragfile: filling fragset (Rosetta fragment file)
        @type fragfile: str
        
        @return: filtered fragset and a list of residue-wise predictions
                 (centroid and torsion angles) 
        @rtype: L{RosettaFragmentMap}, list of L{ResidueEventInfo}
        """
        
        if self._matches is None:
            self.extract_fragments()  
            
        self.log('\n# Building complemented map...')
        
        try:
            filling = rosetta.RosettaFragmentMap.read(fragfile, top=top)
        except IOError as io:
            raise ArgumentIOError(str(io))

        self.log('\n  {0} supplementary fragments loaded'.format(filling.size))
        self.log('    Confidence  Recurrence   Fill?', level=2)               
   
        target = csb.bio.fragments.Target.from_profile(self.query)
        target.assignall(self._matches)
        
        factory = csb.bio.fragments.RosettaFragsetFactory()        
        return factory.make_combined(target, filling, threshold=threshold,
                                     callback=self._merge_event_handler)                


    def build_hybrid_filtered_map(self, fragfile):
        """
        Mix the fragset with the specified (filtered)filling and then filter 
        the mixture. If the filling is a filtered CSfrag library, this will 
        produce a double-filtered map.
        
        @param fragfile: filtered filling (filtered CSfrag fragment file)
        @type fragfile: str        
        
        @rtype: L{RosettaFragmentMap}
        """
        
        if self._matches is None:
            self.extract_fragments()  
            
        self.log('\n# Building hybrid filtered map...')
        
        filling = []
        events = []

        def logger(ri):
            events.append(ri)
            self._filter_event_handler(ri)
        
        try:
            db = csb.bio.io.wwpdb.FileSystemStructureProvider(self.database)
                        
            for f in rosetta.RosettaFragmentMap.read(fragfile):
                filling.append(csb.bio.fragments.Assignment.from_fragment(f, db))
                
        except IOError as io:
            raise ArgumentIOError(str(io))
        except csb.bio.io.wwpdb.StructureNotFoundError, sne:
            msg = "{0} is not a PDBS25-derived fragset (template {1} not found)"
            raise ArgumentIOError(msg.format(fragfile, str(sne)))

        self.log('\n  {0} supplementary fragments loaded'.format(len(filling)))
        self.log('\n    Confidence  Recurrence    Representative       Phi    Psi  Omega', level=2)
        
        if len(filling) > self.query.layers.length:
            msg = "{0} does not look like a filtered fragset (too many centroids)"
            raise ArgumentError(msg.format(fragfile))
                
        target = csb.bio.fragments.Target.from_profile(self.query)
        target.assignall(self._matches)
        target.assignall(filling)

        factory = csb.bio.fragments.RosettaFragsetFactory()
        fragset = factory.make_filtered(target, extend=False, callback=logger)
        
        return fragset, events
        
class SliceContext(hhsearch.Context):
    
    def __init__(self, segment, start, end):
        
        self.start = start
        self.end = end
        
        if not isinstance(segment, csb.core.string):
            segment = segment.to_hmm(convert_scores=True)
        
        super(SliceContext, self).__init__(segment)
        
    @property
    def length(self):
        return self.end - self.start + 1
    
    @property
    def hits(self):
        return self.result
    
    @property
    def recurrence(self):
        return len(self.result)

    def __lt__(self, other):
        
        if self.recurrence == other.recurrence:
            return self.length < other.length
        else:      
            return self.recurrence < other.recurrence


class PredictionBuilder(object):
    
    NULL = '{0:>6}'.format("-")
    HEADER = "rank:int residue:str confidence:float centroid:str phi:float psi:float omega:float" 
    
    @staticmethod
    def format_angle(angle):
        """
        @param angle: torsion angle value
        @type angle: float
        """
        
        if angle is None:
            return '{0:>6}'.format("-")
        else:
            return '{0:6.1f}'.format(angle)  
        
    @staticmethod
    def create(ri):
        """
        @param ri: all predictions
        @type ri: list of L{ResidueEventInfo}
        """
        builder = PredictionBuilder()
        builder.addall(ri)
        return builder
    
    def __init__(self):
        self._tsv = csb.io.tsv.Table(PredictionBuilder.HEADER)
        
    @property
    def product(self):
        """
        @rtype: L{Table}
        """
        return self._tsv
    
    def isnull(self, v):
        if v is None:
            return "-"
        else:
            return v
        
    def add(self, ri):
        """
        @param ri: single residue prediction
        @type ri: L{ResidueEventInfo}
        """
        
        row = [ri.rank, repr(ri.type), ri.confidence]
        
        if ri.rep:
            row.append(ri.rep.id)
            row.append(self.isnull(ri.torsion.phi))
            row.append(self.isnull(ri.torsion.psi))
            row.append(self.isnull(ri.torsion.omega))
        
        else:
            row.extend([PredictionBuilder.NULL] * 4)
        
        self.product.insert(row)
        
    def addall(self, ri):
        """
        @param ri: all predictions
        @type ri: list of L{ResidueEventInfo}
        """
                
        ri = list(ri)
        ri.sort(key=lambda i: i.rank)
        
        for i in ri:
            self.add(i)
            
    
    
if __name__ == '__main__':
    
    AppRunner().run()
    