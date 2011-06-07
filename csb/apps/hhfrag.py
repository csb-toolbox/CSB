"""
HHfrag: build a dynamic variable-length fragment library for protein structure
prediction with Rosetta AbInitio.
"""

import os
import multiprocessing

import csb.apps
import csb.apps.hhsearch as hhsearch

import csb.bio.io
import csb.bio.fragments
import csb.bio.structure
import csb.pyutils


class ExitCodes(csb.apps.ExitCodes):
    
    IO_ERROR = 2
    INVALID_DATA = 3
    HHSEARCH_FAILURE = 4


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

        cmd.add_scalar_option('verbosity', 'v', int, 'verbosity level', default=2)        
        cmd.add_scalar_option('output', 'o', str, 'output directory', default='.')
        cmd.add_scalar_option('gap-filling', 'g', str, 'path to a Rosetta 9-mer fragment file, that will be used '
                              'to complement gaps in the fragment map (if specified, a joint fragment file will be produced)')
        cmd.add_boolean_option('filtered-map', 'f', 'make an additional filtered fragment map', default=False)
        
        cmd.add_positional_argument('QUERY', str, 'query profile HMM (e.g. created with csb.apps.buildhmm)')
                        
        return cmd
    
    
class HHfragApp(csb.apps.Application):
    
    def main(self):
        if not os.path.isdir(self.args.output):
            HHfragApp.exit('Output directory does not exist', code=ExitCodes.INVALID_DATA, usage=True)
            
        try:
            hhf = HHfrag(self.args.QUERY, self.args.hhsearch, self.args.database, logger=self)
            output = os.path.join(self.args.output, hhf.query.id)                
                    
            hhf.slice_query(self.args.min, self.args.max, self.args.step, self.args.cpu)
            hhf.extract_fragments()
            
            fragmap = hhf.build_fragment_map()
            fragmap.dump(output + '.hhfrags.09')
            
            if self.args.filtered_map:
                fragmap = hhf.build_filtered_map()
                fragmap.dump(output + '.filtered.09')
                
            if self.args.gap_filling:
                fragmap = hhf.build_combined_map(self.args.gap_filling)
                fragmap.dump(output + '.complemented.09')
                
            self.log('\nDONE.')

        except ArgumentIOError as ae:
            HHfragApp.exit(str(ae), code=ExitCodes.IO_ERROR)
                    
        except ArgumentError as ae:
            HHfragApp.exit(str(ae), code=ExitCodes.INVALID_DATA, usage=True)

        except csb.pyutils.InvalidCommandError as ose:
            msg = '{0!s}: {0.program}'.format(ose)
            HHfragApp.exit(msg, ExitCodes.IO_ERROR)   
                              
        except csb.pyutils.ProcessError as pe:
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
        self._output = None
        self._aligner = None
        
        self.database = database
        self.aligner = hhsearch.HHsearch(binary, self.pdbs25, cpu=2)

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

        if not 0 < min <= max:
            raise ArgumentError('min and max must be positive numbers, with max >= min')
        if not 0 < step:
            raise ArgumentError('step must be positive number')        

        self.log('\n# Processing query {0}...'.format(self.query.id))
        self.log('', level=2)                
        qp = self.query
        hsqs = []
              
        if not cpu:
            cpu = max - min + 1
        
        for start in range(1, qp.layers.length + 1, step):
            
            self.log('{0:3}.      '.format(start), ending='', level=2)
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
                self.log('{0.start:3} {0.end:3} ({0.length:2} aa)    {0.recurrence:3} hits'.format(rep), level=2)
            else:
                self.log('                     no hits', level=2)
        
        self._hsqs = hsqs
        return tuple(hsqs)
    
    def extract_fragments(self):

        self.log('\n# Extracting fragments...')
                        
        if self._hsqs is None:
            raise InvalidOperationError('The query has to be sliced first')
        
        fragments = []
        
        for si in self._hsqs:
            self.log('\nSEGMENT:   {0.start:3} {0.end:3}   ({0.recurrence})'.format(si), level=2)
            
            for hit in si.hits:
                cn = 0
                for chunk in hit.alignment.segments:
                    chunk.qstart = chunk.qstart + si.start - 1
                    chunk.qend = chunk.qend + si.start - 1
                    cn += 1
                    self.log(' {0.id:5}   L:{0.qstart:3} {0.qend:3}  {0.length:2}aa   P:{0.probability:5.3f}'.format(chunk), ending='', level=2)
                                    
                    sourcefile = os.path.join(self.database, hit.id + '.pdb')
                    source = csb.bio.io.StructureParser(sourcefile).parse().first_chain
                    assert hit.id[-1] == source.id
                    
                    source.compute_torsion()
                    try:
                        fragment = csb.bio.fragments.Assignment(source, 
                                                                chunk.start, chunk.end, hit.id, 
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
            bar = '{0:3} |{1} {2:5.1f}%'.format(length, '\xe2\x8a\x90' * int(percent), percent)
            self.log(bar, level=2)
    
    def build_fragment_map(self):

        self.log('\n# Building dynamic fragment map...')
                
        if self._matches is None:
            raise InvalidOperationError('You need to extract some fragments first')
        
        self._plot_lengths()
        
        target = csb.bio.fragments.Target.from_profile(self.query)
        target.assignall(self._matches)
             
        factory = csb.bio.fragments.RosettaFragsetFactory()
        return factory.make_fragset(target)

    def _filter_event_handler(self, rei):
        self.log('{0.rank:3}.    {0.confidence:5.3f}    {0.count:3}'.format(rei), level=2)
            
    def build_filtered_map(self):

        self.log('\n# Building filtered map...')
        self.log('', level=2)
                
        target = csb.bio.fragments.Target.from_profile(self.query)
        target.assignall(self._matches)

        factory = csb.bio.fragments.RosettaFragsetFactory()
        return factory.make_filtered(target, extend=True,
                                     callback=self._filter_event_handler)

    def _merge_event_handler(self, rei):
        if rei.confidence is None:
            self.log('{0.rank:3}.        -    {0.count:3}'.format(rei), level=2)
        else:
            self.log('{0.rank:3}.    {0.confidence:5.3f}    {0.count:3}'.format(rei), level=2)        
                    
    def build_combined_map(self, fragfile, top=25):

        self.log('\n# Building complemented map...')
        
        try:
            from csb.bio.fragments.rosetta import RosettaFragmentMap
            filling = RosettaFragmentMap.read(fragfile, top=top)
        except IOError as io:
            raise ArgumentIOError(str(io))

        self.log('\n  {0} rosetta fragments loaded'.format(filling.size))
        self.log('', level=2)               
   
        target = csb.bio.fragments.Target.from_profile(self.query)
        target.assignall(self._matches)
        
        factory = csb.bio.fragments.RosettaFragsetFactory()        
        return factory.make_combined(target, filling, threshold=0.5,
                                     callback=self._merge_event_handler)                


class SliceContext(hhsearch.Context):
    
    def __init__(self, segment, start, end):
        
        self.start = start
        self.end = end
        
        if not isinstance(segment, basestring):
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

    def __cmp__(self, other):
        
        if self.recurrence == other.recurrence:
            return cmp(self.length, other.length)
        else:      
            return cmp(self.recurrence, other.recurrence)




if __name__ == '__main__':
    
    AppRunner().run()
    