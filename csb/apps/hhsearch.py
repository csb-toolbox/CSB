"""
Python bindings for the HHsearch program. Capable of executing multiple
HHsearch jobs in parallel.
"""

import multiprocessing as mp

import csb.apps
import csb.io
import csb.bio.io
import csb.pyutils


class ExitCodes(csb.apps.ExitCodes):
    
    IO_ERROR = 2
    INVALID_DATA = 3
    EXT_TOOL_FAILURE = 4
    

class AppRunner(csb.apps.AppRunner):
    
    @property
    def target(self):
        return HHsearchApp
    
    def command_line(self):
        
        cpu = mp.cpu_count()
        cmd = csb.apps.ArgHandler(self.program, __doc__)
          
        cmd.add_scalar_option('binary', 'b', str, 'full path to the HHsearch binary ', default='hhsearch')
        cmd.add_scalar_option('cpu', 'c', int, 'maximum degree of parallelism', default=cpu)            
        cmd.add_scalar_option('database', 'd', str, 'the subject (database) HMM file', required=True)     
        cmd.add_array_argument('query', str, 'query HMM file(s)')      
        
        return cmd
    

class HHsearchApp(csb.apps.Application):
    
    def main(self):
        
        queries = list(self.args.query)
        exe = HHsearch(self.args.binary, self.args.database)
        
        try:
            if len(queries) == 1:
                exe.cpu = self.args.cpu
                context = HHTask(queries[0])
                results = [ exe.run(context) ]
            else:
                context = [ HHTask(q) for q in queries ]
                results = exe.runmany(context, workers=self.args.cpu)

        except IOError as io:
            HHsearchApp.exit(str(io), ExitCodes.IO_ERROR)

        except csb.pyutils.InvalidCommandError as ose:
            msg = '{0!s}: {0.program}'.format(ose)
            HHsearchApp.exit(msg, ExitCodes.IO_ERROR)   
                              
        except csb.pyutils.ProcessError as pe:
            message = 'Bad exit code from HHsearch: #{0.code}.\nSTDERR: {0.stderr}\nSTDOUT: {0.stdout}'.format(pe.context)
            HHsearchApp.exit(message, ExitCodes.EXT_TOOL_FAILURE)

        print '\nRank Hit   Prob  St  End Qst Qend'
        print '-------------------------------------'
                            
        for c in results:
            print '\n\n# QUERY:{0}\n'.format(c.queryfile)
            if c.result:
                for hit in c.result:
                    print '{0.rank:3}. {0.id:5} {0.probability:5.3f} {0.start:3} {0.end:3} {0.qstart:3} {0.qend:3}'.format(hit)


class Context(object):
    
    def __init__(self, query):
        
        self.__query = query
        self.__result = None
    
    @property
    def query(self):
        return self.__query
    
    @property
    def result(self):
        return self.__result
    @result.setter
    def result(self, result):
        self.__result = result
    

class HHTask(Context):
    
    def __init__(self, queryfile):
        
        self.queryfile = queryfile
        query = open(queryfile).read()
        
        super(HHTask, self).__init__(query)
    
        
def _task(args):
    
    try:
        binary, db, cpu, context = args
        return HHsearch(binary, db, cpu=cpu).run(context)
    except (KeyboardInterrupt, SystemExit):
        print '\nTerminating...'
        return

class SecStructureScoring(object):
    
    OFF = 0
    AFTER = 1
    DURING = 2
    AFTER_PREDICTED = 3
    DURING_PREDICTED = 4
    
class HHsearch(object):
    
    class Options(object):
        
        CPU = 'cpu'
        SS = 'ssm'
        MACT = 'mact'
        MAX_HITS = 'Z'
        MAX_ALI = 'B'
        MAX_E = 'E'
        MIN_P = 'p'
    
    def __init__(self, binary, db, cpu=None):
        
        self._program = binary
        self._db = db
        self._opt = {}
        self._parser = csb.bio.io.HHOutputParser()
        
        self.cpu = cpu
        self.ss = None
        self.mac_threshold = None
        self.max_hits = None
        self.max_alignments = None
        self.max_evalue = None
        self.min_probability = None

    @property
    def program(self):
        return self._program
    @program.setter
    def program(self, value):
        self._program = value
      
    @property
    def db(self):
        return self._db
    @db.setter
    def db(self, value):
        self._db = value
    
    @property
    def parser(self):
        return self._parser
    @parser.setter
    def parser(self, value):
        self._parser = value
    
    @property
    def cpu(self):
        return self._get(HHsearch.Options.CPU)
    @cpu.setter
    def cpu(self, value):
        self._opt[HHsearch.Options.CPU] = value
                        
    @property
    def ss(self):
        return self._get(HHsearch.Options.SS)
    @ss.setter
    def ss(self, value):
        self._opt[HHsearch.Options.SS] = value
    
    @property
    def mac_threshold(self):
        return self._get(HHsearch.Options.MACT)
    @mac_threshold.setter
    def mac_threshold(self, value):
        self._opt[HHsearch.Options.MACT] = value
    
    @property
    def max_hits(self):
        return self._get(HHsearch.Options.MAX_HITS)
    @max_hits.setter
    def max_hits(self, value):
        self._opt[HHsearch.Options.MAX_HITS] = value
    
    @property
    def max_alignments(self):
        return self._get(HHsearch.Options.MAX_ALI)
    @max_alignments.setter
    def max_alignments(self, value):
        self._opt[HHsearch.Options.MAX_ALI] = value
    
    @property
    def max_evalue(self):
        return self._get(HHsearch.Options.MAX_E)
    @max_evalue.setter
    def max_evalue(self, value):
        self._opt[HHsearch.Options.MAX_E] = value
    
    @property
    def min_probability(self):
        return self._get(HHsearch.Options.MIN_P)
    @min_probability.setter
    def min_probability(self, value):
        self._opt[HHsearch.Options.MIN_P] = value
    
    def _get(self, option):
        
        if option in self._opt:
            return self._opt[option]
        else:
            return None
        
    def _options(self):
        
        options = []
        
        for option in self._opt:
            value = self._opt[option]
            
            if value is not None and value != '':
                if isinstance(value, bool):
                    options.append('-{0}'.format(option))
                else:    
                    options.append('-{0} {1}'.format(option, value))         
        
        return ' '.join(options)
                        
    def run(self, context):
        
        with csb.io.TempFile() as q:
            
            q.write(context.query)
            q.flush()
            
            with csb.io.TempFile() as o:
                
                cmd = '{0.program} -i {1} -d {0.db} -o {2} {3}'.format(self, q.name, o.name, self._options())                    
                csb.pyutils.Shell.runstrict(cmd)
                
                context.result = self.parser.parse_file(o.name)
                return context
            
    def runmany(self, contexts, workers=mp.cpu_count(), cpu=1):
        
        if workers > len(contexts):
            workers = len(contexts)

        results = []
        taskargs = [(self.program, self.db, cpu, c) for c in contexts]

        pool = mp.Pool(workers)
        
        try:
            for c in pool.map(_task, taskargs):
                results.append(c)
        except KeyboardInterrupt:
            pool.terminate()
        except:
            pool.terminate()
            raise
        
        return results
    
    
    
if __name__ == '__main__':
    
    AppRunner().run()    
