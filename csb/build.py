"""
CSB build related tools and programs. When executed as a program,
this module will run the CSB Build Console.
"""

import os
import sys
import imp
import shutil
import unittest

from csb.pyutils import Shell

ROOT = 'csb'


class Console(object):
    """
    CSB Build Bot. Run with -h for usage.
        
    This is a CLI build program, which can be used to build, test and package
    the entire csb project.
    
    @param input: input source code directory
    @type input: str
    @param output: build output directory
    @type output: str
    @param verbosity: verbosity level
    @type verbosity: int
    """
    
    USAGE = r"""
CSB Build Console. Usage:
 
     python {0} -i input -o output [-v verbosity] [-h]
     
Options:
      -i  input      Source code folder. That should be the output of a fresh
                     'svn checkout'.
      -o  output     Build output directory
      -v  verbosity  Verbosity level, default is 1        
      -h, --help     Display this help
    """    
    
    def __init__(self, input='.', output='.', verbosity=1):
        
        self._program = None
        self._input = None
        self._output = None
        self._temp = None
        self._docs = None         
        self._apidocs = None
        self._root = None
        self._verbosity = None
                
        self.program = sys.argv[0]
        self.input = input
        self.output = output
        self.verbosity = verbosity        
        
    @property
    def program(self):
        return self._program
    @program.setter
    def program(self, value):
        self._program = os.path.basename(value)        

    @property
    def input(self):
        return self._input
    @input.setter
    def input(self, value):
        #value = os.path.dirname(value)
        input = os.path.abspath(value)
        if not os.path.isdir(input):
            raise IOError('Input direcotry not found: {0}'.format(input))
        root = os.path.join(input, ROOT)
        if not os.path.isdir(root):
            raise IOError('Root package directory not found: {0}'.format(root))            
        self._input = input
        
    @property
    def output(self):
        return self._output
    @output.setter
    def output(self, value):
        #value = os.path.dirname(value)        
        self._output = os.path.abspath(value)
        self._temp = os.path.join(self._output, 'build')
        self._docs = os.path.join(self._temp, 'docs')         
        self._apidocs = os.path.join(self._docs, 'api')
        self._root = os.path.join(self._temp, ROOT)    

    @property
    def verbosity(self):
        return self._verbosity
    @verbosity.setter
    def verbosity(self, value):
        self._verbosity = int(value)                     

    def build(self):
        """
        Run the console. 
        
        Since test packages do extensive imports, reload ROOT and
        sub-packages from the output folder. Also modify sys.path temporarily
        so that tests are not executed from the wrong CSB package: the output
        folder is added at the B{beginning} of sys.path.
        """
        self.log('\n# Build started.\n')
        self._init()        
        
        self._revision()    
                
        self.log('\n# Configuring Module Environment...')
        syspath = sys.path
        self._environment(ROOT, [self._temp] + syspath, self._temp)
        
    
        self._test()     
        self._doc()
        version = self._package()

        self.log('\n# Restoring the Normal Module Environment...', level=2)
        self._environment(ROOT, syspath)       
        
        self.log('\n# Done ({0} {1}).\n'.format(ROOT, version))             

    def log(self, message, level=1, ending='\n'):

        if self._verbosity >= level:             
            sys.stdout.write(message)
            sys.stdout.write(ending)
        sys.stdout.flush()
                  
    def _environment(self, root, syspath, prefix=None):
        """
        Reload C{root} and sub-packages using a new C{syspath}.
        
        @param root: root namespace
        @type root: str
        @param syspath: new sys.path
        @type syspath: list
        @param prefix: if defined, make sure that modules have really been reloaded
                       from this path - sanity check
        @type prefix: str
        """
        
        sys.path = list(syspath)
        
        for m in sys.modules.keys():
            if (m == root or m.startswith(root + '.')) and sys.modules[m]:
                
                self.log('{0:70}'.format(sys.modules[m].__file__), level=2, ending='')     
                reload(sys.modules[m])
                
                if prefix and not sys.modules[m].__file__.startswith(prefix):
                    raise ValueError('Module {0} not loaded from {1}'.format(m, prefix))
                
                self.log(' reloaded from {0}'.format(sys.modules[m].__file__), level=2)          

    def _init(self):
        """
        Collect all required stuff in the output folder.
        """
        self.log('# Preparing the file system...')
                
        if not os.path.exists(self._output):
            self.log('Creating output directory {0}'.format(self._output), level=2)
            os.mkdir(self._output)

        if os.path.exists(self._temp):
            self.log('Deleting existing temp directory {0}'.format(self._temp), level=2)         
            shutil.rmtree(self._temp)
                    
        self.log('Copying the source tree to temp directory {0}'.format(self._temp), level=2)            
        shutil.copytree(self._input, self._temp)
                    
        if os.path.exists(self._apidocs):
            self.log('Deleting existing API docs directory {0}'.format(self._apidocs), level=2)            
            shutil.rmtree(self._apidocs)
        if not os.path.isdir(self._docs):
            self.log('Creating docs directory {0}'.format(self._docs), level=2)
            os.mkdir(self._docs)            
        self.log('Creating API docs directory {0}'.format(self._apidocs), level=2)                        
        os.mkdir(self._apidocs)
        
    def _revision(self):
        """
        Write the actual revision number to L{ROOT}.__version__
        """
        self.log('\n# Setting the most recent Revision Number...')        
        root = os.path.join(self._root, '__init__.py')
        
        self.log('Retrieving revision number from {0}'.format(root), level=2)               
        rh = RevisionHandler(root)        
        revision = rh.read().maxrevision
        
        self.log('Writing back revision number {0}'.format(revision), level=2)        
        rh.write(revision, root)        
                
    def _test(self):
        """
        Run tests. Also make sure the current environment loads all modules from
        the output folder (assume a proper self._environment() call has already
        been done).
        """
        import csb.test
        assert csb.test.__file__.startswith(self._temp), 'csb.test not loaded from the output!'     #@UndefinedVariable
        
        self.log('\n# Running the Test Console...')
                
        builder = csb.test.AnyTestBuilder()
        suite = builder.loadTests(ROOT + '.test.cases.*')

        runner = unittest.TextTestRunner(verbosity=self.verbosity)
        result = runner.run(suite)
        if result.wasSuccessful():
            self.log('\n     PASSED all unit tests')
        else:
            self.log('\n     DID NOT PASS. The build might be broken :-(')                             
        
    def _doc(self):
        """
        Build documentation in the output folder.        
        """
        self.log('\n# Emulating ARGV for the Doc Builder...', level=2)        
        argv = sys.argv    
        sys.argv = ['epydoc', '--html', '-o', self._apidocs, '--name', ROOT.upper(), 
                    '--fail-on-error', '--fail-on-warning', '--fail-on-docstring-warning',
                    self._root]
        
        if self._verbosity > 0:                
            sys.argv.append('-v')
            
        import epydoc.cli

        self.log('\n# Running the Doc Console...')
        try:
            epydoc.cli.cli()
            sys.exit(0)
        except SystemExit as ex:
            if ex.code is 0:
                self.log('\n     PASSED all doc tests')
            else:
                if ex.code == 2:
                    self.log('\n     DID NOT PASS. The docs might be broken :-(')                    
                else:
                    self.log('\n     FAIL! Epydoc returned "#{0.code}: {0}"'.format(ex))

        self.log('\n# Restoring the previous ARGV...', level=2)    
        sys.argv = argv    
                        
    def _package(self):
        """
        Make package.
        """
        self.log('\n# Configuring CWD and ARGV for the Setup...', level=2)
        cwd = os.curdir
        os.chdir(self._temp)
                                
        if self._verbosity > 1:
            verbosity = '-v'
        else:
            verbosity = '-q'
        argv = sys.argv            
        sys.argv = ['setup.py', verbosity, 'sdist', '-d', self._output]        
            
        self.log('\n# Building Source Distribution...')
        try:       
            setup = imp.load_source('setupcsb', 'setup.py')
            setup.build()
        except SystemExit as ex:
            if ex.code is not 0:
                self.log('\n     FAIL! Setup returned: \n\n{0}\n'.format(ex))
            
        self.log('\n# Restoring the previous CWD and ARGV...', level=2)
        os.chdir(cwd)
        sys.argv = argv   
        
        return setup.VERSION 
        
    @staticmethod
    def exit(message=None, code=0, usage=True):
        
        if message:
            print message
        if usage:
            print Console.USAGE.format(sys.argv[0])    
                
        sys.exit(code)               

    @staticmethod
    def run(argv=None):
        
        if argv is None:
            argv = sys.argv[1:]
            
        input, output, verb = None, None, 1
            
        import getopt
        
        try:   
            options, dummy = getopt.getopt(argv, 'h:i:o:v:', ['help', 'input=', 'output=', 'verbosity='])
            
            for option, value in options:
                if option in('-h', '--help'):
                    Console.exit(message=None, code=0)
                if option in('-i', '--input'):
                    if not os.path.isdir(value):
                        Console.exit(message='E: Input directory not found "{0}".'.format(value), code=2)
                    input = value
                if option in('-o', '--output'):
                    if not os.path.isdir(value):
                        Console.exit(message='E: Output directory not found "{0}".'.format(value), code=3)
                    output = value
                if option in('-v', '--verbosity'):
                    try:
                        verb = int(value)
                    except ValueError:
                        Console.exit(message='E: Verbosity must be an integer.'.format(value), code=4)                    
        except getopt.GetoptError as oe:
            Console.exit(message='E: ' + str(oe), code=1)        

        if not (input and output):
            Console.exit(code=1, usage=True)
        else:
            # @todo: try..except -> exit
            Console(input, output, verbosity=verb).build()
            
class RevisionError(RuntimeError):
    
    def __init__(self, msg, code, cmd):
        
        super(RevisionError, self).__init__(msg)
        self.code = code
        self.cmd = cmd

class RevisionHandler(object):
    """
    Determines the current SVN revision number of a working copy.
    
    @param path: a local checkout path to be examined
    @type path: str
    @param svn: name of the svn program
    @type svn: str 
    """    
    
    def __init__(self, path, svn='svn'):
        
        self.path = None
        self.svn = None
        
        if os.path.exists(path):
            self.path = path
        else:
            raise IOError('Path not found: {0}'.format(path))
        if Shell.run([svn, 'help']).code is 0:
            self.svn = svn
        else:
            raise RevisionError('SVN probe failed', None, None)   
    
    def read(self):
        """
        Return the current revision information.
        @rtype: L{RevisionInfo}
        
        @todo: we can easily extend the svn output parser to grab more attributes,
               say URL, author, etc. RevisionInfo would also has to be extended
        """
        cmd = '{0.svn} info {0.path}'.format(self)
        revision = None
        maxrevision = None
        
        for line in self._run(cmd):
            if line.startswith('Revision:'):
                revision = int(line[9:] .strip())
                break
            
        for line in self._run(cmd + ' -R'):
            if line.startswith('Revision:'):
                rev = int(line[9:] .strip())
                if rev > maxrevision:
                    maxrevision = rev
        
        return RevisionInfo(self.path, revision, maxrevision)
    
    def write(self, revision, sourcefile):
        """
        Finalize the __version__ = major.minor.micro.{revision} tag.
        Overwrite C{sourcefile} in place by substituting the {revision} macro.
        
        @param revision: revision number to write to the source file.
        @type revision: int
        @param sourcefile: python source file with a __version__ tag, typicaly
                           "csb/__init__.py"
        @type sourcefile: str
        """
        content = open(sourcefile).read()
        content = content.format(revision=revision)
        
        with open(sourcefile, 'w') as src:
            src.write(content)            
    
    def _run(self, cmd):
        
        si = Shell.run(cmd)
        if si.code > 0:
            raise RevisionError('SVN failed ({0.code}): {0.stderr}'.format(si), si.code, si.cmd)
        
        return si.stdout.splitlines()
                 
class RevisionInfo(object):
    
    def __init__(self, item, revision, maxrevision):
        
        self.item = item
        self.revision = revision
        self.maxrevision = maxrevision
        
def main():
    Console.run()
        

if __name__ == '__main__':
    
    main()
