"""
CSB build related tools and programs.

When executed as a program, this module will run the CSB Build Console and
build the source tree it belongs to. The source tree is added at the
B{beginning} of sys.path to make sure that all subsequent imports from the
Test and Doc consoles will import the right thing (think of multiple CSB
packages installed on the same server).

The Console can also be imported and instantiated as a regular Python class.
In this case the Console again builds the source tree it is part of, but
sys.path will remain intact. Therefore, the Console will assume that all
modules currently in memory, as well as those that can be subsequently imported
by the Console itself, belong to the same CSB package.

@note: The CSB build services no longer support the option to build external
       source trees.
@see: [CSB 0000038]
"""
from __future__ import print_function

import os
import sys
import getopt
import traceback
        
if os.path.basename(__file__) == '__init__.py':
    PARENT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
else:
    PARENT = os.path.abspath(os.path.dirname(__file__))

ROOT = 'csb'
SOURCETREE = os.path.abspath(os.path.join(PARENT, ".."))

if __name__ == '__main__':

    # make sure "import io" imports the built in module, not csb.io
    # io is required by tarfile    
    for path in sys.path:
        if path.startswith(SOURCETREE):
            sys.path.remove(path)
            
    import io
    assert hasattr(io, 'BufferedIOBase')   
            
    sys.path = [SOURCETREE] + sys.path


"""
It is now safe to import any modules  
"""
import imp
import shutil
import tarfile

import csb

from csb.io import Shell


class BuildTypes(object):
    """
    Enumeration of build types.
    """
    
    SOURCE = 'source'
    BINARY = 'binary'

    _du = { SOURCE: 'sdist', BINARY: 'bdist' }    
    
    @staticmethod
    def get(key):
        try:
            return BuildTypes._du[key]
        except KeyError:
            raise ValueError('Unhandled build type: {0}'.format(key))
    
        
class Console(object):
    """
    CSB Build Bot. Run with -h for usage.
    
    @param output: build output directory
    @type output: str
    @param verbosity: verbosity level
    @type verbosity: int
    
    @note: The build console automatically detects and builds the csb package
           it belongs to. You cannot build a different source tree with it.
           See the module documentation for more info.
    """
    
    PROGRAM = __file__
    
    USAGE = r"""
CSB Build Console: build, test and package the entire csb project.

Usage:
     python {program} -o output [-v verbosity] [-t type] [-h]
     
Options:
      -o  output     Build output directory
      -v  verbosity  Verbosity level, default is 1
      -t  type       Build type:
                        source - build source code distribution (default)
                        binary - build executable
      -h, --help     Display this help
    """    
    
    def __init__(self, output='.', verbosity=1, buildtype=BuildTypes.SOURCE):
        
        self._input = None
        self._output = None
        self._temp = None
        self._docs = None         
        self._apidocs = None
        self._root = None
        self._verbosity = None
        self._type = buildtype
        self._dist = BuildTypes.get(buildtype) 
            
        if os.path.join(SOURCETREE, ROOT) != PARENT:
            raise IOError('{0} must be a sub-package or sub-module of {1}'.format(__file__, ROOT))
        self._input = SOURCETREE
                
        self.output = output
        self.verbosity = verbosity
        
    @property
    def input(self):
        return self._input        
        
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
        """
        self.log('\n# Building package {0} from {1}\n'.format(ROOT, SOURCETREE))
        
        self._init()        
        v = self._revision()
        self._doc(v)
        self._test()
        vn = self._package()     
        
        self.log('\n# Done ({0}).\n'.format(vn.full))         

    def log(self, message, level=1, ending='\n'):

        if self._verbosity >= level:             
            sys.stdout.write(message)
            sys.stdout.write(ending)
        sys.stdout.flush()
        
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
        version = rh.write(revision, root)

        self.log('  This is {0}.__version__ {1}'.format(ROOT, version), level=1)
        csb.__version__ = version
        
        return version 
                
    def _test(self):
        """
        Run tests. Also make sure the current environment loads all modules from
        the input folder.
        """
        import csb.test
        assert csb.test.__file__.startswith(self._input), 'csb.test not loaded from the input!'     #@UndefinedVariable
        
        from csb.test import unittest
                        
        newdata = os.path.join(self._temp, ROOT, 'test', 'data')
        csb.test.Config.setDefaultDataRoot(newdata)

        self.log('\n# Updating all test pickles in {0} if necessary...'.format(newdata), level=2)        
        csb.test.Config().ensureDataConsistency()

        self.log('\n# Running the Test Console...')
                        
        builder = csb.test.AnyTestBuilder()
        suite = builder.loadTests(ROOT + '.test.cases.*')

        runner = unittest.TextTestRunner(stream=sys.stderr, verbosity=self.verbosity)
        result = runner.run(suite)
        if result.wasSuccessful():
            self.log('\n  Passed all unit tests')
        else:
            self.log('\n  DID NOT PASS: The build might be broken')                             
        
    def _doc(self, version):
        """
        Build documentation in the output folder.        
        """
        self.log('\n# Generating API documentation...')
        try:
            import epydoc.cli
        except ImportError:
            self.log('\n  Skipped: epydoc is missing')
            return
                
        self.log('\n# Emulating ARGV for the Doc Builder...', level=2)        
        argv = sys.argv    
        sys.argv = ['epydoc', '--html', '-o', self._apidocs,
                    '--name', '{0} v{1}'.format(ROOT.upper(), version),
                    '--no-private', '--introspect-only', '--exclude', 'csb.test.cases',
                    '--css', os.path.join(self._temp, 'epydoc.css'),
                    '--fail-on-error', '--fail-on-warning', '--fail-on-docstring-warning',
                    self._root]
        
        if self._verbosity > 0:                
            sys.argv.append('-v')
        
        try:
            epydoc.cli.cli()
            sys.exit(0)
        except SystemExit as ex:
            if ex.code is 0:
                self.log('\n  Passed all doc tests')
            else:
                if ex.code == 2:
                    self.log('\n  DID NOT PASS: The docs might be broken')                    
                else:
                    self.log('\n  FAIL: Epydoc returned "#{0.code}: {0}"'.format(ex))

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
        sys.argv = ['setup.py', verbosity, self._dist, '-d', self._output]        
            
        self.log('\n# Building {0} distribution...'.format(self._type))
        try:       
            setup = imp.load_source('setupcsb', 'setup.py')
            d = setup.build()
            version = setup.VERSION
            package = d.dist_files[0][2]
            
            if self._type == BuildTypes.BINARY:
                self._strip_source(package)
            
        except SystemExit as ex:
            if ex.code is not 0:
                package = 'FAIL'
                self.log('\n  FAIL: Setup returned: \n\n{0}\n'.format(ex))
            
        self.log('\n# Restoring the previous CWD and ARGV...', level=2)
        os.chdir(cwd)
        sys.argv = argv   

        self.log('  Packaged ' + package)   
        return version
    
    def _strip_source(self, package, source='*.py'):
        """
        Delete plain text source code files from the package.
        """    
        cwd = os.getcwd()
        
        try:  
            tmp = os.path.join(self.output, 'tmp')
            os.mkdir(tmp)
        
            self.log('\n# Entering {1} in order to delete .py files from {0}...'.format(package, tmp), level=2)        
            os.chdir(tmp)
                
            oldtar = tarfile.open(package, mode='r:gz')
            oldtar.extractall(tmp)
            oldtar.close()
            
            newtar = tarfile.open(package, mode='w:gz')            
    
            try:
                for i in os.walk('.'):
                    for fn in i[2]:
                        if fn.endswith('.py'):
                            module = os.path.join(i[0], fn);
                            if not os.path.isfile(module.replace('.py', '.pyc')):
                                raise ValueError('Missing bytecode for module {0}'.format(module))
                            else:                                          
                                os.remove(os.path.join(i[0], fn))
                
                for i in os.listdir('.'):
                    newtar.add(i)        
            finally:
                newtar.close()
                
        finally:
            self.log('\n# Restoring the previous CWD...', level=2)            
            os.chdir(cwd)
            if os.path.exists(tmp):
                shutil.rmtree(tmp)    
        
    @staticmethod
    def exit(message=None, code=0, usage=True):
        
        if message:
            print(message)
        if usage:
            print(Console.USAGE.format(program=Console.PROGRAM))    
                
        sys.exit(code)               

    @staticmethod
    def run(argv=None):
        
        if argv is None:
            argv = sys.argv[1:]
            
        output = None
        verb = 1
        buildtype = BuildTypes.SOURCE
                
        try:   
            options, dummy = getopt.getopt(argv, 'o:v:t:h', ['output=', 'verbosity=', 'type=', 'help'])
            
            for option, value in options:
                if option in('-h', '--help'):
                    Console.exit(message=None, code=0)
                if option in('-o', '--output'):
                    if not os.path.isdir(value):
                        Console.exit(message='E: Output directory not found "{0}".'.format(value), code=3)
                    output = value
                if option in('-v', '--verbosity'):
                    try:
                        verb = int(value)
                    except ValueError:
                        Console.exit(message='E: Verbosity must be an integer.', code=4)
                if option in('-t', '--type'):
                    if value not in [BuildTypes.SOURCE, BuildTypes.BINARY]:
                        Console.exit(message='E: Invalid build type "{0}".'.format(value), code=5)
                    buildtype = value                                         
        except getopt.GetoptError as oe:
            Console.exit(message='E: ' + str(oe), code=1)        

        if not output:
            Console.exit(code=1, usage=True)
        else:
            try:
                Console(output, verbosity=verb, buildtype=buildtype).build()
            except Exception as ex:
                msg = 'Unexpected Error: {0}\n\n{1}'.format(ex, traceback.format_exc())
                Console.exit(message=msg, code=99, usage=False)
                
            
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
        revision = 0
        maxrevision = 0
        
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
        @param sourcefile: python source file with a __version__ tag, typically
                           "csb/__init__.py"
        @type sourcefile: str
        
        @return: sourcefile.__version__
        """
        content = open(sourcefile).readlines()
        
        with open(sourcefile, 'w') as src:
            for line in content:
                if line.startswith('__version__'):
                    src.write(line.format(revision=revision))
                else:
                    src.write(line)

        self._delcache(sourcefile)
        return imp.load_source('____source', sourcefile).__version__      
    
    def _run(self, cmd):
        
        si = Shell.run(cmd)
        if si.code > 0:
            raise RevisionError('SVN failed ({0.code}): {0.stderr}'.format(si), si.code, si.cmd)
        
        return si.stdout.splitlines()
    
    def _delcache(self, sourcefile):
        
        compiled = os.path.splitext(sourcefile)[0] + '.pyc'
        if os.path.isfile(compiled):
            os.remove(compiled)
            
        pycache = os.path.join(os.path.dirname(compiled), '__pycache__')
        if os.path.isdir(pycache): 
            shutil.rmtree(pycache)
                
class RevisionInfo(object):
    
    def __init__(self, item, revision, maxrevision):
        
        self.item = item
        self.revision = revision
        self.maxrevision = maxrevision

        
def main():
    Console.run()
        

if __name__ == '__main__':
    
    main()
