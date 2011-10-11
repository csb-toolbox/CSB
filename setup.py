import os
import imp

from distutils.core import setup
from csb.build import ROOT 

try:
    __doc__ = open('README.txt').read()
except IOError:
    pass


LOCATION        = os.path.abspath(os.path.dirname(__file__))
JUNK            = ['.svn']

NAME            = ROOT
VERSION         = imp.load_source('____csb', os.path.join(LOCATION, ROOT, '__init__.py')).__version__
AUTHOR          = "Michael Habeck et al."
EMAIL           = "michael.habeck@tuebingen.mpg.de"
URL             = "http://www.eb.tuebingen.mpg.de/departments/1-protein-evolution/michael-habeck/computational-structural-biology"
SUMMARY         = "Computational Structural Biology Toolbox"
DESCRIPTION     = __doc__
LICENSE         = 'BSD'
REQUIRES        = ['numpy']


def sourcetree(root='csb', junk=('.svn',)):
    """
    Since dustutils requires to HARD CODE the entire package hierarchy here,
    we need this function to load the source tree structure dynamically.
    
    @param root: name of the root package; this is 'csb'. Must be relative!
    @type root: str
    @param junk: skip those directories
    @type junk: tuple
    
    @return: a list of "package" names
    @rtype: list  
    """
    junk = set(junk)
    items = []

    curdir = os.path.abspath(os.curdir)
    cwd = os.path.dirname(__file__) or '.'
    os.chdir(cwd)
    
    if root.startswith(os.path.sep):
        raise ValueError('root must be a relative path')
    elif not os.path.isdir(os.path.join('.', root)):
        raise ValueError('root package "{0}" not found in {1}.'.format(root, cwd))        
    
    for entry in os.walk(root):
            
        directory = entry[0]
        parts = set(directory.split(os.path.sep))
            
        init = os.path.join(directory, '__init__.py')
        # take all packages: directories with __init__, except if a junk 
        # directory is found at any level in the tree
        if os.path.isfile(init) and junk.isdisjoint(parts):
            items.append(directory)

    os.chdir(curdir)
        
    return items


def datatree(package, dataroot, junk=('.svn',), mask='*.*'):
    """
    Since dustutils will crash if the data root folder contains any subfolders,
    we need this function to retrieve the data tree.

    @param package: root "package", containing a data folder. This is a 
                    relative path, e.g. "csb/test"
    @type package: str
    @param dataroot: name of the data root directory for C{package},
                     relative to C{package}
    @type dataroot: str
    @param junk: skip those directories
    @type junk: tuple
    
    @return: a list of all glob patterns with all subdirectories of the data
             root, including the root itself. The paths  returned are relative
             to C{package}  
    @rtype: list      
    """
    junk = set(junk)
    items = []

    curdir = os.path.abspath(os.curdir)
    cwd = os.path.dirname(__file__) or '.'
    os.chdir(cwd)
    
    if package.startswith(os.path.sep):
        raise ValueError('package must be a relative path')
    elif not os.path.isdir(os.path.join('.', package)):
        raise ValueError('package "{0}" not found in {1}.'.format(package, cwd))        

    os.chdir(package)
        
    for entry in os.walk(dataroot):
        
        directory = entry[0]
        parts = set(directory.split(os.path.sep))
        
        # take all directories, except if a junk dir is found at any level in the tree
        if junk.isdisjoint(parts):
            item = os.path.join(directory, mask)
            items.append(item)

    os.chdir(curdir)
        
    return items


def build():
    
    test = os.path.join('csb', 'test')
    nmr = os.path.join('csb', 'bio', 'nmr')
    
    return setup(
              name=NAME,
              packages=sourcetree(NAME, JUNK),
              package_data={
                                test: datatree(test, 'data',      JUNK, '*.*'),
                                nmr:  datatree(nmr,  'resources', JUNK, '*.*')
                            },
              version=VERSION,
              author=AUTHOR,
              author_email=EMAIL,
              url=URL,
              description=SUMMARY,
              long_description=DESCRIPTION,
              license=LICENSE,
              requires=REQUIRES )



if __name__ == '__main__':
    
    build()