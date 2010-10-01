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


def sourcetree(root='csb', junk=('.svn',)):
    """
    Since dustutils requires to HARD CODE the entire package hierarchy here,
    we need this function to load the source tree structure dynamically.
    
    @note: even if some directories are not packages, they will be loaded.
           However, distutils will print a warning of a missing __init__.py
           and then deal with it normally. That is, distutils will fix it.
    
    @param root: name of the root package; this is 'csb'. Must be relative!
    @type root: str
    @param junk: skip those directories
    @type junk: tuple
    
    @return: a list of "package" names
    @rtype: tuple  
    """
    junk = set(junk)
    items = []
    
    if root.startswith(os.path.sep):
        raise ValueError('root must be a relative path')
    elif not os.path.isdir(os.path.join('.', root)):
        raise ValueError('root package "{0}" not found in cwd.'.format(root))        
    
    for entry in os.walk(root):    
        item = entry[0]
        parts = set(item.split(os.path.sep))
        if junk.isdisjoint(parts):
            items.append(item)
    
    return tuple(items)

def build():
    setup(
          name=NAME,
          packages=sourcetree(NAME, JUNK),
          version=VERSION,
          author=AUTHOR,
          author_email=EMAIL,
          url=URL,
          description=SUMMARY,
          long_description=DESCRIPTION,
          license=LICENSE )



if __name__ == '__main__':
    
    build()