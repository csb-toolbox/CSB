"""
CSB is a high-level, object-oriented library used to solve problems in the
field of Computational Structural Biology.


Introduction
============

The library is composed of a set of highly branched python packages
(namespaces). Some of the packages are meant to be directly used by
the clients (core library), while others are utility modules and take
part in the development of the library: 

    1. Core class library -- object-oriented, granular, with an emphasis
       on design and clean interfaces. A Sequence is not a string, and a
       Structure is not a dict or list. Naming conventions matter.
       
    2. Application framework -- executable console applications
       ("protocols"), which consume objects from the core library.
       The framework ensures that each CSB application is also reusable
       and can be instantiated as a regular python object without any
       ugly side effects (sys.exit() and friends). See L{csb.apps} for more
       details. 
       
    3. Test framework -- built on top of the standard unittest as a thin
       wrapping layer. Provides some sugar like transparent management of
       test data files, and modular test execution. L{csb.test} will give
       you all the details. 

The core library is roughly composed of:

    - bioinformatics API: L{csb.bio}, which includes stuff like
      L{csb.bio.io}, L{csb.bio.structure}, L{csb.bio.sequence},
      L{csb.bio.hmm}
    
    - statistics API: L{csb.statistics}, L{csb.numeric}
    
    - utilities - L{csb.io}, L{csb.core}


Getting started
===============
    
Perhaps one of the most frequently used parts of the library is the
L{csb.bio.structure} module, which provides the L{Structure}, L{Chain},
L{Residue} and L{Atom} objects. You could easily build a L{Structure}
from scratch, but a far more common scenario is parsing a structure from
a PDB file using one of the L{AbstractStructureParser}s. All bio IO
objects, including the StructureParser factory, are defined in
L{csb.bio.io} and sub-packages:

    >>> from csb.bio.io.wwpdb import StructureParser
    >>> p = StructureParser("/some/file/pdb1x80.ent")
    >>> s = p.parse_structure()
    >>> print(s)
    <Structure: 1x80, 2 chains>
    
The last statement will return a L{csb.bio.structure.Structure} instance,
which is a composite hierarchical object:

    >>> for chain_id in s.chains:
            chain = s.chains[chain_id]
            for residue in chain.residues:
                for atom_id in residue.atoms:
                    atom = residue.atoms[atom_id]
                    print(atom.vector)

Some of the inner objects in this hierarchy behave just like dictionaries
(but are not):

    >>> s.chains['A']        # access chain A by ID
    <Chain A: Protein>
    >>> s['A']               # the same
    <Chain A: Protein>
    
Others behave like collections:

    >>> chain.residues[10]               # 1-based access to the residues in the chain
    <ProteinResidue [10]: PRO 10>
    >>> chain[10]                        # 0-based, list-like access
    <ProteinResidue [11]: GLY 11>
    
But all entities are iterable because they inherit the C{items} iterator
from L{AbstractEntity}. The above loop can be shortened:

    >>> for chain in s.items:
            for residue in chain.items:
                for atom in residue.items:
                    print(atom.vector)
                    
or even more:

    >>> from csb.bio.structure import Atom
    >>> for atom in s.components(klass=Atom):
            print(atom.vector)

You may also be interested in extracting a sub-chain from this structure:

    >>> s.chains['B'].subregion(3, 20)    # from positions 3 to 20, inclusive
    <Chain B: Protein>
    
or modifying it in some way, for example, in order to append a new residue,
try:

    >>> from csb.bio.structure import ProteinResidue
    >>> from csb.bio.sequence import ProteinAlphabet
    >>> residue = ProteinResidue(401, ProteinAlphabet.ALA)
    >>> s.chains['A'].residues.append(residue)
    
Finally, you would probably want to save your structure back to a PDB file:

    >>> s.to_pdb('/some/file/name.pdb')    


Where to go from here
=====================

If you want to dive into statistics, you could peek inside L{csb.statistics}
and its sub-packages. For example, L{csb.statistics.pdf} contains a collection
of L{probability density objects<csb.statistics.pdf.AbstractDensity>},
like L{Gaussian<csb.statistics.pdf.Normal>} or L{Gamma<csb.statistics.pdf.Gamma>}.

But chances are you would first like to try reading some files, so you could
start exploring L{csb.bio.io} right now. As we have already seen,
L{csb.bio.io.wwpdb} provides PDB L{Structure<csb.bio.structure.Structure>}
parsers, for example L{csb.bio.io.wwpdb.RegularStructureParser} and
L{csb.bio.io.wwpdb.LegacyStructureParser}.

L{csb.bio.io.fasta} is all about reading FASTA
L{Sequence<csb.bio.sequence.AbstractSequence>}s and
L{SequenceAlignment<csb.bio.sequence.AbstractAlignment>}s. Be sure to check out 
L{csb.bio.io.fasta.SequenceParser}, L{csb.bio.io.fasta.SequenceAlignmentReader}
and L{csb.bio.io.fasta.StructureAlignmentFactory}.

If you are working with HHpred (L{ProfileHMM<csb.bio.hmm.ProfileHMM>}s,
L{HHpredHit<csb.bio.hmm.HHpredHit>}s), then L{csb.bio.io.hhpred} is for you.
This package provides L{csb.bio.io.hhpred.HHProfileParser} and
L{csb.bio.io.hhpred.HHOutputParser}, which are used to read *.hhm and *.hhr
files.

Finally, if you want to make some nice plots with matplotlib, you may like the
clean object-oriented interface of our L{Chart<csb.io.plots.Chart>}. See
L{csb.io.plots} and maybe also L{csb.io.tsv} to get started.


Development
===========

When contributing code to CSB, please take into account the following:

    1. New features or bug fixes should always be accompanied by test cases.
       Also, always run the complete test suite before committing. For more
       details on this topic, see L{csb.test}.
       
    2. The source code of CSB must be cross-platform and cross-interpreter
       compatible. L{csb.core} and L{csb.io} will give you all necessary
       details on how to use the CSB compatibility layer.


License
=======

CSB is open source and distributed under OSI-approved MIT license::

    Copyright (c) 2012 Michael Habeck
    
    Permission is hereby granted, free of charge, to any person obtaining
    a copy of this software and associated documentation files (the
    "Software"), to deal in the Software without restriction, including
    without limitation the rights to use, copy, modify, merge, publish,
    distribute, sublicense, and/or sell copies of the Software, and to
    permit persons to whom the Software is furnished to do so, subject to
    the following conditions:
    
    The above copyright notice and this permission notice shall be
    included in all copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
    IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
    CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
    TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    
"""

__version__ = '1.0.0.{revision}'


class Version(object):
    """
    CSB version number.
    """
    
    def __init__(self):
        
        version = __version__.split('.')
        
        if not len(version) in (3, 4):
            raise ValueError(version)
        
        self._package = __name__ 
        
        self._major = version[0]
        self._minor = version[1]
        self._micro = version[2]
        self._revision = None
        
        if len(version) == 4:
            self._revision = version[3]

    def __str__(self):
        return self.short
    
    def __repr__(self):
        return '{0.package} {0.full}'.format(self)
    
    @property
    def major(self):
        """
        Major version (huge, incompatible changes)
        @rtype: int
        """
        return int(self._major)  
     
    @property
    def minor(self):
        """
        Minor version (significant, but compatible changes)
        @rtype: int
        """        
        return int(self._minor)
    
    @property
    def micro(self):
        """
        Micro version (bug fixes and small enhancements)
        @rtype: int
        """        
        return int(self._micro)  
    
    @property
    def revision(self):
        """
        Build number (exact repository revision number)
        @rtype: int
        """          
        try:
            return int(self._revision)
        except:
            return self._revision
    
    @property
    def short(self):
        """
        Canonical three-part version number.
        """
        return '{0.major}.{0.minor}.{0.micro}'.format(self)
    
    @property
    def full(self):
        """
        Full version, including the repository revision number.
        """        
        return '{0.major}.{0.minor}.{0.micro}.{0.revision}'.format(self)
    
    @property
    def package(self):
        return self._package

