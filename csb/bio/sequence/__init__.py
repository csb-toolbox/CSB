"""
Sequence and sequence alignment APIs.

This module defines the base interfaces for biological sequences and alignments:
L{AbstractSequence} and L{AbstractAlignment}. These are the central abstractions
here. This module provides also a number of useful enumerations, like L{SequenceTypes}
and L{SequenceAlphabets}.

Sequences
=========
L{AbstractSequence} has a number of implementations. These are of course interchangeable,
but have different intents and may differ significantly in performance. The standard
L{Sequence} implementation is what you are after if all you need is high performance
and efficient storage (e.g. when you are parsing big files). L{Sequence} objects store
their underlying sequences as strings. L{RichSequence}s on the other hand will store
their residues as L{ResidueInfo} objects, which have the same basic interface as the 
L{csb.bio.structure.Residue} objects. This of course comes at the expense of degraded
performance. A L{ChainSequence} is a special case of a rich sequence, whose residue
objects are I{actually} real L{csb.bio.structure.Residue}s.

Basic usage:

    >>> seq = RichSequence('id', 'desc', 'sequence', SequenceTypes.Protein)
    >>> seq.residues[1]
    <ResidueInfo [1]: SER>
    >>> seq.dump(sys.stdout)
    >desc
    SEQUENCE

See L{AbstractSequence} for details.    

Alignments
==========
L{AbstractAlignment} defines a table-like interface to access the data in an
alignment:

    >>> ali = SequenceAlignment.parse(">a\\nABC\\n>b\\nA-C")
    >>> ali[0, 0]
    <SequenceAlignment>   # a new alignment, constructed from row #1, column #1
    >>> ali[0, 1:3]
    <SequenceAlignment>   # a new alignment, constructed from row #1, columns #2..#3

which is just a shorthand for using the standard 1-based interface:

    >>> ali.rows[1]
    <AlignedSequenceAdapter: a, 3>                        # row #1 (first sequence)
    >>> ali.columns[1]
    (<ColumnInfo a [1]: ALA>, <ColumnInfo b [1]: ALA>)    # residues at column #1

See L{AbstractAlignment} for all details and more examples.

There are a number of L{AbstractAlignment} implementations defined here.
L{SequenceAlignment} is the default one, nothing surprising. L{A3MAlignment}
is a more special one: the first sequence in the alignment is a master sequence.
This alignment is usually used in the context of HHpred. More important is the
L{StructureAlignment}, which is an alignment of L{csb.bio.structure.Chain} objects.
The residues in every aligned sequence are really the L{csb.bio.structure.Residue}
objects taken from those chains.
"""

import re
import csb.core
import csb.io

from abc import ABCMeta, abstractmethod, abstractproperty


class AlignmentFormats(csb.core.enum):
    """
    Enumeration of multiple sequence alignment formats
    """
    A3M='a3m'; FASTA='fa'; PIR='pir'

class SequenceTypes(csb.core.enum):
    """
    Enumeration of sequence types
    """
    NucleicAcid='NA'; DNA='DNA'; RNA='RNA'; Protein='Protein'; Unknown='Unknown'    

class AlignmentTypes(csb.core.enum):
    """
    Enumeration of alignment strategies
    """
    Global='global'; Local='local'

class NucleicAlphabet(csb.core.enum):
    """
    Nucleic sequence alphabet
    """
    Adenine='A'; Cytosine='C'; Guanine='G'; Thymine='T'; Uracil='U'; Purine='R'; Pyrimidine='Y'; Ketone='K';
    Amino='M'; Strong='S'; Weak='W'; NotA='B'; NotC='D'; NotG='H'; NotT='V'; Any='N'; Masked='X'; GAP='-'; INSERTION='.';
    
class ProteinAlphabet(csb.core.enum):
    """
    Protein sequence alphabet
    """
    ALA='A'; ASX='B'; CYS='C'; ASP='D'; GLU='E'; PHE='F'; GLY='G'; HIS='H'; ILE='I'; LYS='K'; LEU='L'; MET='M'; ASN='N';
    PYL='O'; PRO='P'; GLN='Q'; ARG='R'; SER='S'; THR='T'; SEC='U'; VAL='V'; TRP='W'; TYR='Y'; GLX='Z'; UNK='X'; GAP='-';
    INSERTION='.'; STOP='*'
                                    
class StdProteinAlphabet(csb.core.enum):
    """
    Standard protein sequence alphabet
    """      
    ALA='A'; CYS='C'; ASP='D'; GLU='E'; PHE='F'; GLY='G'; HIS='H'; ILE='I'; LYS='K'; LEU='L'; MET='M'; ASN='N';
    PRO='P'; GLN='Q'; ARG='R'; SER='S'; THR='T';  VAL='V'; TRP='W'; TYR='Y'
    
class UnknownAlphabet(csb.core.enum):
    """
    Unknown sequence alphabet
    """  
    UNK='X'; GAP='-'; INSERTION='.'
   
class SequenceAlphabets(object):
    """
    Sequence alphabet enumerations.

    @note: This class is kept for backwards compatibility. The individual
           alphabet classes must be defined in the top level namespace,
           otherwise the new enum types cannot be pickled properly. 
    """
    Nucleic = NucleicAlphabet
    Protein = ProteinAlphabet
    StdProtein = StdProteinAlphabet
    Unknown = UnknownAlphabet
    
    MAP = { SequenceTypes.Protein: ProteinAlphabet,
            SequenceTypes.NucleicAcid: NucleicAlphabet,
            SequenceTypes.DNA: NucleicAlphabet,
            SequenceTypes.RNA: NucleicAlphabet,
            SequenceTypes.Unknown: UnknownAlphabet }
    
    ALL_ALPHABETS = set([ProteinAlphabet, NucleicAlphabet, UnknownAlphabet])

    assert set(MAP) == csb.core.Enum.members(SequenceTypes)
    
    @staticmethod
    def get(type):
        """
        Get the alphabet corresponding to the specified sequence C{type}
        @param type: a member of L{SequenceTypes}
        @type type: L{csb.core.EnumItem}
        @rtype: L{csb.core.enum} 
        """
        return SequenceAlphabets.MAP[type]  
    
    @staticmethod
    def contains(alphabet):
        """
        Return True if C{alphabet} is a sequence alphabet
        @type alphabet: L{csb.core.enum}
        @rtype: bool
        """
        return alphabet in SequenceAlphabets.ALL_ALPHABETS        


class SequenceError(ValueError):
    pass

class PositionError(IndexError):
    
    def __init__(self, index=None, start=1, end=None):
        
        if end == 0:
            start = 0
            
        self.index = index
        self.start = start
        self.end = end
        
        super(PositionError, self).__init__(index, start, end)
        
    def __str__(self):
        
        if self.index is not None:
            s = 'Position {0.index} is out of range [{0.start}, {0.end}]'
        else:
            s = 'Out of range [{0.start}, {0.end}]'
            
        return s.format(self)            
        
class SequencePositionError(PositionError):
    pass

class ColumnPositionError(PositionError):
    pass
     
class SequenceNotFoundError(KeyError):
    pass

class DuplicateSequenceError(KeyError):
    pass           
                
class ResidueInfo(object):
        
    def __init__(self, rank, type):
        
        self._type = None    
        self._rank = rank
        
        self.type = type
                
    @property
    def type(self):
        """
        Residue type - a member of any sequence alphabet
        @rtype: enum item
        """
        return self._type
    @type.setter
    def type(self, type):
        if not SequenceAlphabets.contains(type.enum):
            raise TypeError(type)
        self._type = type
        
    @property
    def rank(self):
        """
        Residue position (1-based)
        @rtype: int
        """
        return self._rank
    
    def __repr__(self):
        return '<{1} [{0.rank}]: {0.type!r}>'.format(self, self.__class__.__name__)

    
class ColumnInfo(object):
    
    def __init__(self, column, id, rank, residue):
        
        self.column = column
        self.id = id
        self.rank = rank
        self.residue = residue

    def __repr__(self):
        return '<{0.__class__.__name__} {0.id} [{0.column}]: {0.residue.type!r}>'.format(self)                
    
class SequenceIndexer(object):
    
    def __init__(self, container):
        self._container = container

    def __getitem__(self, rank):
        
        if not 1 <= rank <= self._container.length:
            raise SequencePositionError(rank, 1, self._container.length)
              
        return self._container._get(rank)
    
    def __iter__(self):
        return iter(self._container)
            
class UngappedSequenceIndexer(SequenceIndexer):

    def __getitem__(self, rank):
        try: 
            return self._container._get_ungapped(rank)
        except SequencePositionError:
            raise SequencePositionError(rank, 1)
    
    def __iter__(self):
        for c in self._container:
            if c.residue.type not in (self._container.alphabet.GAP, self._container.alphabet.INSERTION):
                yield c.residue

class ColumnIndexer(SequenceIndexer):
    
    def __getitem__(self, column):
        
        if not 1 <= column <= self._container.length:
            raise ColumnPositionError(column, 1, self._container.length)
                
        return self._container._get_column(column)
    

class SequenceCollection(csb.core.ReadOnlyCollectionContainer):
    """
    Represents a list of L{AbstractSequence}s.
    """

    def __init__(self, sequences):
        super(SequenceCollection, self).__init__(items=sequences, type=AbstractSequence)

    def to_fasta(self, output_file):
        """
        Dump the whole collection in mFASTA format.
        
        @param output_file: write the output to this file or stream
        @type output_file: str or stream
        """
        from csb.bio.io.fasta import FASTAOutputBuilder
            
        with csb.io.EntryWriter(output_file, close=False) as out:
            builder = FASTAOutputBuilder(out.stream, headers=True)
            
            for s in self:
                builder.add_sequence(s)        

        
class AbstractSequence(object):
    """
    Base abstract class for all Sequence objects.
    
    Provides 1-based access to the residues in the sequence via the
    sequence.residues property. The sequence object itself also behaves like
    a collection and provides 0-based access to its elements (residues).   
        
    @param id: FASTA ID of this sequence (e.g. accession number)
    @type id: str
    @param header: FASTA sequence header
    @type header: str
    @param residues: sequence residues
    @type residues: str or collection of L{ResidueInfo}
    @param type: a L{SequenceTypes} member (defaults to protein)
    @type type: L{EnumItem}
    """
     
    __metaclass__ = ABCMeta
    
    DELIMITER = '>'

    def __init__(self, id, header, residues, type=SequenceTypes.Unknown):

        self._id = None
        self._header = None
        self._residues = []
        self._type = None
          
        self.id = id
        self.header = header
        self.type = type
        
        for residue in residues:
            self._add(residue)
    
    def __getitem__(self, spec):
        
        if isinstance(spec, slice):
            spec = SliceHelper(spec, 0, self.length)
            return self.subregion(spec.start + 1, spec.stop)
        else:
            if not 0 <= spec < self.length:
                raise IndexError(spec)            
            return self._get(spec + 1)
    
    def __iter__(self):
        for index in range(self.length):
            yield self[index]
        
    @abstractmethod
    def _add(self, residue):
        """
        Append a C{residue} to the sequence.
        
        This is a hook method invoked internally for each residue during object
        construction. By implementing this method, sub-classes define how
        residues are attached to the sequence object.   
        """
        pass

    @abstractmethod
    def _get(self, rank):
        """
        Retrieve the sequence residue at the specified position (1-based, positive).
        
        This is a hook method which defines the actual behavior of the sequence
        residue indexer.
          
        @rtype: L{ResidueInfo}
        @raise SequencePositionError: when the supplied rank is out of range
        """
        pass
    
    def _factory(self, *a, **k):
        """
        Return a new sequence of the current L{AbstractSequence} sub-class.
        """
        return self.__class__(*a, **k)    

    def strip(self):
        """
        Remove all gaps and insertions from the sequence.
        
        @return: a new sequence instance, containing no gaps
        @rtype: L{AbstractSequence}
        """
        residues = [r for r in self._residues 
                    if r.type not in (self.alphabet.GAP, self.alphabet.INSERTION)]
        
        return self._factory(self.id, self.header, residues, self.type)
            
    def subregion(self, start, end):
        """
        Extract a subsequence, defined by [start, end]. The start and end
        positions are 1-based, inclusive.
        
        @param start: start position
        @type start: int
        @param end: end position
        @type end: int
        
        @return: a new sequence
        @rtype: L{AbstractSequence}
        
        @raise SequencePositionError: if start/end positions are out of range
        """
        positions = range(start, end + 1)
        return self.extract(positions)
    
    def extract(self, positions):
        """
        Extract a subsequence, defined by a list of 1-based positions.
        
        @param positions: positions to extract
        @type positions: tuple of int
        
        @return: a new sequence
        @rtype: L{AbstractSequence}
        
        @raise SequencePositionError: if any position is out of range
        """

        end = self.length
        residues = []
        
        for rank in sorted(set(positions)):
            if 1 <= rank <= end:
                residues.append(self._get(rank))
            else:
                raise SequencePositionError(rank, 1, end)
            
        return self._factory(self.id, self.header, residues, self.type)
    
    def dump(self, output_file):
        """
        Dump the sequence in FASTA format.
        
        @param output_file: write the output to this file or stream
        @type output_file: str or stream
        """
        from csb.bio.io.fasta import FASTAOutputBuilder
            
        with csb.io.EntryWriter(output_file, close=False) as out:
            FASTAOutputBuilder(out.stream, headers=True).add_sequence(self)
        
    @property
    def length(self):
        """
        Number of residues
        @rtype: int
        """
        return len(self._residues)
    
    @property
    def id(self):
        """
        Sequence identifier
        @rtype: str
        """        
        return self._id
    @id.setter
    def id(self, value):
        if value is not None:
            value = str(value).strip()
        self._id = value
            
    @property
    def header(self):
        """
        Sequence description
        @rtype: str
        """        
        return self._header
    @header.setter
    def header(self, value):
        if not value:
            value = 'sequence'       
        else:
            value = value.strip().lstrip(AbstractSequence.DELIMITER)
        self._header = value
    
    @property  
    def type(self):
        """
        Sequence type - a member of L{SequenceTypes}
        @rtype: enum item
        """
        return self._type
    @type.setter
    def type(self, value):
        if isinstance(value, csb.core.string):
            value = csb.core.Enum.parse(SequenceTypes, value)
        if value.enum is not SequenceTypes:
            raise TypeError(value) 
        self._type = value 

    @property
    def sequence(self): 
        """
        The actual sequence
        @rtype: str
        """
        return ''.join([str(r.type) for r in self._residues])
    
    @property
    def alphabet(self):
        """
        The sequence alphabet corresponding to the current sequence type
        @rtype: L{csb.core.enum}
        """
        return SequenceAlphabets.get(self._type)
    
    @property
    def residues(self):
        """
        Rank-based access to the underlying L{residues<csb.bio.sequence.ResidueInfo>}
        @rtype: L{SequenceIndexer}
        """
        return SequenceIndexer(self)

    def __len__(self):
        return self.length
    
    def __repr__(self):
        return '<{0.__class__.__name__}: {0.id}, {0.length} residues>'.format(self)
    
    def __str__(self):
        return '{0}{1.header}\n{1.sequence}'.format(AbstractSequence.DELIMITER, self)
        
class Sequence(AbstractSequence):
    """
    High-performance sequence object. The actual sequence is stored internally
    as a string. The indexer acts as a residue factory, which creates a new
    L{ResidueInfo} instance each time. 
    
    @note: This class was created with parsing large volumes of data in mind. This
           comes at the expense of degraded performance of the sequence indexer.
    
    @param id: FASTA ID of this sequence (e.g. accession number)
    @type id: str
    @param header: FASTA sequence header
    @type header: str
    @param residues: sequence string
    @type residues: str
    @param type: a L{SequenceTypes} member (defaults to protein)
    @type type: L{EnumItem}
    """
    
    def __init__(self, id, header, residues, type=SequenceTypes.Unknown):

        self._id = None
        self._header = None
        self._residues = ''
        self._type = None
          
        self.id = id
        self.header = header
        self.type = type        

        self._append(residues)
    
    def _append(self, string):
        # this seems to be the fastest method for sanitization and storage        
        self._residues += re.sub('([^\w\-\.])+', '', string)
        
    def _add(self, char):
        self._append(char)
            
    def _get(self, rank):        
        
        type = csb.core.Enum.parse(self.alphabet, self._residues[rank - 1])
        return ResidueInfo(rank, type)
    
    def strip(self):
        residues = self._residues.replace(
                        str(self.alphabet.GAP), '').replace(
                                        str(self.alphabet.INSERTION), '')
        return self._factory(self.id, self.header, residues, self.type)        

    def subregion(self, start, end):

        if not 1 <= start <= end <= self.length:
            raise SequencePositionError(None, 1, self.length)
                       
        residues = self._residues[start - 1 : end]
        return self._factory(self.id, self.header, residues, self.type)                

    def extract(self, positions):

        end = self.length
        residues = []
        
        for rank in sorted(set(positions)):
            if 1 <= rank <= end:
                residues.append(self._residues[rank - 1])
            else:
                raise SequencePositionError(rank, 1, end)
            
        return self._factory(self.id, self.header, ''.join(residues), self.type)
            
    @property
    def sequence(self):
        return self._residues

class RichSequence(AbstractSequence):
    """
    Sequence implementation, which converts the sequence into a list of
    L{ResidueInfo} objects. See L{AbstractSequence} for details.
    """
        
    def _add(self, residue):
        
        if hasattr(residue, 'rank') and hasattr(residue, 'type'):            
            self._residues.append(residue)
            
        else:
            if residue.isalpha() or residue in (self.alphabet.GAP, self.alphabet.INSERTION):
                
                type = csb.core.Enum.parse(self.alphabet, residue)
                rank = len(self._residues) + 1
                self._residues.append(ResidueInfo(rank, type))
            
    def _get(self, rank):
        return self._residues[rank - 1]

    @staticmethod
    def create(sequence):
        """
        Create a new L{RichSequence} from existing L{AbstractSequence}.
        
        @type sequence: L{AbstractSequence}
        @rtype: L{RichSequence}
        """
        return RichSequence(
                sequence.id, sequence.header, sequence.sequence, sequence.type)    

class ChainSequence(AbstractSequence):
    """
    Sequence view for L{csb.bio.structure.Chain} objects.
    See L{AbstractSequence} for details.
    """
            
    def _add(self, residue):
        
        if not (hasattr(residue, 'rank') and hasattr(residue, 'type')):
            raise TypeError(residue)
        else:
            self._residues.append(residue)
            
    def _get(self, rank):
        return self._residues[rank - 1]
    
    @staticmethod
    def create(chain):
        """
        Create a new L{ChainSequence} from existing L{Chain} instance.
        
        @type chain: L{csb.bio.structure.Chain}
        @rtype: L{ChainSequence}
        """        
        return ChainSequence(
                chain.entry_id, chain.header, chain.residues, chain.type)

    
class SequenceAdapter(object):
    """
    Base wrapper class for L{AbstractSequence} objects.
    Needs to be sub-classed (does not do anything special on its own).
    
    @param sequence: adaptee
    @type sequence: L{AbstractSequence}
    """
    
    def __init__(self, sequence):
        
        if not isinstance(sequence, AbstractSequence):
            raise TypeError(sequence)
        
        self._subject = sequence

    def __getitem__(self, i):
        return self._subject[i]
    
    def __iter__(self):
        return iter(self._subject)
                
    def __repr__(self):
        return '<{0.__class__.__name__}: {0.id}, {0.length}>'.format(self)        
    
    def __str__(self):
        return str(self._subject)
    
    def _add(self):
        raise NotImplementedError()
    
    def _get(self, rank):
        return self._subject._get(rank)
    
    def _factory(self, *a, **k):        
        return self.__class__(self._subject._factory(*a, **k))
    
    def strip(self):
        return self._subject.strip()
            
    def subregion(self, start, end):
        return self._subject.subregion(start, end)
    
    def extract(self, positions):
        return self._subject.extract(positions)    

    @property
    def id(self):
        return self._subject.id

    @property
    def length(self):
        return self._subject.length

    @property
    def type(self):
        return self._subject.type

    @property
    def header(self):
        return self._subject.header
    
    @property
    def sequence(self):
        return self._subject.sequence
                    
    @property
    def alphabet(self):
        return self._subject.alphabet

class AlignedSequenceAdapter(SequenceAdapter):
    """
    Adapter, which wraps a gapped L{AbstractSequence} object and makes it
    compatible with the MSA row/entry interface, expected by L{AbstractAlignment}.
    
    The C{adapter.residues} property operates with an L{UngappedSequenceIndexer},
    which provides a gap-free view of the underlying sequence.
    
    The C{adapter.columns} property operates with a standard L{ColumnIndexer},
    the same indexer which is used to provide the column view in multiple 
    alignments. Adapted sequences therefore act as alignment rows and allow for
    MSA-column-oriented indexing.
    
    @param sequence: adaptee
    @type sequence: L{AbstractSequence}    
    """

    def __init__(self, sequence):

        super(AlignedSequenceAdapter, self).__init__(sequence)
        
        self._fmap = {}
        self._rmap = {}
        rank = 0
        
        for column, residue in enumerate(sequence, start=1):
            
            if residue.type not in (self.alphabet.GAP, self.alphabet.INSERTION):
                rank += 1
                self._fmap[column] = rank                
                self._rmap[rank] = column
            else:
                self._fmap[column] = None

    def __getitem__(self, index):
        if not 0 <= index < self.length:
            raise IndexError(index)
        return self._get_column(index + 1)
    
    def __iter__(self):
        for c in sorted(self._fmap):
            yield self._get_column(c)
                    
    @property
    def columns(self):
        """
        Provides 1-based access to the respective columns in the MSA.
        @rtype: L{ColumnIndexer}
        """        
        return ColumnIndexer(self)

    @property
    def residues(self):
        """
        Provides 1-based access to the residues of the unaligned (ungapped)
        sequence.
        @rtype: L{UngappedSequenceIndexer} 
        """
        return UngappedSequenceIndexer(self)

    def _get_column(self, column):
        return ColumnInfo(
                column, self.id, self._fmap[column], self._subject.residues[column])
            
    def _get_ungapped(self, rank):
        return self._subject.residues[self._rmap[rank]]
    
    def map_residue(self, rank):
        """
        Return the MSA column number corresponding to the specified ungapped
        sequence C{rank}.
        
        @param rank: 1-based residue rank
        @type rank: int
        @rtype: int
        """
        return self._rmap[rank]
    
    def map_column(self, column):
        """
        Return the ungapped sequence rank corresponding to the specified MSA
        C{column} number.
        
        @param column: 1-based alignment column number
        @type column: int
        @rtype: int
        """        
        return self._fmap[column]    
    
class SliceHelper(object):
    
    def __init__(self, slice, start=0, stop=0):
        
        s, e, t = slice.start, slice.stop, slice.step
        
        if s is None:
            s = start
        if e is None:
            e = stop
        if t is None:
            t = 1
            
        for value in [s, e, t]:
            if value < 0:
                raise IndexError(value)
            
        self.start = s
        self.stop = e
        self.step = t            

class AlignmentRowsTable(csb.core.BaseDictionaryContainer):
    
    def __init__(self, container):
        
        super(AlignmentRowsTable, self).__init__()
        
        self._container = container
        self._map = {}
        
    def __getitem__(self, item):
        
        try:
            if isinstance(item, int):
                key = self._map[item]
            else:
                key = item
                
            return super(AlignmentRowsTable, self).__getitem__(key)
        
        except KeyError:
            raise SequenceNotFoundError(item)

    def _append(self, sequence):

        n = 0
        sequence_id = sequence.id
        
        while sequence_id in self:
            n += 1
            sequence_id = '{0}:A{1}'.format(sequence.id, n)

        super(AlignmentRowsTable, self)._append_item(sequence_id, sequence)
        self._map[self.length] = sequence_id
    
    def __iter__(self):
        for id in super(AlignmentRowsTable, self).__iter__():
            yield self[id]
        
    
class AbstractAlignment(object):
    """
    Base class for all alignment objects.
    
    Provides 1-based access to the alignment.rows and alignment.columns.
    Alignment rows can also be accessed by sequence ID. In addition, all
    alignments support 0-based slicing:
    
        >>> alignment[rows, columns]
        AbstractAlignment (sub-alignment)
        
    where
        - C{rows} can be a slice, tuple of row indexes or tuple of sequence IDs
        - columns can be a slice or tuple of column indexes
        
    For example:
    
        >>> alignment[:, 2:]
        AbstractAlignment     # all rows, columns [3, alignment.length]
        >>> alignment[(0, 'seqx'), (3, 5)]
        AbstractAlignment     # rows #1 and 'seq3', columns #4 and #5
        
    @param sequences: alignment entries (must have equal length)
    @type sequences: list of L{AbstractSequence}s
    @param strict: if True, raise {DuplicateSequenceError} when a duplicate ID
                   is found (default=True)
    @type strict: bool
    
    @note: if C{strict} is False and there are C{sequences} with redundant identifiers,
           those sequences will be added to the C{rows} collection with :An suffix,
           where n is a serial number. Therefore, rows['ID'] will return only one sequence,
           the first sequence with id=ID. All remaining sequences can be retrieved
           with C{rows['ID:A1']}, {rows['ID:A2']}, etc. However, the sequence objects will
           remain intact, e.g. {rows['ID:A1'].id} still returns 'ID' and not 'ID:A1'. 
    """
    
    __metaclass__ = ABCMeta
    
    def __init__(self, sequences, strict=True):
        
        self._length = None
        self._msa = AlignmentRowsTable(self)
        self._colview = ColumnIndexer(self)
        self._map = {}
        self._strict = bool(strict)
        
        self._construct(sequences)
            
    def __getitem__(self, spec):
        
        # The following code can hardly get more readable than that, sorry.
        # Don't even think of modifying this before there is a 100% unit test coverage 
        
        # 0. expand the input tuple: (rows/, columns/) => (rows, columns)
        if not isinstance(spec, tuple) or len(spec) not in (1, 2):
            raise TypeError('Invalid alignment slice expression')
        
        if len(spec) == 2:
            rowspec, colspec = spec
        else:
            rowspec, colspec = [spec, slice(None)]

        # 1. interpret the row slice: int, iter(int), iter(str) or slice(int) => list(int, 1-based)
        if isinstance(rowspec, slice):
            if isinstance(rowspec.start, csb.core.string) or isinstance(rowspec.stop, csb.core.string):
                raise TypeError("Invalid row slice: only indexes are supported")
            rowspec = SliceHelper(rowspec, 0, self.size)
            rows = range(rowspec.start + 1, rowspec.stop + 1)
        elif isinstance(rowspec, int):
            rows = [rowspec + 1]     
        elif csb.core.iterable(rowspec):
            try:
                rows = []
                for r in rowspec:
                    if isinstance(r, int):
                        rows.append(r + 1)
                    else:
                        rows.append(self._map[r])
            except KeyError as ke:
                raise KeyError('No such Sequence ID: {0!s}'.format(ke))
        else:
                raise TypeError('Unsupported row expression')            

        # 2. interpret the column slice: int, iter(int) or slice(int) => list(int, 1-based)            
        if isinstance(colspec, slice):
            colspec = SliceHelper(colspec, 0, self._length or 0)
            cols = range(colspec.start + 1, colspec.stop + 1)
        elif isinstance(colspec, int):
            cols = [colspec + 1] 
        elif csb.core.iterable(colspec):
            try:
                cols = [ c + 1 for c in colspec ]
            except:            
                raise TypeError('Unsupported column expression')
        else:
            raise TypeError('Unsupported column expression')
        
        # 3. some more checks
        if len(rows) == 0:
            raise ValueError("The expression returns zero rows")
        if len(cols) == 0:
            raise ValueError("The expression returns zero columns")
                
        # 4. we are done
        return self._extract(rows, cols)
    
    def _range(self, slice, start, end):
        
        s, e, t = slice.start, slice.end, slice.step
        
        if s is None:
            s = start
        if e is None:
            e = end
        if t is None:
            t = 1
            
        return range(s, e, t)
    
    def __iter__(self):
        for cn in range(1, self.length + 1):
            yield self._get_column(cn)
            
    @abstractmethod
    def _construct(self, sequences):
        """
        Hook method, called internally upon object construction. Subclasses
        define how the source alignment sequences are handled during alignment
        construction.
        
        @param sequences: alignment entries
        @type sequences: list of L{AbstractSequence}s
        """
        pass
    
    def _initialize(self, rep_sequence):
        """
        Hook method, which is used to initialize various alignment properties
        (such as length) from the first alignned sequence.
        """
        if rep_sequence.length == 0:
            raise SequenceError("Sequence '{0}' is empty".format(rep_sequence.id))
                
        assert self._length is None
        self._length = rep_sequence.length 
        
    def _factory(self, *a, **k):
        """
        Return a new sequence of the current L{AbstractAlignment} sub-class.
        """ 
        return self.__class__(*a, **k)
      
    def add(self, sequence):
        """
        Append a new sequence to the alignment.
        
        @type sequence: L{AbstractSequence}
        @raise SequenceError: if the new sequence is too short/long
        @raise DuplicateSequenceError: if a sequence with same ID already exists  
        """
        
        if self._msa.length == 0:
            self._initialize(sequence)

        if sequence.length != self._length:
            raise SequenceError('{0!r} is not of the expected length'.format(sequence))
        
        if self._strict and sequence.id in self._msa:
            raise DuplicateSequenceError(sequence.id)
        
        self._msa._append(AlignedSequenceAdapter(sequence))
        self._map[sequence.id] = self._msa.length

    @property
    def length(self):
        """
        Number of columns in the alignment
        @rtype: int
        """
        return self._length or 0

    @property
    def size(self):
        """
        Number of rows (sequences) in the alignment
        @rtype: int
        """        
        return self._msa.length    
            
    @property
    def rows(self):
        """
        1-based access to the alignment entries (sequences)
        @rtype: L{AlignmentRowsTable}
        """
        return self._msa 
        
    @property
    def columns(self):
        """
        1-based access to the alignment columns
        @rtype: L{ColumnIndexer}
        """
        return self._colview
    
    def gap_at(self, column):
        """
        Return True of C{column} contains at least one gap.
        @param column: column number, 1-based
        @type column: int
        
        @rtype: bool
        """
        
        for row in self._msa:
            if row.columns[column].residue.type == row.alphabet.GAP:
                return True
            
        return False        
    
    def _get_column(self, column):
        return tuple(row._get_column(column) for row in self.rows)
    
    def _extract(self, rows, cols):
        
        rows = set(rows)
        cols = set(cols)
                
        if not 1 <= min(rows) <= max(rows) <= self.size:
            raise IndexError('Row specification out of range')
                
        if not 1 <= min(cols) <= max(cols) <= self.length:
            raise IndexError('Column specification out of range')
        
        sequences = []
        
        for rn, row in enumerate(self.rows, start=1):
            if rn in rows:
                sequences.append(row.extract(cols))
                
        return self._factory(sequences, strict=self._strict)
    
    def subregion(self, start, end):
        """
        Extract a sub-alignment, ranging from C{start} to C{end} columns.
        
        @param start: starting column, 1-based
        @type start: int
        @param end: ending column, 1-based
        @type end: int
        
        @return: a new alignment of the current type
        @rtype: L{AbstractAlignment}
        
        @raise ColumnPositionError: if start/end is out of range 
        """
        if not 1 <= start <= end <= self.length:
            raise ColumnPositionError(None, 1, self.length)
        
        sequences = []
        
        for row in self.rows:
            sequences.append(row.subregion(start, end))
                
        return self._factory(sequences, strict=self._strict)
    
    def format(self, format=AlignmentFormats.FASTA, headers=True):
        """
        Format the alignment as a string.
        
        @param format: alignment format type, member of L{AlignmentFormats}
        @type format: L{EnumItem}
        @param headers: if False, omit headers
        @type headers: bool
        
        @rtype: str 
        """
        from csb.bio.io.fasta import OutputBuilder

        temp = csb.io.MemoryStream()
                
        try:            
            builder = OutputBuilder.create(format, temp, headers=headers)
            builder.add_alignment(self)
            
            return temp.getvalue()
        
        finally:
            temp.close()          

class SequenceAlignment(AbstractAlignment):
    """
    Multiple sequence alignment. See L{AbstractAlignment} for details.
    """
        
    def _construct(self, sequences):
        
        for sequence in sequences:
            self.add(sequence)
            
    @staticmethod
    def parse(string, strict=True):
        """
        Create a new L{SequenceAlignment} from an mFASTA string.
        
        @param string: MSA-formatted string
        @type string: str
        @param strict: see L{AbstractAlignment}
        @type strict: bool        
        
        @rtype: L{SequenceAlignment}
        """
        from csb.bio.io.fasta import SequenceAlignmentReader
        return SequenceAlignmentReader(strict=strict).read_fasta(string)
        
class StructureAlignment(AbstractAlignment):
    """
    Multiple structure alignment. Similar to a L{SequenceAlignment}, but
    the alignment holds the actual L{csb.bio.structure.ProteinResidue} objects,
    taken from the corresponding source L{csb.bio.structure.Chain}s.
    
    See L{AbstractAlignment} for details.
    """
        
    def _construct(self, sequences):
        
        for sequence in sequences:
            self.add(sequence)
            
    @staticmethod
    def parse(string, provider, id_factory=None, strict=True):
        """
        Create a new L{StructureAlignment} from an mFASTA string. See 
        L{csb.bio.io.fasta.StructureAlignmentFactory} for details. 
        
        @param string: MSA-formatted string
        @type string: str
        @param provider: data source for all structures found in the alignment
        @type provider: L{csb.bio.io.wwpdb.StructureProvider}
        @param strict: see L{AbstractAlignment}
        @type strict: bool
        @param id_factory: callable factory, which transforms a sequence ID into
                           a L{csb.bio.io.wwpdb.EntryID} object. By default
                           this is L{csb.bio.io.wwpdb.EntryID.create}. 
        @type id_factory: callable        
        @rtype: L{StructureAlignment}
        """
        from csb.bio.io.fasta import StructureAlignmentFactory
        
        factory = StructureAlignmentFactory(
                        provider, id_factory=id_factory, strict=strict)
        return factory.make_alignment(string)
    
class A3MAlignment(AbstractAlignment):
    """
    A specific type of multiple alignment, which provides some operations
    relative to a master sequence (the first entry in the alignment). 
    """
    
    def __init__(self, sequences, strict=True):

        self._master = None
        self._matches = 0
        self._insertions = set()
                
        super(A3MAlignment, self).__init__(sequences, strict=strict)

    def _initialize(self, rep_sequence):
        
        super(A3MAlignment, self)._initialize(rep_sequence)
        self._alphabet = rep_sequence.alphabet        
            
    def _construct(self, sequences):
        
        for sequence in sequences:
            
            self.add(sequence)
                    
            for rank, residue in enumerate(sequence, start=1):
                if residue.type == self._alphabet.INSERTION:
                    self._insertions.add(rank)

        if self.size == 0:
            raise SequenceError("At least one sequence is required") 
        
        self._master = list(self._msa)[0]
        self._matches = self._master.strip().length
    
    @property
    def master(self):
        """
        The master sequence
        @rtype: L{AbstractSequence}
        """
        return self._master
    
    def insertion_at(self, column):
        """
        Return True of C{column} contains at least one insertion.
        
        @param column: column number, 1-based
        @type column: int
        @rtype: bool
        """        
        return column in self._insertions
    
    def hmm_subregion(self, match_start, match_end):
        """
        Same as L{AbstractAlignment.subregion}, but start/end positions are
        ranks in the ungapped master sequence.
        """

        if not 1 <= match_start <= match_end <= self.matches:
            raise ColumnPositionError(None, 1, self.matches)
                
        start = self._master.map_residue(match_start)
        end = self._master.map_residue(match_end)
        
        return self.subregion(start, end)

    def format(self, format=AlignmentFormats.A3M, headers=True):
        return super(A3MAlignment, self).format(format, headers) 
    
    @property
    def matches(self):
        """
        Number of match states (residues in the ungapped master).
        @rtype: int
        """
        return self._matches
    
    @staticmethod
    def parse(string, strict=True):
        """
        Create a new L{A3MAlignment} from an A3M string.
        
        @param string: MSA-formatted string
        @type string: str
        @param strict: see L{AbstractAlignment}
        @type strict: bool
        
        @rtype: L{A3MAlignment}
        """        
        from csb.bio.io.fasta import SequenceAlignmentReader
        return SequenceAlignmentReader(strict=strict).read_a3m(string)

    