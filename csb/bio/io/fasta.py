"""
FASTA format sequence parsers.
"""

import cStringIO

import csb.io
import csb.pyutils

from abc import ABCMeta, abstractmethod

from csb.bio.sequence import SequenceTypes, SequenceAlphabets, AlignmentFormats, SequenceError
from csb.bio.sequence import SequenceAlignment, StructureAlignment, A3MAlignment
from csb.bio.sequence import SequenceCollection, AbstractSequence, Sequence, RichSequence, ChainSequence


class BaseSequenceParser(object):
    """
    FASTA parser template. Subclasses must implement the way FASTA strings are
    handled by overriding C{BaseSequenceParser.read_sequence}.
    
    @param product: sequence product factory (an L{AbstractSequence} subclass)
    @type product: type
    @param product_type: default L{SequenceTypes} member for the products
    @type product_type: L{EnumItem}
    """
    __metaclass__ = ABCMeta

    def __init__(self, product=Sequence, product_type=SequenceTypes.Protein):
        
        self._product = None
        
        if not issubclass(product, AbstractSequence):
            raise TypeError(product)
        
        if not product_type.enum is SequenceTypes:
            raise TypeError(product_type)
        
        self._product = product
        self._type = product_type
        
    @property
    def product_factory(self):
        return self._product
    
    @property
    def product_type(self):
        return self._type
        
    @abstractmethod
    def read_sequence(self, string):
        """
        Parse a single FASTA string
        
        @return: a new sequence, created with L{BaseSequenceParser.product_factory}
        @rtype: L{AbstractSequence}
        """
        pass

    def parse_string(self, fasta_string):
        """
        Read FASTA sequences from an (m)FASTA-formatted string

        @param fasta_string: FASTA string to parse
        @type fasta_string: str

        @return: a list of L{Sequence}s
        @rtype: L{SequenceCollection}
        """

        stream = cStringIO.StringIO()
        stream.write(fasta_string)

        return self.parse_file(stream)

    def parse_file(self, fasta_file):
        """
        Read FASTA sequences from a (m)FASTA file

        @param fasta_file: input FASTA file name or opened stream
        @type fasta_file: str, file

        @return: a list of L{Sequence}s
        @rtype: L{SequenceCollection}
        """
        if isinstance(fasta_file, basestring):
            stream = open(fasta_file)
        else:
            stream = fasta_file

        seqs = []
        reader = csb.io.EntryReader(stream, AbstractSequence.DELIMITER, None)

        for entry in reader.entries():
            seqs.append(self.read_sequence(entry))

        return SequenceCollection(seqs)

    def read(self, fasta_file):
        """
        Read FASTA sequences from an (m)FASTA file.

        @param fasta_file: input FASTA file name or opened stream
        @type fasta_file: str, file

        @return: efficient cursor over all L{Sequence}s (parse on demand)
        @rtype: iterator
        """
        if isinstance(fasta_file, basestring):
            stream = open(fasta_file)
        else:
            stream = fasta_file

        reader = csb.io.EntryReader(stream, AbstractSequence.DELIMITER, None)

        for entry in reader.entries():
            yield self.read_sequence(entry)
    
class SequenceParser(BaseSequenceParser):
    """
    Standard FASTA parser. See L{BaseSequenceParser} for details.
    """
    
    def read_sequence(self, string):

        lines = string.strip().splitlines()
        
        if not lines[0].startswith(AbstractSequence.DELIMITER):
            lines = [''] + lines
        if len(lines) < 2:
            raise ValueError('Empty FASTA entry')

        header = lines[0]
        id = header[1:].split()[0]
        sequence = ''.join(lines[1:])

        return self.product_factory(id, header, sequence, self.product_type)
    
class PDBSequenceParser(SequenceParser):
    """
    PDB FASTA parser. Reads the PDB ID and sequence type from the header.
    See L{BaseSequenceParser} for more details.
    """
        
    def read_sequence(self, string):

        seq = super(PDBSequenceParser, self).read_sequence(string)

        if not (seq.header and seq.id) or not (len(seq.id) in(5, 6) and seq.header.find('mol:') != -1):
            raise ValueError('Does not look like a PDB header: {0}'.format(seq.header))

        seq.id = seq.id.replace('_', '')
        stype = seq.header.partition('mol:')[2].partition(' ')[0]
        try:
            seq.type = csb.pyutils.Enum.parse(SequenceTypes, stype)
        except csb.pyutils.EnumValueError:
            seq.type = SequenceTypes.Unknown

        return seq
    

class A3MSequenceIterator(object):
    
    def __init__(self, sequences, insertion):
        
        self._temp = []
        self._insertion = insertion
        
        for sequence in sequences:
            
            row = list(sequence.sequence)
            row.reverse()
            self._temp.append(row)  
    
    def next(self):
        
        sequences = self._temp

        column = [ ]
        has_insertion = False

        for sequence in sequences:
            if len(sequence) > 0 and sequence[-1].islower():
                has_insertion = True

        for sequence in sequences:
            try:
                if has_insertion and not sequence[-1].islower():
                    column.append(self._insertion)
                else:
                    column.append(sequence.pop())
            except IndexError:
                column.append('')

        if not any(column):
            raise StopIteration()

        return column  
    
    def __iter__(self):
        return self    
    
class SequenceAlignmentReader(object):
    """
    Sequence alignment parser.
    
    @param product_type: default L{SequenceTypes} member for the sequence products
    @type product_type: L{EnumItem}
    @param strict: if True, raise exception on duplicate sequence identifiers.
                   See L{csb.bio.sequence.AbstractAlignment} for details
    @type strict: bool    
    """
    
    def __init__(self, product_type=SequenceTypes.Protein, strict=True):
        
        if not product_type.enum is SequenceTypes:
            raise TypeError(product_type)
        
        self._type = product_type
        self._strict = bool(strict)

    @property
    def product_type(self):
        return self._type

    @property
    def strict(self):
        return self._strict
            
    def read_fasta(self, string):
        """
        Parse an alignment in multi-FASTA format.
        
        @param string: alignment string
        @type string: str
        
        @rtype: L{SequenceAlignment}
        """
        
        parser = SequenceParser(RichSequence, self.product_type)
        sequences = parser.parse_string(string)
        
        return SequenceAlignment(sequences, strict=self.strict) 
    
    def read_a3m(self, string):
        """
        Parse an alignment in A3M format.
        
        @param string: alignment string
        @type string: str
        
        @rtype: L{A3MAlignment}
        """
        alphabet = SequenceAlphabets.get(self.product_type)

        # parse all "mis-aligned" sequences as case-sensitive strings        
        parser = SequenceParser(Sequence, self.product_type)
        sequences = parser.parse_string(string)
        
        # storage for expanded sequences
        s = []
        
        for dummy in sequences:
            s.append([])
             
        # expand all sequences with insertion characters and make them equal length
        for column in A3MSequenceIterator(sequences, str(alphabet.INSERTION)):
            for sn, char in enumerate(column):
                s[sn].append(char)
        
        # build normal sequence objects from the equalized sequence strings
        aligned_seqs = []
        
        for sn, seq in enumerate(sequences):
        
            sequence = RichSequence(seq.id, seq.header, s[sn], self.product_type)
            aligned_seqs.append(sequence)
            
        return A3MAlignment(aligned_seqs, strict=self.strict)
    
class StructureAlignmentFactory(object):
    """
    Protein structure alignment parser.
    
    In order to construct the structural alignment, this factory needs a PDB
    structure provider: an object, whose C{provider.get} method returns a 
    L{csb.bio.structute.Structure} for a given sequence identifier. Sequence
    identifiers on the other hand need to be split into 'accession number'
    and 'chain ID'. By default this is done using a standard PDB Entry ID
    factory, but clients are free to provide custom factories. An C{id_factory}
    must be a callable, which accepts a single string identifier and returns
    an EntryID object.
    
    @param provider: data source for all structures found in the alignment
    @type provider: L{csb.bio.io.wwpdb.StructureProvider}
    @param id_factory: callable factory, which transforms a sequence ID into
                       a L{csb.bio.io.wwpdb.EntryID} object. By default
                       this is L{csb.bio.io.wwpdb.EntryID.create}.
    @type id_factory: callable
    @param strict: if True, raise exception on duplicate sequence identifiers.
                   See L{csb.bio.sequence.AbstractAlignment} for details
    @type strict: bool    
    """
    
    def __init__(self, provider, id_factory=None, strict=True):

        from csb.bio.io.wwpdb import EntryID
        
        if id_factory is None:
            id_factory = EntryID.create
        if not hasattr(id_factory, '__call__'):
            raise TypeError(id_factory)
        
        if not hasattr(provider, 'get'):
            raise TypeError(provider)        
                
        self._type = SequenceTypes.Protein
        self._strict = bool(strict)
        self._provider = provider
        self._id_factory = id_factory

    @property
    def product_type(self):
        return self._type
    
    @property
    def provider(self):
        return self._provider  
    
    @property
    def id_factory(self):
        return self._id_factory

    @property
    def strict(self):
        return self._strict
            
    def make_alignment(self, string):
        """
        Build a protein structure alignment given a multi-FASTA string
        and the current structure C{provider}.
        
        @param string: alignment string
        @type string: str
        
        @rtype: L{SequenceAlignment}
        
        @raise SequenceError: when an aligned sequence is not a proper
                              subsequence of its respective source PDB chain
        @raise StructureNotFoundError: if C{provider} can't provide a structure
                                       for a given sequence ID
        @raise InvalidEntryIDError: if a given sequence ID cannot be parsed        
        """

        entries = []
        parser = SequenceParser(Sequence, self.product_type)
        sequences = parser.parse_string(string)
        
        for row in sequences:
            id = self.id_factory(row.id)
            chain = self.provider.get(id.accession).chains[id.chain]
            
            entry = self.make_entry(row, chain)
            entries.append(entry)
            
        return StructureAlignment(entries, strict=self.strict)
    
    def make_entry(self, row, chain):
        """
        Build a protein structure alignment entry, given a sequence alignment
        entry and its corresponding source PDB chain.
        
        @param row: sequence alignment entry (sequence with gaps)
        @type row: L{AbstractSequence}, L{SequenceAdapter}
        @param chain: source PDB chain
        @type chain: L{csb.bio.structure.Chain}
        
        @return: gapped chain sequence, containing cloned residues from the
                 source chain (except for the gaps)
        @rtype: L{ChainSequence}
        @raise SequenceError: when C{row} is not a proper subsequence of C{chain}        
        """
        offset = 1        
        residues = []
        sequence = row.strip().sequence.upper()
        
        try:
            start = chain.sequence.index(sequence) + 1
        except ValueError:
            raise SequenceError('{0}: not a subsequence of {1}'.format(row.id, chain.entry_id))
        
        for rinfo in row.residues:
            
            if rinfo.type == row.alphabet.GAP:
                residues.append(rinfo)
                continue
            else:
                rank = start + offset - 1
                assert chain.residues[rank].type == rinfo.type
                residues.append(chain.residues[rank].clone())
                offset += 1
                continue
            
        return ChainSequence(row.id, row.header, residues, chain.type)
        
    
class OutputBuilder(object):
    """
    Base sequence/alignment string format builder.
    
    @param output: destination stream, where the product is written.
    @type output: file
    @param headers: if False, omit headers
    @type headers: bool
    
    @note: File builders are not guaranteed to check the correctness of the
           product. It is assumed that the client of the builder knows what
           it is doing.   
    """
    __metaclass__ = ABCMeta
    _registry = {}    
        
    def __init__(self, output, headers=True):

        if not hasattr(output, 'write'):
            raise TypeError(output)     

        self._out = output
        self._headers = bool(headers)
        
    @staticmethod
    def create(format, *a, **k):
        """
        Create an output builder, which knows how to handle the specified
        sequence/alignment C{format}. Additional arguments are passed to the
        builder's constructor.
        
        @param format: L{AlignmentFormats} member
        @type format: L{EnumItem}
        @rtype: L{OutputBuilder}
        """
        if format not in OutputBuilder._registry:
            raise ValueError('Unhandled format: {0}'.format(format))
        
        klass = OutputBuilder._registry[format]
        return klass(*a, **k)

    @staticmethod
    def register(format, klass):
        """
        Register a new output builder.
        
        @param format: L{AlignmentFormats} member
        @type format: L{EnumItem}
        @param klass: builder class (L{OutputBuilder} sub-class)
        @type klass: type
        """    
        assert format not in OutputBuilder._registry
        assert issubclass(klass, OutputBuilder)
        
        OutputBuilder._registry[format] = klass

    @property
    def output(self):
        return self._out
    
    @property
    def headers(self):
        return self._headers

    def write(self, text):
        """
        Write a chunk of C{text} to the output stream.
        """
        self._out.write(text)   

    def writeline(self, text):
        """
        Write a chunk of C{text}, followed by a newline terminator.
        """
        self._out.write(text)
        self._out.write('\n')            
    
    @abstractmethod
    def add_sequence(self, sequence):
        """
        Format and append a new sequence to the product.
        @type sequence: L{AbstractSequence}    
        """
        pass
    
    def add_many(self, sequences):
        """
        Format and append a collection of L{AbstractSequence}s to the product.
        @type sequences: iterable of L{AbstractSequence}s
        """
        for s in sequences:
            self.add_sequence(s)
            
    @abstractmethod            
    def add_alignment(self, alignment):
        """
        Format and append an alignment to the product.
        @type alignment: L{AbstractAlignment} 
        """
        pass
    
    def add_separator(self, separator=''):
        """
        Append a sequence separator to the product.
        """
        self.writeline(separator)    

    def add_comment(self, text, comment='#', length=120):
        """
        Append a comment line to the product.
        
        @param text: comment text
        @type text: str
        @param comment: comment prefix
        @type comment: str
        @param length: maximal comment length
        @type length: int
        """
        for i in range(0, len(text), length):
            self.write(comment)
            self.write(' ')            
            self.write(text[i : i + length])
            self.write('\n')
            
class FASTAOutputBuilder(OutputBuilder):
    """
    Formats sequences as standard FASTA strings. See L{OutputBuilder}.
    """
    FORMAT = AlignmentFormats.FASTA
    
    def add_sequence(self, sequence):
        
        if self.headers:
            self.write(AbstractSequence.DELIMITER)
            self.writeline(sequence.header)

        insertion = str(sequence.alphabet.INSERTION)
        gap = str(sequence.alphabet.GAP)
        
        self.writeline(sequence.sequence.replace(insertion, gap))
        
    def add_alignment(self, alignment):
        
        self.add_many(alignment.rows)        
        
class A3MOutputBuilder(OutputBuilder):
    """
    Formats sequences as A3M strings. When appending an alignment, this builder
    will write all insertion-containing columns in lower case. Also, gap symbols
    are omitted if the respective columns contain insertions. 
    
    See L{OutputBuilder}.
    """
    FORMAT = AlignmentFormats.A3M
    
    def add_sequence(self, sequence):
        
        if self.headers:
            self.write(AbstractSequence.DELIMITER)
            self.writeline(sequence.header)

        self.writeline(sequence.sequence)
        
    def add_alignment(self, alignment):

        if isinstance(alignment, A3MAlignment):
            self._add_a3m(alignment)
        else:
            self._add_proper(alignment)
            
    def _add_a3m(self, alignment):
        
        for s in alignment.rows:
            
            if self.headers:
                self.write(AbstractSequence.DELIMITER)
                self.writeline(s.header)
                
            sequence = []
            
            for ci in s.columns:
                
                if ci.residue.type != s.alphabet.INSERTION:
                    char = str(ci.residue.type)
                    
                    if alignment.insertion_at(ci.column):
                        sequence.append(char.lower())
                    else:
                        sequence.append(char)
                else:
                    continue
                        
            self.writeline(''.join(sequence))
            
    def _add_proper(self, alignment):
        
        for s in alignment.rows:
            
            if self.headers:
                self.write(AbstractSequence.DELIMITER)
                self.writeline(s.header)
                
            master = alignment.rows[1]
            sequence = []
            
            for ci in s.columns:
                char = str(ci.residue.type)
                
                if master.columns[ci.column].residue.type == master.alphabet.GAP:
                    if ci.residue.type == s.alphabet.GAP:
                        continue
                    else:
                        sequence.append(char.lower())
                else:
                    sequence.append(char)
                        
            self.writeline(''.join(sequence))            
            
class PIROutputBuilder(OutputBuilder):
    """
    Formats sequences as PIR FASTA strings, recognized by Modeller.
    See L{OutputBuilder} for general alignment documentation.
    """
    FORMAT = AlignmentFormats.PIR
        
    def add_sequence(self, sequence):
        self._add(sequence, template=True)
        
    def add_alignment(self, alignment):
        
        for n, sequence in enumerate(alignment.rows, start=1):
            if n == 1:
                self.add_target(sequence)
            else:
                self.add_template(sequence)
            
        
    def add_target(self, sequence):
        self._add(sequence, template=False)
    
    def add_template(self, sequence):
        self._add(sequence, template=True)
    
    def _add(self, sequence, template=True):
                    
        if self.headers:
            
            if template:
                type = 'structure'
            else:
                type = 'sequence'
            
            id = sequence.id
            start = end = '.'

            if hasattr(sequence, 'start') and hasattr(sequence, 'end'):
                start = sequence.start
                end = sequence.end
                            
            header = 'P1;{0}\n{2}:{0}:{3}:{1}:{4}:{1}::::'.format(id[:-1], id[-1], type, start, end)
                            
            self.write(AbstractSequence.DELIMITER)
            self.writeline(header)

        insertion = str(sequence.alphabet.INSERTION)
        gap = str(sequence.alphabet.GAP)
        
        self.write(sequence.sequence.replace(insertion, gap))
        self.writeline('*')


# register builders
for klass in OutputBuilder.__subclasses__():
    OutputBuilder.register(klass.FORMAT, klass)
