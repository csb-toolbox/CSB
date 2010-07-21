import re
import csb.pyutils
import csb.io

AlignmentFormats = csb.pyutils.enum( A3M='a3m', FASTA='fa' )
"""
@var AlignmentFormats: multiple sequence alignment formats
"""

SequenceTypes = csb.pyutils.enum( NucleicAcid='NA', DNA='DNA', RNA='RNA', Protein='Protein', Unknown='Unknown' )
"""
@var SequenceTypes: Sequence types
"""

AlignmentTypes = csb.pyutils.enum( Global='global', Local='local' )
"""
@var AlignmentTypes: Alignment strategies
"""

class SequenceAlphabets(object):
    """
    Sequence alphabet enumerations
    """
    
    Nucleic = csb.pyutils.enum( Adenine='A', Cytosine='C', Guanine='G', Thymine='T', Uracil='U', Purine='R', Pyrimidine='Y', Ketone='K', Amino='M', 
                                Strong='S', Weak='W', NotA='B', NotC='D', NotG='H', NotT='V', Any='N', Masked='X', Gap='-' )
    """
    @cvar Nucleic: Nucleic sequence alphabet
    """
    
    Protein = csb.pyutils.enum( ALA='A', ASX='B', CYS='C', ASP='D', GLU='E', PHE='F', GLY='G', HIS='H', ILE='I', LYS='K', LEU='L', MET='M', ASN='N', 
                                PYL='O', PRO='P', GLN='Q', ARG='R', SER='S', THR='T', SEC='U', VAL='V', TRP='W', TYR='Y', GLX='Z', UNK='X', GAP='-', 
                                STOP='*' )
    """
    @cvar Protein: Protein sequence alphabet
    """                       
       
    Unknown = csb.pyutils.enum( UNK='X' )
    """
    @cvar Unknown: Unknown sequence alphabet
    """   
    
class StopTraversal(StopIteration):
    pass

class Sequence(object):
    """ 
    Simple FASTA sequence entry.
    
    @param id: id of the sequence
    @type id: str
    @param header: FASTA header, without the '>' sign
    @type header: str
    @param sequence: the sequence string itself
    @type sequence: str
    @param type: sequence type (L{SequenceTypes} member)
    @type type: L{csb.pyutils.EnumItem}
    """
    
    def __init__(self, id, header, sequence, type=SequenceTypes.Unknown):
        
        self.id = id
        self.header = header
        self.sequence = sequence
        self.type = type
    
    def __str__(self):
        return self.to_fasta()
    
    def to_fasta(self):
        """
        Return a FASTA formatted sequence
        
        @rtype: str
        """
        return '{0.header}\n{0.sequence}'.format(self)
            
    @staticmethod
    def from_string(fasta_string):
        """
        Parse a single FASTA-formatted or a plain text sequence.
        
        @param fasta_string: the input FASTA-formatted sequence
        @type fasta_string: str
        
        @return: a new sequence object
        @rtype: L{Sequence}
        
        @raise ValueError: if the sequence contains invalid characters
                           according to the protein FASTA format
        """
        lines = fasta_string.strip().splitlines()
        if not lines[0].startswith('>'):
            lines = ['>sequence'] + lines
        
        seq = Sequence(None, None, None)
        
        seq.header = lines[0]
        seq.id = seq.header[1:].split()[0]
        seq.sequence = ''.join(lines[1:])
        seq.sequence = re.sub('\s+', '', seq.sequence) 
        
        seq_residues = seq.sequence.replace('-', '').replace('*', '') 
        if seq_residues and not seq_residues.isalpha():
            raise ValueError(seq.sequence)
        
        return seq     
    
    @staticmethod
    def parse(fasta_file):
        """
        Read a FASTA formatted file.
        
        @param fasta_file: the input FASTA file name
        @type fasta_file: str
        
        @return: a list of L{Sequence}s
        @rtype: list
        
        @raise ValueError: if the sequence contains invalid characters
                           according to the protein FASTA format
        """
        seqs = []
        reader = csb.io.EntryReader(open(fasta_file), '>', None)
        
        for entry in reader.entries():
            seqs.append(Sequence.from_string(entry))
        
        return seqs   
    
class PDBSequence(Sequence):
    """ 
    Simple FASTA sequence entry.
    
    @param id: id of the sequence
    @type id: str
    @param header: FASTA header, without the '>' sign
    @type header: str
    @param sequence: the sequence string itself
    @type sequence: str
    @param type: sequence type (L{SequenceTypes} member)
    @type type: L{csb.pyutils.EnumItem}
    """
    
    def __init__(self, id, header, sequence, type=SequenceTypes.Unknown):
        
        super(PDBSequence, self).__init__(id, header, sequence, type)

    @property
    def accession(self):
        return self.id[:4].lower()
    
    @property
    def chain(self):
        return self.id[4:]
        
    @staticmethod
    def from_string(fasta_string):
        """
        Parse a single PDB FASTA-formatted sequence.
        
        @param fasta_string: the input PDB FASTA-formatted sequence
        @type fasta_string: str
        
        @return: a new PDB sequence object
        @rtype: L{PDBSequence}
        
        @raise ValueError: if the sequence contains invalid characters
                           according to the protein FASTA format
        """
        
        seq = super(PDBSequence, PDBSequence).from_string(fasta_string)
                
        if not (seq.header and seq.id) or not (len(seq.id) in(5, 6) and seq.header.find('mol:') != -1):
            raise ValueError('Does not look like a PDB header: {0}'.format(seq.header))

        id = seq.id.replace('_', '')
        stype =  seq.header.partition('mol:')[2].partition(' ')[0]
        try:
            stype = csb.pyutils.Enum.parse(SequenceTypes, stype, ignore_case=True)
        except csb.pyutils.EnumValueError:
            stype = SequenceTypes.Unknown
        
        return PDBSequence(id=id, header=seq.header, sequence=seq.sequence, type=stype)
    
    @staticmethod
    def parse(fasta_file):
        """
        Read a PDB FASTA formatted file.
        
        @param fasta_file: the input PDB FASTA file name
        @type fasta_file: str
        
        @return: a list of L{PDBSequence}s
        @rtype: list
        
        @raise ValueError: if the sequence contains invalid characters
                           according to the protein FASTA format
        """
        seqs = []
        reader = csb.io.EntryReader(open(fasta_file), '>', None)
        
        for entry in reader.entries():
            seqs.append(PDBSequence.from_string(entry))
        
        return seqs     

class _A3MFormatter(object):
    
    def __init__(self):
        pass
    
    @staticmethod
    def format_a3m(residue):
        if residue is None:
            return A3MAlignment.GAP
        elif residue is A3MAlignment.INSERTION:
            return ''
        return residue

    @staticmethod    
    def format_fasta(residue):
        if residue is None or residue is A3MAlignment.INSERTION:
            return A3MAlignment.GAP
        return residue                        
            
    def __call__(self, format):
        
        if format == AlignmentFormats.FASTA:
            return _A3MFormatter.format_fasta
        elif format == AlignmentFormats.A3M:
            return _A3MFormatter.format_a3m
        else:
            raise ValueError('Unsupported format {0}'.format(format))

class _A3MAlignmentEntries(object):
    
    def __init__(self, container):
        
        self._container = container
    
    def __getitem__(self, rank):
        try:
            entry = self._container._get_sequence(rank)
            columns = self._container._get_row(rank, AlignmentFormats.A3M)
            return tuple([entry] + columns)
        except (IndexError, KeyError):
            raise IndexError('Column rank ' + str(rank) + ' out of range')
                
class A3MAlignment(object):
    """
    Describes a multiple alignment in A3M format.
    
    @param sequences: a list of L{Sequence}s to build the alignment from
    @type param: list
    @param consensus: a L{Sequence} object, holding the consensus sequence
    @type consensus: L{Sequence}
    """
    
    GAP = '-'
    INSERTION = -1

    def __init__(self, sequences, consensus=None):

        self._msa = { }        
        self._sequences = { }
        self._matches = { }
                
        self._consensus = consensus
        self._formatter = _A3MFormatter()
        self._rows = _A3MAlignmentEntries(self)        
        
        temp = [ ]
        
        for i, s in enumerate(sequences, start=1):
            sequence = csb.pyutils.deepcopy(s)
            self._sequences[i] = sequence
            self._msa[i] = [ ]
            row = list(sequence.sequence)
            row.reverse()
            temp.append(row)
        
        col_rank = -1
        match_rank = 0
        
        while True:
            try:
                next_column = self._next_column(temp)
                col_rank += 1
                has_insertion = False
                for i, residue in enumerate(next_column, start=1):
                    self._msa[i].append(residue)
                    if residue == A3MAlignment.INSERTION:
                        has_insertion = True
                if not has_insertion:
                    match_rank += 1
                    self._matches[match_rank] = col_rank
            except StopTraversal:
                break
    
    def __str__(self):
        return self.to_string(headers=True)
    
    def __len__(self):
        return self.cols_count        
    
    @property
    def consensus(self):
        return self._consensus
    
    @property
    def rows(self):
        return self._rows
            
    @property
    def rows_count(self):
        return len(self._msa)
    
    @property
    def cols_count(self):
        if self.rows_count < 1:
            return 0
        return len(self._msa[ self._msa.keys()[0] ])
    
    @property
    def matches_count(self):
        return len(self._matches)

    @staticmethod
    def parse(a3m_string):
        """
        Create a new L{A3MAlignment} from an a3m-formatted string. 
        
        @param a3m_string: the input A3M alignment
        @type a3m_string: str
        
        @return: a new A3M Alignment instance
        @rtype: L{A3MAlignment}
        """
        from StringIO import StringIO
        
        stream = StringIO()
        stream.write(a3m_string)
        reader = csb.io.EntryReader(stream, '>', None)
        sequences, consensus = [ ], None
        
        for s in reader.entries():
            sequence = Sequence.from_string(s)
            if sequence.id == 'Consensus':
                consensus = sequence
            else:
                sequences.append(sequence)
        
        return A3MAlignment(sequences, consensus)
    
    def subregion(self, start_match, end_match):
        """
        Extract a subregion of the multiple alignment defined by the ranks of
        a start and an end match state. 
        
        @param start_match: start rank
        @type start_match: int
        @param end_match: end rank
        @type end_match: int
                
        @return: a new L{A3MAlignment}
        @rtype: L{A3MAlignment}
        """
        
        alignment = [ ]
                
        try:
            start = self._matches[start_match]
            end = self._matches[end_match]
        except KeyError as ke:
            raise IndexError('Match rank out of range: ' + str(ke))
        
        for seq_id in self._msa:
            
            row = self._get_row(seq_id, AlignmentFormats.A3M)[start : end + 1]
            alignment.append(self._sequences[seq_id].header)
            alignment.append(''.join(row))
        
        if self._consensus:
            alignment.append(self._consensus.header)
            alignment.append(self._consensus.sequence[start_match - 1 : end_match])
            
        return A3MAlignment.parse('\n'.join(alignment))

    def to_fasta(self, headers=True):
        """
        Dump the alignment in m_fasta format.
        
        @param headers: if set to False, headers are be omitted
        @type headers: bool
        
        @return: a FASTA-formatted alignment
        @rtype: str 
        """
        alignment = self._dump(AlignmentFormats.FASTA, headers=headers)
        return '\n'.join(alignment) 
    
    def to_string(self, headers=True):
        """
        Dump the alignment in its native A3M format.
        
        @param headers: if set to False, headers are be omitted
        @type headers: bool
        
        @return: an A3M-formatted alignment
        @rtype: str         
        """      
        alignment = self._dump(AlignmentFormats.A3M, headers=headers)
        return '\n'.join(alignment) 
    
    def _format(self, items, format):
        return map(self._formatter(format), items)
    
    def _dump(self, format, headers=True):
        
        alignment = [ ]        
        
        for seq_id in sorted(self._sequences):
            sequence = self._sequences[seq_id]
            if headers:
                alignment.append(sequence.header)
            row = self._format(self._msa[seq_id], format)
            alignment.append(''.join(row))
        
        return alignment             
    
    def _get_row(self, rank, format):
        return self._format(self._msa[rank], format)
    
    def _get_sequence(self, rank):
        return self._sequences[rank]
             
    def _next_column(self, sequences):
        
        column = [ ]
        has_insertion = False
        
        for sequence in sequences:
            if len(sequence) > 0 and sequence[-1].islower():
                has_insertion = True
        
        for sequence in sequences:
            try:
                if has_insertion and not sequence[-1].islower():
                    column.append(A3MAlignment.INSERTION)
                else:
                    column.append(sequence.pop())
            except IndexError:
                column.append(None)
                        
        if not any(column):
            raise StopTraversal()
        
        return column       
