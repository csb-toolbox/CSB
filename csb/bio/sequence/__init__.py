import re
import csb.pyutils
import csb.io

AlignmentFormats = csb.pyutils.enum( A3M='a3m', FASTA='fa', PIR='pir' )
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

    def __repr__(self):
        return self.to_fasta()
    
    def __len__(self):
        return len(self.sequence)

    def to_fasta(self):
        """
        Return a FASTA formatted sequence
        
        @rtype: str
        """
        marker = '>'
        if self.header and self.header.startswith(marker):
            marker = ''
        return '{0}{1.header}\n{1.sequence}'.format(marker, self)

    def expand(self, gap_characters = ('-','.')):
        """
        Returns the positions in this sequence, which are not gaps
        """
        s = []
        for i,char in enumerate(self.sequence):
            if not char in gap_characters:
                s.append(i)
        return s
        
    def is_gap(self, index, gap_characters=('-','.')):
        """
        Returns whether the symbol at index represents a gap
        """
        return self.sequence[index] in gap_characters

    def sequence_index(self, column_index, gap_characters=('-','.')):
        """
        Returns the postion of the colum column_index
        """
        s = self.sequence[:column_index]
        n = [s.count(c) for c in gap_characters]
        n = reduce(lambda x, y: x+y, n, 0)
        return column_index - n

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
        
        if format in(AlignmentFormats.FASTA, AlignmentFormats.PIR):
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
        Dump the alignment in mFASTA format.
        
        @param headers: if set to False, headers are be omitted
        @type headers: bool
        
        @return: a FASTA-formatted alignment
        @rtype: str 
        """
        alignment = self._dump(AlignmentFormats.FASTA, headers=headers)
        return '\n'.join(alignment) 

    def to_pir(self, target=1):
        """
        Dump the alignment in a Modeller-friendly PIR-like format.
        
        @param target: rank of the modelling target (usually 1)
        @type target: int
        
        @return: a PIR-formatted alignment
        @rtype: str 
        """
        alignment = self._dump(AlignmentFormats.PIR, flags=target)
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
    
    def _dump(self, format, headers=True, flags=None):
        
        alignment = [ ]        
        
        for seq_id in sorted(self._sequences):
            sequence = self._sequences[seq_id]
            if headers:
                if format == AlignmentFormats.PIR:
                    type = 'structure'
                    if seq_id == flags:
                        type = 'sequence'
                    id = '{0:5}'.format(sequence.id)
                    header = '>P1;{0}\n{2}:{0}:.:{1}:.:{1}:::: '.format(id[:4], id[4], type)
                    alignment.append(header) 
                else:
                    alignment.append(sequence.header)
            row = self._format(self._msa[seq_id], format)
            if format == AlignmentFormats.PIR:
                row += '*'
                
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




def insert_breaks(seq, length=60):
    s = [seq[i*length:(i+1)*length] for i in range(len(seq)/length+1)]
    return '\n'.join(s)


def get_accession_number(fasta_record):
    import re
    
    name = fasta_record.title.split()[0].split('|')[0]
    return re.sub('[^a-zA-Z0-9_\-]', '', name)


from csb.pyutils import ordered_dict

class Alignment(ordered_dict):

    formats = ('.fasta', '.fsa', '.fas')
    gap_characters = ('-', '.')

    def read(self, filename, numbers_as_keys = False):

        from Bio import Fasta
        from csb.bio.sequence import Sequence
        
        import os

        Sequence.gap_characters = Alignment.gap_characters
        
        format = os.path.splitext(filename)[1]

        if format not in self.formats:
            message = 'Format of file "%s" is "%s"; supported formats: %s'
            raise IOError, message % (filename, format,
                                      ' / '.join(self.formats))

        parser = Fasta.RecordParser()
        file = open(os.path.expanduser(filename))
        iterator = Fasta.Iterator(file, parser)

        record = iterator.next()

        counter = 1

        while record:

            id = record.title
            seq = record.sequence

            if numbers_as_keys:
                id = counter
                counter += 1

            if self.has_key(id):
                raise ValueError, 'Multiple entries found for %s.' % id

            self[id] = Sequence(id = id,header=record.title, sequence=seq)

            record = iterator.next()

    def write(self, filename):
        import os
        file = open(os.path.expanduser(filename), 'w')
        for id, seq in self.items():
            file.write('>%s\n%s\n' %(id, insert_breaks(seq.sequence)))
        file.close()

    def __str__(self):
        s = ['>%s\n%s' %(id, insert_breaks(seq)) for id, seq in self.items()]
        return '\n'.join(s)
    
    def make_equal_length(self):
        max_length = max(map(len, self.values()))
        for id in self.keys():
            new = self[id].sequence
            new = new + (max_length - len(new)) * '-'
            self[id] = Sequence(id = id, header = self[id].header,
                                sequence =new)

    def column(self, index):
        from numpy import array, sum
        
        columns = array([list(v.sequence) for v in self.values()])

        return columns[index].tolist()


    def columns(self, indices = None):
        from numpy import array, take, sum

        columns = array([list(v.sequence) for v in self.values()])

        if indices is not None:
            indices = list(indices)
            indices.sort()
            columns = take(columns, indices, 1)

        return (columns.transpose()).tolist()


    def sequence(self, member, keep_gaps = False):

        if not member in self:
            raise KeyError('sequence "%s" not in alignment' % str(member))
        
        
        return self[member]
        

    def matches(self, member1 = None, member2 = None, exact_matches = 0):

        a = 0
        b = 0

        if member1 is None: member1 = self.keys()[0]
        if member2 is None: member2 = self.keys()[1]

        x = self[member1].sequence
        y = self[member2].sequence

        matches = []
        
        for i in range(len(x)):

            x_gap = 0
            y_gap = 0

            if x[i] in self.gap_characters:
                x_gap = 1
            if y[i] in self.gap_characters:
                y_gap = 1

            if x_gap == 0 and y_gap == 0:
                if not exact_matches or x[i] == y[i]:
                    matches.append((a,b))

            if not x_gap: a += 1
            if not y_gap: b += 1

        return matches

    def consensus(self, no_gap = False):
        from numpy import argmax
        cons = ''
        for col in self.columns():
            if no_gap and '-' in col:
                cons += '-'
            else:
                counts = [col.count(x) for x in col]
                cons += col[argmax(counts)]
        return cons
                
    def identities(self):

        from numpy import argmax

        cons = ''

        for col in self.columns():

            counts = [col.count(x) for x in col]
            if max(counts) <> len(col):
                cons += '-'
            else:
                cons += col[argmax(counts)]

        return cons

    def rename(self, key1, key2):

        if key2 in self.keys():
            raise KeyError('"%s" already in alignment' % key2)

        self._ordered_dict__keys[self.keys().index(key1)] = key2
        self[key2] = self[key1]
        dict.__delitem__(self, key1)

