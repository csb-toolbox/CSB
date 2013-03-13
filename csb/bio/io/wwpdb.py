"""
PDB structure parsers, format builders and database providers.

The most basic usage is:

    >>> parser = StructureParser('structure.pdb')
    >>> parser.parse_structure()
    <Structure>     # a Structure object (model)

or if this is an NMR ensemble:
    
    >>> parser.parse_models()
    <Ensemble>      # an Ensemble object (collection of alternative Structure-s)
    
This module introduces a family of PDB file parsers. The common interface of all
parsers is defined in L{AbstractStructureParser}. This class has several
implementations:

    - L{RegularStructureParser} - handles normal PDB files with SEQRES fields
    - L{LegacyStructureParser} - reads structures from legacy or malformed PDB files,
      which are lacking SEQRES records (initializes all residues from the ATOMs instead)
    - L{PDBHeaderParser} - reads only the headers of the PDB files and produces structures
      without coordinates. Useful for reading metadata (e.g. accession numbers or just
      plain SEQRES sequences) with minimum overhead
      
Unless you have a special reason, you should use the L{StructureParser} factory,
which returns a proper L{AbstractStructureParser} implementation, depending on the
input PDB file. If the input file looks like a regular PDB file, the factory
returns a L{RegularStructureParser}, otherwise it instantiates L{LegacyStructureParser}.
L{StructureParser} is in fact an alias for L{AbstractStructureParser.create_parser}.

Writing your own, customized PDB parser is easy. Suppose that you are trying to
parse a PDB-like file which misuses the charge column to store custom info. This
will certainly crash L{RegularStructureParser} (for good), but you can create your
own parser as a workaround. All you need to to is to override the virtual
C{_read_charge} hook method::

    class CustomParser(RegularStructureParser):
    
        def _read_charge(self, line):
            try:
                return super(CustomParser, self)._read_charge(line)
            except StructureFormatError:
                return None

Another important abstraction in this module is L{StructureProvider}. It has several
implementations which can be used to retrieve PDB L{Structure}s from various sources:
file system directories, remote URLs, etc. You can easily create your own provider as
well. See L{StructureProvider} for details.

Finally, this module gives you some L{FileBuilder}s, used for text serialization 
of L{Structure}s and L{Ensemble}s:

    >>> builder = PDBFileBuilder(stream)
    >>> builder.add_header(structure)
    >>> builder.add_structure(structure)

where stream is any Python stream, e.g. an open file or sys.stdout.

See L{Ensemble} and L{Structure} from L{csb.bio.structure} for details on these
objects.
"""

import re
import os
import numpy
import datetime
import multiprocessing

import csb.bio.structure
import csb.bio.sequence
import csb.bio.sequence.alignment as alignment
import csb.core
import csb.io

from abc import ABCMeta, abstractmethod
from csb.bio.sequence import SequenceTypes, SequenceAlphabets
from csb.bio.structure import ChemElements, SecStructures


PDB_AMINOACIDS = {
    'PAQ': 'TYR', 'AGM': 'ARG', 'ILE': 'ILE', 'PR3': 'CYS', 'GLN': 'GLN',
    'DVA': 'VAL', 'CCS': 'CYS', 'ACL': 'ARG', 'GLX': 'GLX', 'GLY': 'GLY',
    'GLZ': 'GLY', 'DTH': 'THR', 'OAS': 'SER', 'C6C': 'CYS', 'NEM': 'HIS',
    'DLY': 'LYS', 'MIS': 'SER', 'SMC': 'CYS', 'GLU': 'GLU', 'NEP': 'HIS',
    'BCS': 'CYS', 'ASQ': 'ASP', 'ASP': 'ASP', 'SCY': 'CYS', 'SER': 'SER',
    'LYS': 'LYS', 'SAC': 'SER', 'PRO': 'PRO', 'ASX': 'ASX', 'DGN': 'GLN',
    'DGL': 'GLU', 'MHS': 'HIS', 'ASB': 'ASP', 'ASA': 'ASP', 'NLE': 'LEU',
    'DCY': 'CYS', 'ASK': 'ASP', 'GGL': 'GLU', 'STY': 'TYR', 'SEL': 'SER',
    'CGU': 'GLU', 'ASN': 'ASN', 'ASL': 'ASP', 'LTR': 'TRP', 'DAR': 'ARG',
    'VAL': 'VAL', 'CHG': 'ALA', 'TPO': 'THR', 'CLE': 'LEU', 'GMA': 'GLU',
    'HAC': 'ALA', 'AYA': 'ALA', 'THR': 'THR', 'TIH': 'ALA', 'SVA': 'SER',
    'MVA': 'VAL', 'SAR': 'GLY', 'LYZ': 'LYS', 'BNN': 'ALA', '5HP': 'GLU',
    'IIL': 'ILE', 'SHR': 'LYS', 'HAR': 'ARG', 'FME': 'MET', 'ALO': 'THR',
    'PHI': 'PHE', 'ALM': 'ALA', 'PHL': 'PHE', 'MEN': 'ASN', 'TPQ': 'ALA',
    'GSC': 'GLY', 'PHE': 'PHE', 'ALA': 'ALA', 'MAA': 'ALA', 'MET': 'MET',
    'UNK': 'UNK', 'LEU': 'LEU', 'ALY': 'LYS', 'SET': 'SER', 'GL3': 'GLY',
    'TRG': 'LYS', 'CXM': 'MET', 'TYR': 'TYR', 'SCS': 'CYS', 'DIL': 'ILE',
    'TYQ': 'TYR', '3AH': 'HIS', 'DPR': 'PRO', 'PRR': 'ALA', 'CME': 'CYS',
    'IYR': 'TYR', 'CY1': 'CYS', 'TYY': 'TYR', 'HYP': 'PRO', 'DTY': 'TYR',
    '2AS': 'ASP', 'DTR': 'TRP', 'FLA': 'ALA', 'DPN': 'PHE', 'DIV': 'VAL',
    'PCA': 'GLU', 'MSE': 'MET', 'MSA': 'GLY', 'AIB': 'ALA', 'CYS': 'CYS',
    'NLP': 'LEU', 'CYQ': 'CYS', 'HIS': 'HIS', 'DLE': 'LEU', 'CEA': 'CYS',
    'DAL': 'ALA', 'LLP': 'LYS', 'DAH': 'PHE', 'HMR': 'ARG', 'TRO': 'TRP',
    'HIC': 'HIS', 'CYG': 'CYS', 'BMT': 'THR', 'DAS': 'ASP', 'TYB': 'TYR',
    'BUC': 'CYS', 'PEC': 'CYS', 'BUG': 'LEU', 'CYM': 'CYS', 'NLN': 'LEU',
    'CY3': 'CYS', 'HIP': 'HIS', 'CSO': 'CYS', 'TPL': 'TRP', 'LYM': 'LYS',
    'DHI': 'HIS', 'MLE': 'LEU', 'CSD': 'ALA', 'HPQ': 'PHE', 'MPQ': 'GLY',
    'LLY': 'LYS', 'DHA': 'ALA', 'DSN': 'SER', 'SOC': 'CYS', 'CSX': 'CYS',
    'OMT': 'MET', 'DSP': 'ASP', 'PTR': 'TYR', 'TRP': 'TRP', 'CSW': 'CYS',
    'EFC': 'CYS', 'CSP': 'CYS', 'CSS': 'CYS', 'SCH': 'CYS', 'OCS': 'CYS',
    'NMC': 'GLY', 'SEP': 'SER', 'BHD': 'ASP', 'KCX': 'LYS', 'SHC': 'CYS',
    'C5C': 'CYS', 'HTR': 'TRP', 'ARG': 'ARG', 'TYS': 'TYR', 'ARM': 'ARG',
    'DNP': 'ALA'
    }
"""
Dictionary of non-standard amino acids, which could be found in PDB.
"""


PDB_NUCLEOTIDES = {
    'DA': 'Adenine', 'DG': 'Guanine', 'DC': 'Cytosine', 'DT': 'Thymine',
     'A': 'Adenine', 'G': 'Guanine', 'C': 'Cytosine', 'T': 'Thymine',
     'U': 'Uracil', 'DOC': 'Cytosine', 'R': 'Purine', 'Y': 'Pyrimidine',
     'K': 'Ketone', '  M': 'Amino', 'S': 'Strong', 'W': 'Weak',
     'B': 'NotA', 'D'  : 'NotC', 'H': 'NotG', 'V': 'NotT',
     'N': 'Any', 'X'  : 'Masked'
    }
"""
Dictionary of non-standard nucleotides, which could be found in PDB.
"""

class PDBParseError(ValueError):
    pass

class HeaderFormatError(PDBParseError):
    pass

class SecStructureFormatError(PDBParseError):
    pass

class StructureFormatError(PDBParseError):
    pass

class UnknownPDBResidueError(PDBParseError):
    pass

class StructureNotFoundError(KeyError):
    pass

class InvalidEntryIDError(StructureFormatError):
    pass

class ResidueMappingError(StructureFormatError):
    pass


class EntryID(object):
    """
    Represents a PDB Chain identifier. Implementing classes must define
    how the original ID is split into accession number and chain ID.
    
    @param id: identifier
    @type id: str
    """
    __metaclass__ = ABCMeta

    def __init__(self, id):
        
        self._accession = ''
        self._chain = ''
        
        id = id.strip()        
        self._accession, self._chain = self.parse(id)
    
    @staticmethod
    def create(id):
        """
        Guess the format of C{id} and parse it. 
        
        @return: a new PDB ID of the appropriate type
        @rtype: L{EntryID}
        """
        
        if len(id) in (4, 5):
            return StandardID(id)
        elif len(id) == 6 and id[4] == '_':
            return SeqResID(id)
        else:
            return DegenerateID(id)
    
    @abstractmethod
    def parse(self, id):
        """
        Split C{id} into accession number and chain ID.
        
        @param id: PDB identifier
        @type id: str
        @return: (accession, chain)
        @rtype: tuple of str
        
        @raise InvalidEntryIDError: when C{id} is in an inappropriate format
        """
        pass
    
    def format(self):
        """
        @return: the identifier in its original format
        @rtype: str
        """
        return self.entry_id
        
    @property
    def accession(self):
        """
        Accession number part of the Entry ID
        @rtype: str
        """
        return self._accession

    @property
    def chain(self):
        """
        Chain ID part of the Entry ID
        @rtype: str      
        """
        return self._chain
    
    @property
    def entry_id(self):
        """
        Accession number + Chain ID
        @rtype: str
        """        
        return "{0.accession}{0.chain}".format(self)
    
    def __str__(self):
        return self.entry_id
    
class StandardID(EntryID):
    """
    Standard PDB ID in the following form: xxxxY, where xxxx is the accession
    number (lower case) and Y is an optional chain identifier.
    """    
    def parse(self, id):
        
        if len(id) not in (4, 5):
            raise InvalidEntryIDError(id)
        
        return (id[:4].lower(), id[4:])
        
class DegenerateID(EntryID):
    """
    Looks like a L{StandardID}, except that the accession number may have
    arbitrary length.
    """    
    def parse(self, id):

        if len(id) < 2:
            raise InvalidEntryIDError(id)
        
        return (id[:-1].lower(), id[-1])
        
class SeqResID(EntryID):
    """
    Same as a L{StandardID}, but contains an additional underscore between
    te accession number and the chain identifier.
    """      
    def parse(self, id):
        
        if not (len(id) == 6 and id[4] == '_'):
            raise InvalidEntryIDError(id)
        
        return (id[:4].lower(), id[5:])
    
    def format(self):
        return "{0.accession}_{0.chain}".format(self)


class AbstractStructureParser(object):
    """
    A base PDB structure format-aware parser. Subclasses must implement the
    internal abstract method C{_parse_header} in order to complete the 
    implementation.

    @param structure_file: the input PD file to parse
    @type structure_file: str
    @param check_ss: if True, secondary structure errors in the file will cause 
                     L{SecStructureFormatError} exceptions
    @type check_ss: bool
    @param mapper: residue mapper, used to align ATOM records to SEQRES.
                   If None, use the default (L{CombinedResidueMapper})
    @type mapper: L{AbstractResidueMapper}    

    @raise IOError: when the input file cannot be found
    """

    __metaclass__ = ABCMeta

    @staticmethod
    def create_parser(structure_file, check_ss=False, mapper=None):
        """
        A StructureParser factory, which instantiates and returns the proper parser
        object based on the contents of the PDB file.

        If the file contains a SEQRES section, L{RegularStructureParser} is returned,
        otherwise L{LegacyStructureParser} is instantiated. In the latter case
        LegacyStructureParser will read the sequence data directly from the ATOMs.        
        
        @param structure_file: the PDB file to parse
        @type structure_file: str
        @param check_ss: if True, secondary structure errors in the file will cause 
                         L{SecStructureFormatError} exceptions
        @type check_ss: bool
        @param mapper: residue mapper, used to align ATOM records to SEQRES.
                       If None, use the default (L{CombinedResidueMapper})
        @type mapper: L{AbstractResidueMapper}          
        
        @rtype: L{AbstractStructureParser}
        """
        has_seqres = False
        
        for line in open(structure_file):
            if line.startswith('SEQRES'):
                has_seqres = True
                break
        
        if has_seqres:
            return RegularStructureParser(structure_file, check_ss, mapper)
        else:
            return LegacyStructureParser(structure_file, check_ss, mapper)        

    def __init__(self, structure_file, check_ss=False, mapper=None):

        self._file = None
        self._stream = None
        self._mapper = CombinedResidueMapper()
        self._check_ss = bool(check_ss)

        self.filename = structure_file
        if mapper is not None:
            self.mapper = mapper

    def __del__(self):
        try:
            self._stream.close()
        except:
            pass
        
    @property
    def mapper(self):
        """
        Current residue mapping strategy
        @rtype: L{AbstractResidueMapper}
        """
        return self._mapper
    @mapper.setter
    def mapper(self, value):
        if not isinstance(value, AbstractResidueMapper):
            raise TypeError(value)
        self._mapper = value

    @property
    def filename(self):
        """
        Current input PDB file name
        @rtype: str
        """
        return self._file
    @filename.setter
    def filename(self, name):
        try:
            stream = open(name)
        except IOError:
            raise IOError('File not found: {0}'.format(name))
        
        if self._stream:
            try:
                self._stream.close()
            except:
                pass
        self._stream = stream
        self._file = name

    def models(self):
        """
        Find all available model identifiers in the structure.

        @return: a list of model IDs
        @rtype: list
        """
        models = []
        check = set()

        with open(self._file, 'r') as f:
            for line in f:
                if line.startswith('MODEL'):
                    model_id = int(line[10:14])
                    if model_id in check:
                        raise StructureFormatError('Duplicate model identifier: {0}'.format(model_id))
                    models.append(model_id)
                    check.add(model_id)

        if len(models) > 0:
            if not(min(check) == 1 and max(check) == len(models)):
                raise StructureFormatError('Non-consecutive model identifier(s) encountered')                
            return models
        else:
            return []
        
    def guess_chain_type(self, residue_labels):
        """
        Try to guess what is the sequence type of a chunk of PDB
        C{residue_label}s. The list of labels is probed starting from the middle
        first, because PDB chains often contain modified / unknown residues at
        the termini. If none of the probed residues can be used to determine
        chain's type, just give up and return L{SequenceTypes.Unknown}.
        
        @param residue_labels: an iterable of PDB residue labels
        @type residue_labels: iterable 

        @return: a L{SequenceTypes} enum member
        @rtype: L{csb.core.EnumItem}              
        """
        
        labels = list(residue_labels)
        middle = int(len(labels) / 2)
        
        reordered = labels[middle:] + list(reversed(labels[:middle]))
        
        for label in reordered:
            try:
                type = self.guess_sequence_type(label)
                if type != SequenceTypes.Unknown:
                    return type
                
            except UnknownPDBResidueError:
                continue
            
        return SequenceTypes.Unknown

    def guess_sequence_type(self, residue_label):
        """
        Try to guess what is the sequence type of a PDB C{residue_label}.

        @param residue_label: a PDB-conforming name of a residue
        @type residue_label: str

        @return: a L{SequenceTypes} enum member
        @rtype: L{csb.core.EnumItem}

        @raise UnknownPDBResidueError: when there is no such PDB residue name
                                       in the catalog tables
        """
        if residue_label in PDB_AMINOACIDS:
            return SequenceTypes.Protein                                
        elif residue_label in PDB_NUCLEOTIDES:
            return SequenceTypes.NucleicAcid                            
        else:
            raise UnknownPDBResidueError(residue_label)

    def parse_residue(self, residue_label, as_type=None):
        """
        Try to parse a PDB C{residue_label} and return its closest 'normal'
        string representation. If a sequence type (C{as_type}) is defined,
        guess the alphabet based on that information, otherwise try first to
        parse it as a protein residue.

        @param residue_label: a PDB-conforming name of a residue
        @type residue_label: str
        @param as_type: suggest a sequence type (L{SequenceTypes} member)
        @type L{scb.core.EnumItem}

        @return: a normalized residue name
        @rtype: str

        @raise UnknownPDBResidueError: when there is no such PDB residue name
                                       in the catalog table(s)
        """
        if as_type is None:
            as_type = self.guess_sequence_type(residue_label)         

        try:
            if as_type == SequenceTypes.Protein:                            
                return PDB_AMINOACIDS[residue_label]
            elif as_type == SequenceTypes.NucleicAcid:                      
                return PDB_NUCLEOTIDES[residue_label]
            else:
                raise UnknownPDBResidueError(residue_label)
        except KeyError:
            raise UnknownPDBResidueError(residue_label)

    def parse_residue_safe(self, residue_label, as_type):
        """
        Same as C{parse_residue}, but returns UNK/Any instead of raising
        UnknownPDBResidueError.

        @param residue_label: a PDB-conforming name of a residue
        @type residue_label: str
        @param as_type: suggest a sequence type (L{SequenceTypes} member)
        @type L{scb.core.EnumItem}

        @return: a normalized residue name
        @rtype: str
        """
        try:
            return self.parse_residue(residue_label, as_type)
        except UnknownPDBResidueError:
            if as_type == SequenceTypes.Protein:                        
                return repr(SequenceAlphabets.Protein.UNK)              
            elif as_type == SequenceTypes.NucleicAcid:                  
                return repr(SequenceAlphabets.Nucleic.Any)              
            else:
                return repr(SequenceAlphabets.Unknown.UNK)              

    def parse(self, filename=None, model=None):

        if filename:
            self.filename = filename
        return self.parse_structure(model)
    
    def parse_structure(self, model=None):
        """
        Parse and return the L{Structure} with the specified model identifier.
        If no explicit model is specified, parse the first model in the
        structure.

        @param model: parse exactly the model with this ID
        @type model: str

        @return: object representation of the selected model
        @rtype: L{Structure}

        @raise ValueError: When an invalid model ID is specified
        @raise PDBParseError: When the input PDB file suffers from unrecoverable
                              corruption. More specialized exceptions will be
                              raised depending on the context (see L{PDBParseError}'s
                              subclasses). 
        """
        if model is not None:
            model = int(model)
        
        try:            
            structure = self._parse_header(model)
        except PDBParseError:
            raise
        except ValueError as ex:
            raise HeaderFormatError("Malformed header: {0}".format(ex))
        
        self._parse_atoms(structure, model)
        self._parse_ss(structure)

        return structure    
    
    def parse_models(self, models=()):
        """
        Parse the specified models in the file and build an L{Ensemble}.
        
        @param models: an iterable object providing model identifiers.
                       If not specified, all models will be parsed.
        @type models: tuple
        
        @return: an ensemble with all parsed models
        @rtype: L{Ensemble}
        """
        
        if not models:
            models = self.models()
        else:
            models = list(map(int, models))
        
        ensemble = csb.bio.structure.Ensemble()
        
        if len(models) > 0:
            for model_id in models:
                model = self.parse_structure(model_id)
                ensemble.models.append(model)
        else:
            model = self.parse_structure()
            model.model_id = 1
            ensemble.models.append(model)
        
        return ensemble

    def parse_biomolecule(self, number=1, single=False):
        """
        Parse and return the L{Structure} of the biological unit (quaternary
        structure) as annotated by the REMARK 350 BIOMOLECULE record.

        @param number: biomolecule number
        @type number: int

        @param single: if True, assign new single-letter chain
                       identifiers. If False, assign multi-letter chain identifiers whith a
                       number appended to the original identifier, like "A1", "A2", ...
        @type single: bool

        @return: structure of biological unit
        @rtype: L{Structure}
        """
        remarks = self._parse_remarks()
        if 350 not in remarks:
            raise PDBParseError('There is no REMARK 350')

        current = 1
        biomt = {current: {}}
        chains = tuple()

        def split(line):
            return [c.strip() for c in line.split(',') if c.strip() != '']

        for line in remarks[350]:
            if line.startswith('BIOMOLECULE:'):
                current = int(line[12:])
                biomt[current] = {}
            elif line.startswith('APPLY THE FOLLOWING TO CHAINS:'):
                chains = tuple(split(line[30:]))
            elif line.startswith('                   AND CHAINS:'):
                chains += tuple(split(line[30:]))
            elif line.startswith('  BIOMT'):
                num = int(line[8:12])
                vec = line[12:].split()
                vec = list(map(float, vec))
                biomt[current].setdefault(chains, dict()).setdefault(num, []).extend(vec)

        if number not in biomt or len(biomt[number]) == 0:
            raise KeyError('no BIOMOLECULE number {0}'.format(number))

        asu = self.parse_structure()
        structure = csb.bio.structure.Structure('{0}_{1}'.format(asu.accession, number))

        newchainiditer = iter('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789')

        for chains, matrices in biomt[number].items():
            for num in matrices:
                mat = numpy.array(matrices[num][0:12]).reshape((3,4))
                R, t = mat[:3,:3], mat[:3,3]

                for chain in chains:
                    if chain not in asu:
                        raise PDBParseError('chain {0} missing'.format(chain))
                    copy = asu[chain].clone()
                    copy.transform(R, t)
                    if single:
                        if len(structure.chains) == 62:
                            raise ValueError('too many chains for single=True')
                        copy.id = next(newchainiditer)
                    else:
                        copy.id = '{0}{1}'.format(chain, num)
                    structure.chains.append(copy)

        return structure

    @abstractmethod
    def _parse_header(self, model):
        """
        An abstract method, which subclasses must use to customize the way the
        PDB header (or the absence of a such) is handled. The implementation
        must rely on reading character data from the current internal 
        self._stream and must return a new L{csb.bio.structure.Structure} 
        instance with properly initialized header data: accession, model, 
        molecule identifiers, chains and residues. This structure object is
        then internally passed to the C{_parse_atoms} hook, responsible for
        attachment of the atoms to the residues in the structure.
        
        @param model: model ID to parse
        @type model: str
        
        @rtype: L{Structure}
        """
        pass

    def _scroll_model(self, model, stream):
        """
        Scroll the C{stream} to the specified C{model}.
        """
        
        while True:
            try:
                line = next(stream)
            except StopIteration:
                raise ValueError('No such model {0} in the structure.'.format(model))
            
            if line.startswith('MODEL'):
                model_id = self._read_model(line)
                if model == model_id:
                    return model_id      
        
    def _parse_atoms(self, structure, model):
        """
        Parse the ATOMs from the specified C{model} and attach them to the
        C{structure}.
        
        @param structure: L{Structure} being constructed
        @type structure:L{Structure}
        @param model: model ID to parse
        @type model: str        
        """
        
        structure.model_id = None

        chains = set()
        total_chains = len([c for c in structure.items if c.length > 0])
        
        residues = dict( (chain, []) for chain in structure.chains )
        seen_residues = dict( (chain, {}) for chain in structure.chains )
        
        in_ligands = False                
        in_atom = False
        read_model = False

        self._stream.seek(0)
        while True:

            try:
                line = next(self._stream)
            except StopIteration:
                break
            
            if line.startswith('MODEL'):
                if read_model:
                    break
                else:
                    self._parse_model_line(line, structure, model)
                    model = structure.model_id
                    read_model = True

            elif line.startswith('ATOM') or \
                            (line.startswith('HETATM') and not in_ligands):
                in_atom = True
                
                info = self._parse_atom_line(line, structure)                
                chains.add(info.chain)
                
                if info.id not in seen_residues[info.chain]:
                    residues[info.chain].append(info)
                    seen_residues[info.chain][info.id] = info
                else:
                    atom = info.atoms[0]
                    seen_residues[info.chain][info.id].atoms.append(atom)
                    
            elif in_atom and line.startswith('TER'):
                in_atom = False
                if len(chains) == total_chains:
                    in_ligands = True

            elif line.startswith('ENDMDL'):
                break

            elif line.startswith('END'):
                break

        if structure.model_id != model:
            raise ValueError('No such model {0} in the structure.'.format(model))
                
        self._map_residues(structure, residues)
        
    def _parse_model_line(self, line, structure, model):
        """
        Handle a new MODEL line. The default implementation will read the model
        ID and attach it to the C{structure}.
        
        @param line: raw string line to parse
        @type line: str
        @param structure: L{Structure} being constructed
        @type structure:L{Structure}
        
        @note: this method may have side effects: scrolls the current stream
        
        @return: read model ID
        @rtype: int
        """
        if model and model != self._read_model(line):
            self._scroll_model(model, self._stream)
            structure.model_id = model
        else:
            model = self._read_model(line)
            structure.model_id = model
                            
        return model
        
    def _parse_atom_line(self, line, structure):
        """
        Handle a new ATOM or HETATM line. The default implementation will read
        all data fields and create a new L{Atom}.
        
        @param line: raw string line to parse
        @type line: str
        @param structure: L{Structure} being constructed
        @type structure:L{Structure}
        
        @return: newly constructed atom
        @rtype: L{ResidueInfo}
        """

        atom = self._read_atom(line)
        
        rank = self._read_sequence_number(line)            
        sequence_number = rank            
        insertion_code = self._read_insertion_code(line)        
        
        id = str(sequence_number)        
        if insertion_code:
            id += insertion_code

        chain = self._read_chain_id(line)
        if chain not in structure.chains:
            raise StructureFormatError("Chain {0} is undefined".format(chain))

        type = self._read_residue(line, structure.chains[chain])
        label = self._read_residue_raw(line)
        
        atom.alternate = self._read_alternate(line)
        atom.occupancy = self._read_occupancy(line)
        atom.bfactor = self._read_bfactor(line)
        atom.charge = self._read_charge(line)

        info = ResidueInfo(chain, rank, id, sequence_number, insertion_code, type, label)
        info.atoms = [atom]
        
        return info
    
    def _read_model(self, line):
        """
        @return: model identifier
        @rtype: int
        """          
        try:
            return int(line[10:14])
        except ValueError:
            raise StructureFormatError("Invalid model ID: {0}".format(line))  
    
    def _read_atom(self, line):
        """
        @return: a new atom (serial_number, name, element, vector)
        @rtype: L{Atom}
        """
        try:
            serial_number = int(line[6:11])
            name = line[12:16]
            x, y, z = line[30:38], line[38:46], line[46:54]
            vector = numpy.array([float(x), float(y), float(z)])
        except ValueError as ve:
            raise StructureFormatError("Invalid ATOM line: {0}".format(ve))
        
        element = self._read_element(line)
        return csb.bio.structure.Atom(serial_number, name, element, vector)        

    def _read_sequence_number(self, line):
        """
        @return: PDB sequence number
        @rtype: int
        """
        try:
            return int(line[22:26])
        except ValueError:
            raise StructureFormatError("Invalid sequence number")
        
    def _read_insertion_code(self, line):
        """
        @return: PDB insertion code
        @rtype: str or None
        """       
        code = line[26].strip()
        
        if code:
            return code
        else:
            return None

    def _read_chain_id(self, line):
        """
        @return: chain identifier
        @rtype: str
        """          
        return line[21].strip()
        
    def _read_residue(self, line, chain):
        """
        @param chain: owning L{Chain} object 
        @type chain: L{Chain}
        
        @return: a member of any alphabet (e.g. L{SequenceAlphabets.Protein})
        @rtype: L{EnumItem}
        """
        raw = self._read_residue_raw(line)
        residue = self.parse_residue_safe(raw, as_type=chain.type)
        
        try:
            if chain.type == SequenceTypes.NucleicAcid:                               
                return csb.core.Enum.parsename(SequenceAlphabets.Nucleic, residue)
            else:
                return csb.core.Enum.parsename(SequenceAlphabets.Protein, residue)
        except csb.core.EnumMemberError:
            raise StructureFormatError("{0} is not a valid {1} residue".format(raw, chain.type))  
        
    def _read_residue_raw(self, line):
        """
        @rtype: str
        """
        return line[17:20].strip()              
            
    def _read_element(self, line):
        """
        @return: a member of L{ChemElements}
        @rtype: L{EnumItem} or None
        """
        element = line[76:78].strip()
        if element:
            try:
                element = csb.core.Enum.parsename(ChemElements, element)
            except csb.core.EnumMemberError:
                if element in ('D', 'X'):
                    element = ChemElements.x                            
                else:
                    raise StructureFormatError('Unknown chemical element: {0}'.format(element))
        else:
            element = None
            
        return element
    
    def _read_alternate(self, line):
        """
        @return: alt identifier
        @rtype: str or None
        """
        alternate = line[16].strip()
        
        if not alternate:
            return None
        else:
            return alternate
    
    def _read_occupancy(self, line):
        """
        @return: occupancy
        @rtype: float or None
        """        
        try:
            return float(line[54:60].strip() or 0)
        except ValueError:
            raise StructureFormatError("Malformed occupancy field")

    def _read_bfactor(self, line):
        """
        @return: b-factor
        @rtype: float or None
        """
        try:
            return float(line[60:66].strip() or 0)
        except ValueError:
            raise StructureFormatError("Malformed bfactor field")        
            
    def _read_charge(self, line):
        """
        @return: charge
        @rtype: int or None
        """        
        charge = line[78:80].strip()
        
        if charge:
            if charge in ('+', '-'):
                charge += '1'
            if charge[-1] in ('+', '-'):
                charge = charge[::-1]
            try:
                return int(charge)
            except ValueError:
                raise StructureFormatError("Malformed charge field") 
        else:
            return None

    def _map_residues(self, structure, residues):
        """
        Attach each L{Atom} to its corresponding L{Residue}.
        
        So far we have constructed a sparse (fragmented) chain given the information
        we have read from the ATOM/HETATM records. That includes PDB sequence
        identifiers and insertion codes, which cover only residues with XYZ coordinates
        and often do not correspond to our L{Residue} ranks.
        Our job is to find the right correspondence by matching this unreal sequence
        to what we have got from the SEQRES fields (the complete, gap-free protein
        sequence). Here we delegate this task to the current L{AbstractResidueMapper}
        strategy, which does all the magic.
        
        @param structure: L{Structure} being constructed
        @type structure:L{Structure}
        @param residues: all L{Atom}s which have been constructed so far.
                         This must be a map of the form:
                         <chainID: [L{ResidueInfo}1, L{ResidueInfo}2...]>
        @type residues: dict of L{ResidueInfo}        
        """

        if set(structure.chains) != set(residues.keys()):
            raise PDBParseError("Corrupt PDB file")

        for chain in structure.items:            
            if chain.length == 0 or len(residues[chain.id]) == 0:
                continue
            
            reference = SparseChainSequence.create(chain)
            sparse = SparseChainSequence(
                            "atoms", "", residues[chain.id], chain.type)
            
            aligned = self.mapper.map(sparse, reference)
            assert aligned.length == chain.length 
            
            for residue, mapped in zip(chain.residues, aligned.residues):
                if mapped.type != sparse.alphabet.GAP:
                    residue.id = (mapped.sequence_number, mapped.insertion_code)
                    
                    for atom in mapped.atoms:
                        residue.atoms.append(atom)

    def _parse_ss(self, structure):
        """
        Parse and attach secondary structure data.
        
        @bug: Currently the PDB helix types are ignored. Each HELIX line is treated
              as a regular SecStructures.Helix. This is due to incompatibility 
              between DSSP and PDB helix types.
        @todo: Implement a proper workaround for the previous bug (e.g. skip all
               helices types not included in the DSSP enum)
        
        @warning: In this implementation only the start/end positions of the SS 
                  elements are parsed. Additional data like H-bonding is ignored.
                  
        @bug: Currently structure.to_pdb() is not writing any SS data. 
        """        
        elements = {}
        self._stream.seek(0)
        
        while True:
            try:
                line = next(self._stream)
            except StopIteration:
                break

            if line.startswith('HELIX'):
                
                chain = structure.chains[line[19].strip()]
                if chain.id not in elements:
                    elements[chain.id] = []
                if chain.id != line[31].strip():
                    if self._check_ss:
                        raise SecStructureFormatError('Helix {0} spans multiple chains'.format(line[7:10]))
                    else:
                        continue                
                try:
                    startres = chain.find(line[21:25].strip(), line[25].strip())                
                    endres = chain.find(line[33:37].strip(), line[37].strip())
                except csb.core.ItemNotFoundError as ex:
                    if self._check_ss:
                        raise SecStructureFormatError(
                                'Helix {0} refers to an undefined residue ID: {1}'.format(line[7:10], str(ex)))
                    else:
                        continue
                if not startres.rank <= endres.rank:
                    if self._check_ss:
                        raise SecStructureFormatError('Helix {0} is out of range'.format(line[7:10]))
                    else:
                        continue                    
                helix = csb.bio.structure.SecondaryStructureElement(startres.rank, endres.rank, SecStructures.Helix)        
                elements[chain.id].append(helix)
            
            if line.startswith('SHEET'):
                
                chain = structure.chains[line[21].strip()]
                if chain.id not in elements:
                    elements[chain.id] = []                
                if chain.id != line[32].strip():
                    if self._check_ss:                    
                        raise SecStructureFormatError('Sheet {0} spans multiple chains'.format(line[7:10]))
                    else:
                        continue
                try:
                    startres = chain.find(line[22:26].strip(), line[26].strip())                
                    endres = chain.find(line[33:37].strip(), line[37].strip())
                except csb.core.ItemNotFoundError as ex:
                    if self._check_ss:                    
                        raise SecStructureFormatError(
                                'Sheet {0} refers to an undefined residue ID: {1}'.format(line[7:10], str(ex)))
                    else:
                        continue
                if not startres.rank <= endres.rank:
                    if self._check_ss:
                        raise SecStructureFormatError('Sheet {0} is out of range'.format(line[7:10]))
                    else:
                        continue
                strand = csb.bio.structure.SecondaryStructureElement(startres.rank, endres.rank, SecStructures.Strand)      
                elements[chain.id].append(strand)         
            
            elif line.startswith('MODEL') or line.startswith('ATOM'):
                break        
            
        for chain_id in elements:
            ss = csb.bio.structure.SecondaryStructure()
            for e in elements[chain_id]:
                ss.append(e)
            structure.chains[chain_id].secondary_structure = ss
        
    def _parse_remarks(self):
        """
        Read REMARK lines from PDB file.
        
        @return: dictionary with remark numbers as keys, and lists of lines as values.
        @rtype: dict
        """
        self._stream.seek(0)
        
        remarks = {}
        
        for line in self._stream:
            if line.startswith('REMARK'):
                num = int(line[7:10])
                lstring = line[11:]
                remarks.setdefault(num, []).append(lstring)
            elif line.startswith('DBREF') or line.startswith('ATOM'):
                break        
        
        return remarks


class RegularStructureParser(AbstractStructureParser):
    """
    This is the de facto PDB parser, which is designed to read SEQRES and ATOM
    sections separately, and them map them. Intentionally fails to parse
    malformed PDB files, e.g. a PDB file without a HEADER section.
    """
    
    def _parse_header(self, model):
        """
        Parse the HEADER section of a regular PDB file.
        
        @return: a L{csb.bio.structure.Structure} instance with properly 
                 initialized residues from the SEQRES.
        @rtype: L{csb.bio.structure.Structure}
        
        @raise PDBParseError: if the stream has no HEADER at byte 0
        """
        
        self._stream.seek(0)

        header = next(self._stream)
        if not header.startswith('HEADER'):
            raise PDBParseError('Does not look like a regular PDB file.')

        structure = csb.bio.structure.Structure(header.split()[-1])

        while True:

            try:
                line = next(self._stream)
            except StopIteration:
                break

            if line.startswith('COMPND'):
                if line[10:].lstrip().startswith('MOL_ID:'):
                    mol_id = int(line[18:].replace(';', '').strip())
                    chain_name = ''
                    chains = ''
                    while line.startswith('COMPND'):
                        line = next(self._stream)
                        if line.split()[2].startswith('MOLECULE:'):
                            chain_name += line[20:].strip()
                            while not chain_name.endswith(';'):
                                line = next(self._stream)
                                if not line.startswith('COMPND'):
                                        break
                                chain_name += ' ' + line[11:].strip()
                        else:
                            while not line.split()[2].startswith('CHAIN:'):
                                line = next(self._stream)
                                if not line.startswith('COMPND'):
                                    raise HeaderFormatError('Missing chain identifier in COMPND section')
                            chains = line[17:].strip()
                            while not chains.endswith(';'):
                                line = next(self._stream)
                                if not line.startswith('COMPND'):
                                    break
                                chains += ', ' + line[11:].strip()
                            break

                    chain_ids = chains.replace(';', ' ').replace(',', ' ').split() or ['']  # the second part deals with an empty chain id
                    self._add_chains(structure, chain_name, mol_id, *chain_ids)

            elif line.startswith('REMARK   2 RESOLUTION'):
                structure.resolution = self._read_resolution(line)
                                        
            elif line.startswith('SEQRES'):
                chain_id, residues = self._parse_seqres_line(line, structure)
                chain = structure.chains[chain_id]
                
                for residue in residues:
                    chain.residues.append(residue)
                    if chain.residues.last_index != residue.rank:
                        raise HeaderFormatError("Malformed SEQRES")

            elif line.startswith('MODEL') or line.startswith('ATOM'):
                break
                        
        return structure  
    
    def _add_chains(self, structure, name, mol_id, *chain_ids):

        name = name.strip().rstrip(";")
        
        for chain in chain_ids:
            new_chain = csb.bio.structure.Chain(chain, type=SequenceTypes.Unknown,
                                                name=name, accession=structure.accession)
            new_chain.molecule_id = mol_id
            try:
                structure.chains.append(new_chain)
            except csb.bio.structure.DuplicateChainIDError:
                raise HeaderFormatError('Chain {0} is already defined.'.format(new_chain.id))
    
    def _read_resolution(self, line):
        """
        @return: resolution
        @rtype: float or None
        """
        res = re.search("(\d+(?:\.\d+)?)\s+ANGSTROM", line)
        
        if res and res.groups():
            return float(res.group(1))
        else:
            return None         

    def _parse_seqres_line(self, line, structure):
        """
        Parse a SEQRES line, build and return newly constructed residues.
        If the current sequence type of the chain is unknown, try to guess it
        before parsing the residues.
        
        @return: parsed chain_id and L{Residue}s
        @rtype: 2-tuple: (str, iterable of L{Residue})
        """
        residues = []
        
        rownum = int(line[7:10])
        chain_id = line[11].strip()
        labels = line[18:].split()
        
        if chain_id not in structure.chains:
            raise HeaderFormatError('Chain {0} is undefined'.format(chain_id))
        
        chain = structure.chains[chain_id]
        
        if chain.type == SequenceTypes.Unknown:
            chain.type = self.guess_chain_type(labels)

        for rn, label in enumerate(labels):
            rank = rownum * 13 - (13 - (rn + 1))
            rtype = self.parse_residue_safe(label, as_type=chain.type)
            
            residue = csb.bio.structure.Residue.create(chain.type, rank=rank, type=rtype)
            residue.label = label
            residues.append(residue)
            
        return chain_id, residues
    

class PDBHeaderParser(RegularStructureParser):
    """
    Ultra fast PDB HEADER parser. Does not read any structural data.
    """
    
    def _parse_atoms(self, structure, model):
        pass
    
    def _parse_ss(self, structure):
        pass
    
    def _parse_header(self, model):
        return super(PDBHeaderParser, self)._parse_header(model)
        

class LegacyStructureParser(AbstractStructureParser):
    """
    This is a customized PDB parser, which is designed to read both sequence and
    atom data from the ATOM section. This is especially useful when parsing PDB
    files without a header.
    """
    
    def _parse_header(self, model):
        """
        Initialize a structure with residues from the ATOMs section.
        
        @param model: model identifier (e.g. if multiple models exist)
        @type model: str 
        
        @return: a L{csb.bio.structure.Structure} instance with properly 
                 initialized residues from ATOMs under the specified C{model}.
        @rtype: L{csb.bio.structure.Structure}
        """
         
        self._stream.seek(0)
        in_atom = False
        has_atoms = False
        has_model = False
        chains = csb.core.OrderedDict()

        header = next(self._stream)
        if header.startswith('HEADER'):
            structure = csb.bio.structure.Structure(header.split()[-1])
        else:
            self._stream.seek(0)
            structure = csb.bio.structure.Structure('NONE')

        structure.model_id = None
        
        while True:

            try:
                line = next(self._stream)
            except StopIteration:
                break

            if line.startswith('MODEL'):
                if has_model:
                    break
                else:
                    self._parse_model_line(line, structure, model)
                    model = structure.model_id
                    has_model = True

            elif line.startswith('ATOM') \
                     or (in_atom and line.startswith('HETATM')):
                    in_atom = True
                    has_atoms = True
                    
                    seq_number = self._read_sequence_number(line)
                    ins_code = self._read_insertion_code(line)
                    residue_id = (seq_number, ins_code)
                    label = self._read_residue_raw(line)
                    chain_id = self._read_chain_id(line)
                    
                    if chain_id not in chains:
                        chains[chain_id] = csb.core.OrderedDict()
                        self._add_chain(structure, chain_id)
                        
                    if residue_id not in chains[chain_id]:
                        chains[chain_id][residue_id] = label
                        chain = structure.chains[chain_id] 
                        if chain.type == SequenceTypes.Unknown:
                            self._fix_chain(chain, label)

            elif in_atom and line.startswith('TER'):
                in_atom = False
            elif line.startswith('ENDMDL'):
                break
            elif line.startswith('END'):
                break

        if not has_atoms:
            raise HeaderFormatError("Can't parse legacy structure: no ATOMs found")                                     

        for chain in structure.items:    
            self._build_chain(chain, chains[chain.id])
        
        return structure
    
    def _add_chain(self, structure, chain_id):

        new_chain = csb.bio.structure.Chain(chain_id,
                                            type=SequenceTypes.Unknown, 
                                            accession=structure.accession)
        new_chain.molecule_id = '1'
        structure.chains.append(new_chain)      
                     
    def _build_chain(self, chain, residues):
            
        for residue_id, label in residues.items():
            rank = (chain.residues.last_index or 0) + 1
            
            rname = self.parse_residue_safe(label, as_type=chain.type)
            residue = csb.bio.structure.Residue.create(chain.type, rank=rank, type=rname)
            residue.label = label
            residue.id = residue_id
            chain.residues.append(residue)

    def _fix_chain(self, chain, probe):
                
        try:
            chain.type = self.guess_sequence_type(probe)
        except UnknownPDBResidueError:
            pass  
    
    def _map_residues(self, structure, residues):

        for chain in structure.items:
            for residue_info in residues[chain.id]:
                try:
                    residue = chain.find(residue_info.sequence_number, residue_info.insertion_code)
                    for atom in residue_info.atoms:
                        residue.atoms.append(atom)
                    
                except csb.bio.structure.EntityNotFoundError:
                    pass         
                
    
StructureParser = AbstractStructureParser.create_parser
"""
Alias for L{AbstractStructureParser.create_parser}.
"""


class ResidueInfo(object):
    """
    High-performance struct, which functions as a container for unmapped
    L{Atom}s.
    
    @note: This object must implement the L{csb.bio.sequence.ResidueInfo}
           interface. This is not enforced through inheritance solely
           to save some CPU (by exposing public fields and no properties).
           However, on an abstract level this object is_a ResidueInfo
           and is used to build L{AbstractSequence}s.
    """
    __slots__ = ['chain', 'rank', 'id' , 'sequence_number', 'insertion_code', 'type', 'label', 'atoms']
    
    def __init__(self, chain, rank, id, seq_number, ins_code, type, label):
        
        self.chain = chain
        self.rank = rank
        self.id = id
        self.sequence_number = seq_number
        self.insertion_code = ins_code
        self.type = type
        self.label = label
        self.atoms = []
        
    @property
    def is_modified(self):
        
        if self.type.enum is SequenceAlphabets.Nucleic:
            return self.label != str(self.type)
        else:
            return self.label != repr(self.type)

class SparseChainSequence(csb.bio.sequence.ChainSequence):
    """
    Sequence view for reference (SEQRES) or sparse (ATOM) PDB chains.
    The residue instances passed to the constructor must be
    L{csb.bio.structure.Residue} or L{csb.bio.io.wwpdb.ResidueInfo} objects.
    
    See L{csb.bio.sequence.AbstractSequence} for details.
    """
            
    def _add(self, residue):
        
        if not isinstance(residue, (csb.bio.structure.Residue, ResidueInfo)):
            raise TypeError(residue)
        else:
            self._residues.append(residue)
            
    def _get(self, rank):
        return self._residues[rank - 1]
    
    @staticmethod
    def create(chain):
        """
        Create a new L{SparseChainSequence} from existing L{Chain}.
        
        @type chain: L{csb.bio.structure.Chain}
        @rtype: L{SparseChainSequence}
        """        
        return SparseChainSequence(
                chain.entry_id, chain.header, chain.residues, chain.type)    
    
    
class AbstractResidueMapper(object):
    """
    Defines the base interface of all residue mappers, used to align PDB ATOM
    records to the real (SEQRES) sequence of a chain.
    """
    
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def map(self, sparse, reference):
        """
        Map L{sparse}'s residues to L{reference}. Return all C{sparse} residues,
        aligned over C{reference}, with artificial gap residues inserted at
        relevant positions. The resulting sequence of sparse residues will
        always have the same length as the C{reference} sequence.
        
        @note: C{sparse}'s ranks won't be touched because the C{rank} property
        of the underlying L{ResidueInfo} implementation is not necessarily r/w.
        
        @param sparse: sparse sequence (e.g. derived from ATOMS records)
        @type sparse: L{SparseChainSequence}
        @param reference: reference, complete sequence
                          (e.g. derived from SEQRES records)
        @type reference: L{SparseChainSequence}
        
        @return: all C{sparse} residues, optimally aligned over C{reference}
                 (with gaps)
        @rtype: L{SparseChainSequence}
        
        @raise ResidueMappingError: if the specified sequences are not alignable
        """
        pass
    
    def create_gap(self, alphabet=SequenceAlphabets.Protein):
        """
        Create and return a new gap residue.
        
        @param alphabet: sequence alphabet; a member of L{SequenceAlphabets}
                         which has GAP item
        @type alphabet: L{enum}
        
        @rtype: L{ResidueInfo}
        """
        return ResidueInfo(None, -1, None, None, None, alphabet.GAP, "-")
    
    def _build(self, sparse, aligned):
        
        return SparseChainSequence(
                    sparse.id, sparse.header, aligned, sparse.type)        

class FastResidueMapper(AbstractResidueMapper):
    """
    RegExp-based residue mapper. Fails on heavily malformed input (i.e. it cannot
    insert gaps in the C{reference}), but it is very fast (linear) and memory
    efficient. 
    """
    
    MAX_FRAGMENTS = 20
    
    MIN_UNICODE_CHAR = 300
    FORBIDDEN_CHARS = set('^.*?()-')
    
    CODEC = "utf-8"
    DELIMITER = ").*?(".encode(CODEC).decode(CODEC)
    PATTERN = "^.*?({0}).*?$".encode(CODEC).decode(CODEC)
    
    def __init__(self):
        self._charcode = FastResidueMapper.MIN_UNICODE_CHAR
        self._cache = {} 
    
    def map(self, sparse, reference):

        aligned = []
        mapping = {}

        residues = list(sparse.residues)

        pattern = self._build_pattern(residues)
        seqres = self._encode_sequence(reference)

        matches = re.match(pattern, seqres)
                
        if matches:
            unmapped_item = -1
            
            for fn, fragment in enumerate(matches.groups(), start=1):
                assert fragment != ''
                
                for offset in range(1, len(fragment) + 1):
                    unmapped_item += 1
                    rank = matches.start(fn) + offset 
                    
                    mapped_residue = residues[unmapped_item] 
                    real_residue = reference.residues[rank]                    
                    
                    assert real_residue.type == mapped_residue.type
                    mapping[real_residue] = mapped_residue                    
        
        else:
            raise ResidueMappingError("Can't map ATOM records")
        
        for rank, residue in enumerate(reference.residues, start=1):
            if residue in mapping:
                aligned.append(mapping[residue])
            else:
                aligned.append(self.create_gap(sparse.alphabet))
        
        assert len(aligned) == reference.length
        return self._build(sparse, aligned)
    
    def _build_pattern(self, residues):
        """
        Build and return a sparse regular rexpression for C{residues}.
        """

        fragments = []
        
        for rn, r in enumerate(residues):
            res_name = self._encode(r)
            
            if rn == 0:
                # First residue, start a new fragment:
                fragments.append([res_name])
            elif r.insertion_code: # and not residues[rn - 1].insertion_code:
                # If residue i has an insertion code, initiate a new fragment:
                fragments.append([res_name])
            elif r.sequence_number - residues[rn - 1].sequence_number in (0, 1, -1):
                # If the seq numbers of residues [i-1, i] are consecutive, extend the last fragment:
                fragments[-1].append(res_name)
            else:
                # They are not consecutive, so we better start a new fragment:
                fragments.append([res_name])
        
        for i, frag in enumerate(fragments):
            fragments[i] = ''.join(frag)
        if len(fragments) > FastResidueMapper.MAX_FRAGMENTS:
            # Wow, that's a lot of fragments. Better use a different mapper
            raise ResidueMappingError("Can't map chain with large number of fragments")
        
        blocks = FastResidueMapper.DELIMITER.join(fragments)
        pattern = FastResidueMapper.PATTERN.format(blocks)
        
        return pattern
    
    def _encode(self, r):
        """
        Return a unique single-letter representation of C{r.type}. 
        """
        if not r.is_modified:
            return str(r.type)
        else:             
            return self._register_label(r.label)
        
    def _encode_sequence(self, s):
        return ''.join(map(self._encode, s.residues))
            
    def _register_label(self, label):
        """
        Assign a new unicode character to C{label} and cache it.

        @return: cached single-letter representation of label.
        @rtype: unicode char
        """

        if label not in self._cache:
            if set(label).intersection(FastResidueMapper.FORBIDDEN_CHARS):
                raise ResidueMappingError("Invalid residue label")
            
            self._charcode += 1
            code = self._charcode
            self._cache[label] = csb.io.unichr(code)
            
        return self._cache[label]    
                
        
class RobustResidueMapper(AbstractResidueMapper):
    """
    Exhaustive residue mapper, which uses Needleman-Wunsch global alignment.
    Much slower (quadratic), but fail-proof even with incompatible sequences
    (can insert gaps in both the C{sparse} and the C{reference} sequence).
    
    @param match: score for a match
    @type match: float
    @param mismatch: score for a mismatch (by default mismatches are heavily
                     penalized, while gaps are allowed)
    @type mismatch: float
    @param gap: gap penalty
    @type gap: float  
    """
    
    class GlobalAligner(alignment.GlobalAlignmentAlgorithm):
        
        def _sequence(self, s):
            return [r.label for r in s.residues]
           
    
    def __init__(self, match=1, mismatch=-10, gap=0):
        
        scoring = alignment.IdentityMatrix(match=match, mismatch=mismatch)
        aligner = RobustResidueMapper.GlobalAligner(scoring=scoring, gap=gap)
        
        self._aligner = aligner
        
    def map(self, sparse, reference):
        
        aligned = []
        ali = self._aligner.align(sparse, reference)
        
        if ali.is_empty:
            raise ResidueMappingError("Global alignment failed")
        
        for mapped, residue in zip(ali.query, ali.subject):
            
            if residue.type == reference.alphabet.GAP:
                continue
            elif mapped.type == sparse.alphabet.GAP:
                aligned.append(self.create_gap(sparse.alphabet))
            else:
                aligned.append(mapped)
                
        return self._build(sparse, aligned)

class CombinedResidueMapper(AbstractResidueMapper):
    """
    The best of both worlds: attempts to map the residues using
    L{FastResidueMapper}, but upon failure secures success by switching to
    L{RobustResidueMapper}. 
    """
    FAST = FastResidueMapper()
    ROBUST = RobustResidueMapper()
    
    def map(self, sparse, reference):
        
        try:
            return CombinedResidueMapper.FAST.map(sparse, reference)        
        except ResidueMappingError:
            return CombinedResidueMapper.ROBUST.map(sparse, reference)
                

class FileBuilder(object):
    """
    Base abstract files for all structure file formatters.
    Defines a common step-wise interface according to the Builder pattern.
    
    @param output: output stream (this is where the product is constructed)
    @type output: stream
    """
    
    __metaclass__ = ABCMeta

    def __init__(self, output):

        if not hasattr(output, 'write'):
            raise TypeError(output)
        
        def isnull(this, that, null=None):
            if this is null:
                return that
            else:
                return this        

        self._out = output
        self._isnull = isnull
        
    @property
    def output(self):
        """
        Destination stream
        @rtype: stream
        """
        return self._out
    
    @property
    def isnull(self):
        """
        ISNULL(X, Y) function
        @rtype: callable
        """
        return self._isnull
    
    def write(self, text):
        """
        Write a chunk of text
        """
        self._out.write(text)   

    def writeline(self, text):
        """
        Write a chunk of text and append a new line terminator
        """        
        self._out.write(text)
        self._out.write('\n')
                    
    @abstractmethod
    def add_header(self, master_structure):
        pass
    
    @abstractmethod
    def add_structure(self, structure):
        pass
    
    def finalize(self):
        pass

class PDBFileBuilder(FileBuilder):
    """
    PDB file format builder.
    """
        
    def writeline(self, text):
        self.write('{0:80}\n'.format(text))
        
    def add_header(self, master):
        """
        Write the HEADER of the file using C{master}
        
        @type master: L{Structure}
        """

        isnull = self.isnull
        
        header = 'HEADER    {0:40}{1:%d-%b-%y}   {2:4}'
        self.writeline(header.format('.', datetime.datetime.now(), master.accession.upper()))
        
        molecules = { }
        
        for chain_id in master.chains:
            chain = master.chains[chain_id]
            if chain.molecule_id not in molecules:
                molecules[chain.molecule_id] = [ ]
            molecules[chain.molecule_id].append(chain_id)
        
        k = 0
        for mol_id in sorted(molecules):
            
            chains = molecules[mol_id]
            first_chain = master.chains[ chains[0] ]            
            
            self.writeline('COMPND {0:3} MOL_ID: {1};'.format(k + 1, isnull(mol_id, '0')))
            self.writeline('COMPND {0:3} MOLECULE: {1};'.format(k + 2, isnull(first_chain.name, '')))
            self.writeline('COMPND {0:3} CHAIN: {1};'.format(k + 3, ', '.join(chains)))
            k += 3
            
        for chain_id in master.chains:
            
            chain = master.chains[chain_id]
            res = [ r.label for r in chain.residues ]

            rn = 0
            for j in range(0, chain.length, 13):
                rn += 1
                residues = [ '{0:>3}'.format(r) for r in res[j : j + 13] ]
                self.writeline('SEQRES {0:>3} {1} {2:>4}  {3}'.format(
                                            rn, chain.id, chain.length, ' '.join(residues) ))
                
    def add_structure(self, structure):
        """
        Append a new model to the file
        
        @type structure: L{Structure}
        """

        isnull = self.isnull
         
        for chain_id in structure.chains:
        
            chain = structure.chains[chain_id]
            for residue in chain.residues:
        
                atoms = [ ]
                for an in residue.atoms:
                    atom = residue.atoms[an]
                    if isinstance(atom, csb.bio.structure.DisorderedAtom):
                        for dis_atom in atom: atoms.append(dis_atom)
                    else:
                        atoms.append(atom)
                atoms.sort()
                
                for atom in atoms:

                    alt = atom.alternate
                    if alt is True:
                        alt = 'A'
                    elif alt is False:
                        alt = ' '
                    
                    if atom.element:
                        element = repr(atom.element)
                    else:
                        element = ' '
                    self.writeline('ATOM  {0:>5} {1:>4}{2}{3:>3} {4}{5:>4}{6}   {7:>8.3f}{8:>8.3f}{9:>8.3f}{10:>6.2f}{11:>6.2f}{12:>12}{13:2}'.format(
                                        atom.serial_number, atom._full_name, isnull(alt, ' '), 
                                        residue.label, chain.id, 
                                        isnull(residue.sequence_number, residue.rank), isnull(residue.insertion_code, ' '), 
                                        atom.vector[0], atom.vector[1], atom.vector[2], isnull(atom.occupancy, 0.0), isnull(atom.bfactor, 0.0), 
                                        element, isnull(atom.charge, ' ') ))        

            self.writeline('TER')
        
    def finalize(self):
        """
        Add the END marker
        """
        self.writeline('END')
        self._out.flush()     

class PDBEnsembleFileBuilder(PDBFileBuilder):
    """
    Supports serialization of NMR ensembles.
    
    Functions as a simple decorator, which wraps C{add_structure} with
    MODEL/ENDMDL records.
    """

    def add_structure(self, structure):
        
        model_id = self.isnull(structure.model_id, 1)
        
        self.writeline('MODEL     {0:>4}'.format(model_id))       
        super(PDBEnsembleFileBuilder, self).add_structure(structure)
        self.writeline('ENDMDL')
        

class StructureProvider(object):
    """
    Base class for all PDB data source providers.
    
    Concrete classes need to implement the C{find} method, which abstracts the
    retrieval of a PDB structure file by a structure identifier. This is a hook
    method called internally by C{get}, but subclasses can safely override both
    C{find} and {get} to in order to achieve completely custom behavior. 
    """
    
    __metaclass__ = ABCMeta
    
    def __getitem__(self, id):
        return self.get(id)
    
    @abstractmethod
    def find(self, id):
        """
        Attempt to discover a PDB file, given a specific PDB C{id}.
            
        @param id: structure identifier (e.g. 1x80)
        @type id: str
        @return: path and file name on success, None otherwise
        @rtype: str or None
        """
        pass
    
    def get(self, id, model=None):
        """
        Discover, parse and return the PDB structure, corresponding to the
        specified C{id}.
        
        @param id: structure identifier (e.g. 1x80)
        @type id: str
        @param model: optional model identifier
        @type model: str
                
        @rtype: L{csb.bio.Structure}
        
        @raise StructureNotFoundError: when C{id} could not be found
        """
        pdb = self.find(id)
        
        if pdb is None:
            raise StructureNotFoundError(id)
        else:
            return StructureParser(pdb).parse_structure(model=model)
    
class FileSystemStructureProvider(StructureProvider):
    """
    Simple file system based PDB data source. Scans a list of local directories
    using pre-defined file name templates.
    
    @param paths: a list of paths
    @type paths: iterable or str 
    """
    
    def __init__(self, paths=None):
        
        self._templates = ['pdb{id}.ent', 'pdb{id}.pdb', '{id}.pdb', '{id}.ent']
        self._paths = csb.core.OrderedDict()
        
        if paths is not None:
            if isinstance(paths, csb.core.string):
                paths = [paths]
                            
            for path in paths:
                self.add(path)
                
    @property
    def paths(self):
        """
        Current search paths
        @rtype: tuple
        """
        return tuple(self._paths)
    
    @property
    def templates(self):
        """
        Current file name match templates
        @rtype: tuple
        """
        return tuple(self._templates)
        
    def add(self, path):
        """
        Register a new local C{path}.
        
        @param path: directory name
        @type path: str
        
        @raise IOError: if C{path} is not a valid directory
        """
        if os.path.isdir(path):
            self._paths[path] = path
        else:
            raise IOError(path)
        
    def add_template(self, template):
        """
        Register a custom file name name C{template}. The template must contain
        an E{lb}idE{rb} macro, e.g. pdbE{lb}idE{rb}.ent 
        
        @param template: pattern
        @type template: str 
        """
        
        if '{id}' not in template:
            raise ValueError('Template does not contain an "{id}" macro')
        
        if template not in self._templates:
            self._templates.append(template)
            
    def remove(self, path):
        """
        Unregister an existing local C{path}.
        
        @param path: directory name
        @type path: str
        
        @raise ValueError: if C{path} had not been registered
        """        
        if path not in self._paths:
            raise ValueError('path not found: {0}'.format(path))
        
        del self._paths[path]
            
    def find(self, id):
        
        for path in self._paths:
            for token in self.templates:
                fn = os.path.join(path, token.format(id=id))
                if os.path.exists(fn):
                    return fn
            
        return None
    
class RemoteStructureProvider(StructureProvider):
    """
    Retrieves PDB structures from a specified remote URL.
    The URL requested from remote server takes the form: <prefix>/<ID><suffix>
    
    @param prefix: URL prefix, including protocol
    @type prefix: str
    @param suffix: optional URL suffix (.ent by default)
    @type suffix: str
    """
    
    def __init__(self, prefix='http://www.rcsb.org/pdb/files/pdb', suffix='.ent'):
        
        self._prefix = None
        self._suffix = None
        
        self.prefix = prefix
        self.suffix = suffix
        
    @property
    def prefix(self):
        """
        Current URL prefix
        @rtype: str
        """
        return self._prefix
    @prefix.setter
    def prefix(self, value):
        self._prefix = value
        
    @property
    def suffix(self):
        """
        Current URL suffix
        @rtype: str
        """        
        return self._suffix
    @suffix.setter
    def suffix(self, value):
        self._suffix = value        
    
    def _find(self, id):
        
        try:
            return csb.io.urllib.urlopen(self.prefix + id + self.suffix)
        except:
            raise StructureNotFoundError(id)
        
    def find(self, id):
        
        stream = self._find(id)
        
        try:
            tmp = csb.io.TempFile(dispose=False)
            tmp.write(stream.read().decode('utf-8'))
            tmp.flush()
            return tmp.name
                
        except StructureNotFoundError:
            return None            
        finally:
            stream.close()         
        
    def get(self, id, model=None):
        
        stream = self._find(id)
        
        try:
            with csb.io.TempFile() as tmp:
                tmp.write(stream.read().decode('utf-8'))
                tmp.flush()
                return StructureParser(tmp.name).parse_structure(model=model)
        finally:
            stream.close()
        
class CustomStructureProvider(StructureProvider): 
    """
    A custom PDB data source. Functions as a user-defined map of structure
    identifiers and their corresponding local file names. 
    
    @param files: initialization dictionary of id:file pairs
    @type files: dict-like
    """
        
    def __init__(self, files={}):
        
        self._files = {}        
        for id in files:
            self.add(id, files[id])
            
    @property
    def paths(self):
        """
        List of currently registered file names
        @rtype: tuple
        """        
        return tuple(self._files.values())            

    @property
    def identifiers(self):
        """
        List of currently registered structure identifiers
        @rtype: tuple
        """          
        return tuple(self._files)
            
    def add(self, id, path):
        """
        Register a new local C{id}:C{path} pair.

        @param id: structure identifier
        @type id: str        
        @param path: path and file name
        @type path: str
        
        @raise IOError: if C{path} is not a valid file name
        """
        if os.path.isfile(path):
            self._files[id] = path
        else:
            raise IOError(path)
        
    def remove(self, id):
        """
        Unregister an existing structure C{id}.
        
        @param id: structure identifier
        @type id: str
        
        @raise ValueError: if C{id} had not been registered
        """
        if id not in self._files:
            raise ValueError(id)
        else:
            del self._files[id]
        
    def find(self, id):
        
        if id in self._files:
            return self._files[id]
        else:
            return None
    
def get(accession, model=None, prefix='http://www.rcsb.org/pdb/files/pdb'):
    """
    Download and parse a PDB entry.

    @param accession: accession number of the entry
    @type accession: str
    @param model: model identifier
    @type model: str
    @param prefix: download URL prefix
    @type prefix: str

    @return: object representation of the selected model
    @rtype: L{Structure}
    """
    return RemoteStructureProvider(prefix).get(accession, model=model)

def find(id, paths):
    """
    Try to discover a PDB file for PDB C{id} in C{paths}.
    
    @param id: PDB ID of the entry
    @type id: str
    @param paths: a list of directories to scan
    @type paths: list of str
    
    @return: path and file name on success, None otherwise
    @rtype: str
    """
    return FileSystemStructureProvider(paths).find(id)
    
    
class AsyncParseResult(object):
    
    def __init__(self, result, exception):
        
        self.result = result
        self.exception = exception
    
    def __repr__(self):
        return '<AsyncParseResult: result={0.result}, error={0.exception.__class__.__name__}>'.format(self)


def _parse_async(parser, file, model):
    p = parser(file)
    return p.parse_structure(model)


class AsyncStructureParser(object):
    """
    Wraps StructureParser in an asynchronous call. Since a new process is
    started by Python internally (as opposed to only starting a new thread),
    this makes the parser slower, but provides a way to set a parse timeout
    limit.
    
    If initialized with more than one worker, supports parallel parsing
    through the C{self.parse_async} method. 
    
    @param workers: number of worker threads (1 by default)
    @type workers: int
    """

    def __init__(self, workers=1):

        self._pool = None
        self._workers = 1
        
        if int(workers) > 0:
            self._workers = int(workers)
        else:
            raise ValueError(workers)
        
        self._recycle()
        
    def _recycle(self):

        if self._pool:
            self._pool.terminate()         
               
        self._pool = multiprocessing.Pool(processes=self._workers)

    def parse_structure(self, structure_file, timeout, model=None,
                        parser=RegularStructureParser):
        """
        Call StructureParser.parse_structure() in a separate process and return
        the output. Raise TimeoutError if the parser does not respond within
        C{timeout} seconds.
        
        @param structure_file: structure file to parse
        @type structure_file: str
        @param timeout: raise multiprocessing.TimeoutError if C{timeout} seconds
                        elapse before the parser completes its job
        @type timeout: int
        @param parser: any implementing L{AbstractStructureParser} class
        @type parser: type  
        
        @return: parsed structure
        @rtype: L{csb.structure.Structure}    
        """

        r = self.parse_async([structure_file], timeout, model, parser)
        if len(r) > 0:
            if r[0].exception is not None:
                raise r[0].exception
            else:
                return r[0].result
        return None
    
    def parse_async(self, structure_files, timeout, model=None,
                        parser=RegularStructureParser):
        """
        Call C{self.parse_structure} for a list of structure files
        simultaneously. The actual degree of parallelism will depend on the
        number of workers specified while constructing the parser object.
        
        @note: Don't be tempted to pass a large list of structures to this 
               method. Every time a C{TimeoutError} is encountered, the 
               corresponding worker process in the pool will hang until the
               process terminates on its own. During that time, this worker is
               unusable. If a sufficiently high number of timeouts occur, the 
               whole pool of workers will be unsable. At the end of the method
               however a pool cleanup is performed and any unusable workers
               are 'reactivated'. However, that only happens at B{the end} of
               C{parse_async}.
        
        @param structure_files: a list of structure files
        @type structure_files: tuple of str
        @param timeout: raise multiprocessing.TimeoutError if C{timeout} seconds
                        elapse before the parser completes its job
        @type timeout: int
        @param parser: any implementing L{AbstractStructureParser} class
        @type parser: type
        
        @return: a list of L{AsyncParseResult} objects
        @rtype: list     
        """

        pool =  self._pool
        workers = []
        results = []
        
        for file in list(structure_files):
            result = pool.apply_async(_parse_async, [parser, file, model])
            workers.append(result)
        
        hanging = False
        for w in workers:
            result = AsyncParseResult(None, None)
            try:
                result.result = w.get(timeout=timeout)
            except KeyboardInterrupt as ki:
                pool.terminate()
                raise ki
            except Exception as ex:
                result.exception = ex
                if isinstance(ex, multiprocessing.TimeoutError):
                    hanging = True                    
            results.append(result)
        
        if hanging:
            self._recycle()
            
        return results
