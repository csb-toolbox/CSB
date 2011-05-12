"""
PDB structure parsers.
"""

import re
import multiprocessing
import csb.bio.structure
import csb.pyutils
import csb.io

from abc import ABCMeta, abstractmethod
from csb.bio.sequence import SequenceTypes, SequenceAlphabets
from csb.bio.structure import ChemElements, SecStructures
from numpy import array


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


class UnknownPDBResidueError(PDBParseError):
    pass


class AbstractStructureParser(object):
    """
    A base PDB structure format-aware parser. Subclasses must implement the
    internal abstract method C{_parse_header} in order to complete the 
    implementation.

    @param structure_file: the input PD file to parse
    @type structure_file: str
    @param check_ss: if True, secondary structure errors will result in 
                     a L{PDBParseError} exception
    @type check_ss: bool    

    @raise IOError: when the input file cannot be found
    """

    __metaclass__ = ABCMeta

    """
    List of PDB structures which are known to cause parser hanging.
    """

    @staticmethod
    def create_parser(structure_file, check_ss=False):
        """
        A StructureParser factory, which instantiates and returns the proper parser
        object based on the contents of the PDB file.
        
        @param structure_file: the PDB file to parse. If the file contains a SEQRES
                               section, L{RegularStructureParser} is returned,
                               otherwise L{LegacyStructureParser} is instantiated. 
                               In the latter case LegacyStructureParser will read 
                               the sequence data directly from the ATOMs
        @type structure_file: str
        
        @return: a proper parser instance
        @rtype: L{AbstractStructureParser} subclass
        """
        has_seqres = False
        
        for line in open(structure_file):
            if line.startswith('SEQRES'):
                has_seqres = True
                break
        
        if has_seqres:
            return RegularStructureParser(structure_file)
        else:
            return LegacyStructureParser(structure_file)        

    def __init__(self, structure_file, check_ss=False):

        self._file = None
        self._stream = None
        self._check_ss = bool(check_ss)        

        self.filename = structure_file

    def __del__(self):
        try:
            self._stream.close()
        except:
            pass

    @property
    def filename(self):
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
        Find a list of available model identifiers in the structure.

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
                        raise csb.bio.structure.Broken3DStructureError('Duplicate model identifier: {0}'.format(model_id))
                    models.append(model_id)
                    check.add(model_id)

        if len(models) > 0:
            if not(min(check) == 1 and max(check) == len(models)):
                raise csb.bio.structure.Broken3DStructureError('Non-consecutive model identifier(s) encountered')                
            return models
        else:
            return [None]

    def guess_sequence_type(self, residue_name):
        """
        Try to guess what is the sequence type of a PDB C{residue_name}.

        @param residue_name: a PDB-conforming name of a residue
        @type residue_name: str

        @return: a L{SequenceTypes} enum member
        @rtype: L{csb.pyutils.EnumItem}

        @raise UnknownPDBResidueError: when there is no such PDB residue name
                                       into the catalog tables
        """
        if residue_name in PDB_AMINOACIDS:
            return SequenceTypes.Protein                                #@UndefinedVariable
        elif residue_name in PDB_NUCLEOTIDES:
            return SequenceTypes.NucleicAcid                            #@UndefinedVariable
        else:
            raise UnknownPDBResidueError(residue_name)

    def parse_residue(self, residue_name, as_type=None):
        """
        Try to parse a PDB C{residue_name} and return its closest 'normal'
        string representation. If a sequence type (C{as_type}) is defined,
        guess the alphabet based on that information, otherwise try first to
        parse it as a protein residue.

        @param residue_name: a PDB-conforming name of a residue
        @type residue_name: str
        @param as_type: suggest a sequence type (L{SequenceTypes} member)
        @type L{scb.pyutils.EnumItem}

        @return: a normalized residue name
        @rtype: str

        @raise UnknownPDBResidueError: when there is no such PDB residue name
                                       into the catalog table(s)
        """
        if as_type is None:
            as_type = self.guess_sequence_type(residue_name)         

        try:
            if as_type == SequenceTypes.Protein:                            #@UndefinedVariable
                return PDB_AMINOACIDS[residue_name]
            elif as_type == SequenceTypes.NucleicAcid:                      #@UndefinedVariable
                return PDB_NUCLEOTIDES[residue_name]
            else:
                raise UnknownPDBResidueError(residue_name)
        except KeyError:
            raise UnknownPDBResidueError(residue_name)

    def parse_residue_safe(self, residue_name, as_type):
        """
        Same as C{parse_residue}, but returns UNK/Any instead of raising
        UnknownPDBResidueError.

        @param residue_name: a PDB-conforming name of a residue
        @type residue_name: str
        @param as_type: suggest a sequence type (L{SequenceTypes} member)
        @type L{scb.pyutils.EnumItem}

        @return: a normalized residue name
        @rtype: str
        """
        try:
            return self.parse_residue(residue_name, as_type)
        except UnknownPDBResidueError:
            if as_type == SequenceTypes.Protein:                        #@UndefinedVariable
                return repr(SequenceAlphabets.Protein.UNK)              #@UndefinedVariable
            elif as_type == SequenceTypes.NucleicAcid:                  #@UndefinedVariable
                return repr(SequenceAlphabets.Nucleic.Any)              #@UndefinedVariable
            else:
                return repr(SequenceAlphabets.Unknown.UNK)              #@UndefinedVariable

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

        @raise ValueError: when an invalid model ID is specified
        @raise Broken3DStructureError: on parse error (corrput file)
        """
        if model is not None:
            model = int(model)
                    
        structure = self._parse_header(model)
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
            models = map(int, models)
        
        ensemble = csb.bio.structure.Ensemble()
        for model_id in models:
            model = self.parse_structure(model_id)
            ensemble.models.append(model)
        
        return ensemble

    @abstractmethod
    def _parse_header(self, model):
        """
        An abstract method, which subclasses must use to customize the way the
        PDB header (or the absence of a such) is handled. The implementation
        must rely on reading character data from the current internal 
        self._stream and must return a new L{csb.bio.structure.Structure} 
        instance with properly initialized header data: accession, model, 
        molecule identifiers, chains and residues. This structure object is
        then automatically passed to C{_parse_atoms}, which parses and attaches
        the atom data to the residues in the structure.
        """
        pass

    def _scroll_model(self, model, stream):
        """
        Scroll the C{stream} to the specified C{model}.
        """
        
        while True:
            try:
                line = stream.next()
            except StopIteration:
                raise ValueError('No such model {0} in the structure.'.format(model))
            
            if line.startswith('MODEL'):
                model_id = int(line[10:14])
                if model == model_id:
                    return model_id      
        
    def _parse_atoms(self, structure, model):
        """
        Parse the ATOMs from the specified C{model} and attach them to the
        C{structure}.
        """
        
        structure.model_id = None

        atoms = dict( (chain, []) for chain in structure.chains )
        chains = set()
        in_ligands = False                
        in_atom = False

        self._stream.seek(0)
        while True:

            try:
                line = self._stream.next()
            except StopIteration:
                break

            if line.startswith('MODEL'):
                if model and model != int(line[10:14]):
                    self._scroll_model(model, self._stream)
                    structure.model_id = model
                else:
                    model = int(line[10:14])
                    structure.model_id = model

            elif line.startswith('ATOM') \
                     or (line.startswith('HETATM') and not in_ligands):
                    in_atom = True

                    rank = int(line[22:26])
                    serial_number = int(line[6:11])
                    name = line[12:16]
                    x, y, z = line[30:38], line[38:46], line[46:54]
                    vector = array([float(x), float(y), float(z)])
                    element = line[76:78].strip()
                    if element:
                        try:
                            element = csb.pyutils.Enum.parsename(ChemElements, element)
                        except csb.pyutils.EnumMemberError as ee:
                            if element in ('D', 'X'):
                                element = ChemElements.x                            #@UndefinedVariable
                            else:
                                raise ee
                    else:
                        element = None

                    atom = csb.bio.structure.Atom(serial_number, name, element,
                                                  vector)

                    atom._het = line.startswith('HETATM')
                    atom._rank = rank
                    atom._sequence_number = int(line[22:26].strip())
                    atom._residue_id = str(atom._sequence_number)
                    atom._insertion_code = line[26].strip()
                    if not atom._insertion_code:
                        atom._insertion_code = None
                    else:
                        atom._residue_id += atom._insertion_code

                    atom.alternate = line[16].strip()
                    if not atom.alternate:
                        atom.alternate = None

                    atom._chain = line[21].strip()
                    if atom._chain not in structure.chains:
                        raise csb.bio.structure.Broken3DStructureError('No such chain: {0}'.format(atom._chain))
                    chains.add(atom._chain)
                    residue_name = line[17:20].strip()
                    residue_name = self.parse_residue_safe(residue_name, as_type=structure.chains[atom._chain].type)
                    if structure.chains[atom._chain].type == SequenceTypes.NucleicAcid:                               #@UndefinedVariable
                        atom._residue_name = csb.pyutils.Enum.parsename(SequenceAlphabets.Nucleic, residue_name)
                    else:
                        atom._residue_name = csb.pyutils.Enum.parsename(SequenceAlphabets.Protein, residue_name)

                    atom.occupancy = float(line[54:60])
                    atom.temperature = float(line[60:66])

                    charge = line[78:80].strip()
                    if charge:
                        if charge in ('+', '-'):
                            charge += '1'
                        if charge[-1] in ('+', '-'):
                            charge = charge[::-1]
                        atom.charge = int(charge)

                    atoms[atom._chain].append(atom)

            elif in_atom and line.startswith('TER'):
                in_atom = False
                if len(chains) == len(structure.chains):
                    in_ligands = True

            elif line.startswith('ENDMDL'):
                break

            elif line.startswith('END'):
                break

        if structure.model_id != model:
            raise ValueError('No such model {0} in the structure.'.format(model))

        self._map_residues(structure, atoms)        

    def _map_residues(self, structure, atoms):

        assert set(structure.chains) == set(atoms.keys())

        for chain in structure.chains:

            subject = structure.chains[chain].sequence
            query = ''
            fragments = []

            seq_numbers = []
            lookup = {}

            i = -1
            for a in atoms[chain]:
                if a._residue_id not in lookup:
                    i += 1
                    lookup[a._residue_id] = [a._sequence_number, a._insertion_code]
                    seq_numbers.append(a._residue_id)
                    res_name = a._residue_name.value
                    if a._het:
                        # if it is a HET atom, initiate an optional fragment
                        fragments.append([res_name, '?'])                        
                    elif i == 0 or a._sequence_number - lookup[seq_numbers[i - 1]][0] not in (0, 1, -1):
                        # if residues [i, i-1] are not consecutive or 'overlapping', initiate a new fragment:
                        fragments.append([res_name])
                    else:
                        # then they are consecutive
                        if fragments[-1][-1] == '?':
                            # but the prev was optional (maybe HET), all optionals *MUST* reside in 
                            # singleton fragments, so start a new fragment
                            fragments.append([res_name])
                        else:
                            # append the new residue to the end of the last fragment                            
                            fragments[-1].append(res_name)
            
            for i, frag in enumerate(fragments):
                fragments[i] = ''.join(frag)
            if len(fragments) > 100:
                # Python's regex engine will crash with more than 100 groups
                raise PDBParseError("Can't map structure with more than 100 fragments in ATOMS") 
            query = '^.*?({0}).*?$'.format(').*?('.join(fragments))
            
            matches = re.match(query, subject)
            
            if matches:
                seq_numitem = -1
                for frag_number, frag in enumerate(matches.groups(), start=1):
                    if frag is '':
                        # optional fragment, which did not match (NOTE: all optionals must occupy 
                        # their own fragments of length 1 residue)
                        seq_numitem += 1    # this 1 implies that all optional fragments are 1 residue long
                    else:
                        for i, dummy in enumerate(frag, start=1):
                            seq_numitem += 1
                            # lookup[res_id] is finally set to give the final rank of residue under id res_id:
                            try:
                                lookup[ seq_numbers[seq_numitem] ] = matches.start(frag_number) + i
                            except:
                                raise 

                fixed_residue = None
                for atom in atoms[chain]:
                    if not isinstance(lookup[atom._residue_id], int):
                        continue                                    # this atom was not mapped (e.g. HET)
                    atom._rank = lookup[atom._residue_id]
                    residue = structure.chains[chain].residues[atom._rank]
                    if residue is not fixed_residue:
                        residue.id = atom._sequence_number, atom._insertion_code
                        fixed_residue = residue

                    assert str(residue.type) == subject[atom._rank - 1]
                    residue.atoms.append(atom)
                    
                    del atom._rank
                    del atom._insertion_code
                    del atom._sequence_number
                    del atom._chain
                    del atom._residue_id
                    del atom._residue_name
            else:
                raise csb.bio.structure.Broken3DStructureError('Could not map structure to sequence.')

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
                line = self._stream.next()
            except StopIteration:
                break

            if line.startswith('HELIX'):
                
                chain = structure.chains[line[19].strip()]
                if chain.id not in elements:
                    elements[chain.id] = []
                if chain.id != line[31].strip():
                    if self._check_ss:
                        raise PDBParseError('Helix {0} spans multiple chains'.format(line[7:10]))
                    else:
                        continue                
                try:
                    startres = chain.find(line[21:25].strip(), line[25].strip())                
                    endres = chain.find(line[33:37].strip(), line[37].strip())
                except csb.pyutils.ItemNotFoundError as ex:
                    if self._check_ss:
                        raise PDBParseError('Helix {0} refers to an undefined residue ID: {1}'.format(line[7:10], str(ex)))
                    else:
                        continue
                if not startres.rank <= endres.rank:
                    if self._check_ss:
                        raise PDBParseError('Helix {0} is out of range'.format(line[7:10]))
                    else:
                        continue                    
                helix = csb.bio.structure.SecondaryStructureElement(startres.rank, endres.rank, SecStructures.Helix)        #@UndefinedVariable
                elements[chain.id].append(helix)
            
            if line.startswith('SHEET'):
                
                chain = structure.chains[line[21].strip()]
                if chain.id not in elements:
                    elements[chain.id] = []                
                if chain.id != line[32].strip():
                    if self._check_ss:                    
                        raise PDBParseError('Sheet {0} spans multiple chains'.format(line[7:10]))
                    else:
                        continue
                try:
                    startres = chain.find(line[22:26].strip(), line[26].strip())                
                    endres = chain.find(line[33:37].strip(), line[37].strip())
                except csb.pyutils.ItemNotFoundError as ex:
                    if self._check_ss:                    
                        raise PDBParseError('Sheet {0} refers to an undefined residue ID: {1}'.format(line[7:10], str(ex)))
                    else:
                        continue
                if not startres.rank <= endres.rank:
                    if self._check_ss:
                        raise PDBParseError('Sheet {0} is out of range'.format(line[7:10]))
                    else:
                        continue
                strand = csb.bio.structure.SecondaryStructureElement(startres.rank, endres.rank, SecStructures.Strand)      #@UndefinedVariable
                elements[chain.id].append(strand)         
            
            elif line.startswith('MODEL') or line.startswith('ATOM'):
                break        
            
        for chain_id in elements:
            ss = csb.bio.structure.SecondaryStructure()
            for e in elements[chain_id]:
                ss.append(e)
            structure.chains[chain_id].secondary_structure = ss
        

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
        
        @raise ValueError: if the stream has no HEADER at byte 0
        """
        
        self._stream.seek(0)

        header = self._stream.next()
        if not header.startswith('HEADER'):
            raise ValueError('Does not look like a regular PDB file.')

        structure = csb.bio.structure.Structure(header.split()[-1])

        while True:

            try:
                line = self._stream.next()
            except StopIteration:
                break

            if line.startswith('COMPND'):
                if line[10:].lstrip().startswith('MOL_ID:'):
                    mol_id = int(line[18:].replace(';', '').strip())
                    chain_name = ''
                    chains = ''
                    while line.startswith('COMPND'):
                        line = self._stream.next()
                        if line.split()[2].startswith('MOLECULE:'):
                            chain_name += line[20:].strip()
                            while not chain_name.endswith(';'):
                                line = self._stream.next()
                                if not line.startswith('COMPND'):
                                        break
                                chain_name += ' ' + line[11:].strip()
                        else:
                            while not line.split()[2].startswith('CHAIN:'):
                                line = self._stream.next()
                                if not line.startswith('COMPND'):
                                    raise csb.bio.structure.Broken3DStructureError('COMPND section is not parsable: missing chain identifier.')
                            chains = line[17:].strip()
                            while not chains.endswith(';'):
                                line = self._stream.next()
                                if not line.startswith('COMPND'):
                                    break
                                chains += ', ' + line[11:].strip()
                            break

                    chain_name = chain_name.strip()[:-1]
                    for chain in chains.replace(';', ' ').replace(',', ' ').split() or ['']:  # the second part deals with an empty chain id
                        new_chain = csb.bio.structure.Chain(chain, type=SequenceTypes.Unknown, name=chain_name, accession=structure.accession)
                        new_chain.molecule_id = mol_id
                        try:
                            structure.chains.append(new_chain)
                        except csb.pyutils.DuplicateKeyError:
                            raise csb.bio.structure.Broken3DStructureError('Chain {0} is already defined.'.format(new_chain.id))
                        
            elif line.startswith('SEQRES'):
                # Correct handling of empty chain id
                seq_fields = [line[7:10], line[11], line[13:17] ]
                seq_fields.extend(line[18:].split())

                rank_base = int(seq_fields[0])
                chain_id = seq_fields[1].strip()

                if structure.chains[chain_id].type == SequenceTypes.Unknown:
                    inner_residuerank = int(len(seq_fields[3:]) / 2) + 3
                    for i in [inner_residuerank, 3, -1]:
                        try:
                            structure.chains[chain_id].type = self.guess_sequence_type(seq_fields[i])
                            break
                        except UnknownPDBResidueError:
                            pass

                for i, residue_name in enumerate(seq_fields[3:]):
                    rank = rank_base * 13 - (13 - (i + 1))
                    residue_name = self.parse_residue_safe(residue_name, as_type=structure.chains[chain_id].type)
                    residue = csb.bio.structure.Residue.create(structure.chains[chain_id].type, rank=rank, type=residue_name)
                    structure.chains[chain_id].residues.append(residue)
                    assert structure.chains[chain_id].residues.last_index == rank        

            elif line.startswith('MODEL') or line.startswith('ATOM'):
                break
                        
        return structure   


class LegacyStructureParser(AbstractStructureParser):
    """
    This is a customized PDB parser, which is designed to read both sequence and
    atom data from the ATOM section. This is especially useful when parsing PDB
    files without a HEADER section.
    """
    
    def _parse_header(self, model):
        """
        Initialize a structure with residues from the ATOM section.
        
        @param model: model identifier (e.g. if multiple models exist)
        @type model: str 
        
        @return: a L{csb.bio.structure.Structure} instance with properly 
                 initialized residues from ATOMs under the specified C{model}.
        @rtype: L{csb.bio.structure.Structure}
        """
         
        self._stream.seek(0)
        in_atom = False
        chains = csb.pyutils.OrderedDict()

        header = self._stream.next()
        if header.startswith('HEADER'):
            structure = csb.bio.structure.Structure(header.split()[-1])
        else:
            self._stream.seek(0)
            structure = csb.bio.structure.Structure('NONE')

        structure.model_id = None
        
        while True:

            try:
                line = self._stream.next()
            except StopIteration:
                break

            if line.startswith('MODEL'):
                if model and model != int(line[10:14]):
                    self._scroll_model(model, self._stream)
                    structure.model_id = model
                else:
                    model = int(line[10:14])
                    structure.model_id = model

            elif line.startswith('ATOM') \
                     or (in_atom and line.startswith('HETATM')):
                    in_atom = True
                    
                    residue_id = line[22:27].strip()
                    residue_name = line[17:20].strip()
                    chain_id = line[21].strip()
                    
                    if chain_id not in chains:
                        chains[chain_id] = csb.pyutils.OrderedDict()
                        
                        new_chain = csb.bio.structure.Chain(
                                            chain_id, 
                                            type=SequenceTypes.Unknown, 
                                            accession=structure.accession)
                        new_chain.molecule_id = '1'
                        structure.chains.append(new_chain)                        
                        
                    if residue_id not in chains[chain_id]:
                        chains[chain_id][residue_id] = residue_name
                        
                        if structure.chains[chain_id].type == SequenceTypes.Unknown:
                            try:
                                structure.chains[chain_id].type = self.guess_sequence_type(residue_name)
                            except UnknownPDBResidueError:
                                pass

            elif in_atom and line.startswith('TER'):
                in_atom = False

            elif line.startswith('ENDMDL'):
                break

            elif line.startswith('END'):
                break                            
        
        for chain_id in structure.chains:
            for residue_id in chains[chain_id]:
                residue_name = chains[chain_id][residue_id]
                rank = (structure.chains[chain_id].residues.last_index or 0) + 1
                
                residue_name = self.parse_residue_safe(residue_name, as_type=structure.chains[chain_id].type)
                residue = csb.bio.structure.Residue.create(structure.chains[chain_id].type, rank=rank, type=residue_name)
                structure.chains[chain_id].residues.append(residue)                                  
                        
        return structure    

    
StructureParser = AbstractStructureParser.create_parser


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
    import urllib2

    pdb = csb.io.TempFile()

    browser = urllib2.urlopen(prefix + accession.lower() + '.ent')
    pdb.write(browser.read())
    pdb.flush()

    return StructureParser(pdb.name).parse_structure(model)


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
    import os
    
    if isinstance(paths, basestring):
        paths = [paths]
        
    for path in paths:
        for token in ['pdb{id}.ent', 'pdb{id}.pdb', '{id}.pdb', '{id}.ent']:
            fn = os.path.join(path, token.format(id=id))
            if os.path.exists(fn):
                return fn
        
    return None 


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
