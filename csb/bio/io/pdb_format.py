import os
import re
import multiprocessing
import csb.bio.structure

from csb.bio.sequence import SequenceTypes, SequenceAlphabets
from csb.bio.structure import ChemElements


class PDBParseError(ValueError):
    pass


class UnknownPDBResidueError(PDBParseError):
    pass


class StructureParser(object):
    """
    A PDB structure format-aware parser.

    @param structure_file: the input PD file to parse
    @type structure_file: str

    @raise IOError: when the input file cannot be found
    """

    _pdb_aminoacids = {
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

    _pdb_nucleotides = {
        'DA': 'Adenine', 'DG': 'Guanine', 'DC': 'Cytosine', 'DT': 'Thymine',
         'A': 'Adenine',  'G': 'Guanine',  'C': 'Cytosine',  'T': 'Thymine',
         'U': 'Uracil', 'DOC': 'Cytosine', 'R': 'Purine',    'Y': 'Pyrimidine',
         'K': 'Ketone', '  M': 'Amino',    'S': 'Strong',    'W': 'Weak',
         'B': 'NotA',   'D'  : 'NotC',     'H': 'NotG',      'V': 'NotT',
         'N': 'Any',    'X'  : 'Masked'
        }

    def __init__(self, structure_file):

        if not os.path.exists(structure_file):
            raise IOError('File not found: {0}'.format(structure_file))

        self._file = structure_file
        self._stream = open(structure_file, 'r')

    def __del__(self):
        try:
            self._stream.close()
        except:
            pass

    def models(self):
        """
        Find a list of available model identifiers in the structure.

        @return: a list of model IDs
        @rtype: list
        """
        models = []
        check = {}

        with open(self._file, 'r') as f:
            for line in f:
                if line.startswith('MODEL'):
                    model_id = str(line[10:14]).strip()
                    assert model_id not in check
                    models.append(model_id)
                    check[model_id] = True

        return models

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
        if residue_name in StructureParser._pdb_aminoacids:
            return SequenceTypes.Protein
        elif residue_name in StructureParser._pdb_nucleotides:
            return SequenceTypes.NucleicAcid
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
            if as_type == SequenceTypes.Protein:
                return StructureParser._pdb_aminoacids[residue_name]
            elif as_type == SequenceTypes.NucleicAcid:
                return StructureParser._pdb_nucleotides[residue_name]
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
            if as_type == SequenceTypes.Protein:
                return repr(SequenceAlphabets.Protein.UNK)
            elif as_type == SequenceTypes.NucleicAcid:
                return repr(SequenceAlphabets.Nucleic.Any)
            else:
                return repr(SequenceAlphabets.Unknown.UNK)

    def parse(self, filename, model=None):
        if filename <> self._file:
            self._file = structure_file
            try:
                self._stream.close()
            except:
                pass
            self._stream = open(structure_file, 'r')

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
        self._stream.seek(0)

        header = self._stream.next()
        if not header.startswith('HEADER'):
            raise ValueError('Does not look like a PDB file.')

        structure = csb.bio.structure.Structure(header.split()[-1])
        structure.model_id = None
        if model is not None:
            model = str(model)

        atoms = {}
        sequence = {}
        in_atom = False

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
                    for chain in chains.replace(';', ' ').replace(',', ' ').split():
                        new_chain = csb.bio.structure.Chain(chain, type=SequenceTypes.Unknown, name=chain_name, accession=structure.accession)
                        new_chain.molecule_id = mol_id
                        try:
                            structure.chains.append(new_chain)
                        except csb.pyutils.DuplicateKeyError:
                            raise csb.bio.structure.Broken3DStructureError('Chain {0} is already defined.'.format(new_chain.id))
                        sequence[chain] = []
                        atoms[chain] = []

            elif line.startswith('SEQRES'):

                seq_fields = line[6:].split()
                rank_base = int(seq_fields[0])
                chain_id = seq_fields[1]

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
                    sequence[chain_id].append(str(residue.type))

            elif line.startswith('MODEL'):
                if model and model != line[10:14].strip():
                    while True:
                        try:
                            line = self._stream.next()
                        except StopIteration:
                            raise ValueError('No such model {0} in the structure.'.format(model))

                        if line.startswith('MODEL') \
                               and model == line[10:14].strip():
                            structure.model_id = model
                            break
                else:
                    model = line[10:14].strip()
                    structure.model_id = model

            elif line.startswith('ATOM') \
                     or (in_atom and line.startswith('HETATM')):
                    in_atom = True

                    rank = int(line[22:26])
                    serial_number = int(line[6:11])
                    name = line[12:16]
                    x, y, z = line[30:38], line[38:46], line[46:54]
                    vector = csb.bio.structure.Vector(float(x), float(y), float(z))
                    element = line[76:78].strip()
                    try:
                        element = csb.pyutils.Enum.parsename(ChemElements, element)
                    except csb.pyutils.EnumMemberError as ee:
                        if element in ('D', 'X'):
                            element = ChemElements.x
                        else:
                            raise ee

                    atom = csb.bio.structure.Atom(serial_number, name, element,
                                                  vector)

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
                    residue_name = line[17:20].strip()
                    residue_name = self.parse_residue_safe(residue_name, as_type=structure.chains[atom._chain].type)
                    if structure.chains[atom._chain].type == SequenceTypes.NucleicAcid:
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

            elif line.startswith('ENDMDL'):
                break

            elif line.startswith('END'):
                break

        if structure.model_id != model:
            raise ValueError('No such model {0} in the structure.'.format(model))

        self._map_residues(structure, sequence, atoms)

        return structure

    def _map_residues(self, structure, sequence, atoms):

        assert set(sequence.keys()) == set(atoms.keys())

        for chain in sequence:

            subject = ''.join(sequence[chain])
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
                    if i == 0 or a._sequence_number - lookup[seq_numbers[i-1]][0] not in (0, 1, -1):
                        # if residues [i, i-1] are not consecutive or 'overlapping', initiate a new fragment:
                        fragments.append([a._residue_name.value])
                    else:
                        # else (they are consecutive) append the new residue to the end of the last fragment
                        fragments[-1].append(a._residue_name.value)

            for i, frag in enumerate(fragments):
                fragments[i] = ''.join(frag)
            query = '^.*({0}).*$'.format(').*('.join(fragments))

            matches = re.match(query, subject)
            if matches:
                seq_numitem = -1
                for frag_number, frag in enumerate(matches.groups(), start=1):
                    for i, dummy in enumerate(frag, start=1):
                        seq_numitem += 1
                        # lookup[res_id] is finally set to give the final rank of residue under id res_id:
                        lookup[ seq_numbers[seq_numitem] ] = matches.start(frag_number) + i

                fixed_residue = None
                for atom in atoms[chain]:

                    atom._rank = lookup[atom._residue_id]
                    residue = structure.chains[chain].residues[atom._rank]
                    if residue is not fixed_residue:
                        residue.id = atom._sequence_number, atom._insertion_code
                        fixed_residue = residue
                    residue.structure.append(atom)
                    #assert atom.residue == residue.type

                    del atom._rank
                    del atom._insertion_code
                    del atom._sequence_number
                    del atom._chain
                    del atom._residue_id
                    del atom._residue_name
            else:
                raise csb.bio.structure.Broken3DStructureError('Could not map structure to sequence.')


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
    from tempfile import NamedTemporaryFile

    pdb = NamedTemporaryFile()
    browser = urllib2.urlopen(prefix + accession.lower() + '.ent')
    pdb.write(browser.read())

    return StructureParser(pdb.name).parse_structure(model)


def _parse_async(parser, file, model):
    p = parser(file)
    return p.parse_structure(model)


class AsyncStructureParser(object):
    """
    Wraps StructureParser in an asynchronous call. Since a new process is
    started by Python internally (as opposed to only starting a new thread),
    this makes the parser slower, but provides a way to set a parse timeout
    limit.
    """

    def __init__(self):

        self._worker = None
        self._recycle_worker()

    def _recycle_worker(self):

        if self._worker:
            self._worker.terminate()
        self._worker = multiprocessing.Pool(processes=1)

    def parse_structure(self, structure_file, timeout, model=None,
                        parser=StructureParser):
        """
        Call StructureParser.parse_structure() in a separate process and return
        the output. Raise TimeoutError if the parser does not respond within
        timeout seconds.
        """

        worker = multiprocessing.Pool(processes=1)
        try:
            async_result = worker.apply_async(_parse_async, [parser, structure_file, model])
        except multiprocessing.TimeoutError as ex:
            self._recycle_worker()
            ex.args = (timeout,)
            raise ex

        return async_result.get(timeout=timeout)
