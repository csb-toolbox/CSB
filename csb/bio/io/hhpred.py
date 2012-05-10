"""
HHpred-related format parsers.
"""

import os
import re
import numpy

import csb.io
import csb.bio.io
import csb.bio.structure as structure

from csb.pyutils import Enum, CollectionIndexError, ItemNotFoundError, EnumMemberError

from csb.bio.sequence import Sequence, SequenceAlphabets, A3MAlignment

from csb.bio.hmm import State, Transition, ProfileHMM, HMMLayer, ProfileLength
from csb.bio.hmm import HHpredHitList, HHpredHit, ScoreUnits, States, EVDParameters


class HHProfileFormatError(ValueError):
    pass


class HHOutputFormatError(HHProfileFormatError):
    pass


class HHProfileParser(object):
    """
    A class that is HHpred HMM format aware.
    Produces a L{ProfileHMM} object representation of a given HHpred profile
    HMM.

    @param hhm_file: *.hhm file name to parse
    @type hhm_file: str
    @param pdb_file: an optional hhm structure file, containing the structural
                     data, associated with the HMM.
    @type pdb_file: str
    @note: This is NOT a regular PDB file! It must be specifically formatted
           for use with HHpred.

    @raise IOError: if any of the files does not exist
    """

    def __init__(self, hhm_file, pdb_file=None):

        if not os.path.exists(hhm_file):
            raise IOError("Could not read HHM file {0}".format(hhm_file))

        if pdb_file:
            if not os.path.exists(pdb_file):
                raise IOError("Could not read structure file {0}".format(pdb_file))

        self._file = hhm_file
        self._pdb = pdb_file
        self._properties = None
        self._sequences = None
        self._profile = None

        self._chopped = False
        self._chop()

    def format_structure(self, input_pdb, chain_id, output_file):
        """
        Format input_pdb as a HHpred pre-parsed structure file for the current
        HMM. Formatting means: only chain_id left in the PDB file, residue
        sequence numbers strictly corresponding to the HMM states' ranks.

        @param input_pdb: source PDB file name
        @type input_pdb: str
        @param chain_id: source chain ID (which chain to pick from the
                         structure)
        @type chain_id: str
        @param output_file: save the output PDB file here
        @type output_file: str

        @raise ValuError: if the specified source chain has a different length
                          or residues, compared to the profile HMM; or when a
                          residue's rank is out of range again with reference
                          to the HMM.
        """

        hmm = self.parse()
        s = csb.bio.io.StructureParser(input_pdb).parse_structure()
        chain = s.chains[chain_id]

        if hmm.layers.length != chain.length:
            raise ValueError(
                'Chain {0.id} ({0.length} residues) does not have the expected length ({1})'.format(
                    chain, hmm.layers.length))

        for residue in chain.residues:
            try:
                if hmm.layers[residue.rank].residue.type != residue.type:
                    if residue.type == 'X':
                        residue.type = hmm.layers[residue.rank].residue.type
                    else:
                        raise ValueError(
                            'Residue {0} differs from the residue under the layer with the same rank: {1}'.format(
                                residue, hmm.layers[residue.rank].residue.type))
            except ItemNotFoundError:
                raise ValueError(
                    "Residue {0} is out of the HMM's range {1}..{2}".format(
                        residue, hmm.layers.start_index, hmm.layers.last_index))
            residue.id = str(residue.rank), '$'

        for residue in chain.residues:
            residue.id = str(residue.rank), None

        new_structure = structure.Structure(hmm.id)
        new_structure.chains.append(chain.clone())

        new_structure.to_pdb(output_file)

    def parse(self, units=ScoreUnits.LogScales):
        """
        Parse the HMM profile.

        @param units: also convert the profile score units to the specified
                      L{csb.bio.hmm.ScoreUnits} kind
        @type units: L{csb.pyutils.EnumItem}

        @return: a L{ProfileHMM} instance
        @rtype: L{ProfileHMM}

        @raise HHProfileFormatError: when the hhm file is invalid/corrupt
        """
        assert self._chopped

        hmm = ProfileHMM(units=units, scale=-1000.0, logbase=2)
        if self._profile:
            self._parse_profile(hmm, units)
        else:
            raise HHProfileFormatError('Missing HMM profile table.')

        if self._properties:
            self._parse_properties(hmm)
        else:
            raise HHProfileFormatError('No profile properties found.')

        if self._sequences:
            self._parse_sequences(hmm)
        else:
            raise HHProfileFormatError('No profile MSA and secondary structure found.')

        if hmm.dssp:
            hmm._assign_secstructure()

        if self._pdb:
            self._parse_structure(hmm)

        return hmm

    def _chop(self):
        """
        Chop the HHM file into pieces - HMM properties, secondary structure +
        MSA, HMM.
        """

        try:
            with open(self._file, 'r') as f:
                pr = csb.io.EntryReader(f, 'HHsearch', 'SEQ')
                self._properties = pr.readall()[0].splitlines()

            with open(self._file, 'r') as f:
                sr = csb.io.EntryReader(f, '>', '#')
                self._sequences = sr.readall()

            with open(self._file, 'r') as f:
                hr = csb.io.EntryReader(f, '#', '//')
                self._profile = hr.readall()[0].splitlines()

            self._chopped = True
            
        except IndexError:
            raise HHProfileFormatError('Corrupt HHM file.')

    def _issue(self, hmm, issue):
        """
        Log any known parse error encountered

        @param hmm: the faulty profile HMM
        @type hmm: L{ProfileHMM}
        @param issue: description of the error
        @type issue: str
        """
        hmm._issues.append(issue)

    def _parse_properties(self, hmm):
        """
        Parse the profile properties

        @param hmm: the hmm object being constructed
        @type hmm: L{ProfileHMM}
        @return: the updated hmm
        @rtype: hmm

        @raise NotImplementedError: if an unexpected data field is encountered
        """
        assert self._chopped

        for line in self._properties:

            if line.startswith('NAME '):
                hmm.name = line[6:].strip()

            elif line.startswith('FAM  '):
                hmm.family = line[6:].strip()

            elif line.startswith('LENG '):
                m = re.search('(\d+)\D+(\d+)', line).groups()
                m = map(int, m)
                hmm.length = ProfileLength(m[0], m[1])

            elif line.startswith('FILE '):
                hmm.id = line[6:].strip()

            elif line.startswith('HHsearch'):
                hmm.version = float(line[9:])

            elif line.startswith('NEFF '):
                hmm.effective_matches = float(line[6:])

            elif line.startswith('PCT  '):
                if line[6:].strip().lower() == 'true':
                    hmm.pseudocounts = True

            elif line.startswith('EVD  '):
                lamda, mu = map(float, line[6:].split())
                hmm.evd = EVDParameters(lamda, mu)
                    
            elif line[:5] in ('COM  ', 'DATE ', 'FILT '):
                pass
            else:
                raise NotImplementedError(
                    'Unexpected line. Don\'t know how to parse field {0}.'.format(line[0:5]))

        if not hmm.id and hmm.name:
            hmm.id = hmm.name.split()[0]
        return hmm

    def _parse_sequences(self, hmm):
        """
        Parse secondary structure and multiple alignment sections.

        @param hmm: the hmm object being constructed
        @type hmm: L{ProfileHMM}
        @return: the updated hmm
        @rtype: hmm
        """
        assert self._chopped

        psipred = None
        msa_entries = []

        for entry in self._sequences:

            header_token = entry[:8]
            if header_token in ['>ss_dssp', '>sa_dssp', '>ss_pred', '>ss_conf', '>Consens']:

                lines = entry.strip().splitlines()
                seq = re.sub('\s+', '', ''.join(lines[1:]))

                if header_token == '>ss_dssp':
                    hmm.dssp = structure.SecondaryStructure(seq)
                elif header_token == '>sa_dssp':
                    hmm.dssp_solvent = seq
                elif header_token == '>ss_pred':
                    psipred = seq
                elif header_token == '>ss_conf':
                    conf = seq
                    hmm.psipred = structure.SecondaryStructure(psipred, conf)
                elif header_token == '>Consens':
                    hmm.consensus = Sequence('Consensus', 'Consensus', seq)                    
            else:
                msa_entries.append(entry)

        if msa_entries:
            msa = '\n'.join(msa_entries)
            hmm.alignment = A3MAlignment.parse(msa, strict=False)

        return hmm

    def _parse_profile(self, hmm, units=ScoreUnits.LogScales):
        """
        Parse the HMM profile.

        @param hmm: the hmm object being constructed
        @type hmm: L{ProfileHMM}
        @return: the updated hmm
        @rtype: hmm

        @raise NotImplementedError: when an unknown transition string is
                                    encountered
        """
        assert self._chopped

        # 0. Prepare start and end states
        hmm.start = State(States.Start)
        hmm.end = State(States.End)

        residues = None
        background = {}
        tran_types = None
        tran_lines = []
        start_probs = None

        lines = iter(self._profile)
        pattern = re.compile('^[A-Z\-]\s[0-9]+\s+')

        if units == ScoreUnits.LogScales:

            def parse_probability(v):
                if v.strip() == '*':
                    return None
                else:
                    return float(v)
        else:

            def parse_probability(v):
                if v.strip() == '*':
                    return None
                else:
                    return hmm._convert(units, float(v),
                                        hmm.scale, hmm.logbase)

        # 1. Create all layers (profile columns), create and attach their match states
        
        while True:
            try:
                line = lines.next()
            except StopIteration:
                break

            if line.startswith('NULL'):
                try:
                    backprobs = map(parse_probability, line.split()[1:])

                    line = lines.next()
                    residues = line.split()[1:]
                    residues = [Enum.parse(SequenceAlphabets.Protein, aa) for aa in residues]

                    for pos, aa in enumerate(residues):
                        background[aa] = backprobs[pos]

                    line = lines.next()
                    tran_types = line.split()

                    line = lines.next()
                    start_probs = map(parse_probability, line.split())
                except StopIteration:
                    break

            elif re.match(pattern, line):
                emrow = line
                try:
                    tran_lines.append(lines.next())
                    #junkrow = lines.next()
                except StopIteration:
                    break

                emprobs = emrow.split()
                assert len(emprobs) == 23

                rank = int(emprobs[1])
                residue = structure.ProteinResidue(
                    rank=rank, type=emprobs[0], sequence_number=rank, insertion_code=None)
                if residue.type == SequenceAlphabets.Protein.GAP:                  
                    raise HHProfileFormatError("Layer {0} can't be represented by a gap".format(rank))

                new_layer = hmm.layers.append(HMMLayer(rank, residue))
                assert new_layer == rank, '{0} vs {1}'.format(rank, new_layer)

                match = State(States.Match, emit=Enum.members(SequenceAlphabets.Protein))

                match.rank = rank
                match.background.set(background)

                for col, aa in enumerate(residues):
                    prob = parse_probability(emprobs[col + 2])
                    match.emission.append(aa, prob)

                hmm.layers[new_layer].append(match)
                assert hmm.layers.last_index == match.rank

        # 2. Append starting transitions: S -> M[1] and optionally S -> D[1] and S -> I[0].
        #    States D[1] and I[0] will be created if needed
        #    Note that [0] is not a real layer, I[0] is simply an insertion at the level of Start
        if len(hmm.layers) > 0:

            first_match = hmm.layers[hmm.layers.start_index]

            assert start_probs[0] is not None          # Start -> M[1]
            start_tran = Transition(hmm.start, first_match[States.Match], start_probs[0])
            hmm.start.transitions.append(start_tran)

            if start_probs[1] is not None and start_probs[3] is not None:  # Start -> I[0] -> M[1]
                start_ins = State(States.Insertion, emit=Enum.members(SequenceAlphabets.Protein))
                start_ins.rank = 0
                start_ins.background.set(background)
                start_ins.emission = start_ins.background

                hmm.start_insertion = start_ins
                # Start -> I[0]
                hmm.start.transitions.append(
                        Transition(hmm.start, hmm.start_insertion, start_probs[1]))
                # I[0] -> M[1]
                hmm.start_insertion.transitions.append(
                        Transition(hmm.start_insertion, first_match[States.Match], start_probs[3]))
                # I[0] -> I[0]
                if start_probs[4]:
                    hmm.start_insertion.transitions.append(
                        Transition(hmm.start_insertion, hmm.start_insertion, start_probs[4]))


            if start_probs[2] is None and start_probs[6] is not None:
                # M->D is corrupt (*) at the Start layer, using D->D instead
                start_probs[2] = start_probs[6]

            if start_probs[2] is not None:  # Start -> D[1]
                start_del = State(States.Deletion)
                start_del.rank = 1
                hmm.layers[1].append(start_del)
                start_tran = Transition(hmm.start, first_match[States.Deletion], start_probs[2])
                hmm.start.transitions.append(start_tran)
        else:
            start_tran = Transition(hmm.start, hmm.end, start_probs[0])
            hmm.start.transitions.append(start_tran)


        # 3. Append remaining transitions. I and D states will be created on demand.

        for rank, fields in enumerate(tran_lines, start=hmm.layers.start_index):
            assert hmm.layers[rank][States.Match].rank == rank

            ofields = fields.split()
            fields = map(parse_probability, ofields)
            
            # 3a. Parse all Neff values and create I[i] and D[i] states if NeffX[i] is not None
            for col, neff in enumerate(tran_types[7:10], start=7):

                if fields[col] is not None:
                    neff_value = float(ofields[col]) / abs(hmm.scale)

                    if neff == 'Neff':
                        hmm.layers[rank].effective_matches = neff_value
                        
                    elif neff == 'Neff_I':
                        hmm.layers[rank].effective_insertions = neff_value
                        
                        if States.Insertion not in hmm.layers[rank]:
                            insertion = State(States.Insertion, emit=Enum.members(SequenceAlphabets.Protein))
                            insertion.background.set(background)
                            insertion.emission.set(background)
                            insertion.rank = rank                        
                            hmm.layers[rank].append(insertion)
                        
                    elif neff == 'Neff_D':
                        hmm.layers[rank].effective_deletions = neff_value
    
                        if States.Deletion not in hmm.layers[rank] and neff_value > 0:
                            deletion = State(States.Deletion)
                            deletion.rank = rank
                            hmm.layers[rank].append(deletion)                         
            
            # 3b. Starting from the first layer, parse all transitions and build the HMM graph stepwise
            for col, tran in enumerate(tran_types):

                probability = fields[col]
                if probability is None:
                    continue

                match = hmm.layers[rank][States.Match]
                nextmatch = None
                if rank < hmm.layers.last_index:
                    nextmatch = hmm.layers[rank + 1][States.Match]
                else:
                    nextmatch = hmm.end

                if tran == 'M->M':
                    transition = Transition(match, nextmatch, probability)
                    match.transitions.append(transition)
                                        
                elif tran == 'M->I':
                    insertion = hmm.layers[rank][States.Insertion]
                    transition = Transition(match, insertion, probability)
                    match.transitions.append(transition)
                    
                elif tran == 'M->D':
                    deletion = State(States.Deletion)
                    deletion.rank = rank + 1
                    hmm.layers[rank + 1].append(deletion)
                    transition = Transition(match, deletion, probability)
                    match.transitions.append(transition)
                                        
                elif tran == 'I->M':
                    insertion = hmm.layers[rank][States.Insertion]
                    transition = Transition(insertion, nextmatch, probability)
                    insertion.transitions.append(transition)
                        
                elif tran == 'I->I':
                    insertion = hmm.layers[rank][States.Insertion]
                    selfloop = Transition(insertion, insertion, probability)
                    insertion.transitions.append(selfloop)
                        
                elif tran == 'D->M':
                    deletion = hmm.layers[rank][States.Deletion]
                    transition = Transition(deletion, nextmatch, probability)
                    deletion.transitions.append(transition)
                    
                elif tran == 'D->D':
                    deletion = hmm.layers[rank][States.Deletion]

                    if States.Deletion not in hmm.layers[rank + 1]:
                        nextdeletion = State(States.Deletion)
                        nextdeletion.rank = rank + 1
                        hmm.layers[rank + 1].append(nextdeletion)
                        
                    else:
                        nextdeletion = hmm.layers[rank + 1][States.Deletion]
                        assert match.transitions[States.Deletion].successor == nextdeletion
                        
                    transition = Transition(deletion, nextdeletion, probability)
                    deletion.transitions.append(transition)

                else:
                    if not tran.startswith('Neff'):
                        raise NotImplementedError('Unknown transition "{0}"'.format(tran))

        return hmm

    def _parse_structure(self, hmm):
        """
        Parse and add structural information to an existing HMM.

        @param hmm: the hmm object being constructed
        @type hmm: L{ProfileHMM}
        @return: the updated hmm
        @rtype: hmm

        @raise ValueError: if the structure file does not refer to the HMM
        @raise HHProfileFormatError: when any residue cannot be assigned
                                         to layer from the HMM
        @raise NotImplementedError: when an unknown data field is encountered
                                    in the PDB file
        """

        assert self._pdb
        assert hmm is not None

        with open(self._pdb, 'r') as pdb:

            offset_fix = 0
            for line in pdb:

                if line.startswith('HEADER'):
                    m = re.search('SCOP domain (\w+)\s+', line)
                    if m:
                        fam = m.groups()[0]
                    else:
                        fam = line.split()[-1]
                    if fam.strip().lower() != hmm.id.lower():
                        raise ValueError('This coordinate file refers to domain {0}, while the HMM refers to {1}.'.format(fam, hmm.id))

                elif line.startswith('ATOM'):

                    rank = int(line[22:26])
                    try:
                        while repr(hmm.layers[rank + offset_fix].residue.type) != line[17:20].strip():
                            issue = 'ProteinResidue {0} [{1}] in structure does not match residue {2!r} at layer {3}.'.format(
                                line[17:20],
                                rank,
                                hmm.layers[rank].residue.type,
                                rank + offset_fix)
                            self._issue(hmm, issue)
                            offset_fix += 1
                            
                    except CollectionIndexError:
                        raise HHProfileFormatError(
                            'ProteinResidue {0} at fixed position {1} (former {2}) is out of range.'.format(
                                line[17:20], rank + offset_fix, rank))

                    serial_number = int(line[6:11])
                    name = line[12:16]
                    element = line[76:78].strip()
                    try:
                        element = Enum.parsename(structure.ChemElements, element)
                        
                    except EnumMemberError as ee:
                        if element in ('D', 'X'):
                            element = structure.ChemElements.x              
                        else:
                            raise ee
                        
                    x, y, z = line[30:38], line[38:46], line[46:54]
                    vector = numpy.array([float(x), float(y), float(z)])

                    atom = structure.Atom(serial_number, name, element, vector)
                    atom.rank = rank + offset_fix

                    atom.alternate = line[16].strip()
                    if not atom.alternate:
                        atom.alternate = None
                    atom.occupancy = float(line[54:60])
                    atom.temperature = float(line[60:66])
                    atom.charge = line[78:80].strip()
                    if atom.charge:
                        atom.charge = int(atom.charge)

                    assert repr(hmm.layers[atom.rank].residue.type) == line[17:20].strip(), atom.rank
                    #atom.residue = hmm.layers[atom.rank].residue
                    hmm.layers[atom.rank].residue.atoms.append(atom)

                elif line[:6].strip() in ('AUTHOR', 'REMARK', 'COMPND', 'SEQRES', 'TER', 'END'):
                    pass

                else:
                    raise NotImplementedError(line[:6])

        # compute torsion angles
        for layer in hmm.layers:

            prev_residue, next_residue = None, None

            if layer.rank > 1:
                prev_residue = hmm.layers[layer.rank - 1].residue
            if layer.rank < hmm.layers.last_index:
                next_residue = hmm.layers[layer.rank + 1].residue

            layer.residue.torsion = layer.residue.compute_torsion(
                prev_residue, next_residue, strict=False)

        return hmm


HHpredProfileParser = HHProfileParser


class HHOutputParser(object):
    """
    A parser that is HHsearch's hhr ('HHsearch Result') format aware.

    @param alignments: if set to False, the parser will not read the
                       alignments section in the file at all
    @type alignments: bool
    """

    def __init__(self, alignments=True):

        self.alignments = alignments

    def __repr__(self):
        return "<HHsearch Result Parser>"

    __str__ = __repr__

    def parse_file(self, hhr_file, header_only=False):
        """
        Parse all hits from this HHpred result file.

        @param hhr_file: input HHR file name
        @type hhr_file: str

        @return: parsed hits
        @rtype: HHpredHitList

        @raise HHOutputFormatError: if the output is corrupt
        """

        with open(os.path.expanduser(hhr_file)) as stream:
            return self._parse(stream, header_only)

    def parse_string(self, output, header_only=False):
        """
        Get all hits from an C{output} string.

        @param output: HHpred standard output
        @type output: str

        @return: parsed hits
        @rtype: HHpredHitList

        @raise HHOutputFormatError: if the output is corrupt
        """
        import cStringIO

        stream = cStringIO.StringIO()
        stream.write(output)
        stream.seek(0)

        return self._parse(stream, header_only)

    def _parse(self, stream, header_only):

        qlen = None
        in_hits = False
        in_alis = False
        has_alis = False
        c_rank = 0
        header = {}
        hits = {}
        alis = {}

        for line in stream:

            if not in_hits and not in_alis:

                if line.replace(' ', '').startswith('NoHitProbE-value'):
                    in_hits = True
                    continue
                elif line.strip() == '':
                    continue
                else:  # parse header data (stuff above the hits table)
                    columns = line.strip().split(None, 1)
                    if len(columns) == 2:

                        identifier, data = columns
                        if identifier in ('Query', 'Command'):
                            data = data.strip()
                        elif identifier == 'Neff':
                            data = float(data)
                        elif identifier in ('Searched_HMMs', 'Match_columns'):
                            data = int(data)

                        header[identifier] = data

                        if identifier == 'Match_columns':
                            qlen = data

            if in_hits and not header_only:
                if not line.strip():  # suboptimal way to handle block switch
                    in_hits = False
                    in_alis = True
                    if self.alignments:
                        continue
                    else:
                        break
                elif line.strip() == 'Done':
                    in_hits = False
                    in_alis = False
                    break

                description = line[:34].split()
                rank = int(description[0]) 
                id = description[1]

                pos = line[85:94].strip()
                start, end = map(int, pos.split('-'))

                qpos = line[75:84].strip()
                qstart, qend = map(int, qpos.split('-'))

                probability = float(line[35:40]) / 100.0

                hit = HHpredHit(rank, id, start, end, qstart, qend, probability, qlen)

                hit.evalue = float(line[41:48])
                hit.pvalue = float(line[49:56])
                hit.score = float(line[57:63])
                hit.ss_score = float(line[64:69])

                hit.slength = int(line[94:].replace('(', '').replace(')', ''))

                hits[hit.rank] = hit
                alis[hit.rank] = {'q': [], 's': []}

            elif in_alis and not header_only:
                if line.startswith('Done'):
                    in_alis = False
                    break
                
                elif line.startswith('No '):
                    c_rank = int(line[3:])
                    if c_rank not in hits:
                        raise HHOutputFormatError('Alignment {0}. refers to a non-existing hit'.format(c_rank))
                    
                elif line.startswith('>'):
                    hits[c_rank].name = line[1:].strip()
                    
                elif line.startswith('Probab='):
                    for pair in line.split():
                        key, value = pair.split('=')
                        if key == 'Identities':
                            hits[c_rank].identity = float(
                                value.replace('%', ''))
                        elif key == 'Similarity':
                            hits[c_rank].similarity = float(value)
                        elif key == 'Sum_probs':
                            hits[c_rank].prob_sum = float(value)
                            
                elif line.startswith('Q ') and not line[:11].rstrip() in ('Q Consensus', 'Q ss_pred','Q ss_conf', 'Q ss_dssp'):
                    for residue in line[22:]:
                        if residue.isspace() or residue.isdigit():
                            break
                        else:
                            alis[c_rank]['q'].append(residue)
                            has_alis = True
                            
                elif line.startswith('T ') and not line[:11].rstrip() in ('T Consensus', 'T ss_pred','T ss_conf', 'T ss_dssp'):
                    for residue in line[22:]:
                        if residue.isspace() or residue.isdigit():
                            break
                        else:
                            alis[c_rank]['s'].append(residue)

        if self.alignments and has_alis:
            for rank in alis:
                try:
                    hits[rank].add_alignment(alis[rank]['q'], alis[rank]['s'])
                    
                except (KeyError, ValueError) as er:
                    raise HHOutputFormatError('Corrupt alignment at hit No {0}.\n {1}'.format(rank, er))

        del alis

        hits = HHpredHitList(hits.values())

        hits.sort()

        ## add data obtained from the header to the HHpredHitList
        for identifier, data in header.items():
            if identifier == 'Query':
                hits.query_name = data
            elif identifier == 'Match_columns':
                hits.match_columns = data
            elif identifier == 'No_of_seqs':
                hits.no_of_seqs = data
            elif identifier == 'Neff':
                hits.neff = data
            elif identifier == 'Searched_HMMs':
                hits.searched_hmms = data
            elif identifier == 'Date':
                hits.date = data
            elif identifier == 'Command':
                hits.command = data

        return hits


HHpredOutputParser = HHOutputParser
