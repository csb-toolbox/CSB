"""
HHpred-related format parsers.
"""

import os
import re

import csb.io
import csb.bio.io
import csb.bio.structure as structure

from csb.pyutils import Enum, CollectionIndexError, ItemNotFoundError
from csb.bio.sequence import Sequence, ProteinAlphabet, A3MAlignment

from csb.bio.hmm import State, Transition, ProfileHMM, HMMLayer, ProfileLength
from csb.bio.hmm import HHpredHitList, HHpredHit, ScoreUnits, States, EVDParameters
from csb.bio.hmm import StateNotFoundError, TransitionNotFoundError


class HHProfileFormatError(ValueError):
    pass

class StructureFormatError(HHProfileFormatError):
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
        Format input_pdb as an HHpred pre-parsed structure file for the current
        HMM. Formatting means: only chain_id left in the PDB file, residue
        sequence numbers strictly corresponding to the HMM states' ranks.

        @param input_pdb: source PDB file name
        @type input_pdb: str
        @param chain_id: source chain ID (which chain to pick from the
                         structure)
        @type chain_id: str
        @param output_file: save the output PDB file here
        @type output_file: str

        @raise StructureFormatError: when the specified PDB chain does not
                                     correspond to the HMM.
        """

        hmm = self.parse()
        s = csb.bio.io.StructureParser(input_pdb).parse_structure()
        chain = s.chains[chain_id]

        if s.first_chain.length != hmm.layers.length:
            raise StructureFormatError(
                        "{0}: Incorrect number of residues".format(chain.entry_id)) 

        for layer, residue in zip(hmm.layers, chain.residues):
            
            if residue.type == ProteinAlphabet.UNK:
                residue.type = layer.residue.type
            
            if residue.type != layer.residue.type:
                msg = "Expected residue of type {0} at position {1}"
                raise StructureFormatError(msg.format(layer.residue.type, layer.rank))

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
                m = tuple(map(int, m))
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
                raise NotImplementedError("Unexpected line: {0}.".format(line[0:5]))

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
        @rtype: L{ProfileHMM}

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
                line = next(lines)
            except StopIteration:
                break

            if line.startswith('NULL'):
                try:
                    backprobs = tuple(map(parse_probability, line.split()[1:]))

                    line = next(lines)
                    residues = line.split()[1:]
                    residues = [Enum.parse(ProteinAlphabet, aa) for aa in residues]

                    for pos, aa in enumerate(residues):
                        background[aa] = backprobs[pos]

                    line = next(lines)
                    tran_types = line.split()

                    line = next(lines)
                    start_probs = tuple(map(parse_probability, line.split()))
                except StopIteration:
                    break

            elif re.match(pattern, line):
                emrow = line
                try:
                    tran_lines.append(next(lines))
                    #junkrow = next(lines)
                except StopIteration:
                    break

                emprobs = emrow.split()
                if len(emprobs) != 23:
                    raise HHProfileFormatError(
                            "Unexpected number of data fields: {0}".format(len(emprobs)))

                rank = int(emprobs[1])
                residue = structure.ProteinResidue(
                    rank=rank, type=emprobs[0], sequence_number=rank, insertion_code=None)
                if residue.type == ProteinAlphabet.GAP:                  
                    raise HHProfileFormatError("Layer {0} can't be represented by a gap".format(rank))

                new_layer = hmm.layers.append(HMMLayer(rank, residue))
                if new_layer != rank:
                    raise HHProfileFormatError('Layer {0} defined as {1}'.format(new_layer, rank))

                match = State(States.Match, emit=Enum.members(ProteinAlphabet))

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

            if start_probs[0] is None:
                raise HHProfileFormatError("Transition Start > Match[1] is undefined")
            
            start_tran = Transition(hmm.start, first_match[States.Match], start_probs[0])
            hmm.start.transitions.append(start_tran)

            if start_probs[1] is not None and start_probs[3] is not None:  # Start -> I[0] -> M[1]
                start_ins = State(States.Insertion, emit=Enum.members(ProteinAlphabet))
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
            fields = tuple(map(parse_probability, ofields))
            
            # 3a. Parse all Neff values and create I[i] and D[i] states if NeffX[i] is not None
            for col, neff in enumerate(tran_types[7:10], start=7):

                if fields[col] is not None:
                    neff_value = float(ofields[col]) / abs(hmm.scale)

                    if neff == 'Neff':
                        hmm.layers[rank].effective_matches = neff_value
                        
                    elif neff == 'Neff_I':
                        hmm.layers[rank].effective_insertions = neff_value
                        
                        if States.Insertion not in hmm.layers[rank]:
                            insertion = State(States.Insertion, emit=Enum.members(ProteinAlphabet))
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
                
                if probability is not None:
                    try:
                        self._add_transition(hmm, rank, tran, probability)
                    except (CollectionIndexError, ItemNotFoundError) as ex:
                        msg = "Can't add transition {0} at {1}: {2.__class__.__name__}, {2!s}"
                        raise HHProfileFormatError(msg.format(tran, rank, ex))                

        return hmm
    
    def _add_transition(self, hmm, rank, tran, probability):

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
                raise NotImplementedError('Unknown transition: {0}'.format(tran))        

    def _parse_structure(self, hmm):
        """
        Add structural information to an existing HMM.

        @param hmm: the hmm object being constructed
        @type hmm: L{ProfileHMM}
        @return: the updated hmm
        @rtype: L{ProfileHMM}

        @raise ValueError: if the structure file does not refer to the HMM
        @raise HHProfileFormatError: when a residue cannot be assigned to an HMM layer
        @raise NotImplementedError: when an unknown data field is encountered
                                    in the PDB file
        """

        assert self._pdb
        
        s = csb.bio.io.StructureParser(self._pdb).parse_structure()
                
        if s.chains.length != 1:
            raise StructureFormatError("Must contain exactly one chain")
        if s.first_chain.length != hmm.layers.length:
            raise StructureFormatError("Incorrect number of residues")                   

        chain = s.first_chain
        chain.compute_torsion()
                    
        for layer, residue in zip(hmm.layers, chain.residues):
            
            if residue.type != layer.residue.type:
                msg = "Expected residue of type {0} at position {1}"
                raise StructureFormatError(msg.format(layer.residue.type, layer.rank))
            
            layer.residue.torsion = residue.torsion.copy()
            
            for atom in residue.items:
                layer.residue.atoms.append(atom.clone())

        return hmm

HHpredProfileParser = HHProfileParser


class HHMFileBuilder(object):
    """
    Builder for HHpred's hhm files.
    
    @param output: destination stream
    @type output: file
    """
    
    def __init__(self, output):

        if not hasattr(output, 'write'):
            raise TypeError(output)    

        self._out = output
        
    @property
    def output(self):
        return self._out    

    def write(self, data):
        self._out.write(data)
        
    def writeline(self, data):
        self.write(data)
        self.write('\n')

    def add_hmm(self, hmm):

        if hmm.score_units != ScoreUnits.LogScales:
            raise ValueError('Scores must be converted to LogScales first.')
                
        self.writeline('''HHsearch {0.version}
NAME  {0.name}
FAM   {0.family}
LENG  {0.length.matches} match states, {0.length.layers} columns in multiple alignment
NEFF  {0.effective_matches}
PCT   {0.pseudocounts}'''.format(hmm))
        if hmm.evd:
            self.writeline('EVD   {0.lamda}  {0.mu}'.format(hmm.evd))            

        self.writeline('SEQ')
        if hmm.dssp:
            self.writeline('>ss_dssp')
            self.writeline(hmm.dssp.to_string())
            if hmm.dssp_solvent:
                self.writeline('>sa_dssp')
                self.writeline(hmm.dssp_solvent)                
        if hmm.psipred:
            self.writeline('>ss_pred')
            self.writeline(hmm.psipred.to_string())
            self.writeline('>ss_conf')
            confidence = [''.join(map(str, m.score)) for m in hmm.psipred]
            self.writeline(''.join(confidence))

        if hmm.alignment:
            if hmm.consensus:
                self.writeline(str(hmm.consensus))
            self.writeline(hmm.alignment.format().rstrip('\r\n'))

        self.writeline('#')

        first_match = hmm.layers[1][States.Match]
        null = [int(first_match.background[aa])
                for aa in sorted(map(str, first_match.background))]
        self.writeline('NULL   {0}'.format('\t'.join(map(str, null))))
        self.writeline('HMM    {0}'.format(
            '\t'.join(sorted(map(str, first_match.emission)))))

        tran_types = 'M->M    M->I    M->D    I->M    I->I    D->M    D->D'.split()
        self.writeline('       {0}'.format(
            '\t'.join(tran_types + 'Neff    Neff_I    Neff_D'.split())))

        self.write("       ")
        for tran_type in tran_types:
            source_statekind = Enum.parse(States, tran_type[0])
            target_statekind = Enum.parse(States, tran_type[3])
            if source_statekind == States.Match:
                try:
                    self.write("{0:<7}\t".format(
                        int(hmm.start.transitions[target_statekind].probability)))
                except TransitionNotFoundError:
                    self.write("*\t")
            else:
                self.write("*\t")
        self.writeline('*\t' * 3)

        for layer in hmm.layers:

            self.write("{0} {1:<5}".format(layer.residue.type, layer.rank))
            for aa in sorted(layer[States.Match].emission):
                emission = layer[States.Match].emission[aa]
                if emission is None:
                    emission = '*'
                else:
                    emission = int(emission)
                self.write("{0:<7}\t".format(emission))
            self.writeline("{0}".format(layer.rank))

            self.write("       ")
            for tran_type in tran_types:
                source_statekind = Enum.parse(States, tran_type[0])
                target_statekind = Enum.parse(States, tran_type[3])

                if target_statekind == States.Match and layer.rank == hmm.layers.last_index:
                    target_statekind = States.End

                try:
                    state = layer[source_statekind]
                    self.write("{0:<7}\t".format(
                        int(state.transitions[target_statekind].probability)))
                except StateNotFoundError:
                    self.write("*\t")

            for data in (layer.effective_matches, layer.effective_insertions,
                         layer.effective_deletions):
                if data is None:
                    data = '*'
                else:
                    data = int(data * abs(hmm.scale))
                self.write("{0:<7}\t".format(data))

            self.writeline("\n")

        self.writeline('//')
        
        
class HHOutputParser(object):
    """
    Parser that is HHsearch's hhr ('HHsearch Result') format aware.
    HHsearch Result

    @param alignments: if set to False, the parser will skip the
                       alignments section of the file
    @type alignments: bool
    """

    def __init__(self, alignments=True):
        
        self._alignments = True
        self.alignments = alignments

    def __repr__(self):
        return "<HHsearch Result Parser>"

    @property
    def alignments(self):
        return self._alignments
    @alignments.setter
    def alignments(self, value):
        self._alignments = bool(value)

    def parse_file(self, hhr_file, header_only=False):
        """
        Parse all hits from this HHpred result file.

        @param hhr_file: input HHR file name
        @type hhr_file: str

        @return: parsed hits
        @rtype: HHpredHitList

        @raise HHOutputFormatError: if the hhr file is corrupt
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
        stream = csb.io.MemoryStream()
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
