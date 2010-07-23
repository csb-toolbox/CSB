import os
import sys
import math
import csb.pyutils
import csb.bio.structure as structure
import csb.bio.fragments.isites as isites
import csb.bio.sequence as sequence

from itertools import izip
from collections import namedtuple


class UnobservableStateError(AttributeError):
    pass


class StateExistsError(KeyError):
    pass


class TransitionExistsError(KeyError):
    pass


class EmissionExistsError(KeyError):
    pass

States = csb.pyutils.enum(Match='M', Insertion='I', Deletion='D',
                          Start='S', End='E')
"""
@var States: HMM state types
"""

ScoreUnits = csb.pyutils.enum(LogScales='LogScales', Probability='Probability')
"""
@var ScoreUnits: HMM emission and transition score units
"""


class ProfileHMM(object):
    """
    Describes a protein profile Hidden Markov Model.
    Optional parameters:

    @param units: defines the units of the transition and emission scores
    @type units: L{ScoreUnits}
    @param scale: the scaling factor used to convert emission/transition
                  probabilities
    @type scale: float
    @param logbase: the base of the logarithm used for scaling the emission and
                    transition probabilities
    @type logbase: float
    """

    def __init__(self, units=ScoreUnits.LogScales, scale=-1000., logbase=2):

        self.name = None
        self.id = None
        self.family = None
        self.length = ProfileLength(None, None)
        self.alignment = None
        self.dssp = None
        self.psipred = None
        self.effective_matches = None
        self.version = None
        self.pseudocounts = False
        self.layers = HMMLayersCollection()
        self.start = State(States.Start)
        self.end = State(States.End)
        if units is None:
            self._score_units = ScoreUnits.LogScales
        else:
            self._score_units = units

        self.scale = scale
        self.logbase = logbase

        self._issues = csb.pyutils.CollectionContainer()

    @property
    def residues(self):
        res = [layer.residue for layer in self.layers]
        return csb.pyutils.CollectionContainer(res, type=structure.Residue,
                                               start_index=1)

    @property
    def has_structure(self):
        has = False
        for layer in self.layers:
            if layer.residue.has_structure:
                return True
        return has

    @property
    def score_units(self):
        return self._score_units

    @property
    def known_issues(self):
        return list(self._issues)

    def serialize(self, file_name):
        """
        Serialize this HMM to a file.

        @param file_name: target file name
        @type file_name: str
        """
        import cPickle
        sys.setrecursionlimit(10000)
        cPickle.dump(self, open(file_name, 'w'))
        sys.setrecursionlimit(1000)

    @staticmethod
    def deserialize(file_name):
        """
        De-serialize an HMM from a file.

        @param file_name: source file name
        @type file_name: str
        """
        import cPickle
        sys.setrecursionlimit(10000)
        return cPickle.load(open(file_name, 'r'))
        sys.setrecursionlimit(1000)

    def _convert(self, units, score, scale, logbase):
        if units == ScoreUnits.Probability:
            return logbase ** (score / scale)
        elif units == ScoreUnits.LogScales:
            if score == 0:
                #score = sys.float_info.min
                return None
            return math.log(score, logbase) * scale
        else:
            raise ValueError('Unknown target unit {0}'.format(units))

    def to_hmm(self, output_file=None, convert_scores=False):
        """
        Dump the profile in HHM format.

        @param output_file: the output file name
        @type output_file: str
        @param convert_scores: if True, forces automatic convertion to
                              L{ScoreUnits}.LogScales, which is required
                              by the output file format
        @type convert_scores: bool
        """
        import StringIO

        class MyStringIO(StringIO.StringIO):

            def writeline(self, data):
                self.write(data)
                self.write(os.linesep)

        stream = MyStringIO()
        hmm = self

        if convert_scores:
            hmm.convert_scores(ScoreUnits.LogScales)
        elif hmm.score_units != ScoreUnits.LogScales:
            raise ValueError('Scores must be converted to LogScales first.')

        stream.writeline('''HHsearch {0.version}
NAME  {0.name}
FAM   {0.family}
LENG  {0.length.matches} match states, {0.length.layers} columns in multiple alignment
NEFF  {0.effective_matches}
PCT   {0.pseudocounts}'''.format(hmm))

        stream.writeline('SEQ')
        if hmm.dssp:
            stream.writeline('>ss_dssp')
            stream.writeline(hmm.dssp.tostring())
        if hmm.psipred:
            stream.writeline('>ss_pred')
            stream.writeline(hmm.psipred.tostring())
            stream.writeline('>ss_conf')
            confidence = [''.join(map(str, m.score)) for m in hmm.psipred]
            stream.writeline(''.join(confidence))

        if hmm.alignment:
            if hmm.alignment.consensus:
                stream.writeline(hmm.alignment.consensus.to_fasta())
            stream.writeline(hmm.alignment.to_string())

        stream.writeline('#')

        first_match = hmm.layers[1][States.Match]
        null = [int(first_match.background[aa])
                for aa in sorted(map(str, first_match.background))]
        stream.writeline('NULL   {0}'.format('\t'.join(map(str, null))))
        stream.writeline('HMM    {0}'.format(
            '\t'.join(sorted(map(str, first_match.emission)))))

        tran_types = 'M->M    M->I    M->D    I->M    I->I    D->M    D->D'.split()
        stream.writeline('       {0}'.format(
            '\t'.join(tran_types + 'Neff    Neff_I    Neff_D'.split())))

        stream.write("       ")
        for tran_type in tran_types:
            source_statekind = csb.pyutils.Enum.parse(States, tran_type[0])
            target_statekind = csb.pyutils.Enum.parse(States, tran_type[3])
            if source_statekind == States.Match:
                try:
                    stream.write("{0:<7}\t".format(
                        int(hmm.start.transitions[target_statekind].probability)))
                except csb.pyutils.ItemNotFoundError:
                    stream.write("*\t")
            else:
                stream.write("*\t")
        stream.writeline('*\t' * 3)

        for layer in hmm.layers:

            stream.write("{0} {1:<5}".format(layer.residue.type, layer.rank))
            for aa in sorted(layer[States.Match].emission):
                emission = layer[States.Match].emission[aa]
                if emission is None:
                    emission = '*'
                else:
                    emission = int(emission)
                stream.write("{0:<7}\t".format(emission))
            stream.writeline("{0}".format(layer.rank))

            stream.write("       ")
            for tran_type in tran_types:
                source_statekind = csb.pyutils.Enum.parse(States, tran_type[0])
                target_statekind = csb.pyutils.Enum.parse(States, tran_type[3])

                if target_statekind == States.Match and \
                       layer.rank == hmm.layers.last_index:
                    target_statekind = States.End

                try:
                    state = layer[source_statekind]
                    stream.write("{0:<7}\t".format(
                        int(state.transitions[target_statekind].probability)))
                except csb.pyutils.ItemNotFoundError:
                    stream.write("*\t")

            for data in (layer.effective_matches, layer.effective_insertions,
                         layer.effective_deletions):
                if data is None:
                    data = '*'
                else:
                    data = int(data)
                stream.write("{0:<7}\t".format(data))

            stream.writeline("\n")

        stream.writeline('//')

        data = stream.getvalue()
        stream.close()

        if not output_file:
            return data
        else:
            with open(output_file, 'w') as out:
                out.write(data)

    def segment(self, start, end):
        """
        Extract a sub-segment of the profile.

        @param start: start layer of the segment (rank)
        @type start: int
        @param end: end layer of the segment (rank)
        @type end: int

        @return: a deepcopy of the extracted HMM segment
        @rtype: L{ProfileHMMSegment}
        """
        return ProfileHMMSegment(self, start, end)

    def structure(self, chain_id=None, accession=None):
        """
        Extract the structural information from the HMM.

        @param accession: defines the accession number of the structure
        @type accession: str
        @param chain_id: defines explicitly the chain identifier
        @type chain_id: str

        @return: a shallow L{Structure} wrapper around the residues in the HMM.
        @rtype: L{structure.Structure}
        """
        struct = structure.Structure(accession or self.id)
        chain = self.chain(chain_id)
        struct.chains.append(chain)

        return struct

    def chain(self, chain_id=None):
        """
        Extract the structural information from the HMM.

        @param chain_id: defines explicitly the chain identifier
        @type chain_id: str

        @return: a shallow L{Chain} wrapper around the residues in the HMM.
        @rtype: L{structure.Chain}
        """
        if chain_id is None:
            if self.id:
                chain_id = self.id.rstrip()[-1]
            else:
                chain_id = '_'
        chain = structure.Chain(chain_id, type=sequence.SequenceTypes.Protein, residues=self.residues)

        return chain

    def emission_profile(self):
        """
        Extract the emission scores of all match states in the profile.
        The metric of the emission scores returned depends on the current
        hmm.score_units setting - you may need to call hmm.convert_scores()
        to adjust the hmm to your particular needs.

        @return: a list of dictionaries; each dict key is a single amino acid
        @rtype: list
        """
        profile = []

        for layer in self.layers:
            emission = {}

            for aa in layer[States.Match].emission:
                emission[str(aa)] = layer[States.Match].emission[aa] or 0.0
            profile.append(emission)

        return profile

    def convert_scores(self, units=ScoreUnits.Probability, method=None):
        """
        Convert emission and transition scores to the specified units.

        @param units: the target units for the conversion (a member of
                      L{ScoreUnits}).
        @type units: L{csb.pyutils.EnumItem}
        @param method: if defined, implements the exact mathematical
                       transformation that will be applied. It must be a
                       function or lambda expression with the following
                       signature::
                             def (target_units, score, scale, logbase)

                       and it has to return the score converted to
                       C{target_units}. If method performs a conversion from
                       probabilities to scaled logs, you should also update
                       C{hmm.scale} and C{hmm.logbase}.
        @type method: function, lambda
        """

        if self._score_units == units:
            return

        if method is not None:
            convert = method
        else:
            convert = self._convert

        for layer in self.layers:
            for state_kind in layer:
                state = layer[state_kind]
                if not state.silent:
                    for residue in state.emission:
                        if state.emission[residue] is not None:
                            state.emission.update(residue, convert(
                            units, state.emission[residue],
                            self.scale, self.logbase))
                    for residue in state.background:
                        if state.background[residue] is not None:
                            state.background.update(residue, convert(
                            units, state.background[residue],
                            self.scale, self.logbase))
                for tran_kind in state.transitions:
                    transition = state.transitions[tran_kind]
                    transition.probability = convert(units, transition.probability, self.scale, self.logbase)
            if layer.effective_matches is not None:
                layer.effective_matches = convert(units, layer.effective_matches, self.scale, self.logbase)
            if layer.effective_insertions is not None:
                layer.effective_insertions = convert(units, layer.effective_insertions, self.scale, self.logbase)
            if layer.effective_deletions is not None:
                layer.effective_deletions = convert(units, layer.effective_deletions, self.scale, self.logbase)
        for tran_kind in self.start.transitions:
            transition = self.start.transitions[tran_kind]
            transition.probability = convert(units, transition.probability, self.scale, self.logbase)

        self._score_units = units

    def to_isite(self, cluster_id, start=None, end=None, rep_accession=None,
                 rep_chain=None, alpha=0.2):
        """
        Convert and return this HMM as it was an I-Sites Library fragment.

        @param cluster_id: ID of the I-Site
        @type cluster_id: int
        @param start: start position (optional)
        @type start: int
        @param end: end position (optional)
        @type end: int
        @param rep_accession: accession number for I{cluster.representative}
        @type rep_accession: str
        @param rep_chain: chain ID for I{cluster.representative}
        @type rep_chain: str
        @param alpha: alpha parameter for I{cluster.profile.pssm}
        @type alpha: float

        @return: an I-Sites cluster wrapper
        @rtype: L{isites.Cluster}

        @raise TypeError: when the HMM's structural data is missing or
                          interrupted
        """

        if not self.has_structure:
            raise TypeError('This HMM has no structural data assigned and ' +
                            'cannot be converted.')

        if start is None:
            start = self.layers.start_index
        if end is None:
            end = self.layers.last_index

        score_units = self.score_units
        self.convert_scores(ScoreUnits.Probability)

        background = dict(isites.ProteinProfile.BackgroundFreqs)
        for aa in self.layers[1][States.Match].emission:
            background[str(aa)] = self.layers[1][States.Match].background[aa]

        cluster = isites.Cluster()
        cluster.id = cluster_id
        cluster.overhang = 0
        cluster.profile = isites.ProteinProfile(background, alpha=alpha)

        residues = []
        for layer in self.layers:

            residues.append(csb.pyutils.deepcopy(layer.residue))

            if start <= layer.rank <= end:
                if not layer.residue.has_structure:
                    raise TypeError(
                        "The HMM's structure is interrupted at layer {0}.".format(
                            layer.rank))

                emission = {}

                for aa in layer[States.Match].emission:
                    emission[str(aa)] = layer[States.Match].emission[aa]

                cluster.profile.add_column(**dict(emission))

        cluster.profilelen = cluster.profile.length
        cluster.motiflen = cluster.profile.length - 2 * cluster.overhang

        if rep_accession is None:
            rep_accession = self.id[1:5].lower()
        if rep_chain is None:
            rep_chain = self.id[5].upper()

        chain = structure.Chain(rep_chain, sequence.SequenceTypes.Protein,
                                None, residues, rep_accession)
        chain.compute_torsion()

        source = isites.ChainSequence()
        source.sequence = chain.sequence
        source.accession = chain.accession
        source.id = chain.id
        source.type = chain.type
        source.torsion = chain.torsion

        cluster.representative = isites.RepStructureFragment(chain.accession,
                                                             chain.id, start)
        cluster.representative.source = source
        cluster.representative.structure = chain.subregion(start, end)
        cluster.representative.angles = cluster.representative.structure.torsion

        assert cluster.representative.angles.length == cluster.motiflen
        assert cluster.profile.length == cluster.representative.structure.length
        assert cluster.representative.angles.length == (cluster.profile.length -
                                                        2 * cluster.overhang)

        self.convert_scores(score_units)
        return cluster

    def emission_similarity(self, other):
        """
        Compute the Log-sum-of-odds score between the emission tables of self
        and other (Soeding 2004). If no observable Match state is found at a
        given layer, the Insertion state is used instead.

        @note: This is not a full implementation of the formula since only
        emission vectors are involved in the computation and any transition
        probabilities are ignored.

        @param other: the subject HMM
        @type other: L{ProfileHMM}

        @return: emission log-sum-of-odds similarity between C{self} and
                 C{other}
        @rtype: float

        @raise ValueError: when self and other differ in their length, when the
                           score_units are not Probability, or when no
                           observable states are present
        """
        score = 1

        if self.layers.length != other.layers.length or self.layers.length < 1:
            raise ValueError(
                'Both HMMs must have the same and positive number of layers')
        if self.score_units != ScoreUnits.Probability or \
               other.score_units != ScoreUnits.Probability:
            raise ValueError(
                'Scores must be converted to probabilities first.')

        for q_layer, s_layer in izip(iter(self.layers), iter(other.layers)):

            try:
                if States.Match in q_layer and \
                       not q_layer[States.Match].silent:
                    q_state = q_layer[States.Match]
                else:
                    q_state = q_layer[States.Insertion]

                if States.Match in s_layer and \
                       not s_layer[States.Match].silent:
                    s_state = s_layer[States.Match]
                else:
                    s_state = s_layer[States.Insertion]
            except csb.pyutils.ItemNotFoundError:
                raise ValueError('Query and subject must contain observable ' +
                                 'states at each layer')

            #assert set(q_state.emission) == set(s_state.emission)
            #assert len(q_state.emission) > 0

            emission_dotproduct = 0

            for aa in q_state.emission:
                #assert q_state.background[aa] == s_state.background[aa], 'Identical background probabilities expected'
                #assert q_state.background[aa] not in (0, None), 'Positive background probabilities expected'

                #assert q_state.emission[aa] not in (0, None)
                #assert s_state.emission[aa] not in (0, None), 'subject emission for {0!r}: {1}'.format(aa, s_state.emission[aa])

                q_emission = q_state.emission[aa] or sys.float_info.min
                s_emission = s_state.emission[aa] or sys.float_info.min
                emission_dotproduct += (q_emission * s_emission /
                                        q_state.background[aa])

            score *= emission_dotproduct

        return math.log(score)

    def _assign_secstructure(self):
        """
        Attach references from each profile layer to the relevant DSSP secondary 
        structure element.
        """        
        assert self.dssp is not None
        
        for motif in self.dssp:
            for i in range(motif.start, motif.end + 1):
                self.layers[i].residue.secondary_structure = motif    

class ProfileHMMSegment(ProfileHMM):
    """
    Represents a segment (fragment) of a ProfileHMM.
    
    @param hmm: source HMM
    @type hmm: ProfileHMM
    @param start: start layer of the segment (rank)
    @type start: int
    @param end: end layer of the segment (rank)
    @type end: int   
    
    @raise ValueError: when start or end positions are out of range
    """
    
    def __init__(self, hmm, start, end):

        if start < hmm.layers.start_index or start > hmm.layers.last_index:
            raise ValueError('Start position {0} is out of range {1.start_index} .. {1.last_index}'.format(start, hmm.layers))
        if end < hmm.layers.start_index or end > hmm.layers.last_index:
            raise ValueError('End position {0} is out of range {1.start_index} .. {1.last_index}'.format(end, hmm.layers))
        
        hmm = csb.pyutils.deepcopy(hmm)       
                
        super(ProfileHMMSegment, self).__init__(units=hmm.score_units, scale=hmm.scale, logbase=hmm.logbase)
        self.id = hmm.id
        self.family = hmm.family 
        self.name = hmm.name
        self.pseudocounts = hmm.pseudocounts
        self.effective_matches = hmm.effective_matches
        self.version = hmm.version       
        self.source = hmm.id
        self.source_start = start
        self.source_end = end
        
        self.alignment = hmm.alignment.subregion(start, end)

        layers = hmm.layers[start : end + 1]
        max_score = 1.0
        if hmm.score_units != ScoreUnits.Probability:
            max_score = hmm._convert(hmm.score_units, max_score, hmm.scale, hmm.logbase)
        self._build_graph(layers, max_score)
                        
        if hmm.dssp:
            self.dssp = hmm.dssp.subregion(start, end)
            self._assign_secstructure()
        if hmm.psipred:
            self.psipred = hmm.psipred.subregion(start, end)            
            
        self.length.layers = self.layers.length
        self.length.matches = self.alignment.matches_count            
        
    def _build_graph(self, source_layers, max_score):

        for rank, layer in enumerate(source_layers, start=1):
            
            for atom_kind in layer.residue.structure:
                layer.residue.structure[atom_kind].rank = rank
            layer.residue._rank = rank
            layer.rank = rank                
            
            self.layers.append(layer)
            
            if rank == 1:
                for state_kind in layer:
                    if state_kind in(States.Match, States.Deletion):
                        start_tran = Transition(self.start, layer[state_kind], max_score)
                        self.start.transitions.append(start_tran)
            elif rank == len(source_layers):
                for state_kind in layer:
                    state = layer[state_kind]
                    if not (States.End in state.transitions or States.Match in state.transitions):
                        state.transitions.set({})
                    else:
                        end_tran = Transition(state, self.end, max_score)
                        state.transitions.set({States.End: end_tran}) # TODO: I->I ?                      

    def slideover(self, other):
        """
        Find instances of self by sliding it over other and computing
        emission vector similarity and RSMD.
        
        @param other: the subject HMM
        @type other: L{ProfileHMM}
        
        @return: a list of L{isites.ProfileMatch}-es
        @rtype: list
        """

        query_residues = self.chain()        
        window = self.layers.length
        matches = []
        
        i = other.layers.start_index
        while True:
            if other.layers.last_index - i + 1 < window:
                break
            
            score, tm, tm_out, rmsd_ca = None, None, None, None
            start = i
            end = i + window - 1
            
            subject = ProfileHMMRegion(other, start, end)
            subject_residues = structure.Chain(other.id[4].upper(), residues=other.residues[start : start + window])  

            score = self.emission_similarity(subject)         
                        
            try:
                rmsd_ca = query_residues.rmsd(subject_residues)
            except structure.Broken3DStructureError:
                rmsd_ca = None                          
                                    
            if score is not None and (rmsd_ca is not None):
                match = isites.ProfileMatch(other.id, start, score, rmsd_ca, None, tm, tm_out)
                matches.append(match)               
            
            i += 1
                        
        return matches

class EmissionProfileSegment(ProfileHMMSegment):
    """
    Represents a segment of the Match state emission probabilities of a L{ProfileHMM}.
    Contains only Match states, connected with equal transition probabilities of 100%.   
    """      
    
    def _build_graph(self, source_layers):
        
        match_factory = State.Factory(States.Match)
        
        for rank, source_layer in enumerate(source_layers, start=1):
            
            emission = source_layer[States.Match].emission         
            background = source_layer[States.Match].background
            
            match = match_factory(emission, background)
            match.rank = rank
            
            layer = HMMLayer(rank, source_layer.residue)            
            layer.append(match)
            self.layers.append(layer)
                        
            if rank == 1:
                self.start.transitions.append(Transition(self.start, match, 1.0))
            elif rank < len(source_layers):
                prev_match = self.layers[rank - 1][States.Match]
                prev_match.transitions.append(Transition(prev_match, match, 1.0))
            elif rank == len(source_layers):
                match.transitions.append(Transition(match, self.end, 1.0))
            else:
                assert False
        
class ProfileHMMRegion(ProfileHMM):
    """
    A shallow proxy referring to a sub-region of a given Profile HMM.
    
    @param hmm: source HMM
    @type hmm: L{ProfileHMM}
    @param start: start layer of the segment (rank)
    @type start: int
    @param end: end layer of the segment (rank)
    @type end: int   
    
    @raise ValueError: when start or end positions are out of range    
    """
    
    def __init__(self, hmm, start, end):        
        
        if start < hmm.layers.start_index or start > hmm.layers.last_index:
            raise ValueError('Start position {0} is out of range {1.start_index} .. {1.last_index}'.format(start, hmm.layers))
        if end < hmm.layers.start_index or end > hmm.layers.last_index:
            raise ValueError('End position {0} is out of range {1.start_index} .. {1.last_index}'.format(end, hmm.layers))
        if hmm.score_units != ScoreUnits.Probability:
            raise ValueError('Scores must be converted to probabilities first.')
                
        self.layers = HMMLayersCollection(hmm.layers[start : end + 1])
        self._score_units = hmm.score_units
        self.id = hmm.id
        self.name = hmm.name
        self.family = hmm.family
        self.start = start
        self.end = end
                
class ProfileLength(object):
    
    def __init__(self, matches, layers):
        self.matches = matches
        self.layers = layers

class EmissionTable(csb.pyutils.DictionaryContainer):        
    """ 
    Represents a lookup table of emission probabilities. Provides dictionary-like access:
    
        >>> state.emission['A']
        emission probability of residue A
    
    @param emission: an initialization dictionary of emission probabilities
    @type emission: dict
    @param restrict: a list of residue types allowed for this emission table. 
                     Defaults to the members of I{SequenceAlphabets.Protein}
    @type restrict: list
    """

    def __init__(self, emission=None, restrict=csb.pyutils.Enum.members(sequence.SequenceAlphabets.Protein)):
        super(EmissionTable, self).__init__(emission, restrict)
            
    def append(self, residue, probability):
        """
        Append a new emission probability to the table.
        
        @param residue: residue name (type) - a member of I{SequenceAlphabets.Protein}
        @type residue: L{csb.pyutils.EnumItem}
        @param probability: emission score
        @type probability: float
        
        @raise EmissionExistsError: if residue is already defined
        """
        if residue in self:
            raise EmissionExistsError('State {0} is already emittable. Try x.update()'.format(residue))

        super(EmissionTable, self).append(residue, probability)
    
    def set(self, table):
        """ 
        Set the emission table using the dictionary provided in the argument.
        
        @param table: the new emission table
        @type table: dict
        """          
        super(EmissionTable, self).set(table)
        
    def update(self, residue, probability):
        """ 
        Update the emission C{probability} of a given emission C{residue}.
        
        @param residue: name (type) of the residue to be updated
        @type residue: L{csb.pyutils.EnumItem}
        @param probability: new emission score
        @type probability: float
        """                
        super(EmissionTable, self).update({residue: probability})
        
class TransitionTable(csb.pyutils.DictionaryContainer):        
    """ 
    Represents a lookup table of transitions that are possible from within a given state. 
    
    Provides dictionary-like access, where dictionary keys are target states. These are members of 
    the L{States} enumeration, e.g.:
    
        >>> state.transitions[States.Match]
        transition info regarding transition from the current state to a Match state
        >>> state.transitions[States.Match].predecessor
        state
        >>> state.transitions[States.Match].successor
        the next match state
     
    @param transitions: an initialization dictionary of target L{State} : L{Transition} pairs
    @type transitions: dict
    @param restrict: a list of target states allowed for this transition table. 
                     Defaults to the L{States} enum members
    @type restrict: list
    """
    
    def __init__(self, transitions=None, restrict=csb.pyutils.Enum.members(States)):
        
        super(TransitionTable, self).__init__(transitions, restrict)
    
    def append(self, transition):
        """
        Append a new C{transition} to the table.

        @param transition: transition info
        @type transition: L{Transition}
        
        @raise TransitionExistsError: when a transition to the same target state
                                      already exists for the current state
        """
        if transition.successor.type in self:
            raise TransitionExistsError('A transition to a state of type {0} is already defined.'.format(transition.successor.type))

        super(TransitionTable, self).append(transition.successor.type, transition)
    
    def set(self, table):
        """ 
        Set the transition table using the dictionary provided in the argument.
        
        @param table: the new transition table
        @type table: dict
        """          
        super(TransitionTable, self).set(table)
        
    def update(self, target_statekind, transition):
        """ 
        Update the transition information of a transition, which points to a target 
        state of the specified L{States} kind.
        
        @param target_statekind: transition table key; the key of the transition to be updated
        @type target_statekind: L{csb.pyutils.EnumItem}
        @param transition: new transition info object
        @type transition: L{Transition}
        
        @raise ValueError: if I{transition.successor.type} differs from C{target_statekind}
        """
        if transition.successor.type != target_statekind:
            raise ValueError('Transition successor\'type differs from the specified target_statekind.')
                
        super(TransitionTable, self).update({target_statekind: transition})  
        
class HMMLayersCollection(csb.pyutils.CollectionContainer):
    """
    Provides consecutive, 1-based access to all of the layers in the profile.
    Each profile layer contains a catalog of available states at that index, e.g.:

        >>> profile.layers[i]
        the catalog at profile layer i
        >>> profile.layers[i][States.Deletion]
        the deletion state at index i
        
    @param layers: initialization list of L{HMMLayer}s
    @type layers: list    
    """        
    def __init__(self, layers=None):
        super(HMMLayersCollection, self).__init__(layers, type=HMMLayer, start_index=1)

class HMMLayer(csb.pyutils.DictionaryContainer):
    """
    Provides a dictionary-like catalog of the available states at specific index(layer).
    Lookup keys are members of the L{States} enumeration, e.g.:
    
        >>> profile.layers[i][States.Deletion]
        the deletion state at index i  
    
    @param rank: layer's index
    @type rank: int
    @param residue: a representative L{structure.ProteinResidue} that is associated with this layer
    @type residue: L{structure.ProteinResidue}        
    @param states: initialization dictionary of L{States} type : L{State} info pairs
    @type states: dict 
    """        
    def __init__(self, rank, residue, states=None):
        
        self.rank = int(rank)
        self.residue = residue        
        self.effective_matches = None
        self.effective_insertions = None
        self.effective_deletions = None
        
        super(HMMLayer, self).__init__(states, restrict=csb.pyutils.Enum.members(States))

    def append(self, state):
        """
        Append a new C{state} to the catalog.
        
        @param state: the new state
        @type state: L{State}
        
        @raise StateExistsError: when a state of the same type is already defined
        """
        if state.type in self:
            raise StateExistsError('State {0} is already defined at this position.'.format(state.type))

        super(HMMLayer, self).append(state.type, state)
        
    def update(self, state_kind, state):
        """ 
        Update the sate of the specified kind under the current layer.
        
        @param state_kind: state type (key) - a member of L{States}
        @type state_kind: L{csb.pyutils.EnumItem}
        @param state: the new state info
        @type state: L{State}
        
        @raise ValueError: if state.type differs from state_kind
        """
        if state.type != state_kind:
            raise ValueError('State\'type differs from the specified state_kind.')
                
        super(HMMLayer, self).update({state_kind: state})        
        
class State(object):
    """ 
    Describes a Hidden Markov Model hidden state.
    
    @param type: one of the L{States} enumeration values, e.g. States.Match
    @type type: L{csb.pyutils.EnumItem}
    @param emit: a collection of emittable state names allowed for the state, 
                 e.g. the members of I{SequenceAlphabets.Protein}. If not defined, 
                 the state will be created as a silent (unobservable).
    @type emit: list
    
    @raise ValueError: if type is not a member of the States enum
    """
    
    def __init__(self, type, emit=None):

        if type not in csb.pyutils.Enum.members(States):
            raise ValueError('Unknown state type {0}. Did you use the fragments.hmm.States enum?'.format(type)) 
        
        self.type = type
        self.rank = None      
        self.transitions = TransitionTable()
        self.emission = None
        self.background = None
        
        if emit is not None:
            self.emission = EmissionTable(restrict=emit)
            self.background = EmissionTable(restrict=emit)    
    
    def __repr__(self):
        return "<HMM {0.type.name} State>".format(self)
        
    def __getattribute__(self, name):
        if name == 'emission' and object.__getattribute__(self, name) is None:
            raise UnobservableStateError('{0.name} state is not exposing any observations.'.format(self.type))
        else:
            return object.__getattribute__(self, name)   
        
    @property
    def silent(self):
        try:
            return self.emission is None
        except UnobservableStateError:
            return True     
        
    @staticmethod
    def Factory(state_kind):
        """
        Return a callable function that automates and simplifies the construction 
        of states of the specified kind.
        
        @param state_kind: type of the states, produced by the factory (a L{States} member)
        @type state_kind: L{csb.pyutils.EnumItem}
        
        @return: a callable factory function
        @rtype: function
        
        @raise ValueError: if state_kind is not a member of the States enum
        """
        
        def match_factory(emission, background):
            state = State(States.Match, emit=csb.pyutils.Enum.members(sequence.SequenceAlphabets.Protein))
            state.emission.set(emission)
            state.background.set(background)
            return state

        def insertion_factory(background):
            state = State(States.Insertion, emit=csb.pyutils.Enum.members(sequence.SequenceAlphabets.Protein))
            state.emission.set(background)
            state.background.set(background)
            return state     
        
        def deletion_factory(): 
            return State(States.Deletion)               
            
        if state_kind == States.Match:
            return match_factory
        elif state_kind == States.Insertion:
            return insertion_factory
        elif state_kind == States.Deletion:
            return deletion_factory
        else:
            raise ValueError(state_kind)               
        
class MatchState(State):
    
    def __init__(self, emission, background):
        
        super(MatchState, self).__init__(emit=csb.pyutils.Enum.members(sequence.SequenceAlphabets.Protein))
        
        self.emission.set(emission)
        self.background.set(background)      
        
class InsertionState(State):
    
    def __init__(self, background):
        
        super(InsertionState, self).__init__(emit=csb.pyutils.Enum.members(sequence.SequenceAlphabets.Protein))
        
        self.emission.set(background)
        self.background.set(background)               

class TransitionType(object):

    def __init__(self, source, target):
        self.source_state = source.type
        self.target_state = target.type     
        
    def __repr__(self):
        return '{0}->{1}'.format(self.source_state, self.target_state)   
            
class Transition(object):
    """
    Describes a Hidden Markov Model transition between two states.
    
    @param predecessor: source state
    @type predecessor: L{State}
    @param successor: target state
    @type successor: L{State}
    @param probability: transition score
    @type probability: float
    
    @raise ValueError: if probability is not >= 0
    """
    
    def __init__(self, predecessor, successor, probability):
        
        if not (isinstance(predecessor, State) or isinstance(successor, State)):
            raise TypeError('Transition predecessor and successor states must be State instances.')
        if not (probability >=0):
            raise ValueError('Transition probability must be a positive number.')
         
        assert all([predecessor.type, successor.type])
                
        self.predecessor = predecessor
        self.successor = successor
        self.probability = float(probability)    
        self.type = TransitionType(predecessor, successor)
        
    def __str__(self):
        return '<HMM Transition: {0.type} {0.probability}>'.format(self)
                
_AlignedPair = namedtuple('AlignedPair', 'query subject insertions')


class HHpredHitAlignment(object):
    """
    Represents a query-template alignment in a HHpred result.

    @param query: the query sequence in the alignment region, with gaps
    @type query: str
    @param subject: the subject sequence in the alignment region, with gaps
    @type subject: str

    @raise ValueError: if query and subject differ in their lengths; or when
                       they contain invalid characters
    """

    def __init__(self, hit, query, subject):

        try:
            query = list(query)
            subject = list(subject)
        except:
            raise TypeError('query and subject must be iterable')

        if not isinstance(hit, HHpredHit):
            raise TypeError(hit)

        if not (len(query) == len(subject)) or not len(query) > 0:
            raise ValueError(
                'query and subject must be of the same and positive length (got {0} and {1})'.format(
                    len(query), len(subject)))

        self._query = []
        self._subject = []
        self._hit = hit

        for q, s in izip(query, subject):

            q, s = str(q), str(s)

            for i in (q, s):
                if not (i.isalpha() or i == '-'):
                    raise ValueError('Invalid residue: {0}'.format(i))

            if q != '-':
                self._query.append([q])
                self._subject.append([s])
            else:
                self._query[-1].append(q)
                self._subject[-1].append(s.lower())

    def __str__(self):
        return '\n{0.qstart:4} {1} {0.qend:4}\n{0.start:4} {2} {0.end:4}\n'.format(
            self._hit, self.query, self.subject)

    def __getitem__(self, position):

        if (1 <= position < self._hit.qstart) or \
               (self._hit.qend < position <= self._hit.qlength):
            return _AlignedPair(' ', ' ', [])

        i = position - self._hit.qstart
        try:
            if i < 0:
                raise IndexError(position)
            return _AlignedPair(self._query[i][0], self._subject[i][0],
                                self._subject[i][1:])
        except IndexError:
            raise IndexError(position)

    def __iter__(self):
        for i in range(1, self._hit.qlength + 1):
            yield self[i]

    @property
    def query(self):
        return ''.join([''.join(pos) for pos in self._query])

    @property
    def subject(self):
        return ''.join([''.join(pos) for pos in self._subject])

    def to_a3m(self):
        """
        Format the sequences in the alignment as A3M strings.

        @return: an object containing the sequences of the query and the
                 subject, formatted as A3M strings.
        @rtype: AlignedPair
        """

        query = []
        subject = []

        for i in range(1, self._hit.qlength + 1):

            pos = self[i]

            if pos.query != '-':
                query.append(pos.query.upper().replace(' ', '-'))

            subject.append(pos.subject.upper().replace(' ', '-'))

            if pos.insertions:
                subject.append(''.join(pos.insertions).lower())

        query = '>query\n{0}'.format(''.join(query))
        subject = '>{0}\n{1}'.format(self._hit.id, ''.join(subject))

        return _AlignedPair(query, subject, '')


class HHpredHit(object):
    """
    Represents a single HHsearch hit.

    @param rank: rank of the hit
    @type rank: int
    @param id: id of the hit
    @type id: str
    @param start: subject start
    @type start: int
    @param end: subject end
    @type end: int
    @param qstart: query start
    @type qstart: int
    @param qend: query end
    @type qend: int
    @param probability: probability of the hit
    @type probability: float
    @param qlength: length of the query
    @type qlength: int
    """

    def __init__(self, rank, id, start, end, qstart, qend, probability,
                 qlength):

        self._rank = None
        self._id = None
        self._start = None
        self._end = None
        self._qstart = None
        self._qend = None
        self._probability = None
        self._qlength = None
        self._alignment = None

        self.rank = rank
        self.id = id
        self.start = start
        self.end = end
        self.qstart = qstart
        self.qend = qend
        self.probability = probability
        self.qlength = qlength

        self.slength = None
        self.evalue = None
        self.pvalue = None
        self.score = None
        self.identity = None
        self.similarity = None
        self.prob_sum = None

    def __str__(self):
        return '{0.id} {0.probability} {0.start}-{0.end}'.format(self)

    def __repr__(self):
        return '<HHpredHit: {0!s}>'.format(self)

    def __cmp__(self, other):
        '''TODO: __cmp__ and the rich comparisons (__eq__ and __gt__ in this
        case) should compare the same fields (e.g. rank), otherwise this is
        confusing for users. Ivan wants to check his code for uses of these
        functions, hence they are retained for the moment.'''
        return cmp(self.rank, other.rank)

    def __eq__(self, other):
        '''TODO: __cmp__ and the rich comparisons (__eq__ and __gt__ in this
        case) should compare the same fields (e.g. rank), otherwise this is
        confusing for users. Ivan wants to check his code for uses of these
        functions, hence they are retained for the moment.'''
        return (self.id == other.id and self.start == other.start and
                self.end == other.end)

    def __gt__(self, other):
        '''TODO: __cmp__ and the rich comparisons (__eq__ and __gt__ in this
        case) should compare the same fields (e.g. rank), otherwise this is
        confusing for users. Ivan wants to check his code for uses of these
        functions, hence they are retained for the moment.'''
        if self.length > other.length:
            return True
        elif self.length == other.length:
            if self.probability > other.probability:
                return True
            elif self.probability == other.probability:
                if id(self) > id(other):
                    return True
        return False

    def add_alignment(self, query, subject):
        """
        Add query/subject alignment to the hit.

        @param query: the query sequence within the alignment region, with gaps
        @type query: str
        @param subject: the subject sequence within the alignment region, with
                        gaps
        @type subject: str
        """
        self._alignment = HHpredHitAlignment(self, query, subject)

    @property
    def rank(self):
        return self._rank
    @rank.setter
    def rank(self, value):
        try:
            value = int(value)
        except:
            raise TypeError('rank must be int, not {1}'.format(type(value)))
        self._rank = value

    @property
    def id(self):
        return self._id
    @id.setter
    def id(self, value):
        try:
            value = str(value)
        except:
            raise TypeError('id must be string, not {0}'.format(type(value)))
        self._id = value

    @property
    def start(self):
        return self._start
    @start.setter
    def start(self, value):
        try:
            value = int(value)
        except:
            raise TypeError('start must be int, not {0}'.format(type(value)))
        self._start = value

    @property
    def end(self):
        return self._end
    @end.setter
    def end(self, value):
        try:
            value = int(value)
        except:
            raise TypeError('end must be int, not {0}'.format(type(value)))
        self._end = value

    @property
    def qstart(self):
        return self._qstart
    @qstart.setter
    def qstart(self, value):
        try:
            value = int(value)
        except:
            raise TypeError('qstart must be int, not {0}'.format(type(value)))
        self._qstart = value

    @property
    def qend(self):
        return self._qend
    @qend.setter
    def qend(self, value):
        try:
            value = int(value)
        except:
            raise TypeError('qend must be int, not {0}'.format(type(value)))
        self._qend = value

    @property
    def qlength(self):
        return self._qlength
    @qlength.setter
    def qlength(self, value):
        try:
            value = int(value)
        except:
            raise TypeError('qlength must be int, not {0}'.format(type(value)))
        self._qlength = value

    @property
    def probability(self):
        return self._probability
    @probability.setter
    def probability(self, value):
        try:
            value = float(value)
        except:
            raise TypeError(
                'probability must be float, not {0}'.format(type(value)))
        self._probability = value

    @property
    def alignment(self):
        return self._alignment

    @property
    def length(self):
        try:
            return self.end - self.start + 1
        except:
            return 0

    def includes(self, other, tolerance=1):
        """
        Return True if other overlaps with self

        @type other: HHpredHit
        @param tolerance: allow partial overlaps for that number of residues
                          at either end
        @type tolerance: int

        @rtype: bool
        """
        if self.id == other.id:
            if other.start >= self.start:
                if other.end <= self.end or \
                       other.start <= (self.end - tolerance):
                    return True
            elif other.end <= self.end:
                if other.end >= (self.start + tolerance):
                    return True

        return False


class HHpredHitList(object):

    def __init__(self, hits, query_name='', match_columns=-1, no_of_seqs='',
                 neff=-1., searched_hmms=-1, date='', command=''):
        self._hits = list(hits)

        self.query_name = query_name
        self.match_columns = match_columns
        self.no_of_seqs = no_of_seqs
        self.neff = neff
        self.searched_hmms = searched_hmms
        self.date = date
        self.command = command

    def __getitem__(self, index):
        return self._hits[index]

    def __iter__(self):
        return iter(self._hits)

    def __len__(self):
        return len(self._hits)

    def sort(self):
        from operator import attrgetter
        self._hits.sort(key=attrgetter('rank'))

    def align(self, query_sequence):
        """
        Build a pseudo-multiple HMM alignment, assuming all hits in the
        hitlist have been aligned against the same C{query_sequence}.

        @param query_sequence: must be the complete sequence of the
                               representative in the query HHM file
        @type query_sequence: str

        @return: a new A3M alignment object
        @rtype: L{sequence.A3MAlignment}

        @raise ValueError: if any hit does not contain an alignment with the
                           query; or when the query sequence has a mismatching
                           length
        """
        qseq = ''.join(query_sequence.split()).upper()
        query = '>query\n{0}'.format(qseq)
        mha = [query]

        for hit in self:

            if not hit.alignment:
                raise ValueError(
                    'Hit {0.id} does not contain an alignment'.format(hit))
            elif hit.qlength != len(qseq):
                raise ValueError(
                    'The query sequence is {0} long, while the query in hit {1.rank}.{1.id} was {2}'.format(
                        len(qseq), hit, hit.qlength))
            a3m = hit.alignment.to_a3m()
            mha.append(a3m.subject)

        return sequence.A3MAlignment.parse('\n'.join(mha))
