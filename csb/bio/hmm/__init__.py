"""
HHpred and Hidden Markov Model APIs. 
"""

import sys
import math
import cStringIO

import csb.pyutils
import csb.io
import csb.bio.structure as structure
import csb.bio.sequence as sequence


class UnobservableStateError(AttributeError):
    pass

class StateExistsError(KeyError):
    pass

class TransitionExistsError(KeyError):
    pass

class EmissionExistsError(KeyError):
    pass

class HMMArgumentError(ValueError):
    pass


class States(csb.pyutils.enum):
    """
    Enumeration of HMM state types
    """
    Match='M'; Insertion='I'; Deletion='D'; Start='S'; End='E'

class ScoreUnits(csb.pyutils.enum):
    """
    Enumeration of HMM emission and transition score units
    """
    LogScales='LogScales'; Probability='Probability'
    
BACKGROUND = [ 0.076627178753322270, 0.018866884241976509, 0.053996136712517316,
               0.059788009880742142, 0.034939432842683173, 0.075415244982547675, 
               0.036829356494115069, 0.050485048600600511, 0.059581159080509941, 
               0.099925728794059046, 0.021959667190729986, 0.040107059298840765, 
               0.045310838527464106, 0.032644867589507229, 0.051296350550656143, 
               0.046617000834108295, 0.071051060827250878, 0.072644631719882335, 
               0.012473412286822654, 0.039418044025976547 ]
"""
Background amino acid probabilities
"""

RELATIVE_SA = { 'A': 0.02, 'B': 0.14, 'C': 0.33, 'D': 0.55, 'E': 1.00 }
"""
Relative solvent accessibility codes (upper bounds)
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

        self._name = None
        self._id = None
        self._family = None
        self._length = ProfileLength(0, 0)
        self._alignment = None
        self._consensus = None
        self._dssp = None
        self._dssp_solvent = None
        self._psipred = None
        self._effective_matches = None
        self._evd = EVDParameters(None, None) 
        self._version = None
        self._pseudocounts = False
        self._emission_pseudocounts = False
        self._transition_pseudocounts = False
        self._layers = HMMLayersCollection()
        self._start = State(States.Start)
        self._start_insertion = None
        self._end = State(States.End)
        self._scale = scale
        self._logbase = logbase
        if units is None:
            self._score_units = ScoreUnits.LogScales
        else:
            self._score_units = units        

    @property
    def name(self):
        return self._name
    @name.setter
    def name(self, value):
        self._name = str(value)
    
    @property
    def id(self):
        return self._id
    @id.setter
    def id(self, value):
        self._id = str(value)
    
    @property
    def family(self):
        return self._family
    @family.setter
    def family(self, value):
        self._family = str(value)
    
    @property
    def length(self):
        return self._length
    @length.setter
    def length(self, value):
        if not isinstance(value, ProfileLength):
            raise TypeError(value)
        self._length = value
    
    @property
    def alignment(self):
        return self._alignment
    @alignment.setter
    def alignment(self, value):
        if not isinstance(value, sequence.A3MAlignment):
            raise TypeError(value)        
        self._alignment = value
    
    @property
    def consensus(self):
        return self._consensus
    @consensus.setter
    def consensus(self, value):
        if not isinstance(value, sequence.AbstractSequence):
            raise TypeError(value)         
        self._consensus = value
    
    @property
    def dssp(self):
        return self._dssp
    @dssp.setter
    def dssp(self, value):
        if not isinstance(value, structure.SecondaryStructure):
            raise TypeError(value) 
        self._dssp = value
    
    @property
    def dssp_solvent(self):
        return self._dssp_solvent
    @dssp_solvent.setter
    def dssp_solvent(self, value):
        self._dssp_solvent = str(value)
    
    @property
    def psipred(self):
        return self._psipred
    @psipred.setter
    def psipred(self, value):
        if not isinstance(value, structure.SecondaryStructure):
            raise TypeError(value)         
        self._psipred = value
    
    @property
    def effective_matches(self):
        return self._effective_matches
    @effective_matches.setter
    def effective_matches(self, value):
        self._effective_matches = value
    
    @property
    def evd(self):
        return self._evd
    @evd.setter
    def evd(self, value):
        if not isinstance(value, EVDParameters):
            raise TypeError(value) 
        self._evd = value
    
    @property
    def version(self):
        return self._version
    @version.setter
    def version(self, value):
        self._version = str(value)
    
    @property
    def pseudocounts(self):
        return self._pseudocounts
    @pseudocounts.setter
    def pseudocounts(self, value):
        self._pseudocounts = bool(value)
    
    @property
    def emission_pseudocounts(self):
        return self._emission_pseudocounts
    @emission_pseudocounts.setter
    def emission_pseudocounts(self, value):
        self._emission_pseudocounts = bool(value)
    
    @property
    def transition_pseudocounts(self):
        return self._transition_pseudocounts
    @transition_pseudocounts.setter
    def transition_pseudocounts(self, value):
        self._transition_pseudocounts = bool(value)
    
    @property
    def layers(self):
        return self._layers
    
    @property
    def start(self):
        return self._start
    @start.setter
    def start(self, value):
        if value is None or (isinstance(value, State) and value.type == States.Start):
            self._start = value
        else:
            raise TypeError(value)
    
    @property
    def start_insertion(self):
        return self._start_insertion
    @start_insertion.setter
    def start_insertion(self, value):
        if value is None or (isinstance(value, State) and value.type == States.Insertion):        
            self._start_insertion = value
        else:
            raise TypeError(value)    
    
    @property
    def end(self):
        return self._end
    @end.setter
    def end(self, value):
        if value is None or (isinstance(value, State) and value.type == States.End):        
            self._end = value
        else:
            raise TypeError(value)      
    
    @property
    def scale(self):
        return self._scale
    
    @property
    def logbase(self):
        return self._logbase
    
    @property
    def score_units(self):
        return self._score_units
    
    @property
    def residues(self):
        res = [layer.residue for layer in self.layers]
        return csb.pyutils.ReadOnlyCollectionContainer(
                            res, type=structure.Residue, start_index=1)

    @property
    def all_layers(self):
        """
        A list of layers including start and start_insertion
        """
        complete_layers = []
        first_layer = HMMLayer(rank=0, residue=None)
        first_layer.append(self.start)
        if self.start_insertion:
            first_layer.append(self.start_insertion)
        complete_layers.append(first_layer)
        for layer in self.layers:
            complete_layers.append(layer)
        
        return complete_layers

    @property
    def has_structure(self):
        has = False
        for layer in self.layers:
            if layer.residue.has_structure:
                return True
        return has

    def serialize(self, file_name):
        """
        Serialize this HMM to a file.

        @param file_name: target file name
        @type file_name: str
        """
        import cPickle
        rec = sys.getrecursionlimit()        
        sys.setrecursionlimit(10000)
        cPickle.dump(self, open(file_name, 'w'))
        sys.setrecursionlimit(rec)

    @staticmethod
    def deserialize(file_name):
        """
        De-serialize an HMM from a file.

        @param file_name: source file name (pickle)
        @type file_name: str
        """
        import cPickle
        rec = sys.getrecursionlimit()
        sys.setrecursionlimit(10000)
        return cPickle.load(open(file_name, 'r'))
        sys.setrecursionlimit(rec)

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
        if convert_scores:
            self.convert_scores(ScoreUnits.LogScales)
                        
        temp = cStringIO.StringIO()
        
        builder = HHMFileBuilder(temp)
        builder.add_hmm(self)
        
        data = temp.getvalue()        
        temp.close()
        
        if not output_file:
            return data
        else:
            with csb.io.EntryWriter(output_file, close=False) as out:
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
    
    def subregion(self, start, end):
        
        return ProfileHMMRegion(self, start, end) 

    def add_emission_pseudocounts(self, *a, **k):
        """
        See L{csb.bio.hmm.pseudocounts.PseudocountBuilder}
        """
        from csb.bio.hmm.pseudocounts import PseudocountBuilder
        PseudocountBuilder(self).add_emission_pseudocounts(*a, **k)

    def add_transition_pseudocounts(self, *a, **k):
        """
        See L{csb.bio.hmm.pseudocounts.PseudocountBuilder}
        """        
        from csb.bio.hmm.pseudocounts import PseudocountBuilder
        PseudocountBuilder(self).add_transition_pseudocounts(*a, **k)
                                  
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
            # The Neff-s are interger numbers and should not be transformed
            # (except when writing the profile to a hhm file)

        if self.start_insertion:
            for t_it in self.start_insertion.transitions:
                transition = self.start_insertion.transitions[t_it]
                transition.probability = convert(units, transition.probability, self.scale, self.logbase)

            for residue in self.start_insertion.emission:
                state = self.start_insertion
                if state.emission[residue] is not None:
                    state.emission.update(residue, convert(units, state.emission[residue],self.scale, self.logbase))
                    state.background.update(residue, convert(units, state.background[residue],self.scale, self.logbase))

        for tran_kind in self.start.transitions:
            transition = self.start.transitions[tran_kind]
            transition.probability = convert(units, transition.probability, self.scale, self.logbase)


        self._score_units = units

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

        for q_layer, s_layer in zip(self.layers, other.layers):

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
            raise IndexError('Start position {0} is out of range {1.start_index} .. {1.last_index}'.format(start, hmm.layers))
        if end < hmm.layers.start_index or end > hmm.layers.last_index:
            raise IndexError('End position {0} is out of range {1.start_index} .. {1.last_index}'.format(end, hmm.layers))
        
        #hmm = csb.pyutils.deepcopy(hmm)       
                
        super(ProfileHMMSegment, self).__init__(units=hmm.score_units, scale=hmm.scale, logbase=hmm.logbase)
        self.id = hmm.id
        self.family = hmm.family 
        self.name = hmm.name
        self.pseudocounts = hmm.pseudocounts
        self.evd = hmm.evd
        self.version = hmm.version       
        self.source = hmm.id
        self._source_start = start
        self._source_end = end
        
        if hmm.alignment:
            self.alignment = hmm.alignment.hmm_subregion(start, end)
            self.consensus = hmm.consensus.subregion(start, end)

        layers = csb.pyutils.deepcopy(hmm.layers[start : end + 1])
        max_score = 1.0
        if hmm.score_units != ScoreUnits.Probability:
            max_score = hmm._convert_scores(hmm.score_units, max_score, hmm.scale, hmm.logbase)
        self._build_graph(layers, max_score)
                        
        if hmm.dssp:
            self.dssp = hmm.dssp.subregion(start, end)
            self._assign_secstructure()
        if hmm.psipred:
            self.psipred = hmm.psipred.subregion(start, end)            
            
        self.length.layers = self.layers.length
        self.length.matches = self.layers.length
        self.effective_matches = sum([(l.effective_matches or 0.0) for l in self.layers]) / self.layers.length   

    @property
    def source_start(self):
        return self._source_start

    @property
    def source_end(self):
        return self._source_end
            
    def _build_graph(self, source_layers, max_score):

        for rank, layer in enumerate(source_layers, start=1):
            
            for atom_kind in layer.residue.atoms:
                layer.residue.atoms[atom_kind].rank = rank
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


class EmissionProfileSegment(ProfileHMMSegment):
    """
    Represents a segment of the Match state emission probabilities of a L{ProfileHMM}.
    Contains only Match states, connected with equal transition probabilities of 100%.   
    """      
    
    def _build_graph(self, source_layers):
        
        factory = StateFactory()
        
        for rank, source_layer in enumerate(source_layers, start=1):
            
            emission = source_layer[States.Match].emission         
            background = source_layer[States.Match].background
            
            match = factory.create_match(emission, background)
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
            raise IndexError('Start position {0} is out of range {1.start_index} .. {1.last_index}'.format(start, hmm.layers))
        if end < hmm.layers.start_index or end > hmm.layers.last_index:
            raise IndexError('End position {0} is out of range {1.start_index} .. {1.last_index}'.format(end, hmm.layers))
        if hmm.score_units != ScoreUnits.Probability:
            raise ValueError('Scores must be converted to probabilities first.')
                
        self._layers = HMMLayersCollection(hmm.layers[start : end + 1])
        self._score_units = hmm.score_units
        self.id = hmm.id
        self.name = hmm.name
        self.family = hmm.family
        self._source_start = start
        self._source_end = end

    @property
    def source_start(self):
        return self._source_start
    
    @property
    def source_end(self):
        return self._source_end

            
class ProfileLength(object):
    
    def __init__(self, matches, layers):
        self.matches = matches
        self.layers = layers


class EVDParameters(object):
    
    def __init__(self, lamda, mu):
        self.lamda = lamda
        self.mu = mu
    
    def __nonzero__(self):
        return (self.lamda is not None or self.mu is not None)


class EmissionTable(csb.pyutils.DictionaryContainer):        
    """ 
    Represents a lookup table of emission probabilities. Provides dictionary-like access:
    
        >>> state.emission['A']
        emission probability of residue A
    
    @param emission: an initialization dictionary of emission probabilities
    @type emission: dict
    @param restrict: a list of residue types allowed for this emission table. 
                     Defaults to the members of L{csb.bio.sequence.ProteinAlphabet}
    @type restrict: list
    """

    def __init__(self, emission=None, restrict=csb.pyutils.Enum.members(sequence.ProteinAlphabet)):
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
        super(EmissionTable, self)._set(table)
        
    def update(self, residue, probability):
        """ 
        Update the emission C{probability} of a given emission C{residue}.
        
        @param residue: name (type) of the residue to be updated
        @type residue: L{csb.pyutils.EnumItem}
        @param probability: new emission score
        @type probability: float
        """                
        super(EmissionTable, self)._update({residue: probability})


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
            raise TransitionExistsError('Transition to a {0} state is already defined.'.format(
                                                                        transition.successor.type))
        super(TransitionTable, self).append(transition.successor.type, transition)
    
    def set(self, table):
        """ 
        Set the transition table using the dictionary provided in the argument.
        
        @param table: the new transition table
        @type table: dict
        """          
        super(TransitionTable, self)._set(table)
        
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
            raise ValueError("Transition successor'type differs from the specified target_statekind.")
                
        super(TransitionTable, self)._update({target_statekind: transition})  


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

        super(HMMLayer, self).__init__(states, restrict=csb.pyutils.Enum.members(States))
                
        self._rank = int(rank)
        self._residue = None   
        self._effective_matches = None
        self._effective_insertions = None
        self._effective_deletions = None
        
        self.residue = residue

    @property
    def rank(self):
        return self._rank
    @rank.setter
    def rank(self, value):
        self._rank = int(value)

    @property
    def residue(self):
        return self._residue
    @residue.setter
    def residue(self, residue):
        if residue and residue.type == sequence.SequenceAlphabets.Protein.GAP:          
            raise HMMArgumentError('HMM match states cannot be gaps')
        self._residue = residue
        
    @property
    def effective_matches(self):
        return self._effective_matches
    @effective_matches.setter
    def effective_matches(self, value):
        self._effective_matches = value
    
    @property
    def effective_insertions(self):
        return self._effective_insertions
    @effective_insertions.setter
    def effective_insertions(self, value):
        self._effective_insertions = value
    
    @property
    def effective_deletions(self):
        return self._effective_deletions
    @effective_deletions.setter
    def effective_deletions(self, value):
        self._effective_deletions = value
    
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
            raise ValueError("State's type differs from the specified state_kind")
                
        super(HMMLayer, self)._update({state_kind: state})        

   
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
        
        self._type = None
        self._rank = None      
        self._transitions = TransitionTable()
        self._emission = None
        self._background = None
        
        self.type = type
        
        if emit is not None:
            self._emission = EmissionTable(restrict=emit)
            self._background = EmissionTable(restrict=emit)    
    
    def __repr__(self):
        return "<HMM {0.type!r} State>".format(self)

    @property
    def type(self):
        return self._type
    @type.setter
    def type(self, value):
        if value.enum is not States:
            raise TypeError(value)         
        self._type = value
    
    @property
    def rank(self):
        return self._rank
    @rank.setter
    def rank(self, value):
        self._rank = int(value)
    
    @property
    def transitions(self):
        return self._transitions
    
    @property
    def emission(self):
        if self._emission is None:
            raise UnobservableStateError('Silent {0!r} state'.format(self.type))
        return self._emission
    
    @property
    def background(self):
        return self._background
        
    @property
    def silent(self):
        try:
            return self.emission is None
        except UnobservableStateError:
            return True     
                    
class StateFactory(object):
    """
    Simplifies the construction protein profile HMM states.
    """
            
    def __init__(self):
        self._aa = csb.pyutils.Enum.members(sequence.ProteinAlphabet)        
    
    def create_match(self, emission, background):
        
        state = State(States.Match, emit=self._aa)
        state.emission.set(emission)
        state.background.set(background)
        return state

    def create_insertion(self, background):
        
        state = State(States.Insertion, emit=self._aa)
        state.emission.set(background)
        state.background.set(background)
        return state     
    
    def create_deletion(self): 
        return State(States.Deletion)
                

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
    """
    
    def __init__(self, predecessor, successor, probability):
        
        if not (isinstance(predecessor, State) or isinstance(successor, State)):
            raise TypeError('Transition predecessor and successor states must be State instances.')
                
        self._predecessor = predecessor
        self._successor = successor
        self._probability = None    
        self._type = TransitionType(predecessor, successor)
        
        self.probability = probability
        
    def __str__(self):
        return '<HMM Transition: {0.type} {0.probability}>'.format(self)
    
    @property
    def predecessor(self):
        return self._predecessor
    
    @property
    def successor(self):
        return self._successor
    
    @property
    def probability(self):
        return self._probability
    @probability.setter
    def probability(self, value):
        if not (value >=0):
            raise ValueError('Transition probability must be a positive number.')        
        self._probability = float(value)

    @property
    def type(self):
        return self._type


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
            source_statekind = csb.pyutils.Enum.parse(States, tran_type[0])
            target_statekind = csb.pyutils.Enum.parse(States, tran_type[3])
            if source_statekind == States.Match:
                try:
                    self.write("{0:<7}\t".format(
                        int(hmm.start.transitions[target_statekind].probability)))
                except csb.pyutils.ItemNotFoundError:
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
                source_statekind = csb.pyutils.Enum.parse(States, tran_type[0])
                target_statekind = csb.pyutils.Enum.parse(States, tran_type[3])

                if target_statekind == States.Match and layer.rank == hmm.layers.last_index:
                    target_statekind = States.End

                try:
                    state = layer[source_statekind]
                    self.write("{0:<7}\t".format(
                        int(state.transitions[target_statekind].probability)))
                except csb.pyutils.ItemNotFoundError:
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
    
                
class HHpredHitAlignment(sequence.SequenceAlignment):
    """
    Represents a query-template alignment in an HHpred result.
    
    @param hit: relevant hit object
    @type param: L{HHpredHit}
    @param query: the query sequence in the alignment region, with gaps
    @type query: str
    @param subject: the subject sequence in the alignment region, with gaps
    @type subject: str
    """

    GAP = sequence.ProteinAlphabet.GAP
    
    def __init__(self, hit, query, subject):
        
        if not isinstance(hit, HHpredHit):
            raise TypeError(hit)
        
        self._hit = hit
        
        q = sequence.Sequence('query', '', ''.join(query), type=sequence.SequenceTypes.Protein)
        s = sequence.Sequence(hit.id, '', ''.join(subject), type=sequence.SequenceTypes.Protein)
        
        super(HHpredHitAlignment, self).__init__((q, s))

    @property
    def query(self):
        return self.rows[1].sequence

    @property
    def subject(self):
        return self.rows[2].sequence
    
    @property
    def segments(self):
        """
        Find all ungapped query-subject segments in the alignment.
        Return a generator over all ungapped alignment segments, represented
        by L{HHpredHit} objects
        
        @rtype: generator
        """
        
        def make_segment(sstart, send, qstart, qend):
            
            seg = HHpredHit(self._hit.rank, self._hit.id, sstart, send, 
                            qstart, qend, self._hit.probability, self._hit.qlength)

            seg.slength = self._hit.slength
            seg.evalue = self._hit.evalue
            seg.pvalue = self._hit.pvalue
            seg.score = self._hit.score
            seg.ss_score = self._hit.ss_score
            seg.identity = self._hit.identity
            seg.similarity = self._hit.similarity
            seg.prob_sum = self._hit.prob_sum
            
            return seg
        
        in_segment = False
        qs = self._hit.qstart - 1
        ss = self._hit.start - 1 
        qi, si = qs, ss
        qe, se = qs, ss
        
        for q, s in zip(self.query, self.subject):

            if q != HHpredHitAlignment.GAP:
                qi += 1
            if s != HHpredHitAlignment.GAP:
                si += 1
                
            if HHpredHitAlignment.GAP in (q, s):
                if in_segment:
                    yield make_segment(ss, se, qs, qe)
                    in_segment = False
                    qs, ss = 0, 0
                    qe, se = 0, 0
            else: 
                if not in_segment:
                    in_segment = True
                    qs, ss = qi, si
            
            qe, se = qi, si                        
        
        if in_segment:
            yield make_segment(ss, se, qs, qe)
                    
    def to_a3m(self):
        """
        @return: a query-centric A3M alignment.
        @rtype: L{csb.bio.sequence.A3MAlignment}
        """
        a3m = self.format(sequence.AlignmentFormats.A3M)
        return sequence.A3MAlignment.parse(a3m, strict=False)

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

        # unmanaged fields, should be encapsulated in the future
        self.slength = None
        self.evalue = None
        self.pvalue = None
        self.score = None
        self.ss_score = None
        self.identity = None
        self.similarity = None
        self.prob_sum = None

    def __str__(self):
        return "{0.id} {0.probability} {0.start}-{0.end}".format(self)
    
    def __repr__(self):
        return "<HHpredHit: {0!s}>".format(self)
    
    def __cmp__(self, other):
        return cmp(self.rank, other.rank)

    def equals(self, other):
        """
        Return True if C{self} is completely identical to C{other} (same id, same start
        and end positions).
        
        @param other: right-hand-term
        @type other: HHpredHit
        
        @rtype: bool        
        """
        return (self.id == other.id and self.start == other.start and self.end == other.end)
        
    def surpasses(self, other):
        """
        Return True if C{self} is a superior to C{other} in terms of length 
        and probability. These criteria are applied in the following order:
        
            1. Length (the longer hit is better)
            2. Probability (if they have the same length, the one with the higher
               probability is better)
            3. Address (if they have the same length and probability, the one with
               higher memory ID wins; for purely practical reasons) 
        
        @param other: right-hand-term
        @type other: HHpredHit
        
        @rtype: bool        
        """
        if self.length > other.length:
            return True
        elif self.length == other.length:
            if self.probability > other.probability:
                return True      
            elif self.probability == other.probability:
                if id(self) > id(other):
                    return True
        return False
            
    def includes(self, other, tolerance=1):
        """
        Return True if C{other} overlaps with C{self}, that means C{other}
        is fully or partially included in C{self} when aligned over the query.
        
        @param other: right-hand-term
        @type other: HHpredHit
        @param tolerance: allow partial overlaps for that number of residues at
                          either end
        @type tolerance: int
        
        @rtype: bool 
        """
        if self.id == other.id:
            if other.start >= self.start:
                if (other.end - self.end) <= tolerance:
                    return True
            elif other.end <= self.end:
                if (self.start - other.start) <= tolerance:
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
            raise TypeError('probability must be float, not {0}'.format(type(value)))
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

    def __repr__(self):
        return "HHpredHitList\n\tquery={0.query_name}\n\tmatch_columns={0.match_columns}\n\tno_of_seqs={0.no_of_seqs}\n\tneff={0.neff}\n\tsearched_hmms={0.searched_hmms}\n\tdate={0.date}\n\tcommand={0.command}".format(self)

    __str__ = __repr__

    def __getitem__(self, index):
        return self._hits[index]

    def __iter__(self):
        return iter(self._hits)

    def __len__(self):
        return len(self._hits)

    def sort(self):
        self._hits.sort(key=lambda i: i.rank)

