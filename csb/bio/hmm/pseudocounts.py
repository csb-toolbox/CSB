import sys
import csb.core
import csb.bio.sequence as sequence

from csb.bio.hmm import States, ScoreUnits, Transition, State


GONNET = [10227, 3430, 2875, 3869, 1625, 2393, 4590, 6500, 2352, 3225, 5819, 4172, 1435,
          1579, 3728, 4610, 6264, 418, 1824, 5709, 3430, 7780, 2209, 2589, 584, 2369,
          3368, 3080, 2173, 1493, 3093, 5701, 763, 859, 1893, 2287, 3487, 444, 1338, 2356,
          2875, 2209, 3868, 3601, 501, 1541, 2956, 3325, 1951, 1065, 2012, 2879, 532, 688,
          1480, 2304, 3204, 219, 1148, 1759, 3869, 2589, 3601, 8618, 488, 2172, 6021, 4176,
          2184, 1139, 2151, 3616, 595, 670, 2086, 2828, 3843, 204, 1119, 2015, 1625, 584,
          501, 488, 5034, 355, 566, 900, 516, 741, 1336, 591, 337, 549, 419, 901,
          1197, 187, 664, 1373, 2393, 2369, 1541, 2172, 355, 1987, 2891, 1959, 1587, 1066,
          2260, 2751, 570, 628, 1415, 1595, 2323, 219, 871, 1682, 4590, 3368, 2956, 6021,
          566, 2891, 8201, 3758, 2418, 1624, 3140, 4704, 830, 852, 2418, 2923, 4159, 278, 1268, 2809,
          6500, 3080, 3325, 4176, 900, 1959, 3758, 26066, 2016, 1354, 2741, 3496, 741, 797, 2369,
          3863, 4169, 375, 1186, 2569, 2352, 2173, 1951, 2184, 516, 1587, 2418, 2016, 5409,
          1123, 2380, 2524, 600, 1259, 1298, 1642, 2446, 383, 876, 1691, 3225, 1493, 1065,
          1139, 741, 1066, 1624, 1354, 1123, 6417, 9630, 1858, 1975, 2225, 1260, 1558, 3131,
          417, 1697, 7504, 5819, 3093, 2012, 2151, 1336, 2260, 3140, 2741, 2380, 9630, 25113,
          3677, 4187, 5540, 2670, 2876, 5272, 1063, 3945, 11005, 4172, 5701, 2879, 3616, 591,
          2751, 4704, 3496, 2524, 1858, 3677, 7430, 949, 975, 2355, 2847, 4340, 333, 1451, 2932,
          1435, 763, 532, 595, 337, 570, 830, 741, 600, 1975, 4187, 949, 1300, 1111, 573,
          743, 1361, 218, 828, 2310, 1579, 859, 688, 670, 549, 628, 852, 797, 1259, 2225,
          5540, 975, 1111, 6126, 661, 856, 1498, 1000, 4464, 2602, 3728, 1893, 1480, 2086, 419,
          1415, 2418, 2369, 1298, 1260, 2670, 2355, 573, 661, 11834, 2320, 3300, 179, 876, 2179,
          4610, 2287, 2304, 2828, 901, 1595, 2923, 3863, 1642, 1558, 2876, 2847, 743, 856, 2320,
          3611, 4686, 272, 1188, 2695, 6264, 3487, 3204, 3843, 1197, 2323, 4159, 4169, 2446, 3131,
          5272, 4340, 1361, 1498, 3300, 4686, 8995, 397, 1812, 5172, 418, 444, 219, 204, 187,
          219, 278, 375, 383, 417, 1063, 333, 218, 1000, 179, 272, 397, 4101, 1266, 499,
          1824, 1338, 1148, 1119, 664, 871, 1268, 1186, 876, 1697, 3945, 1451, 828, 4464, 876,
          1188, 1812, 1266, 9380, 2227, 5709, 2356, 1759, 2015, 1373, 1682, 2809, 2569, 1691, 7504,
          11005, 2932, 2310, 2602, 2179, 2695, 5172, 499, 2227.0, 11569.0]
"""
Gonnet matrix frequencies taken from HHpred
"""

        
class PseudocountBuilder(object):
    """
    Constructs profile HMMs with pseudocounts. 
    """
    
    def __init__(self, hmm):
        self._hmm = hmm
        
    @property
    def hmm(self):
        return self._hmm
    
    def add_emission_pseudocounts(self, tau=0.1, pca=2.5, pcb=0.5, pcc=1.0):
        """
        Port from HHpred, it uses the conditional background probabilities,
        inferred from the Gonnet matrix.

        @param tau: admission weight, i.e how much of the final score is
                    determined by the background probabilities.
                    0.0=no pseudocounts. 
        @type tau: float
        """
        from numpy import array, dot, transpose, clip

        if self.hmm.pseudocounts or self.hmm.emission_pseudocounts:
            return
        if abs(tau) < 1e-6:
            return
        
        # Assume probabilities
        if not self.hmm.score_units == ScoreUnits.Probability:
            self.hmm.convert_scores(units=ScoreUnits.Probability)

        alphabet = csb.core.Enum.values(sequence.StdProteinAlphabet)
        
        ## S = SubstitutionMatrix(substitution_matrix)
        s_mat = array(GONNET)
        #Normalize
        s_mat /= s_mat.sum()
        s_mat = s_mat.reshape((len(alphabet), len(alphabet)))
        # Marginalize matrix
        s_marginal = s_mat.sum(-1)
        s_conditional = s_mat / s_marginal
        # Get data and info from hmm 
        em = array([ [layer[States.Match].emission[aa] or 0.0 for aa in alphabet]
                    for layer in self.hmm.layers])

        em = clip(em, sys.float_info.min, 1.)
        
        neff_m = array([l.effective_matches for l in self.hmm.layers])

        g = dot(em, transpose(s_conditional))

        if neff_m is not None:
            tau = clip(pca / (1. + (neff_m / pcb) ** pcc), 0.0, pcc)
            e = transpose((1. - tau) * transpose(em) + tau * transpose(g))
        else:
            e = (1. - tau) * em + tau * g
            
        # Renormalize e
        e = transpose(transpose(e) / e.sum(-1))

        for i, layer in enumerate(self.hmm.layers):
            layer[States.Match].emission.set(dict(zip(alphabet, e[i])))

        self.hmm.emission_pseudocounts = True
        return 
        

        
    def add_transition_pseudocounts(self, gapb=1., gapd=0.15, gape=1.0, gapf=0.6, gapg=0.6, gapi=0.6):
        """
        Add pseudocounts to the transitions. A port from hhsearch
        -gapb 1.0 -gapd 0.15 -gape 1.0 -gapf 0.6 -gapg 0.6 -gapi 0.6
        """

        from numpy import array

        if not self.hmm._score_units == ScoreUnits.Probability:
            self.hmm.convert_scores(units=ScoreUnits.Probability) 

        if self.hmm.pseudocounts or self.hmm.transition_pseudocounts:
            return

        # We need a fully populated HMM so first add all missing states
        states = [States.Match, States.Insertion, States.Deletion] 
        background = self.hmm.layers[1][States.Match].background
        for layer in self.hmm.layers:
            rank = layer.rank
            for state in states:
                if state not in layer:

                    if state is States.Deletion:
                        # Add a new Deletion state
                        deletion = State(States.Deletion)
                        deletion.rank = rank 
                        layer.append(deletion)
                        
                    elif state is States.Insertion:
                        # Add a new Deletion state
                        insertion = State(States.Insertion,
                                          emit=csb.core.Enum.members(
                                          sequence.SequenceAlphabets.Protein))
                        insertion.background.set(background) 
                        insertion.emission.set(background)
                        insertion.rank = rank
                        layer.append(insertion)

        if not self.hmm.start_insertion:
            insertion = State(States.Insertion,
                                          emit=csb.core.Enum.members(
                                          sequence.SequenceAlphabets.Protein))
            insertion.background.set(background) 
            insertion.emission.set(background)
            insertion.rank = 0
            self.hmm.start_insertion = insertion

        # make hmm completly connected
        for i in range(1, self.hmm.layers.length):
            layer = self.hmm.layers[i]
            #Start with match state
            state = layer[States.Match]
            if not States.Insertion in state.transitions:
                state.transitions.append(Transition(state,
                                                    self.hmm.layers[i][States.Insertion],
                                                    0.0))
            if not States.Deletion in state.transitions:
                state.transitions.append(Transition(state,
                                                    self.hmm.layers[i + 1][States.Deletion],
                                                    0.0))
            state = layer[States.Insertion]
            if not States.Insertion in state.transitions:
                state.transitions.append(Transition(state,
                                                    self.hmm.layers[i][States.Insertion],
                                                    0.0))
            if not States.Match in state.transitions:
                state.transitions.append(Transition(state,
                                                    self.hmm.layers[i + 1][States.Match],
                                                    0.0))
            state = layer[States.Deletion]
            if not States.Deletion in state.transitions:
                state.transitions.append(Transition(state,
                                                    self.hmm.layers[i + 1][States.Deletion],
                                                    0.0))
            if not States.Match in state.transitions:
                state.transitions.append(Transition(state,
                                                    self.hmm.layers[i + 1][States.Match],
                                                    0.0))
        # start layer
        state = self.hmm.start
        if not States.Insertion in self.hmm.start.transitions:
            state.transitions.append(Transition(self.hmm.start,
                                                self.hmm.start_insertion,
                                                0.0))
        if not States.Deletion in self.hmm.start.transitions:
            state.transitions.append(Transition(self.hmm.start,
                                                self.hmm.layers[1][States.Deletion],
                                                0.0))

        state = self.hmm.start_insertion
        if not States.Insertion in self.hmm.start_insertion.transitions:
            state.transitions.append(Transition(self.hmm.start_insertion,
                                                self.hmm.start_insertion,
                                                0.0))
        if not States.Match in self.hmm.start_insertion.transitions:
            state.transitions.append(Transition(self.hmm.start_insertion,
                                                self.hmm.layers[1][States.Match],
                                                0.0))

        # last layer
        state = self.hmm.layers[-1][States.Match]
        if not States.Insertion in state.transitions:
            state.transitions.append(Transition(state,
                                                self.hmm.layers[-1][States.Insertion],
                                                0.0))
        state = self.hmm.layers[-1][States.Insertion]
        if not States.Insertion in state.transitions:
            state.transitions.append(Transition(state,
                                                self.hmm.layers[-1][States.Insertion],
                                                0.0))

        if not States.End in state.transitions:
            state.transitions.append(Transition(state,
                                                self.hmm.end,
                                                0.0))
        state = self.hmm.layers[-1][States.Deletion]
        if not States.End in state.transitions:
            state.transitions.append(Transition(state,
                                                self.hmm.end,
                                                0.0))
        

        
        # Now we have created a fully connected HMM
        # Lates add pseuod counts
        # Calculate pseudo counts

        # to be honest I really do not know how they came up with this
        pc_MD = pc_MI = 0.0286 * gapd
        pc_MM = 1. - 2 * pc_MD
        pc_DD = pc_II = gape / (gape - 1 + 1 / 0.75)
        pc_DM = pc_IM = 1. - pc_II

        
        # Get current transtion probabilities
        t_mm = self.hmm.start.transitions[States.Match].probability 
        t_mi = self.hmm.start.transitions[States.Insertion].probability 
        t_md = self.hmm.start.transitions[States.Deletion].probability 

        # Transitions from Match state
        n_eff = self.hmm.effective_matches
        
        t = array([(n_eff - 1) * t_mm + gapb * pc_MM,
                   (n_eff - 1) * t_mi + gapb * pc_MI,
                   (n_eff - 1) * t_md + gapb * pc_MD])
        # normalize to one
        t /= t.sum()
        # Set 
        self.hmm.start.transitions[States.Match].probability = t[0]
        self.hmm.start.transitions[States.Insertion].probability = t[1]
        self.hmm.start.transitions[States.Deletion].probability = t[2]
        
        # Rinse and repeat
        t_im = self.hmm.start_insertion.transitions[States.Match].probability 
        t_ii = self.hmm.start_insertion.transitions[States.Insertion].probability  
        
        t = array([t_im + gapb * pc_IM, t_ii + gapb * pc_II])
        t /= t.sum()

        self.hmm.start_insertion.transitions[States.Match].probability = t[0]
        t_ii = self.hmm.start_insertion.transitions[States.Insertion].probability = t[1]

        # And now for all layers
        for layer in self.hmm.layers[:-1]:
            # Get current transtion probabilities
            t_mm = layer[States.Match].transitions[States.Match].probability 
            t_mi = layer[States.Match].transitions[States.Insertion].probability 
            t_md = layer[States.Match].transitions[States.Deletion].probability 
            n_eff = layer.effective_matches
            t = array([(n_eff - 1) * t_mm + gapb * pc_MM,
                       (n_eff - 1) * t_mi + gapb * pc_MI,
                       (n_eff - 1) * t_md + gapb * pc_MD])
            # normalize to one
            t /= t.sum()
            layer[States.Match].transitions[States.Match].probability = t[0]
            layer[States.Match].transitions[States.Insertion].probability = t[1]
            layer[States.Match].transitions[States.Deletion].probability = t[2]
            
            # Transitions from insert state
            t_im = layer[States.Insertion].transitions[States.Match].probability 
            t_ii = layer[States.Insertion].transitions[States.Insertion].probability
            n_eff = layer.effective_insertions
            t = array([t_im * n_eff + gapb * pc_IM,
                       t_im * n_eff + gapb * pc_II])
            # normalize to one
            t /= t.sum()
            layer[States.Insertion].transitions[States.Match].probability = t[0]
            layer[States.Insertion].transitions[States.Insertion].probability = t[1]

            # Transitions form deletion state
            t_dm = layer[States.Deletion].transitions[States.Match].probability 
            t_dd = layer[States.Deletion].transitions[States.Deletion].probability
            n_eff = layer.effective_deletions
            t = array([t_dm * n_eff + gapb * pc_DM,
                       t_dd * n_eff + gapb * pc_DD])
            # normalize to one
            t /= t.sum()
            layer[States.Deletion].transitions[States.Match].probability = t[0]
            layer[States.Deletion].transitions[States.Deletion].probability = t[1]

        #Last layer

        layer = self.hmm.layers[-1]
        t_mm = layer[States.Match].transitions[States.End].probability 
        t_mi = layer[States.Match].transitions[States.Insertion].probability 
        n_eff = layer.effective_matches
        # No deletion
        t = array([(n_eff - 1) * t_mm + gapb * pc_MM,
                   (n_eff - 1) * t_mi + gapb * pc_MI])
        # normalize to one
        t /= t.sum()
        layer[States.Match].transitions[States.End].probability = t[0]
        layer[States.Match].transitions[States.Insertion].probability = t[1]
        
        # Transitions from insert state
        t_im = layer[States.Insertion].transitions[States.End].probability 
        t_ii = layer[States.Insertion].transitions[States.Insertion].probability
        n_eff = layer.effective_insertions
        t = array([t_im * n_eff + gapb * pc_IM,
                   t_im * n_eff + gapb * pc_II])
        # normalize to one
        t /= t.sum()
        layer[States.Insertion].transitions[States.End].probability = t[0]
        layer[States.Insertion].transitions[States.Insertion].probability = t[1]

        layer[States.Deletion].transitions[States.End].probability = 1.

        self.hmm.transition_pseudocounts = True
        return 
    
