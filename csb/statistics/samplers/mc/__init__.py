"""
Abstract Monte Carlo samplers.
"""

import numpy.random

import csb.numeric
import csb.core

from abc import ABCMeta, abstractmethod, abstractproperty
from csb.statistics.samplers import AbstractSampler, AbstractState, State, EnsembleState


class AbstractMC(AbstractSampler):
    """
    Abstract Monte Carlo sampler class. Subclasses implement various
    Monte carlo equilibrium sampling schemes.
    
    @param state: Initial state
    @type state: L{AbstractState}
    """
    
    __metaclass__ = ABCMeta
    
    def __init__(self, state):
        
        self._state = None
        self.state = state
         
    def _checkstate(self, state):
        
        if not isinstance(state, AbstractState):
            raise TypeError(state)
    
    @abstractproperty
    def energy(self):
        """
        Energy of the current state.
        """
        pass

    @property
    def state(self):
        """
        Current state.
        """
        return self._state
    @state.setter
    def state(self, value):
        self._checkstate(value)
        self._state = value

    @abstractmethod
    def sample(self):
        """
        Draw a sample.
        @rtype: L{AbstractState}
        """
        pass

class AbstractProposalCommunicator(object):
    """
    With the exception of the current state of the Markov chain, this
    holds all the information needed to calculate the acceptance
    probability in an AbstractSingleChainMC-derived sampling
    algorithm.

    @param proposal: Proposal state
    @type proposal: L{State}
    """

    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self, proposal):

        self._proposal = proposal

    @property
    def proposal(self):
        return self._proposal

class RWMCProposalCommunicator(AbstractProposalCommunicator):
    """
    With the exception of the current state of the Markov chain, this
    holds all the information needed to calculate the acceptance
    probability of a proposal state in the random walk Metropolis
    Monte Carlo algorithm, that is, only the proposal state.

    @param proposal: Proposal state
    @type proposal: L{State}
    """

    def __init__(self, proposal):

        super(RWMCProposalCommunicator, self).__init__(proposal)

class HMCProposalCommunicator(AbstractProposalCommunicator):
    """
    With the exception of the current state of the Markov chain, this
    holds all the information needed to calculate the acceptance
    probability of a proposal state in the HMC algorithm, that is,
    only the proposal state.

    @param proposal: Proposal state
    @type proposal: L{State}
    """

    def __init__(self, proposal):

        super(HMCProposalCommunicator, self).__init__(proposal)


class AbstractSingleChainMC(AbstractMC):
    """
    Abstract class for Monte Carlo sampling algorithms simulating
    only one ensemble.
    
    @param pdf: probability density function to sample from
    @type pdf: subclass of L{csb.statistics.pdf.AbstractDensity}

    @param state: Initial state
    @type state: L{State}
    """
    
    __metaclass__ = ABCMeta

    def __init__(self, pdf, state):
        
        super(AbstractSingleChainMC, self).__init__(state)
        
        self._pdf = pdf
        self._nmoves = 0
        self._accepted = 0
        self._last_move_accepted = None
        
    def _checkstate(self, state):
        if not isinstance(state, State):
            raise TypeError(state)

    def sample(self):
        """
        Draw a sample.
        @rtype: L{State}
        """
        
        proposal_communicator = self._propose()
        pacc = self._calc_pacc(proposal_communicator)
        self._accept(proposal_communicator.proposal, pacc)

        return self.state

    @abstractmethod
    def _propose(self):
        """
        Calculate a new proposal state and gather additional information
        needed to calculate the acceptance probability.
        
        @rtype: L{AbstractProposalCommunicator}
        """
        pass

    @abstractmethod
    def _calc_pacc(self, proposal_communicator):
        """
        Calculate probability with which to accept the proposal.

        @param proposal_communicator: Contains information about the proposal
                                      and additional information needed to
                                      calculate the acceptance probability
        @type proposal_communicator: L{AbstractProposalCommunicator}
        """
        pass

    def _accept(self, proposal, pacc):
        """
        Accept / reject proposal with given acceptance probability pacc.

        @param proposal: proposal state
        @type proposal: L{State}

        @param pacc: acceptance probability
        @type pacc: float
        """
        
        self._nmoves += 1
        
        if numpy.random.random() < pacc:
            self.state = proposal
            self._accepted += 1
            self._last_move_accepted = True
            return True
        else:
            self._last_move_accepted = False
            return False

    @property
    def energy(self):
        """
        Log-likelihood of the current state.
        @rtype: float
        """
        return self._pdf.log_prob(self.state.position)

    @property
    def acceptance_rate(self):
        """
        Acceptance rate.
        """
        return float(self._accepted) / float(self._nmoves)

    @property
    def last_move_accepted(self):
        """
        Information whether the last MC move was accepted or not.
        """
        return self._last_move_accepted

class MCCollection(csb.core.BaseCollectionContainer):
    """
    Collection of single-chain samplers.

    @param items: samplers
    @type items: list of L{AbstractSingleChainMC}
    """
    
    def __init__(self, items):
        
        super(MCCollection, self).__init__(items, type=AbstractSingleChainMC)

class AbstractEnsembleMC(AbstractMC):
    """
    Abstract class for Monte Carlo sampling algorithms simulating several ensembles.

    @param samplers: samplers which sample from their respective equilibrium distributions
    @type samplers: list of L{AbstractSingleChainMC}    
    """

    __metaclass__ = ABCMeta

    def __init__(self, samplers):
        
        self._samplers = MCCollection(samplers)
        state = EnsembleState([x.state for x in self._samplers])
        
        super(AbstractEnsembleMC, self).__init__(state)        

    def sample(self):
        """
        Draw an ensemble sample.
        
        @rtype: L{EnsembleState}
        """
        
        sample = EnsembleState([sampler.sample() for sampler in self._samplers])
        self.state = sample

        return sample

    @property
    def energy(self):
        """
        Total ensemble energy.
        """        
        return sum([x.energy for x in self._samplers])

class AbstractSwapParameterInfo(object):
    """
    Subclass instances hold all parameters necessary for performing a swap
    between two given samplers.
    """

    __metaclass__ = ABCMeta
    
    def __init__(self, sampler1, sampler2):
        """
        @param sampler1: First sampler
        @type sampler1: L{AbstractSingleChainMC}

        @param sampler2: Second sampler
        @type sampler2: L{AbstractSingleChainMC}
        """
        
        self._sampler1 = sampler1
        self._sampler2 = sampler2

    @property
    def sampler1(self):        
        return self._sampler1

    @property
    def sampler2(self):        
        return self._sampler2        

class RESwapParameterInfo(AbstractSwapParameterInfo):
    """
    Holds parameters for a standard Replica Exchange swap.
    """
    pass

class MDRENSSwapParameterInfo(RESwapParameterInfo):
    """
    Holds parameters for a MDRENS swap.

    @param sampler1: First sampler
    @type sampler1: L{AbstractSingleChainMC}

    @param sampler2: Second sampler
    @type sampler2: L{AbstractSingleChainMC}

    @param timestep: Integration timestep
    @type timestep: float

    @param traj_length: Trajectory length in number of timesteps
    @type traj_length: int

    @param gradient: Gradient which determines the dynamics during a trajectory
    @type gradient: L{AbstractGradient}       
    """
    
    def __init__(self, sampler1, sampler2, timestep, traj_length, gradient):
        
        super(MDRENSSwapParameterInfo, self).__init__(sampler1, sampler2)
        
        self._timestep = timestep
        self._traj_length = traj_length
        self._gradient = gradient

    @property
    def timestep(self):
        """
        Integration timestep.
        """
        return self._timestep
    @timestep.setter
    def timestep(self, value):
        self._timestep = float(value)

    @property
    def traj_length(self):
        """
        Trajectory length in number of integration steps.
        """
        return self._traj_length
    @traj_length.setter
    def traj_length(self, value):
        self._traj_length = int(value)

    @property
    def gradient(self):
        """
        Gradient which governs the equations of motion.
        """
        return self._gradient

class ThermostattedMDRENSSwapParameterInfo(MDRENSSwapParameterInfo):
    """
    @param sampler1: First sampler
    @type sampler1: subclass instance of L{AbstractSingleChainMC}

    @param sampler2: Second sampler
    @type sampler2: subclass instance of L{AbstractSingleChainMC}

    @param timestep: Integration timestep
    @type timestep: float

    @param traj_length: Trajectory length in number of timesteps
    @type traj_length: int

    @param gradient: Gradient which determines the dynamics during a trajectory
    @type gradient: subclass instance of L{AbstractGradient}

    @param temperature: Temperature interpolation function.
    @type temperature: Real-valued function mapping from [0,1] to R.
        T(0) = temperature of the ensemble sampler1 samples from, T(1) = temperature
        of the ensemble sampler2 samples from

    @param collision_probability: Probability for a collision with the heatbath during one timestep
    @type collision_probability: float

    @param collision_interval: Interval during which collision may occur with probability
        collision_probability
    @type collision_interval: int
    """
        
    def __init__(self, sampler1, sampler2, timestep, traj_length, gradient,
                 temperature=lambda l: 1., collision_probability=0.1, collision_interval=1):
        
        super(ThermostattedMDRENSSwapParameterInfo, self).__init__(
                                    sampler1, sampler2, timestep, traj_length, gradient)
        
        self._collision_probability = None
        self._collision_interval = None
        self._temperature = temperature
        
        self.collision_probability = collision_probability
        self.collision_interval = collision_interval

    @property
    def collision_probability(self):
        """
        Probability for a collision with the heatbath during one timestep.
        """
        return self._collision_probability
    @collision_probability.setter
    def collision_probability(self, value):
        self._collision_probability = float(value)

    @property
    def collision_interval(self):
        """
        Interval during which collision may occur with probability
        C{collision_probability}.
        """
        return self._collision_interval
    @collision_interval.setter
    def collision_interval(self, value):
        self._collision_interval = int(value)

    @property
    def temperature(self):
        return self._temperature

class AbstractSwapCommunicator(object):
    """
    Holds all the information which needs to be communicated between
    distinct swap substeps.

    @param param_info: ParameterInfo instance holding swap parameters
    @type param_info: L{AbstractSwapParameterInfo}

    @param proposal1: Proposal state for first sampler
    @type proposal1: L{State}

    @param proposal2: Proposal state for second sampler
    @type proposal2: L{State}
    """

    __metaclass__ = ABCMeta
    
    def __init__(self, param_info, proposal1, proposal2):
        
        self._sampler1 = param_info.sampler1
        self._sampler2 = param_info.sampler2
        
        self._state1 = self.sampler1.state
        self._state2 = self.sampler2.state
        
        self._proposal1 = proposal1
        self._proposal2 = proposal2
        
        self._acceptance_probability = None
        self._accepted = False
        
    @property
    def sampler1(self):
        return self._sampler1

    @property
    def sampler2(self):
        return self._sampler2
    
    @property
    def state1(self):
        return self._state1    

    @property
    def state2(self):
        return self._state2
        
    @property
    def proposal1(self):
        return self._proposal1
    
    @property
    def proposal2(self):
        return self._proposal2

    @property
    def acceptance_probability(self):
        return self._acceptance_probability
    @acceptance_probability.setter
    def acceptance_probability(self, value):
        self._acceptance_probability = value

    @property
    def accepted(self):
        return self._accepted
    @accepted.setter
    def accepted(self, value):
        self._accepted = value

class RESwapCommunicator(AbstractSwapCommunicator):
    """
    Holds all the information which needs to be communicated between distinct
    RE swap substeps.

    See L{AbstractSwapCommunicator} for constructor signature.
    """
    pass

class RENSSwapCommunicator(AbstractSwapCommunicator):
    """
    Holds all the information which needs to be communicated between distinct
    RENS swap substeps.
    
    @param param_info: ParameterInfo instance holding swap parameters
    @type param_info: L{AbstractSwapParameterInfo}

    @param proposal1: Proposal state for first sampler
    @type proposal1: L{State}

    @param proposal2: Proposal state for second sampler
    @type proposal2: L{State}

    @param heat12: Heat generated during the forward trajectory
    @type heat12: float

    @param heat21: Heat generated during the reverse trajectory
    @type heat21: float
    """
        
    def __init__(self, param_info, proposal1, proposal2, heat12, heat21):

        super(RENSSwapCommunicator, self).__init__(param_info, proposal1, proposal2)
        self._heat12 = heat12
        self._heat21 = heat21

    @property
    def heat12(self):
        return self._heat12

    @property
    def heat21(self):
        return self._heat21

class SingleSwapStatistics(object):
    """
    Tracks swap statistics of a single sampler pair.

    @param param_info: ParameterInfo instance holding swap parameters
    @type param_info: L{AbstractSwapParameterInfo}
    """
    
    def __init__(self, param_info):
        self._total_swaps = 0
        self._accepted_swaps = 0

    @property
    def total_swaps(self):
        return self._total_swaps
    
    @property
    def accepted_swaps(self):
        return self._accepted_swaps
    
    @property
    def acceptance_rate(self):
        """
        Acceptance rate of the sampler pair.
        """
        if self.total_swaps > 0:
            return float(self.accepted_swaps) / float(self.total_swaps)
        else:
            return 0.

    def update(self, accepted):
        """
        Updates swap statistics.
        """        
        self._total_swaps += 1
        self._accepted_swaps += int(accepted)

class SwapStatistics(object):
    """
    Tracks swap statistics for an AbstractExchangeMC subclass instance.

    @param param_infos: list of ParameterInfo instances providing information
                        needed for performing swaps
    @type param_infos: list of L{AbstractSwapParameterInfo}
    """
    
    def __init__(self, param_infos):        
        self._stats = [SingleSwapStatistics(x) for x in param_infos]
        
    @property
    def stats(self):
        return tuple(self._stats)

    @property
    def acceptance_rates(self):
        """
        Returns acceptance rates for all swaps.
        """        
        return [x.acceptance_rate for x in self._stats]
        
class AbstractExchangeMC(AbstractEnsembleMC):
    """
    Abstract class for Monte Carlo sampling algorithms employing some replica exchange method.

    @param samplers: samplers which sample from their respective equilibrium distributions
    @type samplers: list of L{AbstractSingleChainMC}

    @param param_infos: list of ParameterInfo instances providing information needed
        for performing swaps
    @type param_infos: list of L{AbstractSwapParameterInfo}
    """

    __metaclass__ = ABCMeta

    def __init__(self, samplers, param_infos):
        super(AbstractExchangeMC, self).__init__(samplers)
        
        self._swaplist1 = []
        self._swaplist2 = []        
        self._currentswaplist = self._swaplist1
        
        self._param_infos = param_infos
        self._statistics = SwapStatistics(self._param_infos)
        
    def _checkstate(self, state):
        if not isinstance(state, EnsembleState):
            raise TypeError(state)        

    def swap(self, index):
        """
        Perform swap between sampler pair described by param_infos[index]
        and return outcome (true = accepted, false = rejected).

        @param index: index of swap pair in param_infos
        @type index: int

        @rtype: boolean
        """
        param_info = self._param_infos[index]
        swapcom = self._propose_swap(param_info)
        self._calc_pacc_swap(swapcom)
        result = self._accept_swap(swapcom)
        
        self.state = EnsembleState([x.state for x in self._samplers])
        self.statistics.stats[index].update(result)
        
        return result

    @abstractmethod
    def _propose_swap(self, param_info):
        """
        Calculate proposal states for a swap between two samplers.
        
        @param param_info: ParameterInfo instance holding swap parameters
        @type param_info: L{AbstractSwapParameterInfo}
        
        @rtype: L{AbstractSwapCommunicator}
        """ 
        pass

    @abstractmethod
    def _calc_pacc_swap(self, swapcom):
        """
        Calculate probability to accept a swap given initial and proposal states.

        @param swapcom: SwapCommunicator instance holding information to be communicated
                        between distinct swap substeps
        @type swapcom: L{AbstractSwapCommunicator} 
        """
        pass

    def _accept_swap(self, swapcom):
        """
        Accept / reject an exchange between two samplers given proposal states and
        the acceptance probability and returns the outcome (true = accepted, false = rejected).

        @param swapcom: SwapCommunicator instance holding information to be communicated
            between distinct swap substeps
        @type swapcom: L{AbstractSwapCommunicator}

        @rtype: boolean
        """
        
        if numpy.random.random() < swapcom.acceptance_probability:
            swapcom.sampler1.state = swapcom.proposal1
            swapcom.sampler2.state = swapcom.proposal2
            return True
        else:
            return False

    @property
    def acceptance_rates(self):
        """
        Return swap acceptance rates.

        @rtype: list of floats
        """        
        return self.statistics.acceptance_rates

    @property
    def param_infos(self):
        """
        List of SwapParameterInfo instances holding all necessary parameters.

        @rtype: list of L{AbstractSwapParameterInfo}
        """
        return self._param_infos
    
    @property
    def statistics(self):
        return self._statistics

    def _update_statistics(self, index, accepted):
        """
        Update statistics of a given swap process.
        
        @param index: position of swap statistics to be updated
        @type index: int
        
        @param accepted: outcome of the swap
        @type accepted: boolean
        """
        
        self._stats[index][0] += 1
        self._stats[index][1] += int(accepted)

class RENSTrajInfo(object):
    """
    Holds information necessary for calculating a RENS trajectory.

    @param param_info: ParameterInfo instance holding swap parameters
    @type param_info: L{AbstractSwapParameterInfo}

    @param init_state: state from which the trajectory is supposed to start
    @type init_state: L{State}

    @param protocol: Protocol to be used to generate nonequilibrium trajectories
    @type protocol: Real-valued function that maps [0, switching time] to [0, 1]      
    """
    
    def __init__(self, param_info, init_state, protocol):
        
        self._param_info = param_info
        self._protocol = protocol
        self._init_state = init_state
        
    @property
    def param_info(self):
        return self._param_info

    @property
    def protocol(self):
        return self._protocol

    @property
    def init_state(self):
        return self._init_state

class AbstractRENS(AbstractExchangeMC):
    """
    Abstract Replica Exchange with Nonequilibrium Switches
    (RENS, Ballard & Jarzynski 2009) class.
    Subclasses implement various ways of generating trajectories
    (deterministic or stochastic).
    """

    __metaclass__ = ABCMeta

    def _propose_swap(self, param_info):
        
        init_state1 = State(param_info.sampler1.state.position,
                            numpy.random.normal(size=param_info.sampler1.state.position.shape))
        init_state2 = State(param_info.sampler2.state.position,
                            numpy.random.normal(size=param_info.sampler2.state.position.shape))
        
        param_info.sampler1.state = init_state1
        param_info.sampler2.state = init_state2
        
        trajinfo12 = RENSTrajInfo(param_info, init_state1, protocol=lambda t, tau: t / tau)
        trajinfo21 = RENSTrajInfo(param_info, init_state2, protocol=lambda t, tau: (tau - t) / tau)
        
        traj12 = self._run_traj_generator(trajinfo12)
        traj21 = self._run_traj_generator(trajinfo21)

        return RENSSwapCommunicator(param_info, traj21.final, traj12.final, traj21.heat, traj12.heat)

    def _calc_pacc_swap(self, swapcom):
        
        heat12 = swapcom.heat12
        heat21 = swapcom.heat21
        
        proposal1 = swapcom.proposal1
        proposal2 = swapcom.proposal2
        
        state1 = swapcom.state1
        state2 = swapcom.state2
        
        E1 = lambda x:-swapcom.sampler1._pdf.log_prob(x)
        E2 = lambda x:-swapcom.sampler2._pdf.log_prob(x)
        
        w12 = 0.5 * sum(proposal2.momentum ** 2) + E2(proposal2.position) - \
              0.5 * sum(state1.momentum ** 2) - E1(state1.position) - heat12
        w21 = 0.5 * sum(proposal1.momentum ** 2) + E1(proposal1.position) - \
              0.5 * sum(state2.momentum ** 2) - E2(state2.position) - heat21
        
        swapcom.acceptance_probability = csb.numeric.exp(-w12 - w21)

    @abstractmethod
    def _run_traj_generator(self, traj_info):
        """
        Run the trajectory generator which generates a trajectory
        of a given length between the states of two samplers.

        @param traj_info: TrajectoryInfo instance holding information
                          needed to generate a nonequilibrium trajectory   
        @type traj_info: L{RENSTrajInfo}
        
        @rtype: L{Trajectory}
        """
        pass

class AbstractSwapScheme(object):
    """
    Provides the interface for classes defining schemes according to which swaps in
    Replica Exchange-like simulations are performed.

    @param algorithm: Exchange algorithm that performs the swaps
    @type algorithm: L{AbstractExchangeMC}
    """

    __metaclass__ = ABCMeta

    def __init__(self, algorithm):

        self._algorithm = algorithm

    @abstractmethod
    def swap_all(self):
        """
        Advises the Replica Exchange-like algorithm to perform swaps according to
        the some schedule defined here.
        """
        
        pass

class AlternatingAdjacentSwapScheme(AbstractSwapScheme):
    """
    Provides a swapping scheme in which tries exchanges between neighbours only
    following the scheme 1 <-> 2, 3 <-> 4, ... and after a sampling period 2 <-> 3, 4 <-> 5, ...

    @param algorithm: Exchange algorithm that performs the swaps
    @type algorithm: L{AbstractExchangeMC}
    """

    def __init__(self, algorithm):

        super(AlternatingAdjacentSwapScheme, self).__init__(algorithm)
        
        self._current_swap_list = None
        self._swap_list1 = []
        self._swap_list2 = []
        self._create_swap_lists()
    
    def _create_swap_lists(self):
        
        i = 0
        while i < len(self._algorithm.param_infos):
            self._swap_list1.append(i)
            i += 2

        i = 1
        while i < len(self._algorithm.param_infos):
            self._swap_list2.append(i)
            i += 2

        self._current_swap_list = self._swap_list1

    def swap_all(self):
        
        for x in self._current_swap_list:
            self._algorithm.swap(x)

        if self._current_swap_list == self._swap_list1:
            self._current_swap_list = self._swap_list2
        else:
            self._current_swap_list = self._swap_list1

