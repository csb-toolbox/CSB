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

class AbstractPropagationResult(object):
    """
    Abstract class providing the interface for the result
    of a deterministic or stochastic propagation of a state.
    """
    
    __metaclass__ = ABCMeta 
    
    @abstractproperty
    def initial(self):
        """
        Initial state
        """
        pass
    
    @abstractproperty
    def final(self):
        """
        Final state
        """
        pass
    
    @abstractproperty
    def heat(self):
        """
        Heat produced during propagation
        @rtype: float
        """        
        pass    

class PropagationResult(AbstractPropagationResult):
    """
    Describes the result of a deterministic or stochastic
    propagation of a state.

    @param initial: Initial state from which the
                    propagation started
    @type initial: L{State}

    @param final: Final state in which the propagation
                  resulted
    @type final: L{State}

    @param heat: Heat produced during propagation
    @type heat: float
    """
    
    
    def __init__(self, initial, final, heat=0.0):
        
        if not isinstance(initial, AbstractState):
            raise TypeError(initial)
        
        if not isinstance(final, AbstractState):
            raise TypeError(final)        
        
        self._initial = initial
        self._final = final
        self._heat = None
        
        self.heat = heat
        
    @property
    def initial(self):
        return self._initial
    
    @property
    def final(self):
        return self._final
    
    @property
    def heat(self):
        return self._heat
    @heat.setter
    def heat(self, value):
        self._heat = float(value)

class Trajectory(csb.core.CollectionContainer, AbstractPropagationResult):
    """
    Ordered collection of states, representing a phase-space trajectory.

    @param items: list of states defining a phase-space trajectory
    @type items: list of L{AbstractState}
    @param heat: heat produced during the trajectory
    @type heat: float
    @param work: work produced during the trajectory
    @type work: float
    """
    
    def __init__(self, items, heat=0.0, work=0.0):
        
        super(Trajectory, self).__init__(items, type=AbstractState)
        
        self._heat = heat    
        self._work = work
    
    @property
    def initial(self):
        return self[0]
    
    @property
    def final(self):
        return self[self.last_index]
    
    @property
    def heat(self):
        return self._heat
    @heat.setter
    def heat(self, value):
        self._heat = float(value)

    @property
    def work(self):
        return self._work
    @work.setter
    def work(self, value):
        self._work = float(value)

class TrajectoryBuilder(object):
    """
    Allows to  build a Trajectory object step by step.

    @param heat: heat produced over the trajectory
    @type heat: float
    @param work: work produced during the trajectory
    @type work: float
    """
    
    def __init__(self, heat=0.0, work=0.0):
        self._heat = heat
        self._work = work
        self._states = []
        
    @staticmethod
    def create(full=True):
        """
        Trajectory builder factory.

        @param full: if True, a TrajectoryBuilder instance designed
                     to build a full trajectory with initial state,
                     intermediate states and a final state. If False,
                     a ShortTrajectoryBuilder instance designed to
                     hold only the initial and the final state is
                     returned
        @type full: boolean
        """
        
        if full:
            return TrajectoryBuilder()
        else:
            return ShortTrajectoryBuilder()
        
    @property
    def product(self):
        """
        The L{Trajectory} instance build by a specific instance of
        this class
        """
        return Trajectory(self._states, heat=self._heat, work=self._work)

    def add_initial_state(self, state):
        """
        Inserts a state at the beginning of the trajectory

        @param state: state to be added
        @type state: L{State}
        """
        self._states.insert(0, state.clone())
        
    def add_intermediate_state(self, state):
        """
        Adds a state to the end of the trajectory

        @param state: state to be added
        @type state: L{State}
        """
        self._states.append(state.clone())
    
    def add_final_state(self, state):
        """
        Adds a state to the end of the trajectory

        @param state: state to be added
        @type state: L{State}
        """
        self._states.append(state.clone())
    
class ShortTrajectoryBuilder(TrajectoryBuilder):    

    def add_intermediate_state(self, state):
        pass

    @property
    def product(self):
        """
        The L{PropagationResult} instance built by a specific instance of
        this class
        """
        
        if len(self._states) != 2:
            raise ValueError("Can't create a product, two states required")
        
        initial, final = self._states
        return PropagationResult(initial, final, heat=self._heat)

class SimpleProposalCommunicator(object):
    """
    With the exception of the current state of the Markov chain, this
    holds all the information needed to calculate the acceptance
    probability in both the L{RWMCSampler} and L{HMCSampler} classes,
    that is, only the proposal state.
    For more advanced algorithms, one may derive classes capable of
    holding the neccessary additional information from this class.

    @param proposal: Proposal state
    @type proposal: L{State}
    """

    __metaclass__ = ABCMeta

    def __init__(self, proposal):

        self._proposal = proposal

    @property
    def proposal(self):
        return self._proposal

class AbstractSingleChainMC(AbstractMC):
    """
    Abstract class for Monte Carlo sampling algorithms simulating
    only one ensemble.
    
    @param pdf: probability density function to sample from
    @type pdf: subclass of L{csb.statistics.pdf.AbstractDensity}

    @param state: Initial state
    @type state: L{State}

    @param temperature: Pseudo-temperature of the Boltzmann ensemble
                        M{p(x) = 1/N * exp(-1/T * E(x))} with the
                        pseudo-energy defined as M{E(x) = -log(p(x))}
                        where M{p(x)} is the PDF under consideration
    @type temperature: float
    """
    
    __metaclass__ = ABCMeta

    def __init__(self, pdf, state, temperature=1.):
        
        super(AbstractSingleChainMC, self).__init__(state)
        
        self._pdf = pdf
        self._temperature = temperature
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
        
        @rtype: L{SimpleProposalCommunicator}
        """
        pass

    @abstractmethod
    def _calc_pacc(self, proposal_communicator):
        """
        Calculate probability with which to accept the proposal.

        @param proposal_communicator: Contains information about the proposal
                                      and additional information needed to
                                      calculate the acceptance probability
        @type proposal_communicator: L{SimpleProposalCommunicator}
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

    @property
    def temperature(self):
        return self._temperature

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

    @param mass_matrix: Mass matrix
    @type mass_matrix: n-dimensional matrix of type L{InvertibleMatrix} with n being the dimension
                               of the configuration space, that is, the dimension of
                               the position / momentum vectors

    @param traj_length: Trajectory length in number of timesteps
    @type traj_length: int

    @param gradient: Gradient which determines the dynamics during a trajectory
    @type gradient: L{AbstractGradient}
    """

    def __init__(self, sampler1, sampler2, timestep, traj_length, gradient, mass_matrix=None):
        
        super(MDRENSSwapParameterInfo, self).__init__(sampler1, sampler2)
        
        self._mass_matrix = mass_matrix
        if self.mass_matrix is None:
            d = len(sampler1.state.position)
            self.mass_matrix = csb.numeric.InvertibleMatrix(numpy.eye(d), numpy.eye(d))

        self._traj_length = traj_length
        self._gradient = gradient
        self._timestep = timestep
    
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

    @property
    def mass_matrix(self):
        return self._mass_matrix
    @mass_matrix.setter
    def mass_matrix(self, value):
        self._mass_matrix = value

class ThermostattedMDRENSSwapParameterInfo(MDRENSSwapParameterInfo):
    """
    @param sampler1: First sampler
    @type sampler1: subclass instance of L{AbstractSingleChainMC}

    @param sampler2: Second sampler
    @type sampler2: subclass instance of L{AbstractSingleChainMC}

    @param timestep: Integration timestep
    @type timestep: float

    @param mass_matrix: Mass matrix
    @type mass_matrix: n-dimensional L{InvertibleMatrix} with n being the dimension
                       of the configuration space, that is, the dimension of
                       the position / momentum vectors

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
        
    def __init__(self, sampler1, sampler2, timestep, traj_length, gradient, mass_matrix=None,
                 temperature=lambda l: 1., collision_probability=0.1, collision_interval=1):
        
        super(ThermostattedMDRENSSwapParameterInfo, self).__init__(sampler1, sampler2, timestep,
																   traj_length, gradient,
																   mass_matrix=mass_matrix)
        
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

    @param traj12: Forward trajectory
    @type traj12: L{Trajectory}

    @param traj21: Reverse trajectory
    @type traj21: L{Trajectory}
    """

    __metaclass__ = ABCMeta
    
    def __init__(self, param_info, traj12, traj21):
        
        self._sampler1 = param_info.sampler1
        self._sampler2 = param_info.sampler2

        self._traj12 = traj12
        self._traj21 = traj21
        
        self._param_info = param_info
        
        self._acceptance_probability = None
        self._accepted = False
        
    @property
    def sampler1(self):
        return self._sampler1

    @property
    def sampler2(self):
        return self._sampler2
    
    @property
    def traj12(self):
        return self._traj12    

    @property
    def traj21(self):
        return self._traj21

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

    @property
    def param_info(self):
        return self._param_info

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

    See L{AbstractSwapCommunicator} for constructor signature.
    """
    
    pass

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
        swapcom = self._calc_pacc_swap(swapcom)
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

        @rtype: L{AbstractSwapCommunicator}
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
            swapcom.sampler1.state = swapcom.traj21.final
            swapcom.sampler2.state = swapcom.traj12.final
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

        T1 = param_info.sampler1.temperature
        T2 = param_info.sampler2.temperature

        momentum_covariance_matrix1 = T1 * param_info.mass_matrix
        momentum_covariance_matrix2 = T2 * param_info.mass_matrix

        d = len(param_info.sampler1.state.position)

        if param_info.mass_matrix.is_unity_multiple:
            momentum1 = numpy.random.normal(scale=numpy.sqrt(T1 * param_info.mass_matrix[0][0]),
                                            size=d)
            momentum2 = numpy.random.normal(scale=numpy.sqrt(T2 * param_info.mass_matrix[0][0]),
                                            size=d)
        else:
            momentum1 = numpy.random.multivariate_normal(mean=numpy.zeros(d),
                                                         cov=momentum_covariance_matrix1)
            momentum2 = numpy.random.multivariate_normal(mean=numpy.zeros(d),
                                                         cov=momentum_covariance_matrix2)        
        
        init_state1 = State(param_info.sampler1.state.position, momentum1)
        init_state2 = State(param_info.sampler2.state.position, momentum2)

        param_info.sampler1.state = init_state1
        param_info.sampler2.state = init_state2
        
        trajinfo12 = RENSTrajInfo(param_info, init_state1, protocol=lambda t, tau: t / tau)
        trajinfo21 = RENSTrajInfo(param_info, init_state2, protocol=lambda t, tau: (tau - t) / tau)
        
        traj12 = self._run_traj_generator(trajinfo12)
        traj21 = self._run_traj_generator(trajinfo21)

        return RENSSwapCommunicator(param_info, traj12, traj21)

    def _calc_pacc_swap(self, swapcom):

        T1 = swapcom.param_info.sampler1.temperature
        T2 = swapcom.param_info.sampler2.temperature
        
        heat12 = swapcom.traj12.heat
        heat21 = swapcom.traj21.heat
        
        proposal1 = swapcom.traj21.final
        proposal2 = swapcom.traj12.final
        
        state1 = swapcom.traj12.initial
        state2 = swapcom.traj21.initial
        
        if swapcom.param_info.mass_matrix.is_unity_multiple:
            inverse_mass_matrix = 1. / swapcom.param_info.mass_matrix[0][0]
        else:
            inverse_mass_matrix = swapcom.param_info.mass_matrix.inverse
        
        E1 = lambda x:-swapcom.sampler1._pdf.log_prob(x)
        E2 = lambda x:-swapcom.sampler2._pdf.log_prob(x)
        K = lambda x: 0.5 * numpy.dot(x.T, numpy.dot(inverse_mass_matrix, x))

        w12 = (K(proposal2.momentum) + E2(proposal2.position)) / T2 - \
              (K(state1.momentum) + E1(state1.position)) / T1 - heat12 
        w21 = (K(proposal1.momentum) + E1(proposal1.position)) / T1 - \
              (K(state2.momentum) + E2(state2.position)) / T2 - heat21

        swapcom.acceptance_probability = csb.numeric.exp(-w12 - w21)

        return swapcom

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

        if len(self._algorithm.param_infos) == 1:
            self._swap_list1.append(0)
            self._swap_list2.append(0)
        else:
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
