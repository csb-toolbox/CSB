"""
Various Monte Carlo equilibrium sampling algorithms, which simulate only one Markov chain.

Here is how to sample from a PDF using the L{HMCSampler} class. In the following
snippet we draw 5000 samples from a 1D normal distribution and plot them:


    >>> import numpy
    >>> from csb.io.plots import Chart
    >>> from csb.statistics.pdf import Normal
    >>> from csb.statistics.samplers import State
    >>> from csb.statistics.samplers.mc.singlechain import HMCSampler
    
    >>> initial_state = State(numpy.array([1.]))
    >>> grad = lambda q, t: q
    >>> timestep = 1.5
    >>> nsteps = 30
    >>> nsamples = 5000
    
    >>> sampler = HMCSampler(Normal(), initial_state, grad, timestep, nsteps)
    
    >>> states = []
    >>> for i in range(nsamples):
            sampler.sample()
            states.append(sampler.state)
    
    >>> print('acceptance rate:', sampler.acceptance_rate)
    0.8
    
    >>> states = [state.position[0]for state in states]
    >>> chart = Chart()
    >>> chart.plot.hist([numpy.random.normal(size=5000), states], bins=20, normed=True)
    >>> chart.plot.legend(['numpy.random.normal', 'HMC'])
    >>> chart.show()


First, several things which are being needed are imported.
As every sampler in this module implements a Markov Chain, an initial state has to be
chosen. In the following lines, several parameters are set:
    - the gradient of the negative log-probability of the PDF under consideration
    - the integration timestep
    - the number of integration steps to be performed in each iteration, that is, the HMC
      trajectory length
    - the number of samples to be drawn

The empty list states is initialized. It will serve to store the samples drawn.
In the loop, C{sampler.sample()} is repeatedly called. After each call of C{sampler.sample()},
the current state of the Markov Chain is stored in sampler.state and this state is appended
to the sample storage list.

Then the acceptance rate is printed, the numeric values are being extracted from the
L{State} objects in states, a histogram is created and finally plotted.
"""

import numpy
import csb.numeric
import csb.core

from abc import ABCMeta, abstractmethod

from csb.statistics.samplers import State
from csb.statistics.samplers.mc import AbstractMC, MCCollection, augment_state
from csb.statistics.samplers.mc.propagators import MDPropagator
from csb.statistics.samplers.mc.neqsteppropagator import NonequilibriumStepPropagator
from csb.numeric.integrators import FastLeapFrog
from csb.numeric import InvertibleMatrix


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

        accepted = None
        if numpy.random.random() < pacc:
            accepted = True
        else:
            accepted = False
            
        if accepted == True:
            self._accept_proposal(proposal_communicator.proposal_state)
            
        self._update_statistics(accepted)
        self._last_move_accepted = accepted
            
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

    def _accept_proposal(self, proposal_state):
        """
        Accept the proposal state by setting it as the current state of the sampler
        object

        @param proposal_state: The proposal state
        @type proposal_state: L{State}
        """

        self.state = proposal_state
        
    def _update_statistics(self, accepted):
        """
        Update the sampling statistics.

        @param accepted: Whether or not the proposal state has been accepted
        @type accepted: boolean
        """

        self._nmoves += 1
        self._accepted += int(accepted)

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


class HMCSampler(AbstractSingleChainMC):
    """
    Hamilton Monte Carlo (HMC, also called Hybrid Monte Carlo by the inventors,
    Duane, Kennedy, Pendleton, Duncan 1987).

    @param pdf: Probability density function to be sampled from
    @type pdf: L{csb.statistics.pdf.AbstractDensity}

    @param state: Inital state
    @type state: L{State}

    @param gradient: Gradient of the negative log-probability
    @type gradient: L{AbstractGradient}

    @param timestep: Timestep used for integration
    @type timestep: float

    @param nsteps: Number of integration steps to be performed in
                   each iteration
    @type nsteps: int

    @param mass_matrix: Mass matrix
    @type mass_matrix: n-dimensional L{InvertibleMatrix} with n being the dimension
                       of the configuration space, that is, the dimension of
                       the position / momentum vectors

    @param integrator: Subclass of L{AbstractIntegrator} to be used for
                       integrating Hamiltionian equations of motion
    @type integrator: L{AbstractIntegrator}

    @param temperature: Pseudo-temperature of the Boltzmann ensemble
                        M{p(x) = 1/N * exp(-1/T * E(x))} with the
                        pseudo-energy defined as M{E(x) = -log(p(x))}
                        where M{p(x)} is the PDF under consideration
    @type temperature: float
    """

    def __init__(self, pdf, state, gradient, timestep, nsteps,
				 mass_matrix=None, integrator=FastLeapFrog, temperature=1.):

        super(HMCSampler, self).__init__(pdf, state, temperature)
        
        self._timestep = None
        self.timestep = timestep
        self._nsteps = None
        self.nsteps = nsteps

        self._d = len(self.state.position)

        self._mass_matrix = mass_matrix
        if self.mass_matrix is None:
            self.mass_matrix = InvertibleMatrix(numpy.eye(self._d), numpy.eye(self._d))
            
        self._momentum_covariance_matrix = self._temperature * self.mass_matrix

        self._integrator = integrator
        self._gradient = gradient

        self._propagator = MDPropagator(self._gradient, self._timestep,
                                        mass_matrix=self._mass_matrix,
                                        integrator=self._integrator)

    def _propose(self):

        current_state = self.state.clone()
        current_state = augment_state(current_state, self.temperature, self.mass_matrix)
        proposal_state = self._propagator.generate(current_state, self._nsteps).final
        
        return SimpleProposalCommunicator(current_state, proposal_state)

    def _calc_pacc(self, proposal_communicator):

        current_state = proposal_communicator.current_state
        proposal_state = proposal_communicator.proposal_state
        
        E = lambda x: -self._pdf.log_prob(x)
        K = lambda x: 0.5 * numpy.dot(x.T, numpy.dot(self.mass_matrix.inverse, x))

        pacc = csb.numeric.exp((-K(proposal_state.momentum) - E(proposal_state.position)
                               + K(current_state.momentum) + E(current_state.position)) / 
                               self.temperature)

        if self.state.momentum is None:
            proposal_communicator.proposal_state.momentum = None
        else:
            proposal_communicator.proposal_state.momentum = self.state.momentum

        return pacc

    @property
    def timestep(self):
        return self._timestep

    @timestep.setter
    def timestep(self, value):
        self._timestep = float(value)
        if "_propagator" in dir(self):
            self._propagator.timestep = self._timestep

    @property
    def nsteps(self):
        return self._nsteps

    @nsteps.setter
    def nsteps(self, value):
        self._nsteps = int(value)

    @property
    def mass_matrix(self):
        return self._mass_matrix
    @mass_matrix.setter
    def mass_matrix(self, value):
        self._mass_matrix = value
        if "_propagator" in dir(self):
           self._propagator.mass_matrix = self._mass_matrix


class RWMCSampler(AbstractSingleChainMC):
    """
    Random Walk Metropolis Monte Carlo implementation
    (Metropolis, Rosenbluth, Teller, Teller 1953; Hastings, 1970).
    
    @param pdf: Probability density function to be sampled from
    @type pdf: L{csb.statistics.pdf.AbstractDensity}

    @param state: Inital state
    @type state: L{State}

    @param stepsize: Serves to set the step size in
                     proposal_density, e.g. for automatic acceptance
                     rate adaption
    @type stepsize: float

    @param proposal_density: The proposal density as a function f(x, s)
                             of the current state x and the stepsize s.
                             By default, the proposal density is uniform,
                             centered around x, and has width s.
    @type proposal_density: callable

    @param temperature: Pseudo-temperature of the Boltzmann ensemble
                        M{p(x) = 1/N * exp(-1/T * E(x))} with the
                        pseudo-energy defined as M{E(x) = -log(p(x))}
                        where M{p(x)} is the PDF under consideration
    @type temperature: float
    """
    
    def __init__(self, pdf, state, stepsize=1., proposal_density=None, temperature=1.):
        
        super(RWMCSampler, self).__init__(pdf, state, temperature)
        self._stepsize = None
        self.stepsize = stepsize
        if proposal_density == None:
            self._proposal_density = lambda x, s: x.position + \
                                     s * numpy.random.uniform(size=x.position.shape, low=-1., high=1.)
        else:
            self._proposal_density = proposal_density

    def _propose(self):

        current_state = self.state.clone()
        proposal_state = self.state.clone()
        proposal_state.position = self._proposal_density(current_state, self.stepsize)
        
        return SimpleProposalCommunicator(current_state, proposal_state)

    def _calc_pacc(self, proposal_communicator):

        current_state = proposal_communicator.current_state
        proposal_state = proposal_communicator.proposal_state
        E = lambda x:-self._pdf.log_prob(x)
        
        pacc = csb.numeric.exp((-(E(proposal_state.position) - E(current_state.position))) /
                               self.temperature)
        return pacc

    @property
    def stepsize(self):
        return self._stepsize

    @stepsize.setter
    def stepsize(self, value):
        self._stepsize = float(value)


class AbstractNCMCSampler(AbstractSingleChainMC):
    """
    Implementation of the NCMC sampling algorithm (Nilmeier et al., "Nonequilibrium candidate Monte
    Carlo is an efficient tool for equilibrium simulation", PNAS 2011) for sampling from one
    ensemble only.
    Subclasses have to specify the acceptance probability, which depends on the kind of
    perturbations and propagations in the protocol.

    @param state: Inital state
    @type state: L{State}

    @param protocol: Nonequilibrium protocol with alternating perturbations and propagations
    @type protocol: L{Protocol}

    @param reverse_protocol: The reversed version of the protocol, that is, the order of
                             perturbations and propagations in each step is reversed
    @type reverse_protocol: L{Protocol}
    """

    __metaclass__ = ABCMeta
    
    def __init__(self, state, protocol, reverse_protocol):

        self._protocol = None
        self.protocol = protocol
        self._reverse_protocol = None
        self.reverse_protocol = reverse_protocol

        pdf = self.protocol.steps[0].perturbation.sys_before.hamiltonian
        temperature = self.protocol.steps[0].perturbation.sys_before.hamiltonian.temperature

        super(AbstractNCMCSampler, self).__init__(pdf, state, temperature)

    def _pick_protocol(self):
        """
        Picks either the protocol or the reversed protocol with equal probability.

        @return: Either the protocol or the reversed protocol
        @rtype: L{Protocol}
        """
        
        if numpy.random.random() < 0.5:
            return self.protocol
        else:
            return self.reverse_protocol

    def _propose(self):

        protocol = self._pick_protocol()
        
        gen = NonequilibriumStepPropagator(protocol)

        traj = gen.generate(self.state)

        return NCMCProposalCommunicator(traj)

    def _accept_proposal(self, proposal_state):

        if self.state.momentum is not None:
            proposal_state.momentum *= -1.0
        else:
            proposal_state.momentum = None
        
        super(AbstractNCMCSampler, self)._accept_proposal(proposal_state)

    @property
    def protocol(self):
        return self._protocol
    @protocol.setter
    def protocol(self, value):
        self._protocol = value

    @property
    def reverse_protocol(self):
        return self._reverse_protocol
    @reverse_protocol.setter
    def reverse_protocol(self, value):
        self._reverse_protocol = value


class SimpleProposalCommunicator(object):
    """
    This holds all the information needed to calculate the acceptance
    probability in both the L{RWMCSampler} and L{HMCSampler} classes,
    that is, only the proposal state.
    For more advanced algorithms, one may derive classes capable of
    holding the neccessary additional information from this class.

    @param current_state: Current state
    @type current_state: L{State}
    
    @param proposal_state: Proposal state
    @type proposal_state: L{State}
    """

    __metaclass__ = ABCMeta

    def __init__(self, current_state, proposal_state):

        self._current_state = current_state
        self._proposal_state = proposal_state

    @property
    def current_state(self):
        return self._current_state
        
    @property
    def proposal_state(self):
        return self._proposal_state


class NCMCProposalCommunicator(SimpleProposalCommunicator):
    """
    Holds all information (that is, the trajectory with heat, work, Hamiltonian difference
    and jacbian) needed to calculate the acceptance probability in the AbstractNCMCSampler class.

    @param traj: Non-equilibrium trajectory stemming from a stepwise protocol
    @type traj: NCMCTrajectory
    """

    def __init__(self, traj):

        self._traj = None
        self.traj = traj

        super(NCMCProposalCommunicator, self).__init__(traj.initial, traj.final)
