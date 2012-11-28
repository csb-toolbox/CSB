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

from csb.statistics.samplers import State
from csb.statistics.samplers.mc import AbstractSingleChainMC
from csb.statistics.samplers.mc import SimpleProposalCommunicator
from csb.statistics.samplers.mc.propagators import MDPropagator
from csb.numeric.integrators import FastLeapFrog
from csb.numeric import InvertibleMatrix

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

        if not self.mass_matrix.is_unity_multiple:
            momenta = numpy.random.multivariate_normal(mean=numpy.zeros(self._d), 
                                                       cov=self._momentum_covariance_matrix)
        else:
            mass = self.mass_matrix[0][0]
            momenta = numpy.random.normal(size=self._d, scale=numpy.sqrt(self.temperature * mass))
            
        self.state = State(self.state.position, momenta)
        proposal = self._propagator.generate(self.state, self._nsteps).final
        
        return SimpleProposalCommunicator(proposal)

    def _calc_pacc(self, proposal_communicator):

        proposal = proposal_communicator.proposal
        E = lambda x: -self._pdf.log_prob(x)
        K = lambda x: 0.5 * numpy.dot(x.T, numpy.dot(self.mass_matrix.inverse, x))

        pacc = csb.numeric.exp((-K(proposal.momentum) - E(proposal.position)
                               + K(self.state.momentum) + E(self.state.position)) / self.temperature)
        
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
            self._proposal_density = lambda x, s: x.position + s * numpy.random.uniform(size=x.position.shape, low=-1., high=1.)
        else:
            self._proposal_density = proposal_density

    def _propose(self):
        
        return SimpleProposalCommunicator(State(self._proposal_density(self._state, self.stepsize)))

    def _calc_pacc(self, proposal_communicator):

        proposal = proposal_communicator.proposal
        E = lambda x:-self._pdf.log_prob(x)
        
        pacc = csb.numeric.exp((-(E(proposal.position) - E(self.state.position))) / self.temperature)
        return pacc

    @property
    def stepsize(self):
        return self._stepsize

    @stepsize.setter
    def stepsize(self, value):
        self._stepsize = float(value)
