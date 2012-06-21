"""
Various Monte Carlo equilibrium sampling algorithms, which simulate only one Markov chain.

In the following, I will shortly explain how sample from a PDF using the HMCSampler class.
The follwing short code sniplet draws 5000 samples from a 1D normal distribution and plots
them using matplotlib:


import numpy
import matplotlib.pyplot as plt

from csb.statistics.pdf import Normal
from csb.statistics.samplers import State
from csb.statistics.samplers.mc.singlechain import HMCSampler

initial_state = State(numpy.array([1.]))
grad = lambda q, t: q
timestep = 1.5
nsteps = 30
nsamples = 5000

sampler = HMCSampler(Normal(), initial_state, grad, timestep, nsteps)

states = []

for i in range(nsamples):
    sampler.sample()
    states.append(sampler.state)

print 'acceptance rate:', sampler.acceptance_rate

states = [state.position[0]for state in states]

plt.hist([numpy.random.normal(size=5000), states], bins=20, normed=True, \
         label=['numpy.random.normal', 'HMC'])
plt.legend()
plt.show()


First, several things which are being needed are imported.
As every sampler in this module implements a Markov Chain, an initial state has to be
chosen. In the following lines, several parameters are set:
- the gradient of the negative log-probability of the PDF under consideration
- the integration timestep
- the number of integration steps to be performed in each iteration, that is, the HMC
  trajectory length
- the number of samples to be drawn

The empty list states is initialized. It will serve to store the samples drawn.
In the loop, sampler.sample() is repeatedly called. After each call of sampler.sample(),
the current state of the Markov Chain is stored in sampler.state and this state is appended
to the sample storage list.

Then the acceptance rate is printed, the numeric values are being extracted from the State
objects in states, a histogram is created and finally plotted.
"""

import numpy
import csb.numeric

from csb.statistics.samplers import State
from csb.statistics.samplers.mc import AbstractSingleChainMC
from csb.statistics.samplers.mc.propagators import MDPropagator
from csb.numeric.integrators import FastLeapFrog


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

    @param integrator: Subclass of L{AbstractIntegrator} to be used for
                       integrating Hamiltionian equations of motion
    @type integrator: type
    """

    def __init__(self, pdf, state, gradient, timestep, nsteps, integrator=FastLeapFrog):
        
        super(HMCSampler, self).__init__(pdf, state)
        
        self._timestep = None
        self.timestep = timestep
        self._nsteps = None
        self.nsteps = nsteps
        self._integrator = integrator
        self._gradient = gradient

    def _propose(self):
        
        gen = MDPropagator(self._gradient, self._timestep, self._integrator)
        momenta = numpy.random.normal(size=self.state.position.shape)
        self.state = State(self.state.position, momenta)
        proposal = gen.generate(self.state, self._nsteps).final
        
        return proposal

    def _calc_pacc(self, proposal):
        
        E = lambda x:-self._pdf.log_prob(x)
        
        pacc = csb.numeric.exp(-0.5 * sum(proposal.momentum ** 2) - E(proposal.position)
                               + 0.5 * sum(self.state.momentum ** 2) + E(self.state.position))

        return pacc

    @property
    def timestep(self):
        return self._timestep

    @timestep.setter
    def timestep(self, value):
        self._timestep = float(value)

    @property
    def nsteps(self):
        return self._nsteps

    @nsteps.setter
    def nsteps(self, value):
        self._nsteps = int(value)

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
    """
    
    def __init__(self, pdf, state, stepsize=1., proposal_density=None):
        
        super(RWMCSampler, self).__init__(pdf, state)
        self._stepsize = None
        self.stepsize = stepsize
        if proposal_density == None:
            self._proposal_density = lambda x, s: x.position + s * numpy.random.uniform(size=x.position.shape, low=-1., high=1.)
        else:
            self._proposal_density = proposal_density

    def _propose(self):
        
        return State(self._proposal_density(self._state, self.stepsize))

    def _calc_pacc(self, proposal):
        
        E = lambda x:-self._pdf.log_prob(x)
        
        pacc = csb.numeric.exp(-(E(proposal.position) - E(self.state.position)))
        return pacc

    @property
    def stepsize(self):
        return self._stepsize

    @stepsize.setter
    def stepsize(self, value):
        self._stepsize = float(value)
