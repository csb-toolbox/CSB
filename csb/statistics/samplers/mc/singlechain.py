"""
Various Monte Carlo equilibrium sampling algorithms, which simulate
only one Markov chain.
"""

import numpy
import csb.numeric

from csb.statistics.samplers import State
from csb.statistics.samplers.mc import AbstractSingleChainMC
from csb.statistics.samplers.mc.propagators import MDPropagator


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

    @param timestep: Timestep used for integrating Hamiltonian equations of motion
    @type timestep: float

    @param nsteps: Number of integration steps to be performed
    @type nsteps: int

    @param integrator: Integrator to be used
    @type integrator: L{AbstractIntegrator}
    """

    def __init__(self, pdf, state, gradient, timestep, nsteps, integrator):
        
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
        proposal = gen.generate(self.state, self._nsteps).last
        
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
    (Metropolis, Rosenbluth, Teller, Teller 1953; Hastings, 1970)
    with uniform proposal density.
    
    @param pdf: Probability density function to be sampled from
    @type pdf: L{csb.statistics.pdf.AbstractDensity}

    @param state: Inital state
    @type state: L{State}

    @param max_proposal_distance: Maximum proposal distance
    @type max_proposal_distance: float
    """
    
    def __init__(self, pdf, state, max_proposal_distance):
        
        super(RWMCSampler, self).__init__(pdf, state)
        self._max_proposal_distance = None
        self.max_proposal_distance = max_proposal_distance

    def _propose(self):
        
        return State(self.state.position
                     + (2 * numpy.random.uniform(size=self.state.position.shape) - 1)
                     * self._max_proposal_distance)

    def _calc_pacc(self, proposal):
        
        E = lambda x:-self._pdf.log_prob(x)
        
        pacc = csb.numeric.exp(-(E(proposal.position) - E(self.state.position)))
        return pacc

    @property
    def max_proposal_distance(self):
        return self._max_proposal_distance

    @max_proposal_distance.setter
    def max_proposal_distance(self, value):
        self._max_proposal_distance = float(value)
