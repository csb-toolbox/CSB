"""
Provides various deterministic and stochastic propagators.
"""

import numpy

import csb.core

from abc import ABCMeta, abstractmethod, abstractproperty
from csb.statistics.samplers import AbstractState, State
from csb.statistics.samplers.mc import Trajectory, AbstractPropagationResult
from csb.statistics.samplers.mc import TrajectoryBuilder ,ShortTrajectoryBuilder
from csb.numeric.integrators import FastLeapFrog, VelocityVerlet

class AbstractPropagator(object):
    """
    Abstract propagator class. Subclasses serve to propagate
    an inital state by some dynamics to a final state.
    """

    __metaclass__ = ABCMeta

    @abstractmethod
    def generate(self, init_state, length, return_trajectory=False):
        """
        Generate a trajectory, starting from an initial state with a certain length.

        @param init_state: Initial state from which to propagate
        @type init_state: L{State}

        @param length: Length of the trajectory (in integration steps or stochastic moves)
        @type length: int

        @param return_trajectory: Return complete L{Trajectory} instead of the initial
                                  and final states only (L{PropagationResult})                                 
        @type return_trajectory: boolean

        @rtype: L{AbstractPropagationResult}
        """
        pass    

class MDPropagator(AbstractPropagator):
    """
    Molecular Dynamics propagator. Generates a trajectory
    by integration of Hamiltionian equations of motion.

    @param gradient: Gradient of potential energy. Guides the dynamics.
    @type gradient: L{AbstractGradient}

    @param timestep: Timestep to be used for integration
    @type timestep: float

    @param integrator: Subclass of L{AbstractIntegrator} to be used to integrate
                       Hamiltonian equations of motion
    @type integrator: type
    """

    def __init__(self, gradient, timestep, integrator=FastLeapFrog):
        
        super(MDPropagator, self).__init__()

        self._gradient = None
        self.gradient = gradient
        
        self._timestep = None
        self.timestep = timestep
        
        self._integrator = integrator

    @property
    def gradient(self):
        return self._gradient
    @gradient.setter
    def gradient(self, value):
        self._gradient = value

    @property
    def timestep(self):
        return self._timestep
    @timestep.setter
    def timestep(self, value):
        self._timestep = float(value)

    def generate(self, init_state, length, return_trajectory=False):
        
        integrator = self._integrator(self.timestep, self.gradient)
        
        result = integrator.integrate(init_state, length, return_trajectory)
        return result

class ThermostattedMDPropagator(MDPropagator):
    """
    Thermostatted Molecular Dynamics propagator. Employs the Andersen thermostat
    which simulates collision with particles of a heat bath at a given temperature.

    @param gradient: Gradient of potential energy. Guides the dynamics.
    @type gradient: L{AbstractGradient}

    @param timestep: Timestep to be used for integration
    @type timestep: float

    @param temperature: Time-dependent temperature
    @type temperature: Real-valued function

    @param collision_probability: collision probability within duration of one timestep
    @type collision_probability: float

    @param update_interval: Interval with which momenta are redrawn
    @type update_interval: int

    @param integrator: Subclass of L{AbstractIntegrator} to be used to perform
                       integration steps between momentum updates
    @type integrator: type
    """

    def __init__(self, gradient, timestep, temperature, collision_probability=0.1,
                 update_interval=1, integrator=VelocityVerlet):
        
        super(ThermostattedMDPropagator, self).__init__(gradient, timestep, integrator)
        
        self._collision_probability = collision_probability
        self._update_interval = update_interval
        self._temperature = temperature

    def _update(self, momentum, T, collision_probability):
        """
        Simulate collision with heat bath particles.

        @param momentum: Momentum
        @type momentum: one-dimensional numpy array of numbers
        
        @param T: Temperature of the heat bath
        @type T: float
        
        @param collision_probability: collision probability within duration of one timestep
        @type collision_probability: float

        @rtype: tuple (updated momentum, heat induced by the update)
        """
        
        heat = 0.
        update_list = []
        for k in range(len(momentum)):
            if numpy.random.random() < collision_probability:
                update_list.append(k)
        if len(update_list) > 0:
            ke_old = 0.5 * sum(momentum ** 2)
            for k in range(len(momentum)):
                if k in update_list: momentum[k] = numpy.random.normal(scale=numpy.sqrt(T))
            heat = 0.5 * sum(momentum ** 2) - ke_old

        return momentum, heat
    
    def _step(self, i, state, heat, integrator):
        """
        Performs one step consisting of an integration step
        and possibly a momentum update

        @param i: integration step count
        @type i: int

        @param state: state to be updated
        @type state: L{State}

        @param heat: heat produced up to the current integration step
        @type heat: float

        @param integrator: integration scheme used to evolve the state deterministically
        @type integrator: L{AbstractIntegrator}
        """

        state = integrator.integrate_once(state, i)
        
        if i % self._update_interval == 0:
            state.momentum, stepheat = self._update(state.momentum,
                                                    self._temperature(i * self.timestep),
                                                    self._collision_probability)
            
            heat += stepheat
            
        return state, heat

    def generate(self, init_state, length, return_trajectory=False):
        
        integrator = self._integrator(self.timestep, self.gradient)
        builder = TrajectoryBuilder.create(full=return_trajectory)

        builder.add_initial_state(init_state)

        heat = 0.
        state = init_state.clone()
        
        for i in range(length - 1):
            state, heat = self._step(i, state, heat, integrator)
            builder.add_intermediate_state(state)

        state, heat = self._step(length - 1, state, heat, integrator)        
        builder.add_final_state(state)            

        traj = builder.product
        traj.heat = heat
        
        return traj

class AbstractMCPropagator(AbstractPropagator):
    """
    Provides the interface for MC trajectory generators. Implementations
    generate a sequence of states according to some implementation of
    L{AbstractSingleChainMC}.

    @param pdf: PDF to sample from
    @type pdf: L{AbstractDensity}
    """

    __metaclass__ = ABCMeta

    def __init__(self, pdf):

        self._pdf = pdf
        self._acceptance_rate = 0.0

    def generate(self, init_state, length, return_trajectory=True):

        self._init_sampler()
        self._sampler.state = init_state

        builder = TrajectoryBuilder.create(full=return_trajectory)

        builder.add_initial_state(init_state)

        for i in range(length):
            self._sampler.sample()
            builder.add_intermediate_state(self._sampler.state)

        self._acceptance_rate = self._sampler.acceptance_rate

        return builder.product

    @abstractmethod
    def _init_sampler(self):
        """
        Initializes the sampler with which to obtain the MC state
        trajectory.
        """
        
        pass

    @property
    def acceptance_rate(self):
        """
        Acceptance rate of the MC sampler that generated the
        trajectory.
        """
        return self._acceptance_rate

class RWMCPropagator(AbstractMCPropagator):
    """
    Draws a number of samples from a PDF using the L{RWMCSampler} and
    returns them as a L{Trajectory}.

    @param pdf: PDF to sample from
    @type pdf: L{AbstractDensity}
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

    def __init__(self, pdf, stepsize=1., proposal_density=None):

        super(RWMCPropagator, self).__init__(pdf)

        self._stepsize = stepsize
        self._proposal_density = proposal_density

        self._init_sampler()

    def _init_sampler(self):

        from csb.statistics.samplers.mc.singlechain import RWMCSampler

        dummy_state = State(numpy.array([0.]))
        self._sampler = RWMCSampler(self._pdf, dummy_state, self._stepsize,
                                    self._proposal_density)

class HMCPropagator(AbstractMCPropagator):
    """
    Draws a number of samples from a PDF using the L{HMCSampler} and
    returns them as a L{Trajectory}.

    @param pdf: PDF to sample from
    @type pdf: L{AbstractDensity}
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
    
    def __init__(self, pdf, gradient, timestep, nsteps, integrator=FastLeapFrog):

        super(HMCPropagator, self).__init__(pdf)

        self._gradient = gradient
        self._timestep = timestep
        self._nsteps = nsteps
        self._integrator = integrator

    def _init_sampler(self):

        from csb.statistics.samplers.mc.singlechain import HMCSampler
        
        dummy_state = State(numpy.array([0.]))        
        self._sampler = HMCSampler(self._pdf, dummy_state, self._gradient,
                                   self._timestep, self._nsteps, self._integrator)
