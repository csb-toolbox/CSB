"""
Provides various deterministic and stochastic propagators.
"""

import numpy

from abc import ABCMeta, abstractmethod

from csb.statistics.samplers.mc import TrajectoryBuilder
from csb.numeric.integrators import FastLeapFrog, VelocityVerlet
from csb.numeric import InvertibleMatrix

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

    @param mass_matrix: Mass matrix
    @type mass_matrix: n-dimensional L{InvertibleMatrix} with n being the dimension
                               of the configuration space, that is, the dimension of
                               the position / momentum vectors

    @param integrator: Subclass of L{AbstractIntegrator} to be used to integrate
                       Hamiltonian equations of motion
    @type integrator: type
    """

    def __init__(self, gradient, timestep, mass_matrix=None, integrator=FastLeapFrog):
        
        super(MDPropagator, self).__init__()

        self._gradient = None
        self.gradient = gradient
        self._mass_matrix = mass_matrix
        self._timestep = None
        self.timestep = timestep
        self._integrator = integrator

        self._first_run = True

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

    @property
    def mass_matrix(self):
        return self._mass_matrix
    @mass_matrix.setter
    def mass_matrix(self, value):
        self._mass_matrix = value

    def generate(self, init_state, length, return_trajectory=False):
        
        integrator = self._integrator(self.timestep, self.gradient)
        
        result = integrator.integrate(init_state, length,
                                      mass_matrix=self.mass_matrix,
                                      return_trajectory=return_trajectory)
        
        return result

class Looper(object):
    """
    Implements an iterable list with a ring-like topology,
    that is, if the iterator points on the last element,
    next() returns the first element.
    """
    
    def __init__(self, items):
        
        self._items = items
        self._n_items = len(self._items)
        self._current = 0
        
    def __iter__(self):
        
        return self
    
    def next(self):
        
        if self._current == self._n_items:
            self._current = 0
            
        self._current += 1
        
        return self._items[self._current - 1]


class ThermostattedMDPropagator(MDPropagator):
    """
    Thermostatted Molecular Dynamics propagator. Employs the Andersen thermostat
    which simulates collision with particles of a heat bath at a given temperature.

    @param gradient: Gradient of potential energy. Guides the dynamics.
    @type gradient: L{AbstractGradient}

    @param timestep: Timestep to be used for integration
    @type timestep: float

    @param mass_matrix: Mass matrix
    @type mass_matrix: n-dimensional L{InvertibleMatrix} with n being the dimension
                               of the configuration space, that is, the dimension of
                               the position / momentum vectors

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

    def __init__(self, gradient, timestep, mass_matrix=None, temperature=lambda t: 1.,
                 collision_probability=0.1, update_interval=1, integrator=VelocityVerlet):
        
        super(ThermostattedMDPropagator, self).__init__(gradient, timestep,
                                                        mass_matrix, integrator)
        
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

        d = len(momentum)
        
        heat = 0.
        update_list = numpy.where(numpy.random.random(d) < collision_probability)[0]

        if len(update_list) > 0:
            K = None
            if self.mass_matrix.is_unity_multiple:
                K = lambda x: 0.5 * sum(x ** 2) / self.mass_matrix[0][0]
            else:
                K = lambda x: 0.5 * numpy.dot(x.T, numpy.dot(self.mass_matrix.inverse, x))

            ke_old = K(momentum)
            
            updated_momentum = [numpy.sqrt(T) * self._random_loopers[i].next() for i in update_list]
            momentum[update_list] = updated_momentum
            heat = (K(momentum) - ke_old) / T

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

        state = integrator.integrate_once(state, i, mass_matrix=self.mass_matrix)

        if i % self._update_interval == 0:
            state.momentum, stepheat = self._update(state.momentum,
                                                    self._temperature(i * self.timestep),
                                                    self._collision_probability)
            
            heat += stepheat
            
        return state, heat

    def generate(self, init_state, length, return_trajectory=False):

        if self._first_run == True and self.mass_matrix is None:
            d = len(init_state.position)
            self.mass_matrix = InvertibleMatrix(numpy.eye(d), numpy.eye(d))

        integrator = self._integrator(self.timestep, self.gradient)
        builder = TrajectoryBuilder.create(full=return_trajectory)

        builder.add_initial_state(init_state)

        heat = 0.
        state = init_state.clone()

        d = len(state.position)

        n_randoms = int(1.5 * length * self._collision_probability / float(self._update_interval))
        
        if n_randoms < 5:
            n_randoms = 5

        if not self.mass_matrix.is_unity_multiple:
            randoms = numpy.random.multivariate_normal(mean=numpy.zeros(d),
                                                       cov=self.mass_matrix,
                                                       size=n_randoms).T
        else:
            randoms = numpy.random.normal(scale=numpy.sqrt(self.mass_matrix[0][0]),
                                          size=(d, n_randoms))
        self._random_loopers = [Looper(x) for x in randoms]
        
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

    @param temperature: See documentation of L{AbstractSingleChainMC}
    @type temperature: float
    """

    __metaclass__ = ABCMeta

    def __init__(self, pdf, temperature=1.):

        self._pdf = pdf
        self._temperature = temperature
        self._acceptance_rate = 0.0

    def generate(self, init_state, length, return_trajectory=True):

        self._init_sampler(init_state)
        self._sampler.state = init_state

        builder = TrajectoryBuilder.create(full=return_trajectory)

        builder.add_initial_state(init_state)

        for i in range(length):
            self._sampler.sample()
            if i != length - 1:
                builder.add_intermediate_state(self._sampler.state)

        builder.add_final_state(self._sampler.state)
                
        self._acceptance_rate = self._sampler.acceptance_rate
        
        return builder.product

    @abstractmethod
    def _init_sampler(self, init_state):
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

    @param temperature: See documentation of L{AbstractSingleChainMC}
    @type temperature: float
    """

    def __init__(self, pdf, stepsize=1., proposal_density=None, temperature=1.):

        super(RWMCPropagator, self).__init__(pdf, temperature)

        self._stepsize = stepsize
        self._proposal_density = proposal_density

    def _init_sampler(self, init_state):

        from csb.statistics.samplers.mc.singlechain import RWMCSampler

        self._sampler = RWMCSampler(self._pdf, init_state, self._stepsize,
                                    self._proposal_density, self._temperature)

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

    @param mass_matrix: Mass matrix
    @type mass_matrix: n-dimensional L{InvertibleMatrix} with n being the dimension
                               of the configuration space, that is, the dimension of
                               the position / momentum vectors

    @param integrator: Subclass of L{AbstractIntegrator} to be used for
                       integrating Hamiltionian equations of motion
    @type integrator: type

    @param temperature: See documentation of L{AbstractSingleChainMC}
    @type temperature: float
    """
    
    def __init__(self, pdf, gradient, timestep, nsteps, mass_matrix=None,
                 integrator=FastLeapFrog, temperature=1.):

        super(HMCPropagator, self).__init__(pdf, temperature)

        self._gradient = gradient
        self._timestep = timestep
        self._nsteps = nsteps
        self._mass_matrix = mass_matrix
        self._integrator = integrator

    def _init_sampler(self, init_state):

        from csb.statistics.samplers.mc.singlechain import HMCSampler

        self._sampler = HMCSampler(self._pdf, init_state, self._gradient,
                                   self._timestep, self._nsteps,
                                   mass_matrix=self.mass_matrix,
                                   integrator=self._integrator, temperature=self._temperature)

    @property
    def mass_matrix(self):
        return self._mass_matrix
    @mass_matrix.setter
    def mass_matrix(self, value):
        self._mass_matrix = value


class AbstractNCMCPropagator(AbstractMCPropagator):
    """
    Draws a number of samples from a PDF using the L{AbstractNCMCSampler}.

    @param protocol: The nonequilibrium protocol specifying a sequence of
                     perturbation and propagation steps
    @type protocol: L{Protocol}

    @param reverse_protocol: The protocol with the order of perturbation and
                             propagation reversed in each step.
    @type reverse_protocol: L{Protocol}
    """

    __metaclass__ = ABCMeta
    
    def __init__(self, protocol, reverse_protocol):

        self._protocol = None
        self.protocol = protocol
        self._reverse_protocol = None
        self.reverse_protocol = reverse_protocol

        pdf = self.protocol.steps[0].perturbation.sys_before.hamiltonian
        temperature = self.protocol.steps[0].perturbation.sys_before.hamiltonian.temperature

        super(AbstractNCMCPropagator, self).__init__(pdf, temperature)

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
