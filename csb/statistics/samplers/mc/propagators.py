"""
Provides various deterministic and stochastic propagators.
"""

import numpy

import csb.pyutils

from abc import ABCMeta, abstractmethod, abstractproperty
from csb.statistics.samplers import AbstractState


class AbstractPropagationResult(object):
    
    __metaclass__ = ABCMeta 
    
    @abstractproperty
    def initial(self):
        pass
    
    @abstractproperty
    def final(self):
        pass

class PropagationResult(AbstractPropagationResult):
    
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
    
class Trajectory(csb.pyutils.CollectionContainer, AbstractPropagationResult):
    """
    Ordered collection of states, representing a phase-space trajectory.

    @param items: list of states
    @type items: list of L{AbstractState}
    @param heat: heat produced during the trajectory
    @type heat: float
    """
    
    def __init__(self, items, heat=0.0):
        
        super(Trajectory, self).__init__(items, type=AbstractState)
        
        self._heat = None    
        self.heat = heat
    
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

class TrajectoryBuilder(object):
    
    def __init__(self, heat=0.0):
        self._heat = heat
        self._states = []
        
    @staticmethod
    def create(full=True):
        
        if full:
            return TrajectoryBuilder()
        else:
            return ShortTrajectoryBuilder()
        
    @property
    def product(self):
        return Trajectory(self._states, heat=self._heat)

    def add_initial_state(self, state):
        self._states.insert(0, state.clone())
        
    def add_intermediate_state(self, state):
        self._states.append(state.clone())
    
    def add_final_state(self, state):
        self._states.append(state.clone())
    
class ShortTrajectoryBuilder(TrajectoryBuilder):    

    def add_intermediate_state(self, state):
        pass

    @property
    def product(self):
        
        if len(self._states) != 2:
            raise ValueError("Can't create a product, two states required")
        
        initial, final = self._states
        return PropagationResult(initial, final, heat=self._heat)

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

    @param integrator: Integrator to be used
    @type integrator: L{AbstractIntegrator}
    """

    def __init__(self, gradient, timestep, integrator):
        
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

    @param integrator: Integrator to be used
    @type integrator: L{AbstractIntegrator}

    @param temperature: Time-dependent temperature
    @type temperature: Real-valued function

    @param collision_probability: collision probability within duration of one timestep
    @type collision_probability: float

    @param update_interval: Interval with which momenta are redrawn
    @type update_interval: int
    """

    def __init__(self, gradient, timestep, integrator, temperature,
                 collision_probability=0.1, update_interval=1):
        
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


