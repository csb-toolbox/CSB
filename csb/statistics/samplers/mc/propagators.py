"""
Provides various deterministic and stochastic propagators.
"""

import copy
import numpy

from abc import ABCMeta, abstractmethod
from csb.statistics.samplers import Trajectory


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

        @param return_trajectory: Switch to determine whether the complete trajectory
                                 is returned or only the final state
        @type return_trajectory: boolean

        @rtype: L{Trajectory}
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
        
        traj = integrator.integrate(copy.deepcopy(init_state), length, return_trajectory)
        return traj

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

    def generate(self, init_state, length, return_trajectory=False):
        
        integrator = self._integrator(self.timestep, self.gradient)

        result = []
        if return_trajectory:
            result.append(init_state)

        heat = 0.
        state = copy.deepcopy(init_state)
        for i in range(length):
            integrator.integrate_once(state, i)
            if i % self._update_interval == 0:
                state.momentum, stepheat = self._update(state.momentum,
                                                        self._temperature(i * self.timestep),
                                                        self._collision_probability)
                
                heat += stepheat
            if return_trajectory:
                result.append(state)

        if not return_trajectory:
            result.append(state)

        traj = Trajectory(result, heat)
        return traj
