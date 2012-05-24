"""
Provides various integration schemes and an abstract gradient class.
"""

from abc import ABCMeta, abstractmethod

from csb.statistics.samplers import Trajectory
from csb.statistics.samplers.mc import State


class AbstractIntegrator(object):
    """
    Abstract integrator class. Subclasses implement different integration
    schemes for solving deterministic equations of motion.

    @param timestep: Integration timestep
    @type timestep: float

    @param gradient: Gradient of potential energy
    @type gradient: L{AbstractGradient}
    """
    
    __metaclass__ = ABCMeta

    def __init__(self, timestep, gradient):

        self._timestep = timestep
        self._gradient = gradient

    def integrate(self, init_state, length, return_trajectory=False):
        """
        Integrates equations of motion starting from an initial state a certain
        number of steps.

        @param init_state: Initial state from which to start integration
        @type init_state: L{State}
        
        @param length: Nubmer of integration steps to be performed
        @type length: int
        
        @param return_trajectory: Return complete trajectory instead of final
                                  state only
        @type return_trajectory: boolean

        @rtype: L{Trajectory}
        """
        
        result = []
        if return_trajectory:
            result.append(init_state)

        state = init_state
        for i in range(length):
            self.integrate_once(state, i)
            if return_trajectory:
                result.append(state)

        if not return_trajectory:
            result.append(state)

        return Trajectory(result)

    @abstractmethod
    def integrate_once(self, state, current_step):
        """
        Integrates one step starting from an initial state and an initial time
        given by the product of the timestep and the current_step parameter.

        @param state: State which to evolve one integration step
        @type state: L{State}
        
        @param current_step: Current integration step
        @type current_step: int
        
        @rtype: L{State}
        """
        pass

class LeapFrog(AbstractIntegrator):
    """
    Leap Frog integration scheme implementation.
    """
    
    def integrate_once(self, state, current_step):
        
        i = current_step
        if i == 0:
            self._oldgrad = self._gradient(state.position, 0.)
            
        momentumhalf = state.momentum - 0.5 * self._timestep * self._oldgrad
        state.position = state.position + self._timestep * momentumhalf
        self._oldgrad = self._gradient(state.position, (i + 1) * self._timestep)
        state.momentum = momentumhalf - 0.5 * self._timestep * self._oldgrad

        return state

class VelocityVerlet(AbstractIntegrator):
    """
    Velocity Verlet integration scheme implementation.
    """

    def integrate_once(self, state, current_step):
        
        i = current_step        
        if i == 0:
            self._oldgrad = self._gradient(state.position, 0.)
            
        state.position = state.position + self._timestep * state.momentum \
                         - 0.5 * self._timestep ** 2 * self._oldgrad
        newgrad = self._gradient(state.position, (i + 1) * self._timestep)
        state.momentum = state.momentum - 0.5 * self._timestep * (self._oldgrad + newgrad)
        self._oldgrad = newgrad

        return state

class AbstractGradient(object):
    """
    Abstract gradient class. Implementations evaluate the gradient of an energy
    function.
    """

    __metaclass__ = ABCMeta

    @abstractmethod
    def evaluate(self, q, t):
        """
        Evaluates the gradient at position q and time t.

        @param q: Position array
        @type q:  One-dimensional numpy array
        
        @param t: Time
        @type t: float
        
        @rtype: numpy array
        """
        pass

    def __call__(self, q, t):
        """
        Evaluates the gradient at position q and time t.

        @param q: Position array
        @type q:  One-dimensional numpy array
        
        @param t: Time
        @type t: float
        
        @rtype: numpy array
        """
        State.check_flat_array(q)
        return self.evaluate(q, t)
