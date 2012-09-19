"""
provides various integration schemes and an abstract gradient class.
"""

import numpy

from abc import ABCMeta, abstractmethod
from csb.statistics.samplers.mc import State, TrajectoryBuilder


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

    def integrate(self, init_state, length, mass_matrix=None, return_trajectory=False):
        """
        Integrates equations of motion starting from an initial state a certain
        number of steps.

        @param init_state: Initial state from which to start integration
        @type init_state: L{State}
        
        @param length: Nubmer of integration steps to be performed
        @type length: int

        @param mass_matrix: Mass matrix
        @type mass_matrix:  n-dimensional numpy array with n being the dimension
                            of the configuration space, that is, the dimension of
                            the position / momentum vectors

        @param return_trajectory: Return complete L{Trajectory} instead of the initial
                                  and final states only (L{PropagationResult}). This reduces
                                  performance.
        @type return_trajectory: boolean

        @rtype: L{AbstractPropagationResult}
        """

        d = len(init_state.position)

        self._mass_matrix = mass_matrix
        
        if self._mass_matrix == None:
            self._mass_matrix = numpy.eye(d)
            self._inverse_mass_matrix = mass_matrix
        else:
            self._mass_matrix = numpy.dot(self._mass_matrix, numpy.eye(d))
            self._inverse_mass_matrix = numpy.linalg.inv(self._mass_matrix)
        
        builder = TrajectoryBuilder.create(full=return_trajectory)
            
        builder.add_initial_state(init_state)
        state = init_state.clone()
        
        for i in range(length - 1):
            state = self.integrate_once(state, i, inverse_mass_matrix=self._inverse_mass_matrix)
            builder.add_intermediate_state(state)

        state = self.integrate_once(state, length - 1, inverse_mass_matrix=self._inverse_mass_matrix)
        builder.add_final_state(state)

        return builder.product

    @abstractmethod
    def integrate_once(self, state, current_step, inverse_mass_matrix=None):
        """
        Integrates one step starting from an initial state and an initial time
        given by the product of the timestep and the current_step parameter.
        The input C{state} is changed in place.

        @param state: State which to evolve one integration step
        @type state: L{State}
        
        @param current_step: Current integration step
        @type current_step: int

        @param inverse_mass_matrix: Inverse mass matrix
        @type inverse_mass_matrix:  n-dimensional numpy array with n being the dimension
                                    of the configuration space, that is, the dimension of
                                    the position / momentum vectors
        @return: the altered state
        @rtype: L{State}
        """
        pass

class LeapFrog(AbstractIntegrator):
    """
    Leap Frog integration scheme implementation that calculates position and
    momenta at equal times. Slower than FastLeapFrog, but intermediate points
    in trajectories obtained using
    LeapFrog.integrate(init_state, length, return_trajectoy=True) are physical.
    """
    
    def integrate_once(self, state, current_step, inverse_mass_matrix=None):

        if inverse_mass_matrix == None:
                inverse_mass_matrix = 1.
        
        i = current_step
        
        if i == 0:
            self._oldgrad = self._gradient(state.position, 0.)
            
        momentumhalf = state.momentum - 0.5 * self._timestep * self._oldgrad
        state.position = state.position + self._timestep * numpy.dot(inverse_mass_matrix, momentumhalf)
        self._oldgrad = self._gradient(state.position, (i + 1) * self._timestep)
        state.momentum = momentumhalf - 0.5 * self._timestep * self._oldgrad

        return state

class FastLeapFrog(LeapFrog):
    """
    Leap Frog integration scheme implementation that calculates position and
    momenta at unequal times by concatenating the momentum updates of two
    successive integration steps.
    WARNING: intermediate points in trajectories obtained by
    FastLeapFrog.integrate(init_state, length, return_trajectories=True)
    are NOT to be interpreted as phase-space trajectories, because
    position and momenta are not given at equal times! In the initial and the
    final state, positions and momenta are given at equal times.
    """

    def integrate(self, init_state, length, mass_matrix=None, return_trajectory=False):
		
        d = len(init_state.position)
        
        if mass_matrix == None:
            mass_matrix = numpy.eye(d)
            inverse_mass_matrix = mass_matrix
        else:
            mass_matrix = numpy.dot(mass_matrix, numpy.eye(d))
            inverse_mass_matrix = numpy.linalg.inv(mass_matrix)

        builder = TrajectoryBuilder.create(full=return_trajectory)
            
        builder.add_initial_state(init_state)
        state = init_state.clone()
        
        state.momentum = state.momentum - 0.5 * self._timestep * self._gradient(state.position, 0.)
        
        for i in range(length-1):
            state.position = state.position + self._timestep * numpy.dot(inverse_mass_matrix, state.momentum)
            state.momentum = state.momentum - self._timestep * \
                             self._gradient(state.position, (i + 1) * self._timestep)
            builder.add_intermediate_state(state)

        state.position = state.position + self._timestep * numpy.dot(inverse_mass_matrix, state.momentum)
        state.momentum = state.momentum - 0.5 * self._timestep * \
                         self._gradient(state.position, length * self._timestep)
        builder.add_final_state(state)
        
        return builder.product

class VelocityVerlet(AbstractIntegrator):
    """
    Velocity Verlet integration scheme implementation.
    """

    def integrate_once(self, state, current_step, inverse_mass_matrix=None):

        if inverse_mass_matrix == None:
            inverse_mass_matrix = 1.

        i = current_step
        
        if i == 0:
            self._oldgrad = self._gradient(state.position, 0.)
            
        state.position = state.position + self._timestep * numpy.dot(inverse_mass_matrix, state.momentum) \
                         - 0.5 * self._timestep ** 2 * numpy.dot(inverse_mass_matrix, self._oldgrad)
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
