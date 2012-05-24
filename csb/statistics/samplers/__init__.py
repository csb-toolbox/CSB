"""
Defines abstract samplers.
"""

import numpy as np
import csb.pyutils

from abc import ABCMeta, abstractmethod, abstractproperty


class DimensionError(TypeError):
    
    pass

class AbstractSampler(object):
    """
    Abstract interface for sampling algorithms.
    """
    
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def sample(self):
        """
        Draw a sample.
        @rtype: L{AbstractState}
        """
        pass

class AbstractState(object):
    """
    Represents a point in phase-space.
    """
    
    __metaclass__ = ABCMeta    
    
    @abstractproperty
    def position(self):
        pass
    
    @abstractproperty
    def momentum(self):
        pass
    
class State(AbstractState):
    """
    Represents a point in phase-space.
    """
    
    @staticmethod
    def check_flat_array(*args):
        """
        Check whether arguments are flat, one-dimensional numpy arrays.
        """
        
        for q in args:
            if not isinstance(q, np.ndarray):
                raise TypeError(q, 'numpy.ndarray expected!')
    
            if not len(q.squeeze().shape) <= 1:
                raise DimensionError(q, '1d numpy.ndarray expected!')
        
    @staticmethod
    def check_equal_length(q, p):
        """
        Check whether arguments have equal length.
        """
        
        if len(q) != len(p):
            raise DimensionError(p, 'momentum needs to have the same dimension as coordinates!')
    
    def __init__(self, position, momentum=None):
        
        self._position = None
        self._momenum = None

        self.position = position
        self.momentum = momentum

    @property
    def position(self):        
        return self._position.copy()
    @position.setter
    def position(self, value):        
        State.check_flat_array(value)        
        self._position = np.array(value)

    @property
    def momentum(self):
        return self._momentum.copy()
    @momentum.setter
    def momentum(self, value):
        if not value == None:
            State.check_flat_array(value)
            State.check_equal_length(value, self.position)
        self._momentum = np.array(value)
        
class EnsembleState(csb.pyutils.BaseCollectionContainer, AbstractState):
    """
    Defines an Ensemble Monte Carlo state; it is a read-only collection
    of State objects.

    @param items: initialization list of states
    @type items: list of L{States}
    """

    def __init__(self, items):   
        super(EnsembleState, self).__init__(items, type=State)
    
    @property
    def position(self):        
        return np.array([s.positions for s in self])

    @property
    def momentum(self):
        return np.array([s.momentum for s in self])

class Trajectory(csb.pyutils.CollectionContainer):
    """
    Ordered collection of states, representing a phase-space trajectory.

    @param items: initialization list
    @type items: list
    @param heat: heat produced during the trajectory
    @type heat: float
    """
    
    def __init__(self, items, heat=0.0):
        
        super(Trajectory, self).__init__(items)
        
        self._heat = None    
        self.heat = heat
    
    @property
    def heat(self):
        return self._heat
    @heat.setter
    def heat(self, value):
        self._heat = float(value)
        
    @property
    def last(self):
        return self._items[-1]

