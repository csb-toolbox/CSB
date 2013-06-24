"""
Defines abstract samplers.
"""

import numpy as np
import csb.core

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
        self._momentum = None

        self.position = position
        self.momentum = momentum

    def __eq__(self, other):

        return self.position == other.position and self.momentum == other.momentum

    @property
    def position(self):        
        return self._position.copy()
    @position.setter
    def position(self, value):        
        State.check_flat_array(value)        
        self._position = np.array(value)

    @property
    def momentum(self):
        if self._momentum is None:
            return None
        else:
            return self._momentum.copy()
    @momentum.setter
    def momentum(self, value):
        if not value is None:
            State.check_flat_array(value)
            State.check_equal_length(value, self.position)
            self._momentum = np.array(value)
        else:
            self._momentum = None
        
    def clone(self):
        if self.momentum is not None:
            return self.__class__(self.position.copy(), self.momentum.copy())
        else:
            return self.__class__(self.position.copy())
        
        
class EnsembleState(csb.core.BaseCollectionContainer, AbstractState):
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
        return np.array([s.position for s in self])

    @property
    def momentum(self):
        return np.array([s.momentum for s in self])
