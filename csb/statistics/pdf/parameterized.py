"""
Probability density functions with support for shared and computed parameters.

This module extends the functionality of L{csb.statistics.pdf} with a specialized 
and more sophisticated L{AbstractDensity} -- the L{ParameterizedDensity}, which 
works with L{AbstractParameter} objects rather than simple floats. 

Each L{AbstractParameter} holds two properties - L{AbstractParameter.name} and
L{AbstractParameter.value}:

    >>> class Sigma(AbstractParameter):
    >>>     def _validate(self, value):
    >>>         return float(value)     
    >>>     def _compute(self, base_value):                
    >>>         return 1.0 / base_value ** 0.5
    >>>                            
    >>> sigma = Sigma(3)
    >>> sigma.name, sigma.value
    sigma, 3

L{AbstractParameter}s holding a single float value are indistinguishable from 
the simple float parameters used in L{csb.statistics.pdf.BaseDensity}. 
However, each L{AbstractParameter} supports the concept of binding:

    >>> sigma.is_virtual
    False
    >>> precision = Precision(1) 
    >>> sigma.bind_to(precision)
    >>> sigma.is_virtual
    True    
    >>> precision.set(2)  # triggers an implicit, lazy update in sigma
    >>> sigma.set(1)
    ParameterizationError: Virtual parameters can't be updated explicitly
    
The instance of Sigma is now a virtual parameter which receives automatic updates 
from another base parameter using a predefined rule (L{AbstractParameter._compute}).
This is a lazy, on demand process. As soon as Sigma's computed value is 
requested (via C{s.value}), Sigma will query the parameter it depends on 
(Precision), which in turn will get recomputed first based on its own base, etc. 
Thus, the L{AbstractParameter} model supports a parameter dependency chain with 
linear or tree-like topologies::
            
               sigma -- y   
              /              
    precision -- sigma2 -- x 
    
In this graph precision is a base (non-virtual) parameter and sigma, sigma2, x, and y
are all virtual (computed). Binding precision to another parameter will immediately
turn it into a virtual one. However, cycles are not allowed (e.g. precision can't 
be bound to sigma2 or x) and each virtual parameter must have exactly one base.   

Computed parameters can then be used to implement custom PDFs with dependent 
parameters within one PDF or spanning multiple PDFs.
"""

import csb.statistics.pdf as pdf

from abc import abstractmethod


class ParameterizationError(ValueError):
    pass

class ParameterValueError(pdf.ParameterValueError):
    pass


class ParameterizedDensity(pdf.AbstractDensity):
    """
    Base abstract class for all PDFs, which operate on simple or computed 
    (chained) parameters. Parameters of type different from L{AbstractParameter}
    will trigger TypeError-s.
    """
    
    def _set(self, param, value):
        
        if not isinstance(value, AbstractParameter):
            raise TypeError(value)
        
        super(ParameterizedDensity, self)._set(param, value)    
    
    
class AbstractParameter(object):
    """
    Abstract parameterization, which can exist independently or be coupled 
    to other parameters upon request. Virtual/coupled/derived parameters cannot
    be overwritten explicitly, but their values will get recomputed once their
    corresponding base parameters get updated. This is a lazy process - parameter
    recalculation happens only when an out of date parameter is requested. 
    
    @param value: initial value (defaults to None / AbstractParameter.NULL)
    @type value: object
    @param name: name of parameter (this is the name of the class by default)
    @type name: str
    @param base: optional base parameter to compute this instance from
    @type base: L{AbstractParameter}
    """
    
    NULL = None
    
    def __init__(self, value=NULL, name=None, base=None):
        
        self._derivatives = set()
        self._base = None
        self._consistent = True
        
        if name is None:
            name = self.__class__.__name__.lower()
        
        self._name = str(name)
        self._value = AbstractParameter.NULL
        
        self._update(value)
        
        if base is not None:
            self.bind_to(base)
            
    #def __repr__(self):
    #    return "<Parameter {0._name}={0._value}>".format(self)
    
    @property
    def name(self):
        """
        Parameter name
        """
        return self._name
    
    @property
    def value(self):
        """
        Parameter value (guaranteed to be up to date)
        """
        self._ensure_consistency()
        return self._value
    
    @property
    def is_virtual(self):
        """
        True if this parameter is virtual (computed)
        """
        return self._base is not None

    def set(self, value):
        """
        Update the value of this parameter. This is not possible for 
        virtual parameters.
        
        @param value: new value
        @type value: object
        
        @raise ParameterizationError: if this is a virtual parameter
        @raise ParameterValueError: on invalid value
        """
        if self.is_virtual:
            raise ParameterizationError(
                            "Virtual parameters can't be updated explicitly")
            
        self._update(value)
        
        self._invalidate()
        self._consistent = True

    def bind_to(self, parameter):
        """
        Bind the current parameter to a base parameter. This converts 
        the current parameter to a virtual one, whose value will get 
        implicitly updated to be consistent with its base.
        
        Note that virtual parameters must have exactly one base; computing a
        parameter from multiple bases is not allowed. Cycles are also not 
        allowed; the topology must always stay a tree with a non-virtual
        parameter at the root.        
        
        @param parameter: base parameter to compute this instance from
        @param parameter: L{AbstractParameter}
        
        @raise ParameterizationError: if this parameter is already virtual
        @raise ParameterizationError: on attempt to create a circular dependency
                                        
        """

        if not isinstance(parameter, AbstractParameter):
            raise TypeError(parameter)
        
        if parameter.find_base_parameter() is self:
            raise ParameterizationError("Circular dependency detected")
        
        if self.is_virtual:
            msg = "Parameter {0.name} is already bound to {1.name}"
            raise ParameterizationError(msg.format(self, self._base))
        
        self._set_base(parameter)
        self._base._add_derived(self)
        
        self._invalidate()
        
    def _set_base(self, parameter):
        self._base = parameter
        
    def _add_derived(self, parameter):
        self._derivatives.add(parameter) 
        
    def _invalidate(self):
        """
        Mark self and its virtual children as inconsistent
        """
        for p in self._derivatives: 
            p._invalidate()
            
        self._consistent = False                  
            
    def _update(self, value):
        """
        Overwrite the current value of the parameter. This triggers
        an abstract (custom) validation hook, but has no side effects 
        (i.e. it doesn't propagate!)
        """
        sanitized = self._validate(value)
        self._value = sanitized
        
    @abstractmethod
    def _validate(self, value):
        """
        Validate and sanitize the specified value before assignment.
        @return: sanitized value
        
        @raise ParameterValueError: on invalid value
        """
        return value
    
    def _compute(self, base_value):
        """
        Compute a new value for the current parameter given the value
        of a base parameter (assuming self.is_virtual). By default this returns
        the value of the base parameter (i.e. self just inherits the value 
        of its base untouched).
        """
        return base_value
            
    def _ensure_consistency(self, force=False):
        """
        Make sure that the current value is up to date. If it isn't,
        trigger a real-time cascaded update from the non-virtual root 
        to all virtual nodes. Also mark all nodes consistent in the course of
        doing this update. 
        """        
        if self._consistent and not force:
            return
        
        root = self.find_base_parameter()
        root._recompute_derivatives()
        
    def _recompute_derivatives(self):
        """
        Recompute all derived parameters starting from self and mark 
        them consistent.
        """
        if self.is_virtual:
            value = self._compute(self._base._value)
            self._update(value)
        
        self._consistent = True
        
        for p in self._derivatives:
            p._recompute_derivatives()        
        
    def find_base_parameter(self):
        """
        Find and return the non-virtual base parameter that is the root
        of the current hierarchy. If self is not virtual, return self.
        
        @return: base parameter
        @rtype: L{AbstractParameter}
        """
        root = self
        
        while root.is_virtual:
            root = root._base
            
        return root 


class Parameter(AbstractParameter):
    """
    Default parameter implementation which accepts float values only.
    """

    def __init__(self, value=0.0, name=None, base=None):
        super(Parameter, self).__init__(value, name, base)
            
    def _validate(self, value):
        
        try:
            return float(value)
        except (ValueError, TypeError):
            raise ParameterValueError(self.name, value)
        
        
class NonVirtualParameter(Parameter):
    """
    A float L{Parameter} that is explicitly non-computed and cannot be 
    bound to another L{Parameter}.
    """
    
    def bind_to(self, parameter):
        raise ParameterizationError(
                            "This parameter is explicitly non-computed")
    
    @property
    def is_virtual(self):
        return False
        






