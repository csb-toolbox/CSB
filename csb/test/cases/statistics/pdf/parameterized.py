import numpy
import csb.test as test

from csb.numeric import log

from csb.statistics.pdf.parameterized import ParameterizedDensity
from csb.statistics.pdf.parameterized import ParameterValueError, ParameterizationError
from csb.statistics.pdf.parameterized import AbstractParameter, Parameter, NonVirtualParameter



class Location(NonVirtualParameter):
    
    def _validate(self, value):
        return float(value)

class Scale(Parameter):    
       
    def _validate(self, value):
        return float(value)
    
    def _compute(self, base_value):
        
        if base_value == 0.0:
            return numpy.inf
        else:
            return 1.0 / base_value ** 0.5
        
    def bind_to(self, base):
        
        if base.name != "precision":
            raise ValueError(base)
                
        super(Scale, self).bind_to(base)
        
class DoubleScale(Parameter):    
       
    def _validate(self, value):
        return float(value)
    
    def _compute(self, base_value):
        return base_value * 2.0     
    
class Precision(Parameter):

    def _validate(self, value):
        
        if value < 0:
            raise ParameterValueError(self.name, value)
        return float(value)
    
                
class FancyGaussian(ParameterizedDensity):        
        
    def __init__(self, mu=0, precision=1):
        
        super(FancyGaussian, self).__init__()
        
        self._register('mu')
        self._register('sigma')
        self._register('precision')
        
        loc = Location(mu)
        prec = Precision(precision) 
        sigma = Scale(0)
        sigma.bind_to(prec)
        
        self.set_params(mu=loc, sigma=sigma, precision=prec)

    @property
    def mu(self):
        return self['mu'].value

    @property
    def sigma(self):
        return self['sigma'].value
        
    @property
    def precision(self):
        return self['precision'].value  
            
    def log_prob(self, x):

        mu = self.mu
        sigma = self.sigma
        
        return log(1.0 / numpy.sqrt(2 * numpy.pi * sigma ** 2)) - (x - mu) ** 2 / (2 * sigma ** 2)


@test.unit
class TestAbstractGenericParameter(test.Case):
    """
    Use AbstractParameter as a generic class which accepts values 
    of any type.
    """
    
    def setUp(self):
        
        class Value(object):
            pass
        
        class Param(AbstractParameter):
            
            def _validate(self, value):
                if not isinstance(value, Value):
                    raise TypeError(value)
                return value
            
        self.value = Value()        
        self.param = Param(self.value)
        
    def testValue(self):
        self.assertIs(self.param.value, self.value)
        
    def testSet(self):
        self.assertRaises(TypeError, self.param.set, 3)

        
@test.unit
class TestParameter(test.Case):
    """
    This is the main test case with complete coverage for AbstractParameter's
    methods and behavior. Covers also Parameter.
    
              computed -- leaf
             /
        base -- computed2
            \
             computed3
    """
    
    def setUp(self):        
        self.base = Precision(1.2)
        self.computed = Scale(100, base=self.base)
        self.computed2 = Scale(200, base=self.base)
        self.computed3 = Scale(300, base=self.base)
        self.leaf = DoubleScale(400, base=self.computed)
        
    def testConstrucor(self):
        # make sure newly constructed parameters are left in a consistent state 
        # to avoid unnecessary consistency updates 
        self.assertTrue(self.base._consistent)
        self.assertTrue(Scale(1)._consistent)

    def testName(self):
        
        self.assertEqual(self.base.name, "precision")
        self.assertEqual(self.computed.name, "scale")
        self.assertEqual(Scale(name="TesT").name, "TesT")
        
    def testValue(self):
        
        self.assertEqual(self.base.value, 1.2)
        self.assertEqual(self.computed.value, 1.0 / numpy.sqrt(self.base.value))
        self.assertEqual(self.computed2.value, 1.0 / numpy.sqrt(self.base.value))
        self.assertEqual(self.leaf.value, self.computed.value * 2)
        
        # turn self.base into a virtual parameter
        self.base.bind_to(Precision(12.2))
        self.assertEqual(self.base.value, 12.2)        
        
    def testIsVirtual(self):
        
        self.assertFalse(self.base.is_virtual)
        self.assertTrue(self.computed.is_virtual)
        
        self.base.bind_to(Precision(12.2))
        self.assertTrue(self.base.is_virtual)
        
    def testSet(self):
        
        base_initial_value = self.base._value
        
        # recompute all derivatives from the initial value of base
        self.assertEqual(self.computed._value, 100)
        self.leaf._ensure_consistency()
        self.computed2._ensure_consistency()
        self.computed3._ensure_consistency()
        
        # set self.base - it should remain consistent because it is not computed
        self.assertTrue(self.base._consistent)
        self.base.set(2.2)
        self.assertTrue(self.base._consistent)        
        self.assertEqual(self.base.value, 2.2)
        
        # self.computed and self.leaf should be inconsistent now that their base is updated
        self.assertFalse(self.computed._consistent)   
        self.assertFalse(self.leaf._consistent)
        self.assertEqual(self.computed._value, 1.0 / numpy.sqrt(base_initial_value))
        self.assertEqual(self.leaf._value, 2.0 / numpy.sqrt(base_initial_value))
        # retrieving self.computed's value should trigger updates up to self.computed
        recomputed = self.computed.value
        self.assertTrue(self.computed._consistent)
        self.assertEqual(recomputed, 1.0 / numpy.sqrt(self.base._value))
        # self.leaf is still inconsistent
        self.assertFalse(self.leaf._consistent)
        self.assertEqual(self.leaf._value, 2.0 / numpy.sqrt(base_initial_value))
        self.assertIs(self.leaf._nearest_consistent_base()[-1], self.computed)
        # until we request its value       
        recomputed = self.leaf.value
        self.assertTrue(self.leaf._consistent)
        self.assertEqual(recomputed, 2.0 / numpy.sqrt(self.base._value))
        self.assertEqual(recomputed, 2.0 * self.computed._value)
        # make sure the other two branches are still inconsistent
        initial_value = 1.0 / numpy.sqrt(base_initial_value)
        self.assertEqual(self.computed2._value, initial_value)
        self.assertEqual(self.computed3._value, initial_value)
        # until they get used
        recomputed = self.computed2.value
        self.assertTrue(self.computed2._consistent)
        self.assertEqual(recomputed, 1.0 / numpy.sqrt(self.base._value))               
        
        # attempt to set self.computed - not allowed
        self.assertRaises(ParameterizationError, self.computed.set, 2)

        # attempt to set a negative Precision        
        self.assertRaises(ParameterValueError, self.base.set, -2)
        # attempt to assigned non-float - not allowed in the Parameter specialization        
        self.assertRaises(ParameterValueError, Parameter().set, object())
        
    def testBindTo(self):
                
        # can't bind self.base to itself
        self.assertRaises(ParameterizationError, self.base.bind_to, self.base)
        # deeper circular dependency 
        self.assertRaises(ParameterizationError, self.base.bind_to, self.computed)        

        # self.base is not virtual and therefore must be consistent        
        self.assertTrue(self.base._consistent)
        
        # make it virtual - should get inconsistent now
        self.base.bind_to(Precision(12.2))
        self.assertFalse(self.base._consistent)
        self.assertTrue(self.base.is_virtual)
        
        # retrieving its value should trigger the consistency cascade
        self.assertEqual(self.base.value, 12.2)
        self.assertTrue(self.base._consistent)
        
    def testFindBaseParameter(self):

        self.assertIs(self.base.find_base_parameter(), self.base)        
        self.assertIs(self.computed.find_base_parameter(), self.base)
    

@test.unit
class TestNonVirtualParameter(test.Case):
    """
    Make sure explicit NonVirtualParameter-s are updatable and 
    refuse binding requests
    """
    
    def setUp(self):
        self.param = Location()
        
    def testConstructor(self):
        base = Parameter()
        self.assertRaises(ParameterizationError, lambda: Location(base=base))
        
    def testIsVirtual(self):
        self.assertFalse(self.param.is_virtual)
        
    def testBindTo(self):
        base = Parameter()
        self.assertRaises(ParameterizationError, self.param.bind_to, base)
        
    def testSet(self):
        self.param.set(22)
        self.assertEqual(self.param.value, 22)
        

@test.unit
class TestParameterizedDensity(test.Case):
    
    def setUp(self):
        self.pdf = FancyGaussian(2, 5)
        
    def testConstructor(self):
        
        class Density(ParameterizedDensity):
            
            def __init__(self, p):
                super(Density, self).__init__()                
                self._register('p')
                self.set_params(p=p)
                
            def log_prob(self, x):
                return x
        
        self.assertRaises(TypeError, Density, 4)
    
    def testProperties(self):
        
        self.assertEqual(self.pdf.mu, 2)
        self.assertEqual(self.pdf.precision, 5)
        self.assertAlmostEqual(self.pdf.sigma, 0.4472, places=3)
        
    def testParameterChaining(self):
        
        self.assertEqual(self.pdf.precision, 5)
        self.assertAlmostEqual(self.pdf.sigma, 0.4472, places=3)
        
        self.pdf['precision'].set(2)
        
        self.assertEqual(self.pdf.precision, 2)
        self.assertAlmostEqual(self.pdf.sigma, 0.7071, places=3)
    
    def testAssignment(self):
        
        self.pdf['sigma'] = Scale(55)
        self.assertEqual(self.pdf.sigma, 55)
        self.assertEqual(self.pdf['sigma'].name, 'scale')
        
        def assign(i):
            self.pdf['sigma'] = i
            
        self.assertRaises(TypeError, assign, 55) 

    

if __name__ == '__main__':
    test.Console()        
    
