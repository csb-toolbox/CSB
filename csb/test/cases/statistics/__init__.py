
import numpy
import numpy.random

import csb.test as test

from csb.statistics import Cumulative
from csb.statistics import kurtosis, skewness,autocorrelation


@test.functional
class TestStatFunction(test.Case):


    def testCumulative(self):
        from scipy.stats import norm
        
        x = numpy.linspace(-5.,5.,200)
        samples = numpy.random.normal(size=100000)
        cumula = Cumulative(samples)
        c = cumula(x)
        
        cx = norm.cdf(x)
        for i in range(199):
            self.assertWithinDelta(cx[i], c[i], delta = 1e-2)
        

    def testKurtosis(self):
        samples = numpy.random.normal(size=100000)
        self.assertWithinDelta(kurtosis(samples), 0., delta = 1e-1)

        samples = numpy.random.uniform(-2.,2.,size=100000)
        self.assertWithinDelta(kurtosis(samples), -1.2, delta = 1e-1)


    def testSkewness(self):
        samples = numpy.random.gamma(2.,0.5, size=100000)
        self.assertWithinDelta(skewness(samples), 2./numpy.sqrt(2.), delta = 1e-1)

    def testAutorcorrelation(self):
        x  = numpy.random.normal(size=1000) + numpy.sin(numpy.linspace(0., 2 * numpy.pi, 1000))
        n = 10
        ac = autocorrelation(x,n)
        self.assertWithinDelta(ac[0],1.)
        
    def testEntropy(self):
        pass

    def testCircvar(self):
        pass

    def testCircmean(self):
        pass

    
        

