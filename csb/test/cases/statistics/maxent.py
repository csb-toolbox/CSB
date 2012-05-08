
import numpy
import numpy.random

import csb.test as test
from scipy.optimize import fmin_powell

import csb
from csb.math import log_sum_exp
from csb.statistics.maxent import MaxentModel, MaxentPosterior

@test.functional
class TestMaxent(test.Case):

    def setUp(self):
        super(TestMaxent, self).setUp()
        self.data_fn = self.config.getTestFile('maxent.pickle')
        

    def testMaxent(self):
        k = 2
        data = csb.io.load(self.data_fn)
        model = MaxentModel(k)
        model.sample_weights()
        posterior = MaxentPosterior(model, data[:100000]/ 180. * numpy.pi  )

        x = model.get() * 1.

        x0 = posterior.model.get().flatten()
        target = lambda w: -posterior(w, n = 50)
        x = fmin_powell(target, x0, disp = False)

        self.assertTrue(x != None)
        self.assertTrue(len(x) == k * k * 4)

        posterior.model.set(x)
        posterior.model.normalize(True)

        xx = numpy.linspace(0 ,2 * numpy.pi, 500)
        fx = posterior.model.log_prob(xx,xx)

        self.assertTrue(fx != None)
        z = numpy.exp(log_sum_exp(numpy.ravel(fx))) 
        self.assertAlmostEqual(z * xx[1]**2, 1., places = 1)


if __name__ == '__main__':
    
    test.Console()
        
