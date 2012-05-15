import numpy
import numpy.random

import csb
import csb.test as test

from csb.statmech.ensembles import BoltzmannEnsemble
from csb.statmech.wham import WHAM, ImplicitWHAM

class FunnyGaussian(object):

    def __init__(self, d, k=100.):

        self.d = int(d)
	self.k = float(k)

    def sample(self, n_samples, inv_T=1):

	from numpy.random import standard_normal
	from numpy import sqrt, sum
	from csb.statistics.rand import truncated_gamma

	x = standard_normal((self.d, n_samples))
	x /= sqrt(sum(x**2,0))

	r = truncated_gamma(n_samples, 0.5 * self.d, self.k * inv_T, 0., 0.5)
	r = (2*r)**0.5

	return (x*r).T

    def energy(self, x):

	from numpy import sum
	x = numpy.array(x)
	return 0.5 * self.k * sum(x**2,-1)

    def log_Z(self, beta=1.):

	from numpy import pi
	from csb.math import log
	from scipy.special import gammainc, gammaln

	return log(0.5*self.d) + log(gammainc(0.5*self.d,0.5*self.k)) + \
	       gammaln(0.5*self.d) + (0.5*self.d) * (log(2) - log(self.k))

	
    def log_g(self, energies):

        from csb.math import log

        return (0.5*self.d - 1) * log(2*energies/self.k) + log(self.d/self.k)

 


@test.functional

class TestWHAM(test.Case):

    def setUp(self):
        self.betas = numpy.linspace(1e-5,1.,10)
	self.n = n = 1000
		
	gaussian = FunnyGaussian(10, 100.)

	self.samples = []
	self.raw_energies = []


	for beta in self.betas:
	    self.samples.append(gaussian.sample(n,beta))
	    self.raw_energies.append(gaussian.energy(self.samples[-1]))

	self.raw_energies = numpy.array(self.raw_energies)
	self.ensembles = [BoltzmannEnsemble(beta = beta) for beta in self.betas]

	self.log_z = gaussian.log_Z()
	self.log_g = gaussian.log_g(numpy.ravel(self.raw_energies))
	
    def testWHAM(self):

	w = WHAM(self.ensembles, 
		 numpy.ravel(self.raw_energies), 
		 numpy.array([self.n] * len(self.betas)))
	w.estimate()

	self.assertWithinDelta(numpy.dot(numpy.array([1,-1]),
					 w.log_z(numpy.array([1.,0.]))),
				self.log_z, delta = 0.5)



    def testImplicitWHAM(self):
	w = ImplicitWHAM(self.ensembles, 
			 numpy.ravel(self.raw_energies), 
			 [self.n] * len(self.betas))
	w.estimate()
	ens = [BoltzmannEnsemble(beta = 1.,),
	       BoltzmannEnsemble(beta = 0.)]
	self.assertWithinDelta(numpy.dot(numpy.array([1,-1]),
					 w.log_z(ensembles = ens)),
				self.log_z, delta = 0.5)


if __name__ == '__main__':

    test.Console()

