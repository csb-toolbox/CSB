
import numpy
import numpy.random

import csb.test as test

from csb.statistics.pdf import Normal, GeneralizedNormal, Laplace, Gamma, InverseGamma
from csb.statistics.pdf import MultivariateGaussian, Dirichlet, InverseGaussian

@test.functional
class TestLogProb(test.Case):

    def testDirichlet(self):

        alpha = numpy.array([2,2])
        pdf = Dirichlet(alpha)

        x =  numpy.array([0.5,0.5])
        self.assertAlmostEqual(pdf(x), 1.5 , places=2)

        x =  numpy.array([0.1,0.9])
        self.assertAlmostEqual(pdf(x), .54 , places=2)

        x =  numpy.array([0.9,0.1])
        self.assertAlmostEqual(pdf(x), .54 , places=2)

        x =  numpy.array([[0.9,0.1],
                          [0.1,0.9],
                          [0.5,0.5]])

        self.assertAlmostEqual(pdf(x)[0], .54 , places=2)
        self.assertAlmostEqual(pdf(x)[1], .54 , places=2)
        self.assertAlmostEqual(pdf(x)[2], 1.5 , places=2)


    def testNormal(self):
        n = Normal(0.,1.)

        self.assertWithinDelta(n(0.), 1./numpy.sqrt( 2 * numpy.pi))
        self.assertWithinDelta(n(1.), numpy.exp(-0.5)/numpy.sqrt( 2 * numpy.pi) )
        

    def testInverseGaussian(self):
        ig = InverseGaussian(1.,1.)
        
        self.assertWithinDelta(ig(1.), numpy.sqrt(1 / (2 * numpy.pi)), delta=1e-5)
        self.assertWithinDelta(ig(2.), (numpy.sqrt(1 / (2 * numpy.pi * 2**3)) *
                                        numpy.exp(- 1. / 4.))
                               , delta=1e-5)


    def testGeneralizedNormal(self):
        from scipy.special import gamma

        gn = GeneralizedNormal(0.,1.,1.)

        self.assertWithinDelta(gn(0.), 1./2./gamma(1) , delta=1e-5)

        self.assertWithinDelta(gn(1.), 1./2./gamma(1) * numpy.exp(-1) , delta=1e-5)


    def testGamma(self):
        g = Gamma(1.,1.)
        self.assertWithinDelta(g(1.), numpy.exp(-1.) , delta=1e-5)
        self.assertWithinDelta(g(2.), numpy.exp(-2.) , delta=1e-5)

        g = Gamma(2.,1.)
        self.assertWithinDelta(g(1.),  1. * numpy.exp(-1) , delta=1e-5)
        self.assertWithinDelta(g(2.),  2. * numpy.exp(-2) , delta=1e-5)


    def testInverseGamma(self):
        ig = InverseGamma(1.,1.)
        self.assertWithinDelta(ig(1.), numpy.exp(-1.) , delta=1e-5)
        self.assertWithinDelta(ig(2.), 2**-2  * numpy.exp(-0.5) , delta=1e-5)

        ig = InverseGamma(2.,1.)
        self.assertWithinDelta(ig(1.), 1. * numpy.exp(-1) , delta=1e-5)
        self.assertWithinDelta(ig(2.), 2.**-3 * numpy.exp(-0.5) , delta=1e-5)
        

    def testMultivariateGaussian(self):
        mg = MultivariateGaussian(numpy.zeros(2),numpy.eye(2))

        self.assertWithinDelta(mg(numpy.zeros(2)), 1./( 2 * numpy.pi))
        self.assertWithinDelta(mg(numpy.ones(2)), 1./( 2 * numpy.pi) * numpy.exp(-0.5) )
    

        
        
@test.functional
class TestRandom(test.Case):

    def testInverseGaussian(self):

        ig = InverseGaussian(1.,1.)
        samples = ig.random(1000000)
        mu = numpy.mean(samples)
        var = numpy.var(samples)
        
        self.assertWithinDelta(ig.mu, mu, delta=1e-1)
        self.assertWithinDelta(ig.mu**3/ig.llambda, var, delta=1e-1)

        ig = InverseGaussian(3.,6.)

        samples = ig.random(1000000)
        mu = numpy.mean(samples)
        var = numpy.var(samples)
        
        self.assertWithinDelta(ig.mu, mu, delta=1e-1)
        self.assertWithinDelta(ig.mu**3/ig.llambda, var, delta=5e-1)

    def testGamma(self):

        gamma = Gamma(1., 1.)
        samples = gamma.random(10000)
        mu = numpy.mean(samples)
        var = numpy.var(samples)
        
        self.assertWithinDelta(gamma.alpha/gamma.beta, mu, delta=1e-1)
        self.assertWithinDelta(gamma.alpha/gamma.beta**2, var, delta=1e-1)
        
@test.functional
class TestParameterEstimation(test.Case):
        
    def testNormal(self):
        
        mu, sigma = 2.2, 0.3
        data = numpy.random.normal(mu, sigma, 100000)
        
        pdf = Normal(1, 1)
        pdf.estimate(data)
        
        self.assertAlmostEqual(pdf.mu, mu, places=2)
        self.assertAlmostEqual(pdf.sigma, sigma, places=2)

    def testGeneralizedNormal(self):
        
        mu, alpha, beta = -0.04, 2.11, 1.90
        data = [-1.1712, -2.5608, -0.7143, 2.6218, -2.0655, 0.7544, 1.208, -0.5289, 0.0045, 1.1746, 
                -1.0766, 1.1198, 1.2785, -0.6051, 2.2913, -3.6672, -0.2525, 0.8782, -0.0617, -0.0239]
        
        pdf = GeneralizedNormal(1, 1, 1)
        pdf.estimate(data)
        
        self.assertAlmostEqual(pdf.mu, mu, places=2)
        self.assertAlmostEqual(pdf.alpha, alpha, places=2)
        self.assertAlmostEqual(pdf.beta, beta, places=1)
        
    def testLaplace(self):
        
        mu, b = 2.2, 2
        data = numpy.random.laplace(mu, b, 100000)
        
        pdf = Laplace(1, 1)
        pdf.estimate(data)
        
        self.assertAlmostEqual(pdf.mu, mu, places=1)
        self.assertAlmostEqual(pdf.b, b, places=1)

    def testGamma(self):

        alpha = 0.1
        beta = 0.1

        data = numpy.random.gamma(alpha, 1. / beta, 10000)
        pdf = Gamma(1, 1)
        pdf.estimate(data)
        
        self.assertAlmostEqual(pdf.alpha, alpha, places=2)
        self.assertAlmostEqual(pdf.beta, beta, places=1)

    def testMultivariateGaussian(self):
        d = 3
        mu = numpy.ones(d)
        sigma = numpy.eye(d)

        pdf = MultivariateGaussian(mu, sigma)
        samples = pdf.random(100000)
        pdf.estimate(samples)

        for i in range(d):
            self.assertWithinDelta(pdf.mu[i], mu[i], delta=1e-1)
            for j in range(d):
                self.assertWithinDelta(pdf.sigma[i,j], sigma[i,j], delta=1e-1)
                

        d = 3
        mu = numpy.array([0., 1., 2.])
        sigma = numpy.random.random((d,d))
        sigma = numpy.dot(sigma, sigma.T)
        pdf = MultivariateGaussian(mu, sigma)
        samples = pdf.random(1000000)
        pdf.estimate(samples)

        for i in range(d):
            self.assertWithinDelta(pdf.mu[i], mu[i], delta=1e-1)
            for j in range(d):
                self.assertWithinDelta(pdf.sigma[i,j], sigma[i,j], delta=1e-1)

        
    def testDirichlet(self):

        alpha = numpy.array([1.,1.,1.])
        pdf = Dirichlet(alpha)

        self.assertWithinDelta(pdf(numpy.array([1.,0.,0])),
                               pdf(numpy.array([0.,1.,0])),
                               delta = 1e-1)
        
        self.assertWithinDelta(pdf(numpy.array([1./3.,1./3.,1./3])),
                               pdf(numpy.array([0.,1.,0])),
                               delta = 1e-1)

        alpha = numpy.array([16.,16.,16.])
        pdf = Dirichlet(alpha)
        

        self.assertWithinDelta(pdf(numpy.array([1.,0.,0])),
                               pdf(numpy.array([0.,1.,0])),
                               delta = 1e-1)

        self.assertWithinDelta(pdf(numpy.array([1.,0.,0])),
                               pdf(numpy.array([0.,0.,1])),
                               delta = 1e-1)
        
        self.assertNotEqual(pdf(numpy.array([1./3.,1./3.,1./3])),
                            pdf(numpy.array([0.,1.,0])))
        
        alpha = numpy.array([1.,1.,1.])
        pdf = Dirichlet(alpha)

        samples = pdf.random(100000)
        pdf.estimate(samples)
 
        for i in range(len(alpha)):
            self.assertWithinDelta(pdf.alpha[i], alpha[i], delta = 1e-1)

 
        alpha = numpy.array([200.,60.,20.])
        pdf = Dirichlet(alpha)

        samples = pdf.random(100000)
        pdf.estimate(samples)
 
        for i in range(len(alpha)):
            self.assertWithinDelta(pdf.alpha[i], alpha[i], delta = 1e1)

            
if __name__ == '__main__':
    
    test.Console()
        
