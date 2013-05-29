
import numpy
import numpy.random

import csb.test as test

from csb.statistics.pdf import ParameterValueError, EstimationFailureError
from csb.statistics.pdf import Normal, GeneralizedNormal, Laplace, Gamma, InverseGamma, GumbelMaximum, GumbelMinimum
from csb.statistics.pdf import MultivariateGaussian, Dirichlet, InverseGaussian, GeneralizedInverseGaussian


@test.regression
class TestDensityRegressions(test.Case):
    
    def testParameterValidation(self):
        """
        @see: [CSB 0000132]
        """
        self.assertRaises(ParameterValueError, InverseGaussian, -1, 1)       
        
        pdf = InverseGaussian()
        data = [ 3.1, 1.15, 6.86, 0.69, -1.73, -2.91, -2.2800744]
        
        def testMu(v):
            pdf['mu'] = v
            pdf.mu = v
            
        self.assertRaises(ParameterValueError, testMu, -1)        
        self.assertRaises(EstimationFailureError, pdf.estimate, data)        


@test.unit
class TestNormal(test.Case):    
    
    def testParameters(self):
    
        pdf = Normal(2.5, 1.5)
        pdf.mu = 33
        pdf.sigma = 44
        
        self.assertEqual(pdf.mu, 33)
        self.assertEqual(pdf.sigma, 44)
        self.assertEqual(pdf.sigma, pdf['sigma'])        
        
    def testLogProb(self):

        n = Normal(0, 1)

        self.assertAlmostEqual(n(0.), 1. / numpy.sqrt(2 * numpy.pi), delta=1e-1)
        self.assertAlmostEqual(n(1.), numpy.exp(-0.5) / numpy.sqrt(2 * numpy.pi), delta=1e-1)
        
    def testParameterEstimation(self):
        
        mu, sigma = 2.2, 0.3
        data = numpy.random.normal(mu, sigma, 100000)
        
        pdf = Normal(1, 1)
        pdf.estimate(data)
        
        self.assertAlmostEqual(pdf.mu, mu, delta=0.05)
        self.assertAlmostEqual(pdf.sigma, sigma, delta=0.05)

    def testRandom(self):
        
        pdf = Normal(4, 2)
        data = pdf.random(100000)
         
        self.assertAlmostEqual(numpy.mean(data), pdf.mu, delta=0.05)
        self.assertAlmostEqual(numpy.std(data), pdf.sigma, delta=0.05)    


@test.unit
class TestLaplace(test.Case):        
    
    def testParameters(self):
    
        pdf = Laplace(1.5, 2.5)
        pdf.mu = 33
        pdf.b = 44
        
        self.assertEqual(pdf.mu, 33)
        self.assertEqual(pdf.b, 44)
        self.assertEqual(pdf.b, pdf['b'])
        
        def propertyAssignment():
            pdf.b = -1
        def directAssignment():
            pdf['b'] = -1
        
        self.assertRaises(ParameterValueError, propertyAssignment)
        self.assertRaises(ParameterValueError, directAssignment)
        
    def testLogProb(self):

        pdf = Laplace(0, 2)
        self.assertEqual(pdf(0), 1.0 / (2 * pdf.b))
        self.assertEqual(pdf(1), 1 / (2 * pdf.b) * numpy.exp(-1 / pdf.b))
        
    def testParameterEstimation(self):
        
        mu, b = 2.2, 2
        data = numpy.random.laplace(mu, b, 100000)
        
        pdf = Laplace(1, 1)
        pdf.estimate(data)
        
        self.assertAlmostEqual(pdf.mu, mu, places=1)
        self.assertAlmostEqual(pdf.b, b, places=1)

    def testRandom(self):
        
        pdf = Laplace(4, 2)
        data = pdf.random(100000)
         
        self.assertAlmostEqual(numpy.median(data), pdf.mu, delta=0.05)   


@test.unit
class TestInverseGaussian(test.Case):        
    
    def testParameters(self):
    
        pdf = InverseGaussian(1.5, 2.5)
        pdf.mu = 3
        pdf.shape = 4
        
        self.assertEqual(pdf.mu, 3)
        self.assertEqual(pdf.shape, 4)
        self.assertEqual(pdf.mu, pdf['mu'])
        
        def propertyAssignment():
            pdf.shape = -1
        def directAssignment():
            pdf['shape'] = -1
        
        self.assertRaises(ParameterValueError, propertyAssignment)
        self.assertRaises(ParameterValueError, directAssignment)
        
    def testLogProb(self):

        ig = InverseGaussian(1., 1.)
        
        self.assertAlmostEqual(ig(1.), numpy.sqrt(1 / (2 * numpy.pi)), delta=1e-5)
        self.assertAlmostEqual(ig(2.), (numpy.sqrt(1 / (2 * numpy.pi * 2 ** 3)) * 
                                        numpy.exp(-1. / 4.)), delta=1e-5)
        
    def testParameterEstimation(self):
        
        mu, shape = 1, 3
        data = numpy.random.wald(mu, shape, 100000)
        
        pdf = InverseGaussian(1, 1)
        pdf.estimate(data)
        
        self.assertAlmostEqual(pdf.mu, mu, places=2)
        self.assertAlmostEqual(pdf.shape, shape, delta=0.05)

    def testRandom(self):
        
        ig = InverseGaussian(1., 1.)
        samples = ig.random(1000000)
        mu = numpy.mean(samples)
        var = numpy.var(samples)
        
        self.assertAlmostEqual(ig.mu, mu, delta=1e-1)
        self.assertAlmostEqual(ig.mu ** 3 / ig.shape, var, delta=1e-1)

        ig = InverseGaussian(3., 6.)

        samples = ig.random(1000000)
        mu = numpy.mean(samples)
        var = numpy.var(samples)
        
        self.assertAlmostEqual(ig.mu, mu, delta=1e-1)
        self.assertAlmostEqual(ig.mu ** 3 / ig.shape, var, delta=5e-1)


@test.unit
class TestGeneralizedNormal(test.Case):        
    
    def testParameters(self):
    
        pdf = GeneralizedNormal(0, 1, 1)
        pdf.mu = 11
        pdf.alpha = 2
        pdf.beta = 3
        
        self.assertEqual(pdf.mu, 11)
        self.assertEqual(pdf.alpha, 2)
        self.assertEqual(pdf.beta, 3)
        self.assertEqual(pdf.beta, pdf['beta'])
        
        def propertyAssignment():
            pdf.beta = -1
        def directAssignment(p):
            pdf[p] = -1
        
        self.assertRaises(ParameterValueError, propertyAssignment)
        self.assertRaises(ParameterValueError, directAssignment, 'alpha')
        self.assertRaises(ParameterValueError, directAssignment, 'beta')        
        
    def testLogProb(self):

        from scipy.special import gamma
        
        pdf = GeneralizedNormal(0, 1, 1)
        a = pdf.alpha        
        b = pdf.beta
        
        self.assertEqual(pdf(0), b / (2 * a * gamma(1 / b)))
        self.assertAlmostEqual(pdf(1), 1. / (2 * gamma(1)) * numpy.exp(-1), delta=1e-5)        
        
    def testParameterEstimation(self):
        
        mu, alpha, beta = -0.04, 2.11, 1.90
        data = [-1.1712, -2.5608, -0.7143, 2.6218, -2.0655, 0.7544, 1.208, -0.5289, 0.0045, 1.1746,
                -1.0766, 1.1198, 1.2785, -0.6051, 2.2913, -3.6672, -0.2525, 0.8782, -0.0617, -0.0239]
        
        pdf = GeneralizedNormal(1, 1, 1)
        pdf.estimate(data)
        
        self.assertAlmostEqual(pdf.mu, mu, places=2)
        self.assertAlmostEqual(pdf.alpha, alpha, places=2)
        self.assertAlmostEqual(pdf.beta, beta, places=1) 


@test.unit
class TestGeneralizedInverseGaussian(test.Case):        
    
    def testParameters(self):
    
        pdf = GeneralizedInverseGaussian(1, 1, 1)
        pdf.a = 2
        pdf.b = 3
        pdf.p = 4
        
        self.assertEqual(pdf.a, 2)
        self.assertEqual(pdf.b, 3)
        self.assertEqual(pdf.p, 4)
        self.assertEqual(pdf.p, pdf['p'])
        
        def propertyAssignment():
            pdf.a = -1
        def directAssignment(p):
            pdf[p] = -1
        
        self.assertRaises(ParameterValueError, propertyAssignment)
        self.assertRaises(ParameterValueError, directAssignment, 'a')
        self.assertRaises(ParameterValueError, directAssignment, 'b')
        self.assertRaises(ParameterValueError, directAssignment, 'p')        
        
    def testLogProb(self):

        fx = [2.7939e-24, 1.5305e-03, 2.3028e-02, 6.6302e-02, 1.1759e-01, 1.6810e-01,
              2.1365e-01, 2.5255e-01, 2.8438e-01, 3.0937e-01, 3.2806e-01, 3.4111e-01,
              3.4922e-01, 3.5307e-01, 3.5329e-01, 3.5048e-01, 3.4515e-01, 3.3778e-01,
              3.2878e-01, 3.1850e-01, 3.0726e-01, 2.9531e-01, 2.8289e-01, 2.7018e-01,
              2.5736e-01, 2.4454e-01, 2.3185e-01, 2.1937e-01, 2.0717e-01, 1.9532e-01,
              1.8385e-01, 1.7280e-01, 1.6220e-01, 1.5205e-01, 1.4236e-01, 1.3315e-01,
              1.2440e-01, 1.1611e-01, 1.0828e-01, 1.0088e-01, 9.3917e-02, 8.7363e-02,
              8.1207e-02, 7.5432e-02, 7.0021e-02, 6.4958e-02, 6.0224e-02, 5.5804e-02]
        
        x = numpy.arange(0.01, 5., 0.1)

        a = 2.
        b = 1.
        p = 2
        gig = GeneralizedInverseGaussian(a, b, p)
        fx2 = gig(x)

        for i in range(len(fx)):
            self.assertAlmostEqual(fx[i], fx2[i], delta=1e-1)

    def testRandom(self):
        
        from scipy.special import kv
        from numpy import sqrt

        a = 2.
        b = 1.
        p = 1
        gig = GeneralizedInverseGaussian(a, b, p)
        samples = gig.random(10000)

        mu_analytical = sqrt(b) * kv(p + 1, sqrt(a * b)) / (sqrt(a) * kv(p, sqrt(a * b)))

        var_analytical = b * kv(p + 2, sqrt(a * b)) / a / kv(p, sqrt(a * b)) - mu_analytical ** 2
        
        mu = numpy.mean(samples)
        var = numpy.var(samples)
        
        self.assertAlmostEqual(mu_analytical, mu, delta=1e-1)
        self.assertAlmostEqual(var_analytical, var, delta=1e-1)
        
        
@test.unit
class TestGamma(test.Case):        
    
    def testParameters(self):
    
        pdf = Gamma(1, 1)
        pdf.alpha = 2
        pdf.beta = 3
        
        self.assertEqual(pdf.alpha, 2)
        self.assertEqual(pdf.beta, 3)
        self.assertEqual(pdf.beta, pdf['beta'])
        
        def propertyAssignment():
            pdf.beta = -1
        def directAssignment(p):
            pdf[p] = -1
        
        self.assertRaises(ParameterValueError, propertyAssignment)
        self.assertRaises(ParameterValueError, directAssignment, 'alpha')
        self.assertRaises(ParameterValueError, directAssignment, 'beta')        
        
    def testLogProb(self):

        g = Gamma(1., 1.)
        self.assertAlmostEqual(g(1.), numpy.exp(-1.), delta=1e-5)
        self.assertAlmostEqual(g(2.), numpy.exp(-2.), delta=1e-5)

        g = Gamma(2., 1.)
        self.assertAlmostEqual(g(1.), 1. * numpy.exp(-1), delta=1e-5)
        self.assertAlmostEqual(g(2.), 2. * numpy.exp(-2), delta=1e-5)
        
    def testParameterEstimation(self):
        
        alpha = 0.1
        beta = 0.1

        data = numpy.random.gamma(alpha, 1. / beta, 10000)
        pdf = Gamma(1, 1)
        pdf.estimate(data)
        
        self.assertAlmostEqual(pdf.alpha, alpha, places=2)
        self.assertAlmostEqual(pdf.beta, beta, places=1)

    def testRandom(self):
        
        gamma = Gamma(1., 1.)
        samples = gamma.random(10000)
        mu = numpy.mean(samples)
        var = numpy.var(samples)
        
        self.assertAlmostEqual(gamma.alpha / gamma.beta, mu, delta=1e-1)
        self.assertAlmostEqual(gamma.alpha / gamma.beta ** 2, var, delta=1e-1)  
        

@test.unit
class TestInverseGamma(test.Case):        
    
    def testParameters(self):
    
        pdf = InverseGamma(1, 1)
        pdf.alpha = 2
        pdf.beta = 3
        
        self.assertEqual(pdf.alpha, 2)
        self.assertEqual(pdf.beta, 3)
        self.assertEqual(pdf.beta, pdf['beta'])
        
        def propertyAssignment():
            pdf.beta = -1
        def directAssignment(p):
            pdf[p] = -1
        
        self.assertRaises(ParameterValueError, propertyAssignment)
        self.assertRaises(ParameterValueError, directAssignment, 'alpha')
        self.assertRaises(ParameterValueError, directAssignment, 'beta')        
        
    def testLogProb(self):

        ig = InverseGamma(1., 1.)
        self.assertAlmostEqual(ig(1.), numpy.exp(-1.), delta=1e-5)
        self.assertAlmostEqual(ig(2.), 2 ** -2 * numpy.exp(-0.5), delta=1e-5)

        ig = InverseGamma(2., 1.)
        self.assertAlmostEqual(ig(1.), 1. * numpy.exp(-1), delta=1e-5)
        self.assertAlmostEqual(ig(2.), 2.** -3 * numpy.exp(-0.5), delta=1e-5)


@test.unit
class TestMultivariateGaussian(test.Case):    
    
    def testParameters(self):
    
        pdf = MultivariateGaussian(numpy.zeros(3), numpy.eye(3))
        d = 2
        pdf.mu = numpy.zeros(d)
        pdf.sigma = numpy.eye(d)
        
        self.assertEqual(list(pdf.mu), [0, 0])
        self.assertEqual(pdf.sigma.tolist(), numpy.eye(d).tolist())    
        
    def testLogProb(self):

        mg = MultivariateGaussian(numpy.zeros(2), numpy.eye(2))

        self.assertAlmostEqual(mg(numpy.zeros(2)), 1. / (2 * numpy.pi), delta=1e-1)
        self.assertAlmostEqual(mg(numpy.ones(2)), 1. / (2 * numpy.pi) * numpy.exp(-0.5), delta=1e-1)
        
    def testParameterEstimation(self):
        
        d = 3
        mu = numpy.ones(d)
        sigma = numpy.eye(d)

        pdf = MultivariateGaussian(mu, sigma)
        samples = pdf.random(100000)
        pdf.estimate(samples)

        for i in range(d):
            self.assertAlmostEqual(pdf.mu[i], mu[i], delta=1e-1)
            for j in range(d):
                self.assertAlmostEqual(pdf.sigma[i, j], sigma[i, j], delta=1e-1)
                

        d = 3
        mu = numpy.array([0., 1., 2.])
        sigma = numpy.random.random((d, d))
        sigma = numpy.dot(sigma, sigma.T)
        pdf = MultivariateGaussian(mu, sigma)
        samples = pdf.random(1000000)
        pdf.estimate(samples)

        for i in range(d):
            self.assertAlmostEqual(pdf.mu[i], mu[i], delta=1e-1)
            for j in range(d):
                self.assertAlmostEqual(pdf.sigma[i, j], sigma[i, j], delta=1e-1)
                
                
@test.unit
class TestDirichlet(test.Case):        
    
    def testParameters(self):
    
        pdf = Dirichlet([1, 1])
        
        alpha = numpy.array([2, 2])        
        pdf.alpha = numpy.array([2, 2])
        
        self.assertEqual(pdf.alpha.tolist(), alpha.tolist())
        self.assertEqual(pdf.alpha.tolist(), pdf['alpha'].tolist())
        
    def testLogProb(self):

        alpha = numpy.array([2, 2])
        pdf = Dirichlet(alpha)

        x = numpy.array([0.5, 0.5])
        self.assertAlmostEqual(pdf(x), 1.5, places=2)

        x = numpy.array([0.1, 0.9])
        self.assertAlmostEqual(pdf(x), .54, places=2)

        x = numpy.array([0.9, 0.1])
        self.assertAlmostEqual(pdf(x), .54, places=2)

        x = numpy.array([[0.9, 0.1],
                          [0.1, 0.9],
                          [0.5, 0.5]])

        self.assertAlmostEqual(pdf(x)[0], .54, places=2)
        self.assertAlmostEqual(pdf(x)[1], .54, places=2)
        self.assertAlmostEqual(pdf(x)[2], 1.5, places=2)
        
    def testParameterEstimation(self):
        
        alpha = numpy.array([1., 1., 1.])
        pdf = Dirichlet(alpha)

        self.assertAlmostEqual(pdf(numpy.array([1., 0., 0])),
                               pdf(numpy.array([0., 1., 0])),
                               delta=1e-1)
        
        self.assertAlmostEqual(pdf(numpy.array([1. / 3., 1. / 3., 1. / 3])),
                               pdf(numpy.array([0., 1., 0])),
                               delta=1e-1)

        alpha = numpy.array([16., 16., 16.])
        pdf = Dirichlet(alpha)
        

        self.assertAlmostEqual(pdf(numpy.array([1., 0., 0])),
                               pdf(numpy.array([0., 1., 0])),
                               delta=1e-1)

        self.assertAlmostEqual(pdf(numpy.array([1., 0., 0])),
                               pdf(numpy.array([0., 0., 1])),
                               delta=1e-1)
        
        self.assertNotEqual(pdf(numpy.array([1. / 3., 1. / 3., 1. / 3])),
                            pdf(numpy.array([0., 1., 0])))
        
        alpha = numpy.array([1., 1., 1.])
        pdf = Dirichlet(alpha)

        samples = pdf.random(100000)
        pdf.estimate(samples)
 
        for i in range(len(alpha)):
            self.assertAlmostEqual(pdf.alpha[i], alpha[i], delta=1e-1)

 
        alpha = numpy.array([200., 60., 20.])
        pdf = Dirichlet(alpha)

        samples = pdf.random(100000)
        pdf.estimate(samples)
 
        for i in range(len(alpha)):
            self.assertAlmostEqual(pdf.alpha[i], alpha[i], delta=1e1)
        
        
@test.unit
class TestGumbelMinimum(test.Case):
    
    @property
    def klass(self):
        return GumbelMinimum 
    
    def testParameters(self):
    
        pdf = self.klass(1, 1)
        pdf.mu = 2
        pdf.beta = 3
        
        self.assertEqual(pdf.mu, 2)
        self.assertEqual(pdf.beta, 3)
        self.assertEqual(pdf.beta, pdf['beta'])
        
        def propertyAssignment():
            pdf.beta = -1
        def directAssignment():
            pdf['beta'] = -1
        
        self.assertRaises(ParameterValueError, propertyAssignment)
        self.assertRaises(ParameterValueError, directAssignment)
        
    def testLogProb(self):

        pdf = self.klass(0, 2)
        self.assertAlmostEqual(pdf(0), 1.0 / (pdf.beta * numpy.e))
        
    def testParameterEstimation(self):
        
        mu, beta = 4, 2
        data = -numpy.random.gumbel(-mu, beta, 1000000)
        
        pdf = GumbelMinimum(1, 1)
        pdf.estimate(data)

        self.assertAlmostEqual(pdf.mu, mu, delta=0.01)        
        self.assertAlmostEqual(pdf.beta, beta, delta=0.01)

    def testRandom(self):
        
        pdf = self.klass(4, 2)
        data = pdf.random(100000)
         
        beta = numpy.sqrt(6 * numpy.var(data)) / numpy.pi
        self.assertAlmostEqual(pdf.beta, beta, delta=0.05)          
        
        
@test.unit
class TestGumbelMaximum(TestGumbelMinimum):
    
    @property
    def klass(self):
        return GumbelMaximum   
            
    def testParameterEstimation(self):
        
        mu, beta = 1, 2
        data = numpy.random.gumbel(mu, beta, 1000000)
        
        pdf = GumbelMaximum(1, 1)
        pdf.estimate(data)
        
        self.assertAlmostEqual(pdf.mu, mu, delta=0.01)
        self.assertAlmostEqual(pdf.beta, beta, delta=0.01)

            
if __name__ == '__main__':
    
    test.Console()
        
