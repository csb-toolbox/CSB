"""
Probability density functions.

This module defines L{AbstractDensity}: a common interface for all PDFs.
Each L{AbstractDensity} describes a specific type of probability distribution,
for example L{Normal} is an implementation of the Gaussian distribution:

    >>> pdf = Normal(mu=10, sigma=1.1)
    >>> pdf.mu, pdf['sigma']
    10.0, 1.1

Every PDF provides an implementation of the L{AbstractDensity.evaluate} 
method, which evaluates the PDF for a list of input data points:

    >>> pdf.evaluate([10, 9, 11, 12])
    array([ 0.3626748 ,  0.2399147 ,  0.2399147 ,  0.06945048])
    
PDF instances also behave like functions:
    
    >>> pdf(data)    # the same as pdf.evaluate(data)
    
Some L{AbstractDensity} implementations may support drawing random numbers from
the distribution (or raise an exception otherwise):

    >>> pdf.random(2)
    array([ 9.86257083,  9.73760515])
    
Each implementation of L{AbstractDensity} may support infinite number of estimators,
used to estimate and re-initialize the PDF parameters from a set of observed data
points:

    >>> pdf.estimate([5, 5, 10, 10])
    >>> pdf.mu, pdf.sigma
    (7.5, 2.5)
    >>> pdf.estimator
    <csb.statistics.pdf.GaussianMLEstimator>
    
Estimators implement the L{AbstractEstimator} interface. They are treated as
pluggable tools, which can be exchanged through the L{AbstractDensity.estimator}
property (you could create, initialize and plug your own estimator as well).
This is a classic Strategy pattern.  
"""

import numpy.random
import scipy.special
import csb.core

from abc import ABCMeta, abstractmethod
from csb.core import OrderedDict

from csb.numeric import log, exp, psi, inv_psi
from scipy.special import gammaln
from numpy import array, fabs, power, sqrt, pi, mean, median, clip


class IncompatibleEstimatorError(TypeError):
    pass

class ParameterNotFoundError(AttributeError):
    pass

class ParameterValueError(ValueError):
    
    def __init__(self, param, value):
                        
        self.param = param
        self.value = value

        super(ParameterValueError, self).__init__(param, value)


class AbstractEstimator(object):
    """
    Density parameter estimation strategy.
    """
    
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def estimate(self, context, data):
        """
        Estimate the parameters of the distribution from same {data}.
        
        @param context: context distribution
        @type context: L{AbstractDensity}
        @param data: sample values
        @type data: array
        
        @return: a new distribution, initialized with the estimated parameters
        @rtype: L{AbstractDensity}
        """
        pass
       
class NullEstimator(AbstractEstimator):
    """
    Does not estimate anything.
    """
    def estimate(self, context, data):
        raise NotImplementedError()

class LaplaceMLEstimator(AbstractEstimator):
    
    def estimate(self, context, data):
         
        x = array(data)
        
        mu = median(x)
        b = mean(fabs(x - mu))
        
        return Laplace(b, mu)
    
class GaussianMLEstimator(AbstractEstimator):
    
    def estimate(self, context, data):
         
        x = array(data)
        
        mu = mean(x)
        sigma = sqrt(mean((x - mu) ** 2))
        
        return Normal(mu, sigma)
    
class InverseGaussianMLEstimator(AbstractEstimator):
    
    def estimate(self, context, data):
         
        x = array(data)
        
        mu = mean(x)
        il = mean((1.0 / x) - (1.0 / mu))
        
        return InverseGaussian(mu, 1.0 / il)

class GammaMLEstimator(AbstractEstimator):

    def __init__(self):
        super(GammaMLEstimator, self).__init__()
        self.n_iter = 1000
        

    def estimate(self, context, data):
        
        mu = mean(data)
        logmean = mean(log(data))

        a = 0.5 / (log(mu) - logmean)

        for dummy in range(self.n_iter):

            a = inv_psi(logmean - log(mu) + log(a))

        return Gamma(a, a / mu)

class GenNormalBruteForceEstimator(AbstractEstimator):
    
    def __init__(self, minbeta=0.5, maxbeta=8.0, step=0.1):
        
        self._minbeta = minbeta
        self._maxbeta = maxbeta
        self._step = step
        
        super(GenNormalBruteForceEstimator, self).__init__()
        
    def estimate(self, context, data):
        
        pdf = GeneralizedNormal(1, 1, 1)
        data = array(data)
        logl = []
        
        for beta in numpy.arange(self._minbeta, self._maxbeta, self._step):
            
            self.update(pdf, data, beta)
            
            l = pdf.log_prob(data).sum()       
            logl.append([beta, l])
            
        logl = numpy.array(logl)
        
        # optimal parameters:
        beta = logl[ numpy.argmax(logl[:, 1]) ][0]
        self.update(pdf, data, beta)
        
        return pdf
    
    def estimate_with_fixed_beta(self, data, beta):
        
        mu = median(data)
        v = mean((data - mu) ** 2)
        alpha = sqrt(v * exp(gammaln(1. / beta) - gammaln(3. / beta)))
    
        return mu, alpha
    
    def update(self, pdf, data, beta):
        
        mu, alpha = self.estimate_with_fixed_beta(data, beta)        
        
        pdf.mu = mu
        pdf.alpha = alpha
        pdf.beta = beta
        
        return pdf

class MultivariateGaussianMLEstimator(AbstractEstimator):

    def __init__(self):
        super(MultivariateGaussianMLEstimator, self).__init__()

    def estimate(self, context, data):
        return MultivariateGaussian(numpy.mean(data, 0), numpy.cov(data.T))
    

class DirichletEstimator(AbstractEstimator):

    def __init__(self):
        super(DirichletEstimator, self).__init__()
        self.n_iter = 1000
        self.tol = 1e-5

    def estimate(self, context, data):

        log_p = numpy.mean(log(data), 0)
        
        e = numpy.mean(data, 0)
        v = numpy.mean(data ** 2, 0)
        q = (e[0] - v[0]) / (v[0] - e[0] ** 2)

        a = e * q
        y = a * 0
        k = 0
        while(sum(abs(y - a)) > self.tol and k < self.n_iter):
            y = psi(sum(a)) + log_p
            a = numpy.array(list(map(inv_psi, y)))
            k += 1 

        return Dirichlet(a)
        
            
class AbstractDensity(object):
    """
    Defines the interface and common operations for all probability density
    functions.
    
    Subclasses must complete the implementation by implementing the
    L{AbstractDensity.log_prob} method. Subclasses could also consider--but
    are not obliged to--override the L{AbstractDensity.random} method. 
    """

    __metaclass__ = ABCMeta
    

    def __init__(self):

        self._params = OrderedDict()
        self._estimator = None
        
        self.estimator = NullEstimator()

    def __getitem__(self, param):
        
        if param in self._params: 
            return self._params[param]
        else:
            raise ParameterNotFoundError(param)
        
    def __setitem__(self, param, value):
        
        if param in self._params: 
            if csb.core.iterable(value):
                value = array(value)
            else:
                value = float(value)
            
            self._validate(param, value)
            self._params[param] = value
        else:
            raise ParameterNotFoundError(param)
        
    @property
    def estimator(self):
        return self._estimator
    @estimator.setter
    def estimator(self, strategy):
        if not isinstance(strategy, AbstractEstimator):
            raise TypeError(strategy)
        self._estimator = strategy

    def __call__(self, x):
        return self.evaluate(x)

    def __str__(self):
        name = self.__class__.__name__
        params = ', '.join([ '{0}={1}'.format(p, v) for p, v in self._params.items() ])
        
        return '{0}({1})'.format(name, params)           
        
    def _register(self, name):
        """
        Register a new parameter name.
        """
        if name not in self._params:
            self._params[name] = None
            
    def _validate(self, param, value):
        """
        Parameter value validation hook.
        """
        pass

    def get_params(self):
        return [self._params[name] for name in  self.parameters]
    
    def set_params(self, *values, **named_params):
        
        for p, v in zip(self.parameters, values):
            self[p] = v
            
        for p in named_params:
            self[p] = named_params[p]
    
    @property
    def parameters(self):
        """
        Get a list of all distribution parameter names.
        """
        return tuple(self._params)

    @abstractmethod
    def log_prob(self, x):
        """
        Evaluate the logarithm of the probability of observing values C{x}.

        @param x: values
        @type x: array
        @rtype: array        
        """
        pass
    
    def evaluate(self, x):
        """
        Evaluate the probability of observing values C{x}.
        
        @param x: values
        @type x: array        
        @rtype: array
        """
        x = numpy.array(x)
        return exp(self.log_prob(x))      
    
    def random(self, size=None):
        """
        Generate random samples from the probability distribution.
        
        @param size: number of values to sample
        @type size: int
        """
        raise NotImplementedError()

    def estimate(self, data):
        """
        Estimate and load the parameters of the distribution from sample C{data}
        using the current L{AbstractEstimator} strategy.
        
        @param data: sample values
        @type data: array
                
        @raise NotImplementedError: when no estimator is available for this
                                    distribution
        """
        pdf = self.estimator.estimate(self, data)

        try:
            for param in pdf.parameters:
                self[param] = pdf[param]
                    
        except ParameterNotFoundError:
            raise IncompatibleEstimatorError(self.estimator)

class Laplace(AbstractDensity):
        
    def __init__(self, b, mu):
        
        super(Laplace, self).__init__()
        
        self._register('b')
        self._register('mu')
        
        self.set_params(b=b, mu=mu)
        self.estimator = LaplaceMLEstimator()
        
    def _validate(self, param, value):
        
        if param == 'b' and value < 0:
            raise ParameterValueError(param, value)
        
    @property
    def b(self):
        return self['b']
    @b.setter
    def b(self, value):
        self['b'] = value

    @property
    def mu(self):
        return self['mu']
    @mu.setter
    def mu(self, value):
        self['mu'] = value
            
    def log_prob(self, x):

        b = self.b
        mu = self.mu
        
        return log(1 / (2. * b)) - fabs(x - mu) / b

    def random(self, size=None):
        
        loc = self.mu
        scale = self.b
        
        return numpy.random.laplace(loc, scale, size)
    
class Normal(AbstractDensity):
    
    def __init__(self, mu=0, sigma=1):
        
        super(Normal, self).__init__()
        
        self._register('mu')
        self._register('sigma')
        
        self.set_params(mu=mu, sigma=sigma)
        self.estimator = GaussianMLEstimator()

    @property
    def mu(self):
        return self['mu']
    @mu.setter
    def mu(self, value):
        self['mu'] = value

    @property
    def sigma(self):
        return self['sigma']
    @sigma.setter
    def sigma(self, value):
        self['sigma'] = value
            
    def log_prob(self, x):

        mu = self.mu
        sigma = self.sigma
        
        return log(1.0 / sqrt(2 * pi * sigma ** 2)) - (x - mu) ** 2 / (2 * sigma ** 2)
    
    def random(self, size=None):
        
        mu = self.mu
        sigma = self.sigma
                
        return numpy.random.normal(mu, sigma, size)

class InverseGaussian(AbstractDensity):

    def __init__(self, mu=1., llambda=1.):

        super(InverseGaussian, self).__init__()

        self._register('mu')
        self._register('llambda')

        self.set_params(mu=mu, llambda=llambda)
        self.estimator = InverseGaussianMLEstimator()

    @property
    def mu(self):
        return self['mu']
    @mu.setter
    def mu(self, value):
        if value <= 0.:
            raise ValueError("Mean mu should be greater than 0")
        self['mu'] = value

    @property
    def llambda(self):
        return self['llambda']
    @llambda.setter
    def llambda(self, value):
        if value <= 0.:
            raise ValueError("Shape Parameter lambda should be greater than 0")
        self['llambda'] = value
            
    def log_prob(self, x):

        mu = self.mu
        _lambda = self.llambda

        y = -0.5 * _lambda * (x - mu) ** 2 / (mu ** 2 * x)
        z = 0.5 * (log(_lambda) - log(2 * pi * x ** 3))
        return  z + y 


    def random(self, size=None):

        mu = self.mu
        _lambda = self.llambda

        mu_2l = mu / _lambda / 2.
        Y = numpy.random.standard_normal(size)
        Y = mu * Y ** 2
        X = mu + mu_2l * (Y - sqrt(4 * _lambda * Y + Y ** 2))
        U = numpy.random.random(size)

        m = numpy.less_equal(U, mu / (mu + X))

        return m * X + (1 - m) * mu ** 2 / X
 
class GeneralizedNormal(AbstractDensity):
    
    def __init__(self, mu, alpha, beta):
        
        super(GeneralizedNormal, self).__init__()
        
        self._register('mu')
        self._register('alpha')
        self._register('beta')
        
        self.set_params(mu=mu, alpha=alpha, beta=beta)
        self.estimator = GenNormalBruteForceEstimator()

    @property
    def mu(self):
        return self['mu']
    @mu.setter
    def mu(self, value):
        self['mu'] = value

    @property
    def alpha(self):
        return self['alpha']
    @alpha.setter
    def alpha(self, value):
        self['alpha'] = value
        
    @property
    def beta(self):
        return self['beta']
    @beta.setter
    def beta(self, value):
        self['beta'] = value
            
    def log_prob(self, x):

        mu = self.mu
        alpha = self.alpha
        beta = self.beta
             
        return log(beta / (2.0 * alpha)) - gammaln(1. / beta) - power(fabs(x - mu) / alpha, beta)

class GeneralizedInverseGaussian(AbstractDensity):

    def __init__(self, a=1., b=1., p=1.):
        super(GeneralizedInverseGaussian, self).__init__()

        self._register('a')
        self._register('b')
        self._register('p')
        self.set_params(a=a, b=b, p=p)

        self.estimator = NullEstimator()

    @property
    def a(self):
        return self['a']

    @a.setter
    def a(self, value):
        if value <= 0:
            raise ValueError("Parameter a is nonnegative")
        else:
            self['a'] = value

    @property
    def b(self):
        return self['b']

    @b.setter
    def b(self, value):
        if value <= 0:
            raise ValueError("Parameter b is nonnegative")
        else:
            self['b'] = value

    @property
    def p(self):
        return self['p']

    @p.setter
    def p(self, value):
        if value <= 0:
            raise ValueError("Parameter p is nonnegative")
        else:
            self['p'] = value

    def log_prob(self, x):

        a = self['a']
        b = self['b']
        p = self['p']

        lz = 0.5 * p * (log(a) - log(b)) - log(2 * scipy.special.kv(p, sqrt(a * b)))

        return lz + (p - 1) * log(x) - 0.5 * (a * x + b / x)
        
    def random(self, size=None):

        from csb.statistics.rand import inv_gaussian
        
        rvs = []
        burnin = 10
        a = self['a']
        b = self['b']
        p = self['p']

        s = a * 0. + 1.

        if p < 0:
            a, b = b, a

        if size == None:
            size = 1
        for i in range(int(size)):
            for j in range(burnin):

                l = b + 2 * s
                m = sqrt(l / a)

                x = inv_gaussian(m, l, shape=m.shape)
                s = numpy.random.gamma(abs(p) + 0.5, x)

            if p >= 0:
                rvs.append(x)
            else:
                rvs.append(1 / x)

        return numpy.array(rvs)
        
class Gamma(AbstractDensity):

    def __init__(self, alpha=1, beta=1):
        super(Gamma, self).__init__()

        self._register('alpha')
        self._register('beta')

        self.set_params(alpha=alpha, beta=beta)
        self.estimator = GammaMLEstimator()

    @property
    def alpha(self):
        return self['alpha']
    @alpha.setter
    def alpha(self, value):
        self['alpha'] = value
        
    @property
    def beta(self):
        return self['beta']

    @beta.setter
    def beta(self, value):
        self['beta'] = value

    def log_prob(self, x):
            
        a, b = self['alpha'], self['beta']

        return a * log(b) - gammaln(clip(a, 1e-308, 1e308)) + \
               (a - 1) * log(clip(x, 1e-308, 1e308)) - b * x

    def random(self, size=None):
        return numpy.random.gamma(self['alpha'], 1 / self['beta'], size)

class InverseGamma(AbstractDensity):

    def __init__(self, alpha=1, beta=1):
        super(InverseGamma, self).__init__()

        self._register('alpha')
        self._register('beta')

        self.set_params(alpha=alpha, beta=beta)
        self.estimator = NullEstimator()

    @property
    def alpha(self):
        return self['alpha']

    @alpha.setter
    def alpha(self, value):
        self['alpha'] = value
        
    @property
    def beta(self):
        return self['beta']

    @beta.setter
    def beta(self, value):
        self['beta'] = value

    def log_prob(self, x):
        a, b = self['alpha'], self['beta']
        return a * log(b) - gammaln(a) - (a + 1) * log(x) - b / x

    def random(self, size=None):
        return 1. / numpy.random.gamma(self['alpha'], 1 / self['beta'], size)
    
class MultivariateGaussian(Normal):

    def __init__(self, mu=numpy.zeros(2), sigma=numpy.eye(2)):
                
        super(MultivariateGaussian, self).__init__(mu, sigma)
        self.estimator = MultivariateGaussianMLEstimator()
        
    def random(self, size=None):
        return numpy.random.multivariate_normal(self.mu, self.sigma, size)

    def log_prob(self, x):

        from numpy.linalg import det
        
        mu = self.mu
        S = self.sigma
        D = len(mu)
        q = self.__q(x)
        return -0.5 * (D * log(2 * pi) + log(abs(det(S)))) - 0.5 * q ** 2

    def __q(self, x):
        from numpy import sum, dot, reshape
        from numpy.linalg import inv

        mu = self.mu
        S = self.sigma
        
        return sqrt(clip(sum(reshape((x - mu) * dot(x - mu, inv(S).T.squeeze()), (-1, len(mu))), -1), 0., 1e308))

    def conditional(self, x, dims):
        """
        Returns the distribution along the dimensions
        dims conditioned on x

        @param x: conditional values
        @param dims: new dimensions
        """
        from numpy import take, dot
        from numpy.linalg import inv

        dims2 = [i for i in range(self['mu'].shape[0]) if not i in dims]

        mu1 = take(self['mu'], dims)
        mu2 = take(self['mu'], dims2)

        # x1 = take(x, dims)
        x2 = take(x, dims2)

        A = take(take(self['Sigma'], dims, 0), dims, 1)
        B = take(take(self['Sigma'], dims2, 0), dims2, 1)
        C = take(take(self['Sigma'], dims, 0), dims2, 1)

        mu = mu1 + dot(C, dot(inv(B), x2 - mu2))
        Sigma = A - dot(C, dot(inv(B), C.T))
        
        return MultivariateGaussian((mu, Sigma))

class Dirichlet(AbstractDensity):

    def __init__(self, alpha):
        super(Dirichlet, self).__init__()

        self._register('alpha')

        self.set_params(alpha=alpha)
        self.estimator = DirichletEstimator()

    @property
    def alpha(self):
        return self['alpha']

    @alpha.setter
    def alpha(self, value):
        self['alpha'] = numpy.ravel(value)

    def log_prob(self, x):
        #TODO check wether x is in the probability simplex
        alpha = self.alpha
        return gammaln(sum(alpha)) - sum(gammaln(alpha)) \
              + numpy.dot((alpha - 1).T, log(x).T) 
        
    def random(self, size=None):
        return numpy.random.mtrand.dirichlet(self.alpha, size)
