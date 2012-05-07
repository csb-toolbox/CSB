
import numpy
import numpy.random


import csb
import csb.statistics.ars
import csb.statistics.rand

from csb.math import log, exp
from csb.statistics import harmonic_mean, geometric_mean
from csb.statistics.pdf import AbstractEstimator, AbstractDensity, Gamma, InverseGamma

class ScaleMixturePrior(object):

    def get_estimator(self):
        """
        Returns an appropriate estimator for the scales of the mixture distribution
        under this prior
        """
        raise NotImplementedError()

        

class ScaleMixture(AbstractDensity):

    def __init__(self, scales = numpy.array([1.]), prior = None, d = 3):
        
        super(ScaleMixture, self).__init__()

        self._register('scales')
       
        self._d = d
        self.set_params(scales=scales)
        self._prior = prior

        if self._prior is not None:
            self.estimator = self._prior.get_estimator()

    @property
    def scales(self):
        return numpy.squeeze(self['scales'])

    @scales.setter
    def scales(self, value):
        if not isinstance(value, basestring) and \
           isinstance(value,(numpy.ndarray, list, tuple)):
            self['scales'] = numpy.array(value)
        else:
            raise ValueError("numpy array expected")

    @property
    def prior(self):
        return self._prior

    @prior.setter
    def prior(self, value):
        if not isinstance(value, ScaleMixturePrior):
            raise TypeError(value)
        self._prior = value



    def log_prob(self, x):
        from csb.math import log_sum_exp

        dim = self._d
        s = self.scales
        
        log_p = numpy.squeeze(-numpy.multiply.outer(x*x, 0.5*s)) + \
                numpy.squeeze(dim * 0.5  * (log(s)-log(2 * numpy.pi)))

        if self._prior is not None:
            log_p +=  numpy.squeeze(self._prior.log_prob(s)) 
        return log_sum_exp(log_p.T, 0)


    def random(self, shape = None):
        
        s = self.scales

        if shape is None:

            return numpy.random.normal() * s[numpy.random.sample(xrange(len(s)))]
        
        else:
            n = s.shape[0]

            nrv = numpy.random.normal(shape)
        
            indices = numpy.random.sample(xrange(n), shape)
             
            return s[indices] * nrv 



class ARSPosteriorAlpha(csb.statistics.ars.LogProb):
    """
    Posterior of alpha using a gamma distribution as
    prior for alpha
    """

    def __init__(self, a, b, n):

        self.a = float(a)
        self.b = float(b)
        self.n = float(n)

    def __call__(self, x):

        from scipy.special import gammaln, psi

        return self.a * x -\
               self.n * gammaln(numpy.clip(x, 1e-308, 1e308)) +\
               self.b * log(x), \
               self.a - self.n * psi(x) + self.b / x

    def initial_values(self, tol=1e-10):
        """
        Generate initial values by doing fixed point
        iterations to solve for alpha
        """
        from csb.math import psi, inv_psi, d_approx_psi

        n, a, b = self.n, self.a, self.b

        if abs(b) < 1e-10:
            alpha = inv_psi(a / n)
        else:
            alpha = 1.

        z = tol+1.

        while abs(z) > tol:

            z = n * psi(alpha) - \
                b / numpy.clip(alpha, 1e-300, 1e300) - a

            alpha -= z / (n * d_approx_psi(alpha) - b
                          / (alpha**2 + 1e-300))
            alpha = numpy.clip(alpha, 1e-100, 1e300)

        return numpy.clip(alpha - 1/(n + 1e-300), 1e-100, 1e300), \
               alpha + 1 / (n + 1e-300), alpha


class GammaPosteriorSampler(AbstractEstimator):

    def __init__(self):
        super(GammaPosteriorSampler, self).__init__()
        self.n_samples = 2


    def estimate(self, context, data):
        """
        Generates samples from the posterior of alpha and beta

        For beta the posterior is a gamma distribution and analytically
        acessible

        The posterior of alpha can not be expressed analytically and is
        aproximated using adaptive rejection sampling 
        """
        from csb.statistics.ars import ARS

        pdf = GammaPrior()

        ## sufficient statistics

        a = numpy.mean(data)
        b = exp(numpy.mean(log(data)))
        v = numpy.std(data)**2
        n = len(data)

        samples = []
        beta  = a / v
        alpha = beta * a

        for i in range(self.n_samples):

            ## sample beta from Gamma distribution
            beta = numpy.random.gamma(n*alpha+context._hyper_beta.alpha,
                                      1 / (n * a + context._hyper_beta.beta))

            ## sample alpha with ARS
            logp = ARSPosteriorAlpha(n * log(beta * b)\
                                     - context.hyper_alpha.beta,
                                     context.hyper_alpha.alpha-1., n)
            ars = ARS(logp)
            ars.initialize(logp.initial_values()[:2], z0 = 0.)
            alpha = ars.sample()

            if alpha is None:
                raise ValueError("ARS failed")

            samples.append((alpha, beta))

        pdf.alpha, pdf.beta = samples[-1]

        return pdf


class InvGammaPosteriorSampler(AbstractEstimator):

    def __init__(self):
        super(InvGammaPosteriorSampler, self).__init__()
        self.n_samples = 2


    def estimate(self, context, data):
        """
        Generates samples from the posterior of alpha and beta

        For beta the posterior is a gamma distribution and analytically
        acessible

        The posterior of alpha can not be expressed analytically and is
        aproximated using adaptive rejection sampling 
        """
        from csb.statistics.ars import ARS

        pdf = GammaPrior()

        ## sufficient statistics
        
        h = harmonic_mean(numpy.clip(data,1e-308,1e308))
        g = geometric_mean(numpy.clip(data,1e-308,1e308))

        n = len(data) 

        samples = []

        a = numpy.mean(1 / data)
        v = numpy.std(1 / data)**2

        beta  = a / v
        alpha = beta * a

        for i in range(self.n_samples):

            ## sample alpha with ARS

            logp = ARSPosteriorAlpha(n*(log(beta)-log(g))-context.hyper_alpha.beta,
                                     context.hyper_alpha.alpha-1., n)
            ars = ARS(logp)
            ars.initialize(logp.initial_values()[:2],z0=0.)

            alpha = ars.sample()

            if alpha is None:
                raise ValueError("Sampling failed")

            ## sample beta from Gamma distribution

            beta = numpy.random.gamma(n * alpha + context.hyper_beta.alpha, \
                                      1 / (n / h + context.hyper_beta.beta))

            samples.append((alpha,beta))

        pdf.alpha, pdf.beta = samples[-1]
        return pdf


class GammaScaleSampler(AbstractEstimator):


    def estimate(self, context, data):
        """
        The posterior distribution of the scales if we assume a gamma prior
        is again a gamma distribution with new parameters
        """
        pdf = ScaleMixture()
        alpha = context.prior.alpha
        beta = context.prior.beta
        d = context._d

        a = alpha + 0.5 * d
        b = beta  + 0.5 * data**2

        s = numpy.clip(numpy.random.gamma(a, 1./b), 1e-20, 1e10)
        pdf.scales = s

        context.prior.estimate(s)
        pdf.prior = context.prior

        return pdf

class InvGammaScaleSampler(AbstractEstimator):


    def estimate(self, context, data):
        """
        The posterior distribution of the scales if we assume a gamma prior
        is again a gamma distribution with new parameters
        """
        pdf = ScaleMixture()
        alpha = context.prior.alpha
        beta = context.prior.beta
        d = context._d

        p = -alpha + 1.5 
        b = 2*  beta 
        a = 1e-5 + data**2 

        s = csb.statistics.rand.gen_inv_gaussian(a, b, p)
        s = numpy.clip(s,1e-300,1e300)

        pdf.scales = s
        context.prior.estimate(s)
        pdf.prior = context.prior

        return pdf
    


class GammaPrior(Gamma,ScaleMixturePrior):
    """
    Gamma prior on mixture weights including a weak gamma prior on its parameter
    """

    def __init__(self, alpha = 1. , beta = 1., hyper_alpha=(1e-3, 1e-3),
                 hyper_beta = (1e-3, 1e-3)):

        super(GammaPrior, self).__init__(alpha,beta)

        self._hyper_alpha = Gamma(hyper_alpha[0],hyper_alpha[0])
        self._hyper_beta = Gamma(hyper_alpha[0],hyper_alpha[0])
        self.estimator = GammaPosteriorSampler()


    @property
    def hyper_beta(self):
        return self._hyper_beta

    @hyper_beta.setter
    def hyper_beta(self, value):
        if isinstance(value, AbstractDensity):
            self._hyper_beta = value
        else:
            raise ValueError(value)


    @property
    def hyper_alpha(self):
        return self._hyper_alpha

    @hyper_alpha.setter
    def hyper_alpha(self, value):
        if isinstance(value, AbstractDensity):
            self._hyper_beta = value
        else:
            raise ValueError(value)



    def log_prob(self, x):

        from scipy.special import gammaln
        
        a, b = self['alpha'], self['beta']

        l_a = self._hyper_alpha(a)
        l_b = self._hyper_beta(b)

        return super(GammaPrior, self).log_prob(x) + l_a + l_b


    def get_estimator(self):
        return GammaScaleSampler()


        
class InvGammaPrior(InverseGamma,ScaleMixturePrior):
    """
    Inverse Gamma prior on mixture weights including a weak gamma prior on its parameter
    """

    def __init__(self, alpha = 1. , beta = 1., hyper_alpha=(4., 1.),
                 hyper_beta = (2.,1.)):

        super(InvGammaPrior, self).__init__(alpha,beta)

        self._hyper_alpha = Gamma(hyper_alpha[0],hyper_alpha[0])
        self._hyper_beta = Gamma(hyper_alpha[0],hyper_alpha[0])
        self.estimator = InvGammaPosteriorSampler()


    @property
    def hyper_beta(self):
        return self._hyper_beta

    @hyper_beta.setter
    def hyper_beta(self, value):
        if isinstance(value, AbstractDensity):
            self._hyper_beta = value
        else:
            raise ValueError(value)


    @property
    def hyper_alpha(self):
        return self._hyper_alpha

    @hyper_alpha.setter
    def hyper_alpha(self, value):
        if isinstance(value, AbstractDensity):
            self._hyper_beta = value
        else:
            raise ValueError(value)



    def log_prob(self, x):

        from scipy.special import gammaln
        
        a, b = self['alpha'], self['beta']

        l_a = self._hyper_alpha(a)
        l_b = self._hyper_beta(b)

        return super(InvGammaPrior, self).log_prob(x) + l_a + l_b


    def get_estimator(self):
        return InvGammaScaleSampler()
        
         
