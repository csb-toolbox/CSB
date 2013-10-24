"""
Approximation of a distribution as a mixture of gaussians with a zero mean but different sigmas
"""

import numpy.random

import csb.core
import csb.statistics.ars
import csb.statistics.rand

from csb.core import typedproperty
from abc import abstractmethod, ABCMeta

from csb.numeric import log, exp, approx_psi, inv_psi, d_approx_psi
from scipy.special import  psi, kve
from csb.statistics import harmonic_mean, geometric_mean
from csb.statistics.pdf import AbstractEstimator, AbstractDensity, BaseDensity
from csb.statistics.pdf import Gamma, InverseGamma, NullEstimator


def inv_digamma_minus_log(y, tol=1e-10, n_iter=100):
    """
    Solve y = psi(alpha) - log(alpha) for alpha by fixed point
    integration.
    """
    if y >= -log(6.):
        x = 1 / (2 * (1 - exp(y)))
    else:
        x = 1.e-10
    for _i in range(n_iter):
        z = approx_psi(x) - log(x) - y
        if abs(z) < tol:
            break
        x -= x * z / (x * d_approx_psi(x) - 1)
        x = abs(x)
    return x

class ScaleMixturePriorEstimator(AbstractEstimator):
    """
    Prior on the scales of a L{ScaleMixture}, which determines how the scales 
    are estimated.
    """
    
    __metaclass__ = ABCMeta
    
    @abstractmethod
    def get_scale_estimator(self):
        """
        Return an appropriate estimator for the scales of the mixture distribution
        under this prior.
        """
        pass
    

class ScaleMixturePrior(object):
    """
    Prior on the scales of a L{ScaleMixture}, which determines how the scales 
    are estimated.
    """

    def __init__(self, *args):
        super(ScaleMixturePrior, self).__init__(*args)
        self._scale = None
        self._scale_estimator = NullEstimator()
        
    @property
    def scale_estimator(self):
        return self._scale_estimator

    @property
    def estimator(self):
        return self._estimator

    @estimator.setter
    def estimator(self, strategy):
        if not isinstance(strategy, AbstractEstimator):
            raise TypeError(strategy)
        self._estimator = strategy
        if isinstance(strategy, ScaleMixturePriorEstimator):
            self._scale_estimator = strategy.get_scale_estimator()
        else:
            self._scale_estimator = NullEstimator()


class ScaleMixture(BaseDensity):
    """
    Robust probabilistic superposition and comparison of protein structures
    Martin Mechelke and Michael Habeck

    Represenation of a distribution as a mixture of gaussians with a mean of 
    zero and different inverse variances/scales. The number of scales equals 
    the number of datapoints.

    The underlying family is determined by a prior L{ScaleMixturePrior} on the
    scales. Choosing a L{GammaPrior} results in Stundent-t posterior, wheras a 
    L{InvGammaPrior} leads to a K-Distribution as posterior.
    """

    def __init__(self, scales=numpy.array([1., 1.]), prior=None, d=3):
        
        super(ScaleMixture, self).__init__()

        self._register('scales')
       
        self._d = d
        self.set_params(scales=scales)
        self._prior = prior

        if self._prior is not None:
            self.estimator = self._prior.scale_estimator

    @property
    def scales(self):
        return numpy.squeeze(self['scales'])

    @scales.setter
    def scales(self, value):
        if not isinstance(value, csb.core.string) and \
           isinstance(value, (numpy.ndarray, list, tuple)):
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
        self.estimator = self._prior.scale_estimator

    def log_prob(self, x):
        from csb.numeric import log_sum_exp

        dim = self._d
        s = self.scales
        
        log_p = numpy.squeeze(-numpy.multiply.outer(x * x, 0.5 * s)) + \
                numpy.squeeze(dim * 0.5 * (log(s) - log(2 * numpy.pi)))

        if self._prior is not None:
            log_p += numpy.squeeze(self._prior.log_prob(s)) 
        return log_sum_exp(log_p.T, 0)


    def random(self, shape=None):
        
        s = self.scales

        if shape is None:
            return numpy.random.normal() * s[numpy.random.randint(len(s))]
        
        else:
            #n = s.shape[0]
            nrv = numpy.random.normal(size=shape)
            indices = numpy.random.randint(len(s), size=shape)
             
            return s[indices] * nrv 


class ARSPosteriorAlpha(csb.statistics.ars.LogProb):
    """
    This class represents the posterior distribution of the alpha parameter
    of the Gamma and Inverse Gamma prior, and allows sampling using adaptive
    rejection sampling L{ARS}.
    """

    def __init__(self, a, b, n):

        self.a = float(a)
        self.b = float(b)
        self.n = float(n)

    def __call__(self, x):

        from scipy.special import gammaln

        return self.a * x - \
               self.n * gammaln(numpy.clip(x, 1e-308, 1e308)) + \
               self.b * log(x), \
               self.a - self.n * psi(x) + self.b / x

    def initial_values(self, tol=1e-10):
        """
        Generate initial values by doing fixed point
        iterations to solve for alpha
        """
        n, a, b = self.n, self.a, self.b

        if abs(b) < 1e-10:
            alpha = inv_psi(a / n)
        else:
            alpha = 1.

        z = tol + 1.

        while abs(z) > tol:

            z = n * psi(alpha) - \
                b / numpy.clip(alpha, 1e-300, 1e300) - a

            alpha -= z / (n * d_approx_psi(alpha) - b
                          / (alpha ** 2 + 1e-300))
            alpha = numpy.clip(alpha, 1e-100, 1e300)

        return numpy.clip(alpha - 1 / (n + 1e-300), 1e-100, 1e300), \
               alpha + 1 / (n + 1e-300), alpha


class GammaPosteriorSampler(ScaleMixturePriorEstimator):

    def __init__(self):
        super(GammaPosteriorSampler, self).__init__()
        self.n_samples = 2

    def get_scale_estimator(self):
        return GammaScaleSampler()
    
    def estimate(self, context, data):
        """
        Generate samples from the posterior of alpha and beta.

        For beta the posterior is a gamma distribution and analytically
        acessible.

        The posterior of alpha can not be expressed analytically and is
        aproximated using adaptive rejection sampling.
        """
        pdf = GammaPrior()

        ## sufficient statistics

        a = numpy.mean(data)
        b = exp(numpy.mean(log(data)))
        v = numpy.std(data) ** 2
        n = len(data)

        beta = a / v
        alpha = beta * a
        samples = []

        for _i in range(self.n_samples):

            ## sample beta from Gamma distribution
            beta = numpy.random.gamma(n * alpha + context._hyper_beta.alpha,
                                      1 / (n * a + context._hyper_beta.beta))

            ## sample alpha with ARS
            logp = ARSPosteriorAlpha(n * log(beta * b)\
                                     - context.hyper_alpha.beta,
                                     context.hyper_alpha.alpha - 1., n)
            ars = csb.statistics.ars.ARS(logp)
            ars.initialize(logp.initial_values()[:2], z0=0.)
            alpha = ars.sample()

            if alpha is None:
                raise ValueError("ARS failed")

            samples.append((alpha, beta))

        pdf.alpha, pdf.beta = samples[-1]

        return pdf

class GammaPosteriorMAP(ScaleMixturePriorEstimator):

    def __init__(self):
        super(GammaPosteriorMAP, self).__init__()

    def get_scale_estimator(self):
        return GammaScaleMAP()

    def estimate(self, context, data):
        """
        Estimate alpha and beta from their posterior
        """

        pdf = GammaPrior()

        s = data[0].mean()
        y = data[1].mean() - log(s) - 1.

        alpha = abs(inv_digamma_minus_log(numpy.clip(y,
                                                     - 1e308,
                                                     - 1e-300)))
        beta = alpha / s

        
        pdf.alpha, pdf.beta = alpha, beta
        return pdf


        
        
class InvGammaPosteriorSampler(ScaleMixturePriorEstimator):

    def __init__(self):
        super(InvGammaPosteriorSampler, self).__init__()
        self.n_samples = 2

    def get_scale_estimator(self):
        return InvGammaScaleSampler()

    def estimate(self, context, data):
        """
        Generate samples from the posterior of alpha and beta.

        For beta the posterior is a gamma distribution and analytically
        acessible.

        The posterior of alpha can not be expressed analytically and is
        aproximated using adaptive rejection sampling. 
        """
        pdf = GammaPrior()

        ## sufficient statistics
        
        h = harmonic_mean(numpy.clip(data, 1e-308, 1e308))
        g = geometric_mean(numpy.clip(data, 1e-308, 1e308))

        n = len(data) 

        samples = []

        a = numpy.mean(1 / data)
        v = numpy.std(1 / data) ** 2

        beta = a / v
        alpha = beta * a

        for i in range(self.n_samples):

            ## sample alpha with ARS

            logp = ARSPosteriorAlpha(n * (log(beta) - log(g)) - context.hyper_alpha.beta,
                                     context.hyper_alpha.alpha - 1., n)
            ars = csb.statistics.ars.ARS(logp)
            ars.initialize(logp.initial_values()[:2], z0=0.)

            alpha = numpy.abs(ars.sample())

            if alpha is None:
                raise ValueError("Sampling failed")


            ## sample beta from Gamma distribution

            beta = numpy.random.gamma(n * alpha + context.hyper_beta.alpha, \
                                      1 / (n / h + context.hyper_beta.beta))

            samples.append((alpha, beta))

        pdf.alpha, pdf.beta = samples[-1]
        return pdf

        
class InvGammaPosteriorMAP(ScaleMixturePriorEstimator):

    def __init__(self):
        super(InvGammaPosteriorMAP, self).__init__()

    def get_scale_estimator(self):
        return InvGammaScaleMAP()

    def estimate(self, context, data):
        """
        Generate samples from the posterior of alpha and beta.

        For beta the posterior is a gamma distribution and analytically
        acessible.

        The posterior of alpha can not be expressed analytically and is
        aproximated using adaptive rejection sampling. 
        """
        pdf = GammaPrior()

        
        y = log(data).mean() - log((data ** -1).mean())

        alpha = inv_digamma_minus_log(numpy.clip(y,
                                                 - 1e308,
                                                 - 1e-300))
        alpha = abs(alpha)
    
        beta = numpy.clip(alpha / 
                           (data ** (-1)).mean(),
                           1e-100, 1e100)

        pdf.alpha, pdf.beta = alpha, beta
        return pdf



class GammaScaleSampler(AbstractEstimator):
    """
    Sample the scalies given the data
    """

    def estimate(self, context, data):
        pdf = ScaleMixture()
        alpha = context.prior.alpha
        beta = context.prior.beta
        d = context._d

        if len(data.shape) == 1:
            data = data[:, numpy.newaxis]
        
        a = alpha + 0.5 * d * len(data.shape)
        b = beta + 0.5 * data.sum(-1) ** 2

        s = numpy.clip(numpy.random.gamma(a, 1. / b), 1e-20, 1e10)
        pdf.scales = s

        context.prior.estimate(s)
        pdf.prior = context.prior
        
        return pdf


class GammaScaleMAP(AbstractEstimator):
    """
    MAP estimator of the scales
    """

    def estimate(self, context, data):
        pdf = ScaleMixture()
        alpha = context.prior.alpha
        beta = context.prior.beta
        d = context._d

        if len(data.shape) == 1:
            data = data[:, numpy.newaxis]
        
        a = alpha + 0.5 * d * len(data.shape)
        b = beta + 0.5 * data.sum(-1) ** 2

        s = a / b
        log_s = psi(a) - log(b)

        pdf.scales = s
        context.prior.estimate([s, log_s])
        pdf.prior = context.prior

        return pdf

        
    
class InvGammaScaleSampler(AbstractEstimator):
    """
    Sample the scales given the data
    """

    def estimate(self, context, data):
        pdf = ScaleMixture()
        alpha = context.prior.alpha
        beta = context.prior.beta
        d = context._d

        if len(data.shape) == 1:
            data = data[:, numpy.newaxis]

        p = -alpha + 0.5 * d 
        b = 2 * beta 
        a = 1e-5 + data.sum(-1) ** 2

        s = csb.statistics.rand.gen_inv_gaussian(a, b, p)
        s = numpy.clip(s, 1e-300, 1e300)

        pdf.scales = s
        context.prior.estimate(s)
        pdf.prior = context.prior

        return pdf


        
class InvGammaScaleMAP(AbstractEstimator):
    """
    MAP estimator of the scales
    """
    def estimate(self, context, data):
        pdf = ScaleMixture()
        alpha = context.prior.alpha
        beta = context.prior.beta
        d = context._d

        if len(data.shape) == 1:
            data = data[:, numpy.newaxis]

        p = -alpha + 0.5 * d 
        b = 2 * beta 
        a = 1e-5 + data.sum(-1) ** 2 

        s = numpy.clip((numpy.sqrt(b) * kve(p + 1, numpy.sqrt(a * b)))
                       / (numpy.sqrt(a) * kve(p, numpy.sqrt(a * b))),
                       1e-10, 1e10)
        pdf.scales = s
        context.prior.estimate(s)
        pdf.prior = context.prior

        return pdf



class GammaPrior(ScaleMixturePrior, Gamma):
    """
    Gamma prior on mixture weights including a weak gamma prior on its parameter.
    """

    def __init__(self, alpha=1., beta=1., hyper_alpha=(4., 1.),
                 hyper_beta=(2., 1.)):

        super(GammaPrior, self).__init__(alpha, beta)

        self._hyper_alpha = Gamma(*hyper_alpha)
        self._hyper_beta = Gamma(*hyper_beta)
        self.estimator = GammaPosteriorSampler()
   
    @typedproperty(AbstractDensity)
    def hyper_beta():
        pass

    @typedproperty(AbstractDensity)
    def hyper_alpha():
        pass

    def log_prob(self, x):

        a, b = self['alpha'], self['beta']

        l_a = self._hyper_alpha(a)
        l_b = self._hyper_beta(b)

        return super(GammaPrior, self).log_prob(x) + l_a + l_b


        
class InvGammaPrior(ScaleMixturePrior, InverseGamma):
    """
    Inverse Gamma prior on mixture weights including a weak gamma
    prior on its parameter.
    """

    def __init__(self, alpha=1., beta=1., hyper_alpha=(4., 1.),
                 hyper_beta=(2., 1.)):

        super(InvGammaPrior, self).__init__(alpha, beta)

        self._hyper_alpha = Gamma(*hyper_alpha)
        self._hyper_beta = Gamma(*hyper_beta)
        self.estimator = InvGammaPosteriorSampler()

    @typedproperty(AbstractDensity)
    def hyper_beta():
        pass

    @typedproperty(AbstractDensity)
    def hyper_alpha():
        pass

    def log_prob(self, x):

        a, b = self['alpha'], self['beta']

        l_a = self._hyper_alpha(a)
        l_b = self._hyper_beta(b)

        return super(InvGammaPrior, self).log_prob(x) + l_a + l_b




