"""
Mixture models for multi-dimensional data.

Reference: Hirsch M, Habeck M. - Bioinformatics. 2008 Oct 1;24(19):2184-92
"""
import numpy

from abc import ABCMeta, abstractmethod


class GaussianMixture(object):
    """
    Gaussian mixture model for multi-dimensional data.
    """
    _axis = None

    # prior for variance (inverse Gamma distribution)
    ALPHA_SIGMA = 0.0001
    BETA_SIGMA = 0.01
    MIN_SIGMA = 0.0

    use_cache = True

    def __init__(self, X, K, train=True, axis=None):
        """
        @param X: multi dimensional input vector with samples along first axis
        @type X: (M,...) numpy array

        @param K: number of components
        @type K: int

        @param train: train model
        @type train: bool

        @param axis: component axis in C{X}
        @type axis: int
        """
        if self._axis is not None:
            if axis is not None and axis != self._axis:
                raise ValueError('axis is fixed for {0}'.format(type(self).__name__))
            axis = self._axis
        elif axis is None:
            axis = 0
        self._axis = axis

        N = X.shape[axis]
        self._X = X
        self._dimension = numpy.prod(X.shape) / N

        c = numpy.linspace(0, K, N, False).astype(int)
        self._scales = numpy.equal.outer(range(K), c).astype(float)
        self._means = numpy.zeros((K,) + X.shape[1:])
        self.del_cache()

        if train:
            self.em()

    @property
    def K(self):
        """
        Number of components
        @rtype: int
        """
        return len(self.means)

    @property
    def N(self):
        """
        Length of component axis
        @rtype: int
        """
        return self._scales.shape[1]

    @property
    def M(self):
        """
        Number of data points
        @rtype: int
        """
        return len(self._X)

    def del_cache(self):
        """Clear model parameter cache (force recalculation)"""
        self._w = None
        self._sigma = None
        self._delta = None

    @property
    def dimension(self):
        """
        Dimensionality of the mixture domain
        @rtype: int
        """
        return self._dimension

    @property
    def means(self):
        """
        @rtype: (K, ...) numpy array
        """
        return self._means

    @means.setter
    def means(self, means):
        if means.shape != self._means.shape:
            raise ValueError('shape mismatch')
        self._means = means
        self.del_cache()

    @property
    def scales(self):
        """
        @rtype: (K, N) numpy array
        """
        return self._scales

    @scales.setter
    def scales(self, scales):
        if scales.shape != self._scales.shape:
            raise ValueError('shape mismatch')
        self._scales = scales
        self.del_cache()

    @property
    def w(self):
        """
        Component weights
        @rtype: (K,) numpy array
        """
        if not self.use_cache or self._w is None:
            self._w = self.scales.mean(1)
        return self._w

    @property
    def sigma(self):
        """
        Component variations
        @rtype: (K,) numpy array
        """
        if not self.use_cache or self._sigma is None:
            alpha = self.dimension * self.scales.sum(1) + self.ALPHA_SIGMA
            beta = (self.delta * self.scales.T).sum(0) + self.BETA_SIGMA
            self._sigma = numpy.sqrt(beta / alpha).clip(self.MIN_SIGMA)
        return self._sigma

    @property
    def delta(self):
        """
        Squared "distances" between data and components
        @rtype: (N, K) numpy array
        """
        if not self.use_cache or self._delta is None:
            self._delta = numpy.transpose([[d.sum()
                for d in numpy.swapaxes([(self.means[k] - self.datapoint(m, k)) ** 2
                    for m in range(self.M)], 0, self._axis)]
                for k in range(self.K)])
        return self._delta

    @property
    def log_likelihood_reduced(self):
        """
        Log-likelihood of the marginalized model (no auxiliary indicator variables)
        @rtype: float
        """
        from csb.numeric import log, log_sum_exp
        s_sq = (self.sigma ** 2).clip(1e-300, 1e300)
        log_p = log(self.w) - 0.5 * \
                (self.delta / s_sq + self.dimension * log(2 * numpy.pi * s_sq))
        return log_sum_exp(log_p.T).sum()

    @property
    def log_likelihood(self):
        """
        Log-likelihood of the extended model (with indicators)
        @rtype: float
        """
        from csb.numeric import log
        from numpy import pi, sum
        n = self.scales.sum(1)
        N = self.dimension
        Z = self.scales.T
        s_sq = (self.sigma ** 2).clip(1e-300, 1e300)
        return sum(n * log(self.w)) - 0.5 * \
                (sum(Z * self.delta / s_sq) + N * sum(n * log(2 * pi * s_sq)) + sum(log(s_sq)))

    def datapoint(self, m, k):
        """
        Training point number C{m} as if it would belong to component C{k}
        @rtype: numpy array
        """
        return self._X[m]

    def estimate_means(self):
        """
        Update means from current model and samples
        """
        n = self.scales.sum(1)
        self.means = numpy.array([numpy.sum([self.scales[k, m] * self.datapoint(m, k)
            for m in range(self.M)], 0) / n[k]
            for k in range(self.K)])

    def estimate_scales(self, beta=1.0):
        """
        Update scales from current model and samples
        @param beta: inverse temperature
        @type beta: float
        """
        from csb.numeric import log, log_sum_exp, exp
        s_sq = (self.sigma ** 2).clip(1e-300, 1e300)
        Z = (log(self.w) - 0.5 * (self.delta / s_sq + self.dimension * log(s_sq))) * beta
        self.scales = exp(Z.T - log_sum_exp(Z.T))

    def randomize_means(self):
        """
        Pick C{K} samples from C{X} as means
        """
        import random
        self.means = numpy.asarray(random.sample(self._X, self.K))
        self.estimate_scales()

    def randomize_scales(self, ordered=True):
        """
        Random C{scales} initialization
        """
        from numpy.random import random, multinomial
        if ordered:
            K, N = self.scales.shape
            Ks = numpy.arange(K)
            w = random(K) + (5. * K / N) # with pseudocounts
            c = numpy.repeat(Ks, multinomial(N, w / w.sum()))
            self.scales = numpy.equal.outer(Ks, c).astype(float)
        else:
            s = random(self.scales.shape)
            self.scales = s / s.sum(0)
        self.estimate_means()

    def e_step(self, beta=1.0):
        """
        Expectation step for EM
        @param beta: inverse temperature
        @type beta: float
        """
        self.estimate_scales(beta)

    def m_step(self):
        """
        Maximization step for EM
        """
        self.estimate_means()

    def em(self, n_iter=100, eps=1e-30):
        """
        Expectation maximization

        @param n_iter: maximum number of iteration steps
        @type n_iter: int

        @param eps: log-likelihood convergence criterion
        @type eps: float
        """
        LL_prev = -numpy.inf
        for i in range(n_iter):
            self.m_step()
            self.e_step()

            if eps is not None:
                LL = self.log_likelihood
                if abs(LL - LL_prev) < eps:
                    break
                LL_prev = LL

    def anneal(self, betas):
        """
        Deterministic annealing

        @param betas: sequence of inverse temperatures
        @type betas: iterable of floats
        """
        for beta in betas:
            self.m_step()
            self.e_step(beta)

    def increment_K(self, train=True):
        """
        Split component with largest sigma

        @returns: new instance of mixture with incremented C{K}
        @rtype: L{GaussianMixture} subclass
        """
        i = self.sigma.argmax()

        # duplicate column
        Z = numpy.vstack([self.scales, self.scales[i]])

        # mask disjoint equal sized parts
        mask = Z[i].cumsum() / Z[i].sum() > 0.5
        Z[i, mask] *= 0.0
        Z[-1, ~mask] *= 0.0

        new = type(self)(self._X, self.K + 1, False, self._axis)
        new.scales = Z
        new.m_step()
        if train:
            new.em()

        return new

    @classmethod
    def series(cls, X, start=1, stop=9):
        """
        Iterator with mixture instances for C{K in range(start, stop)}

        @type X: (M,...) numpy array
        @type start: int
        @type stop: int
        @rtype: generator
        """
        mixture = cls(X, start)
        yield mixture

        for K in range(start + 1, stop):                        #@UnusedVariable
            mixture = mixture.increment_K()
            yield mixture

    @classmethod
    def new(cls, X, K=0):
        """
        Factory method with optional C{K}. If C{K=0}, guess best C{K} according
        to L{BIC<GaussianMixture.BIC>}.

        @param X: multi dimensional input vector with samples along first axis
        @type X: (M,...) numpy array

        @return: Mixture instance
        @rtype: L{GaussianMixture} subclass
        """
        if K > 0:
            return cls(X, K)

        mixture_it = cls.series(X)
        mixture = next(mixture_it)

        # increase K as long as next candidate looks better
        for candidate in mixture_it:
            if candidate.BIC >= mixture.BIC:
                break
            mixture = candidate

        return mixture

    @property
    def BIC(self):
        """
        Bayesian information criterion, calculated as
        BIC = M * ln(sigma_e^2) + K * ln(M)

        @rtype: float
        """
        from numpy import log

        n = self.M
        k = self.K
        error_variance = sum(self.sigma ** 2 * self.w)

        return n * log(error_variance) + k * log(n)

    @property
    def membership(self):
        """
        Membership array
        @rtype: (N,) numpy array
        """
        return self.scales.argmax(0)

    def overlap(self, other):
        """
        Similarity of two mixtures measured in membership overlap

        @param other: Mixture or membership array
        @type other: L{GaussianMixture} or sequence

        @return: segmentation overlap
        @rtype: float in interval [0.0, 1.0]
        """
        if isinstance(other, GaussianMixture):
            other_w = other.membership
            K = min(self.K, other.K)
        elif isinstance(other, (list, tuple, numpy.ndarray)):
            other_w = other
            K = min(self.K, len(set(other)))
        else:
            raise TypeError('other')

        self_w = self.membership
        if len(self_w) != len(other_w):
            raise ValueError('self.N != other.N')

        # position numbers might be permutated, so count equal pairs
        ww = tuple(zip(self_w, other_w))
        same = sum(sorted(ww.count(i) for i in set(ww))[-K:])

        return float(same) / len(ww)

class AbstractStructureMixture(GaussianMixture):
    """
    Abstract mixture model for protein structure ensembles.
    """
    __metaclass__ = ABCMeta

    def __init__(self, X, K, *args, **kwargs):
        if len(X.shape) != 3 or X.shape[-1] != 3:
            raise ValueError('X must be array of shape (M,N,3)')

        self._R = numpy.zeros((len(X), K, 3, 3))
        self._t = numpy.zeros((len(X), K, 3))

        super(AbstractStructureMixture, self).__init__(X, K, *args, **kwargs)

    @property
    def R(self):
        """
        Rotation matrices
        @rtype: (M,K,3,3) numpy array
        """
        return self._R

    @property
    def t(self):
        """
        Translation vectors
        @rtype: (M,K,3) numpy array
        """
        return self._t

    def datapoint(self, m, k):
        return numpy.dot(self._X[m] - self._t[m, k], self._R[m, k])

    def m_step(self):
        self.estimate_means()
        self.estimate_T()

    @abstractmethod
    def estimate_T(self):
        """
        Estimate superpositions
        """
        raise NotImplementedError

class SegmentMixture(AbstractStructureMixture):
    """
    Gaussian mixture model for protein structure ensembles using a set of segments

    If C{X} is the coordinate array of a protein structure ensemble which
    can be decomposed into 2 rigid segments, the segmentation will be found by:

    >>> mixture = SegmentMixture(X, 2)

    The segment membership of each atom is given by:

    >>> mixture.membership
    array([0, 0, 0, ..., 1, 1, 1])
    """
    _axis = 1

    def estimate_T(self):
        from csb.bio.utils import wfit
        for m in range(self.M):
            for k in range(self.K):
                self._R[m, k], self._t[m, k] = wfit(self._X[m], self.means[k], self.scales[k])

    def estimate_means(self):
        # superpositions are weighted, so do unweighted mean here
        self.means = numpy.mean([[self.datapoint(m, k)
            for m in range(self.M)]
            for k in range(self.K)], 1)

class ConformerMixture(AbstractStructureMixture):
    """
    Gaussian mixture model for protein structure ensembles using a set of conformers

    If C{mixture} is a trained model, the ensemble coordinate array of
    structures from C{X} which belong to conformation C{k} is given by:

    >>> indices = numpy.where(mixture.membership == k)[0]
    >>> conformer = [mixture.datapoint(m, k) for m in indices]
    """
    _axis = 0

    def estimate_T(self):
        from csb.bio.utils import fit
        for m in range(self.M):
            for k in range(self.K):
                self._R[m, k], self._t[m, k] = fit(self._X[m], self.means[k])

# vi:expandtab:smarttab:sw=4
