"""
Mixture models for protein structure ensembles

Reference: Hirsch M, Habeck M. - Bioinformatics. 2008 Oct 1;24(19):2184-92
"""

from abc import ABCMeta, abstractmethod

class GaussianMixture(object):
    """
    Abstract mixture model for protein structure ensembles
    """

    __metaclass__ = ABCMeta

    n_inner = 1           ## inner loop in M-step

    alpha_sigma = 0.0001  ## prior for variance (inverse Gamma
    beta_sigma  = 0.01    ## distribution)
    min_sigma   = 0.0     ##

    use_cache = True      ## cache for distance matrix

    def __init__(self, N, M, K=1, beta=1.):
        """
        @param N: number of atoms
        @type N: int

        @param M: number of structures
        @type M: int

        @param K: number of components
        @type K: int

        @param beta: reciprocal temperature (for annealed EM)
        @type beta: float
        """
        from numpy import zeros, ones, identity, array

        self.K = int(K)
        self.N = int(N)
        self.M = int(M)
        self.beta = float(beta)

        self.Y = zeros((self.K, self.N, 3))
        self.R = array([[identity(3) for k in range(self.K)]
                        for m in range(self.M)])
        self.t = zeros((self.M, self.K, 3))
        self.sigma = ones(self.K)
        self.w = 1./self.K + zeros(self.K)

        self.Z = None

        self.del_cache()

    def del_cache(self):
        self.cache = None

    def get_dimension(self, Z):
        """
        dimensionality of the mixture domain
        """
        return 3 * (int(len(Z) == self.M) * self.N + \
                    int(len(Z) == self.N) * self.M)

    def log_likelihood_reduced(self, X):
        """
        log-likelihood of the marginalized model (no auxiliary indicator
        variables)
        """
        from csb.math import log, log_sum_exp
        from numpy import pi, sum, transpose, clip

        D = self.get_delta(X)
        N = self.get_dimension(D)

        log_p = - 0.5 * D / clip(self.sigma**2, 1e-300, 1e300) \
                - 0.5 * N * log(2 * pi * self.sigma**2) \
                + log(self.w)

        return sum(log_sum_exp(transpose(log_p)))

    def log_likelihood(self, X, Z):
        """
        log-likelihood of the extended model (with indicators)
        """
        from csb.math import log
        from numpy import pi, sum, clip

        D = self.get_delta(X)
        n = sum(Z, 0)
        N = self.get_dimension(Z)

        L = - 0.5 * sum(Z * D / clip(self.sigma**2, 1e-300, 1e300)) \
            - 0.5 * N * sum(n * log(2 * pi * self.sigma**2)) \
            + sum(n * log(self.w)) \
            - 0.5 * sum(log(self.sigma**2)) ## Jeffreys' prior

        return L

    def get_delta(self, X, update_cache=False):
        """
        calculate squared 'distance matrix' between data and components
        """
        if update_cache and self.use_cache:
            D = self.calculate_delta(X)
            self.cache = D

        if self.cache is not None:
            return self.cache
        return self.calculate_delta(X)

    @abstractmethod
    def calculate_delta(self, X):
        pass

    def estimate_w(self, Z):
        """
        estimate weights of components from indicators
        """
        from numpy import sum

        n = sum(Z, 0)

        self.w[:] = n / sum(n)

    def estimate_sigma(self, X, Z):
        """
        estimate variances of the components
        """
        from numpy import sum, sqrt

        D = self.get_delta(X, update_cache=True)
        N = self.get_dimension(Z)

        alpha = N * sum(Z, 0) + self.alpha_sigma
        beta  = sum(Z*D, 0) + self.beta_sigma

        self.sigma = sqrt(beta / alpha).clip(self.min_sigma)

    @abstractmethod
    def estimate_Y(self, X, Z):
        """
        estimate conformers/segments
        """
        pass

    @abstractmethod
    def estimate_T(self, X, Z):
        """
        estimate superpositions
        """
        pass

    def e_step(self, X):
        """
        expectation step: calculation of indicators
        """
        from numpy import clip
        from csb.math import log, log_sum_exp, exp

        Z  = -0.5 * self.get_delta(X) / clip(self.sigma**2, 1e-300, 1e300)
        Z -=  0.5 * self.get_dimension(Z) * log(self.sigma**2)
        Z += log(self.w)
        Z  = self.beta * Z.T
        Z -= log_sum_exp(Z)

        self.del_cache()

        return exp(Z).T

    def m_step(self, X, Z):
        """
        maximization step: maximize posterior given indicators
        """
        for i in range(self.n_inner):
            self.estimate_Y(X, Z)
            self.estimate_T(X, Z)

        self.estimate_w(Z)
        self.estimate_sigma(X, Z)

    @abstractmethod
    def initialize(self, X, randomize=True):
        """
        @param X: MxNx3 input vector (M: #structures, N: #atoms)
        @type X: numpy array

        @param randomize: Randomize segment sizes. If False, all components will
        be equal sized and continuous.
        @type randomize: bool
        """
        pass

    def em(self, X, n_iter=100, verbose=0, initialize=True,
            store_loglikelihood=False, eps=1e-30):
        """
        Expectation maximization

        @param X: MxNx3 input vector (M: #structures, N: #atoms)
        @type X: numpy array

        @param n_iter: maximum number of iteration steps
        @type n_iter: int

        @param verbose: if non-zero, print out log likelihood every C{verbose} step.
        @type verbose: int

        @param initialize: if True, call L{GaussianMixture.initialize}
        @type initialize: bool
        """
        if initialize:
            self.initialize(X)

        if store_loglikelihood:
            self.LL = []

        LL_prev = 1e10

        for i in range(n_iter):

            self.Z = self.e_step(X)
            self.m_step(X, self.Z)

            output = verbose and not (i % verbose)

            if not (store_loglikelihood or output) and eps is None:
                continue

            LL = self.log_likelihood(X, self.Z)

            if store_loglikelihood:
                self.LL.append(LL)

            if output:
                print i, LL

            # convergence criteria
            if abs(LL - LL_prev) < eps:
                break

            LL_prev = LL

    def anneal(self, X, betas, verbose=0):
        """
        deterministic annealing
        """
        self.beta = betas[0]
        self.initialize(X)

        if verbose:
            self.LL = []

        i = 0
        for beta in betas:

            self.beta = beta
            self.Z = Z = self.e_step(X)

            self.m_step(X, Z)

            if verbose and not (i % verbose):
                self.LL.append(self.log_likelihood(X, Z))
                print beta, self.LL[-1]

            i += 1

    def increment_K(self, X):
        """
        Split component with largest sigma

        @returns: new instance of mixture with incremented K
        @rtype: L{GaussianMixture} subclass
        """
        from numpy import argmax, hstack

        i = argmax(self.sigma)

        # duplicate column
        Z = hstack([self.Z, self.Z[:,i].reshape((-1, 1))])

        # mask disjoint equal sized parts
        mask = Z[:,i].cumsum() / Z[:,i].sum() > 0.5
        Z[mask,i] *= 0.0
        Z[mask == False,-1] *= 0.0

        new = type(self)(self.N, self.M, self.K + 1, self.beta)
        new.Z = Z
        new.m_step(X, new.Z)

        return new

    @classmethod
    def series(cls, X, n_iter=100, verbose=0, start=1, stop=9, eps=1e-30):
        """
        For param description, see L{GaussianMixture.em}.

        @return: Iterator with mixture instances for C{K in range(start, stop)}
        @rtype: generator
        """
        m = cls(X.shape[1], X.shape[0], start)
        m.initialize(X, False)

        for K in range(start+1, stop):
            m.em(X, n_iter, verbose, False, False, eps)
            yield m

            m = m.increment_K(X)

        yield m

    @classmethod
    def from_coords(cls, X, K=0, n_iter=100, verbose=0, randomize=False, eps=1e-30):
        """
        Factory method with coodinate array.

        @param X: MxNx3 input vector (M: #structures, N: #atoms)
        @type X: numpy array

        @param K: Number of components. If zero, then find best C{K} with
        non-random initialization and sequential splitting of component with
        largest sigma. Best C{K} is judged heuristically as mixture with
        smallest L{BIC<GaussianMixture.BIC>}.
        @type K: int 

        @param n_iter: maximum number of iteration steps
        @type n_iter: int

        @param verbose: if non-zero, print out log likelihood every C{verbose} step.
        @type verbose: int

        @param randomize: if non-zero, do random initialization of components
        C{randomize} times and take mixture with maximum log(likelihood).
        @type randomize: int

        @return: Mixture instance
        @rtype: L{GaussianMixture} subclass
        """
        if isinstance(K, int):
            if K > 0:
                if randomize > 1:
                    return max([cls.from_coords(X, K, n_iter, verbose, 1, eps)
                        for _ in range(randomize)],
                        key = lambda mixture: mixture.log_likelihood_reduced(X))

                mixture = cls(X.shape[1], X.shape[0], K)
                mixture.initialize(X, randomize)
                mixture.em(X, n_iter, verbose, False, False, eps)

                return mixture

            K = (1, 9)

        if randomize:
            raise ValueError("randomize=True and K=0 are exclusive")

        try:
            start, stop = K
        except:
            raise TypeError("K must be integer or sequence of 2 integers")

        mixture_it = cls.series(X, n_iter, verbose, start, stop, eps)
        mixture = mixture_it.next()

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
        """
        from numpy import log

        n = self.M
        k = self.K
        error_variance = sum(self.sigma**2 * self.w)

        return n * log(error_variance) + k * log(n)

    def overlap(self, other):
        """
        Similarity of two mixtures measured in membership overlap

        @param other: Mixture with C{len(self.Z) == len(other.Z)} or a
        segmentation array with C{len(self.Z) == len(other)}.
        @type other: L{GaussianMixture} or sequence

        @return: segmentation overlap
        @rtype: float in interval [0.0, 1.0]
        """
        from numpy import argmax, ndarray

        if isinstance(other, GaussianMixture):
            other_w = argmax(other.Z, 1)
            K = min(self.K, other.K)
        elif isinstance(other, (list, tuple, ndarray)):
            other_w = other
            K = min(self.K, len(set(other)))
        else:
            raise TypeError

        self_w = argmax(self.Z, 1)
        assert len(self_w) == len(other_w)

        # position numbers might be permutated, so count equal pairs
        ww = zip(self_w, other_w)
        same = sum(sorted(ww.count(i) for i in set(ww))[-K:])

        return float(same) / len(ww)


class SegmentMixture(GaussianMixture):
    """
    Gaussian mixture model for protein structure ensembles using a set of
    segments
    
    If C{X} is the coordinate array of a protein structure ensemble which
    can be decomposed into 2 rigid segments, the segmentation will be found by:

    >>> mixture = SegmentMixture.from_coords(X, 2)
    >>> membership = numpy.argmax(mixture.Z, 1)

    Superposition on first segment:

    >>> X_superposed = mixture.clusters(X)[0]

    """

    def calculate_delta(self, X):

        from numpy import dot, array, sum, transpose

        D = array([[sum((self.Y[k] - dot(X[m] - self.t[m,k], self.R[m,k]))**2, 1)
                    for k in range(self.K)]
                   for m in range(self.M)])

        return transpose(sum(D, 0))

    def estimate_Y(self, X, Z):
        """
        (unweighted) average of backtransformed structures
        """
        from numpy import sum, dot

        for k in range(self.K):
            self.Y[k,:,:] = sum([dot(X[m] - self.t[m,k], self.R[m,k])
                                 for m in range(self.M)], 0) / self.M

    def estimate_T(self, X, Z):
        """
        weighted least-squares fit
        """
        from csb.bio.utils import wfit

        for m in range(self.M):
            for k in range(self.K):
                self.R[m,k,:,:], self.t[m,k,:] = wfit(X[m], self.Y[k], Z[:,k])

    def initialize(self, X, randomize=True):

        from numpy.random import random, multinomial
        from numpy import repeat, equal, arange, linspace

        if randomize:
            w = random(self.K) + (5. * self.K / self.N)
            w/= sum(w)

            c = repeat(arange(self.K), multinomial(self.N, w))
        else:
            c = linspace(0, self.K, self.N, False).astype(int)

        self.Z = 1. * equal.outer(c, arange(self.K))

        self.m_step(X, self.Z)

    def clusters(self, X):
        """
        Return C{K} copies of the ensemble coordinate array, each superposed
        on a different segment.
        """

        from numpy import argmax, array, dot

        C = [array([dot(X[m] - self.t[m,k], self.R[m,k])
                    for m in range(self.M)])
             for k in range(self.K)]

        return C

class SegmentMixture2(SegmentMixture):
    """
    Using a single coordinate array
    """

    def __init__(self, *args, **kw):

        super(SegmentMixture2, self).__init__(*args, **kw)
        self.Y = self.Y[0]

    def calculate_delta(self, X):

        from numpy import dot, array, sum, transpose

        D = array([[sum((self.Y - dot(X[m] - self.t[m,k], self.R[m,k]))**2, 1)
                    for k in range(self.K)]
                   for m in range(self.M)])

        return transpose(sum(D, 0))

    def estimate_Y(self, X, Z):
        """
        (unweighted) average of backtransformed structures
        """
        from numpy import sum, dot, clip, transpose

        s = 1. / clip(self.sigma**2, 1e-308, 1e308)
        n = self.M * dot(Z, s)

        self.Y[:,:] = 0.

        for k in range(self.K):

            Y = sum([dot(X[m] - self.t[m,k], self.R[m,k])
                     for m in range(self.M)], 0)
            self.Y += transpose(transpose(Y) * Z[:,k] * s[k])

        self.Y = transpose(transpose(self.Y) / n)

    def estimate_T(self, X, Z):
        """
        weighted least-squares fit
        """
        from csb.bio.utils import wfit

        for m in range(self.M):
            for k in range(self.K):
                self.R[m,k,:,:], self.t[m,k,:] = wfit(X[m], self.Y, Z[:,k])

class ConformerMixture(GaussianMixture):
    """
    Gaussian mixture model for protein structure ensembles using a set of
    conformers
    """

    def calculate_delta(self, X):

        from numpy import dot, array, sum

        return array([[sum((self.Y[k] - dot(X[m] - self.t[m,k], self.R[m,k]))**2)
                       for k in range(self.K)]
                      for m in range(self.M)])

    def estimate_Y(self, X, Z):
        """
        weighted average of backtransformed structures
        """
        from numpy import sum, dot, clip

        n = clip(sum(Z, 0), 1e-300, 1e300)

        for k in range(self.K):
            self.Y[k,:,:] = sum([Z[m,k] * dot(X[m] - self.t[m,k], self.R[m,k])
                                 for m in range(self.M)], 0)[:,:] / n[k]

    def estimate_T(self, X, Z):
        """
        (unweighted) least-squares fit
        """
        from csb.bio.utils import fit

        for m in range(self.M):
            for k in range(self.K):
                self.R[m,k,:,:], self.t[m,k,:] = fit(X[m], self.Y[k])

    def initialize(self, X, randomize=True):

        from numpy.random import permutation
        from numpy import linspace

        if randomize:
            self.Y = X[permutation(self.M)[:self.K]]
        else:
            self.Y = X[linspace(0, self.M-1, self.K).astype(int)]

        self.estimate_T(X, None)
        Z = self.e_step(X)
        self.estimate_w(Z)
        self.estimate_sigma(X, Z)

    def clusters(self, X):
        """
        Return C{K} ensemble coordinate arrays, representing C{K} conformations.
        """

        from numpy import argmax, array, dot, where

        Z = self.e_step(X)
        c = argmax(Z, 1)

        C = [array([dot(X[m] - self.t[m,k], self.R[m,k])
                    for m in where(c == k)[0]])
             for k in range(self.K)]

        return C

# vi:expandtab:smarttab:sw=4
