"""
Adaptive Rejection Sampling (ARS)

The ARS class generates a single random sample from a
univariate distribution specified by an instance of the
LogProb class, implemented by the user. An instance of
LogProb returns the log of the probability density and
its derivative. The log probability function passed must
be concave.

The user must also supply initial guesses.  It is not
essential that these values be very accurate, but performance
will generally depend on their accuracy.
"""

from numpy import exp, log

class Envelope(object):
    """
    Envelope function for adaptive rejection sampling.
    
    The envelope defines a piecewise linear upper and lower
    bounding function of the concave log-probability.
    """
    def __init__(self, x, h, dh):

        from numpy import array, inf

        self.x = array(x)
        self.h = array(h)
        self.dh = array(dh)
        self.z0 = -inf
        self.zk = inf
        
    def z(self):
        """
        Support intervals for upper bounding function.
        """
        from numpy import concatenate

        h = self.h
        dh = self.dh
        x = self.x

        z = (h[1:] - h[:-1] + x[:-1] * dh[:-1] - x[1:] * dh[1:]) / \
            (dh[:-1] - dh[1:])

        return concatenate(([self.z0], z, [self.zk]))

    def u(self, x):
        """
        Piecewise linear upper bounding function.
        """
        z = self.z()[1:-1]
        j = (x > z).sum()

        return self.h[j] + self.dh[j] * (x - self.x[j])

    def u_max(self):

        z = self.z()[1:-1]

        return (self.h + self.dh * (z - self.x)).max()

    def l(self, x):
        """
        Piecewise linear lower bounding function.
        """
        from numpy import inf

        j = (x > self.x).sum()

        if j == 0 or j == len(self.x):
            return -inf
        else:
            j -= 1
            return ((self.x[j + 1] - x) * self.h[j] + (x - self.x[j]) * self.h[j + 1]) / \
                   (self.x[j + 1] - self.x[j])

    def insert(self, x, h, dh):
        """
        Insert new support point for lower bounding function
        (and indirectly for upper bounding function).
        """
        from numpy import concatenate

        j = (x > self.x).sum()

        self.x = concatenate((self.x[:j], [x], self.x[j:]))
        self.h = concatenate((self.h[:j], [h], self.h[j:]))
        self.dh = concatenate((self.dh[:j], [dh], self.dh[j:]))

    def log_masses(self):
        
        from numpy import  abs, putmask

        z = self.z()
        b = self.h - self.x * self.dh
        a = abs(self.dh)
        m = (self.dh > 0)
        q = self.x * 0.        
        putmask(q, m, z[1:])
        putmask(q, 1 - m, z[:-1])
        
        log_M = b - log(a) + log(1 - exp(-a * (z[1:] - z[:-1]))) + \
                self.dh * q

        return log_M

    def masses(self):

        z = self.z()
        b = self.h - self.x * self.dh
        a = self.dh
        
        return exp(b) * (exp(a * z[1:]) - exp(a * z[:-1])) / a

    def sample(self):

        from numpy.random import random
        from numpy import add
        from csb.numeric import log_sum_exp
        
        log_m = self.log_masses()
        log_M = log_sum_exp(log_m)
        c = add.accumulate(exp(log_m - log_M))
        u = random()
        j = (u > c).sum()

        a = self.dh[j]
        z = self.z()
        
        xmin, xmax = z[j], z[j + 1]

        u = random()

        if a > 0:
            return xmax + log(u + (1 - u) * exp(-a * (xmax - xmin))) / a
        else:
            return xmin + log(u + (1 - u) * exp(a * (xmax - xmin))) / a


class LogProb(object):

    def __call__(self, x):
        raise NotImplementedError()

class Gauss(LogProb):

    def __init__(self, mu, sigma=1.):

        self.mu = float(mu)
        self.sigma = float(sigma)

    def __call__(self, x):

        return -0.5 * (x - self.mu) ** 2 / self.sigma ** 2, \
               - (x - self.mu) / self.sigma ** 2


class ARS(object):

    from numpy import inf

    def __init__(self, logp):

        self.logp = logp

    def initialize(self, x, z0=-inf, zmax=inf):

        from numpy import array

        self.hull = Envelope(array(x), *self.logp(array(x)))
        self.hull.z0 = z0
        self.hull.zk = zmax

    def sample(self, maxiter=100):

        from numpy.random import random

        for i in range(maxiter):

            x = self.hull.sample()
            l = self.hull.l(x)
            u = self.hull.u(x)
            w = random()

            if w <= exp(l - u): return x

            h, dh = self.logp(x)

            if w <= exp(h - u): return x

            self.hull.insert(x, h, dh)
