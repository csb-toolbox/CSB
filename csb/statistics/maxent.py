"""
A Maximum-Entropy model for backbone torsion angles.
Reference: Rowicka and Otwinowski 2004
"""

import numpy

from csb.statistics.pdf import AbstractDensity


class MaxentModel(AbstractDensity):
    """
    Fourier expansion of a biangular log-probability density
    """
    def __init__(self, n, beta=1.):
        """

        @param n: order of the fourier expansion
        @type n: int

        @param beta: inverse temperature
        @type beta: float
        """
        super(MaxentModel, self).__init__()
        
        self._n = int(n)

        self._cc = numpy.zeros((self._n, self._n))
        self._ss = numpy.zeros((self._n, self._n))
        self._cs = numpy.zeros((self._n, self._n))
        self._sc = numpy.zeros((self._n, self._n))

        self._beta = float(beta)

    @property
    def beta(self):
        """
        Inverse temperature

        @rtype: float
        """
        return self._beta

    @property
    def n(self):
        """
        Order of the fourier expansion

        @rtype: int
        """
        return self._n

    def load_old(self, aa, f_name):
        """
        Load set of expansion coefficients from isd.

        @param aa: Amino acid type
        @param f_name: File containing ramachandran definition
        """

        import os
        params, _energies = eval(open(os.path.expanduser(f_name)).read())
        params = params[self._n - 1]

        for k, l, x, f, g in params[aa]:

            if f == 'cos' and g == 'cos':
                self._cc[k, l] = -x
            elif f == 'cos' and g == 'sin':
                self._cs[k, l] = -x
            elif f == 'sin' and g == 'cos':
                self._sc[k, l] = -x
            elif f == 'sin' and g == 'sin':
                self._ss[k, l] = -x

    def load(self, aa, f_name):
        """
        Load set of expansion coefficients from isd+.

        @param aa: Amino acid type
        @param f_name: File containing ramachandran definition
        """
        import os
        from numpy import reshape, array
        from csb.io import load

        f_name = os.path.expanduser(f_name)
        params, _energies = load(f_name)
        params = params[self._n]

        a, b, c, d = params[aa]
        a, b, c, d = reshape(array(a), (self._n, self._n)).astype('d'), \
                  reshape(array(b), (self._n, self._n)).astype('d'), \
                  reshape(array(c), (self._n, self._n)).astype('d'), \
                  reshape(array(d), (self._n, self._n)).astype('d')
        # Not a typo, I accidently swichted cos*sin and sin*cos
        self._cc, self._cs, self._sc, self._ss = -a, -c, -b, -d

    def _periodicities(self):
        return numpy.arange(self._n)

    def log_prob(self, x, y):
        """
        Return the energy at positions (x,y).

        @param x: x-coordinates for evaluation
        @type x: array-like

        @param y: y-coordinates for evaluation
        @type y: array-like
        """
        return -self.energy(x, y)
    
    def set(self, coef):
        """
        Set the fourier expansion coefficients and calculations the 
        new partation function.

        @param coef: expansion coefficents
        @type coef: array like, with shape (4,n,n)
        """
        self._cc[:, :], self._ss[:, :], self._cs[:, :], self._sc[:, :] = \
                    numpy.reshape(coef, (4, self._n, self._n))
        self.normalize()

    def get(self):
        """
        Return current expansion coefficients.
        """
        return numpy.array([self._cc, self._ss, self._cs, self._sc])

    def energy(self, x, y=None):
        """
        Return the energy at positions (x,y). 

        @param x: x-coordinates for evaluation
        @type x: array-like

        @param y: y-coordinates for evaluation
        @type y: array-like
        """
        from numpy import sin, cos, dot, multiply

        k = self._periodicities()
        cx, sx = cos(multiply.outer(k, x)), sin(multiply.outer(k, x))
        if y is not None:
            cy, sy = cos(multiply.outer(k, y)), sin(multiply.outer(k, y))
        else:
            cy, sy = cx, sx

        return dot(dot(cx.T, self._cc), cy) + \
               dot(dot(cx.T, self._cs), sy) + \
               dot(dot(sx.T, self._sc), cy) + \
               dot(dot(sx.T, self._ss), sy)

    def sample_weights(self):
        """
        Create a random set of expansion coefficients.
        """
        from numpy import add
        from numpy.random import standard_normal

        k = self._periodicities()
        k = add.outer(k ** 2, k ** 2)
        self.set([standard_normal(k.shape) for i in range(4)])   
        self.normalize(True)

    def prob(self, x, y):
        """
        Return the probability of the configurations x cross y.
        """
        from csb.math import exp
        return exp(-self.beta * self(x, y))

    def z(self):
        """
        Calculate the partion function .
        """
        from scipy.integrate import dblquad
        from numpy import pi

        return dblquad(self.prob, 0., 2 * pi, lambda x: 0., lambda x: 2 * pi)

    def log_z(self, n=500, integration='simpson'):
        """
        Calculate the log partion function.
        """
        from numpy import pi, linspace, max
        from csb.math import log, exp

        if integration == 'simpson':
            from csb.numeric import simpson_2d
            x = linspace(0., 2 * pi, 2 * n + 1)
            dx = x[1] - x[0]

            f = -self.beta * self.energy(x)
            f_max = max(f)
            f -= f_max

            I = simpson_2d(exp(f))
            return log(I) + f_max + 2 * log(dx)

        elif integration == 'trapezoidal':

            from csb.numeric import trapezoidal_2d
            x = linspace(0., 2 * pi, n)
            dx = x[1] - x[0]

            f = -self.beta * self.energy(x)
            f_max = max(f)
            f -= f_max
            I = trapezoidal_2d(exp(f))
            return log(I) + f_max + 2 * log(dx)
        else:
            raise NotImplementedError(
                'Choose from trapezoidal and simpson-rule Integration')
    
    def entropy(self, n=500):
        """
        Calculate the entropy of the model.

        @param n: number of integration points for numerical integration
        @type n: integer
        """
        from csb.numeric import trapezoidal_2d
        from numpy import pi, linspace, max
        from csb.math import log, exp

        x = linspace(0., 2 * pi, n)
        dx = x[1] - x[0]

        f = -self.beta * self.energy(x)
        f_max = max(f)

        log_z = log(trapezoidal_2d(exp(f - f_max))) + f_max + 2 * log(dx)
        average_energy = trapezoidal_2d(f * exp(f - f_max))\
                         * exp(f_max + 2 * log(dx) - log_z)

        return -average_energy + log_z

    def calculate_statistics(self, data):
        """
        Calculate the sufficient statistics for the data.
        """
        from numpy import cos, sin, dot, multiply

        k = self._periodicities()
        cx = cos(multiply.outer(k, data[:, 0]))
        sx = sin(multiply.outer(k, data[:, 0]))
        cy = cos(multiply.outer(k, data[:, 1]))
        sy = sin(multiply.outer(k, data[:, 1]))

        return dot(cx, cy.T), dot(sx, sy.T), dot(cx, sy.T), dot(sx, cy.T)

    def normalize(self, normalize_full=True):
        """
        Remove parameter, which do not have any influence on the model
        and compute the partition function.

        @param normalize_full: compute partition function
        @type normalize_full: boolean
        """
        self._cc[0, 0] = 0.
        self._ss[:, 0] = 0.
        self._ss[0, :] = 0.
        self._cs[:, 0] = 0.
        self._sc[0, :] = 0.

        if normalize_full:
            self._cc[0, 0] = self.log_z()


class MaxentPosterior(object):
    """
    Object to hold and calculate the posterior (log)probability
    given an exponential family model and corresponding data.
    """

    def __init__(self, model, data):
        """
        @param model: MaxentModel
        @param data: two dimensonal data
        """
        
        self._model = model
        self._data = numpy.array(data)
        self._stats = self.model.calculate_statistics(self._data)

        self._log_likelihoods = []

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, value):
        self._model = value
        self._stats = self.model.calculate_statistics(self._data)
        
    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, value):
        self._data = numpy.array(value)
        self._stats = self.model.calculate_statistics(value)

    @property
    def stats(self):
        return self._stats

    def __call__(self, weights=None, n=100):
        """
        Returns the log posterior likelihood

        @param weights: optional expansion coefficients of the model,
                         if none are specified those of the model are used
        @param n: number of integration point for calculating the partition function
        """
        from numpy import sum

        if weights is not None:
            self.model.set(weights)

        a = sum(self._stats[0] * self.model._cc)
        b = sum(self._stats[1] * self.model._ss)
        c = sum(self._stats[2] * self.model._cs)
        d = sum(self._stats[3] * self.model._sc)

        log_z = self.data.shape[0] * self.model.log_z(n=n)

        log_likelihood = -self.model.beta * (a + b + c + d) - log_z

        self._log_likelihoods.append(log_likelihood)

        return log_likelihood
