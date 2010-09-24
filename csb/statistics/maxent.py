"""
A Maximum-Entropy model for backbone torsion angles.
Reference:
Rowicka and  Otwinowski 2004
"""

def load_data(aa, ss='all', path='~/data/dunbrack', isd=True):
    """
    Load raw rotamer data

    @param aa: Three-letter aminoacid code
    @param ss: Secondary structure element (H,E,all)
    @param path: search path for data files
    @paran isd map
    """
    from numpy import array, compress
    from csb.math import degree2radian
    from csb.io import load

    data = array(load(path + '/%s_%s' % (aa, ss)).tolist())
    data = compress(data[:, 0] != 360., data, 0)
    data = compress(data[:, 1] != 360., data, 0)
    if isd:
        data = map_to_isd(data)
    else:
        data = degree2radian(data)

    return data


def map_to_isd(angles):
    """
    map dihedrals angles in torsion degrees to the ISD definition
    """
    from csb.math import degree2radian
    from numpy import array, pi, fmod

    angles = fmod(array(angles), 180.)

    return 2 * pi - array(degree2radian(angles))



def torsion_to_radian(phi, psi, radian=0):
    """
    Comverts a torsian angle to radian
    """
    from csb.math import degree2radian

    if radian == 1:
        phi = degree2radian(phi)
        psi = degree2radian(psi)

    return phi, psi



class MaxentModel:
    """
    Fourier expansion of a biangular log-probability density
    """
    def __init__(self, n, beta=1.):
        from numpy import zeros

        self.n = int(n)

        self.cc = zeros((self.n, self.n))
        self.ss = zeros((self.n, self.n))
        self.cs = zeros((self.n, self.n))
        self.sc = zeros((self.n, self.n))
        self.beta = float(beta)


    def load_old(self, aa, f_name):
        """
        Loads set of expansion coefficents from isd

        param aa: Amino acid type
        param f_name: File containing ramachandran definition
        """

        import os
        params, _energies = eval(open(os.path.expanduser(f_name)).read())
        params = params[self.n - 1]

        for k, l, x, f, g in params[aa]:

            if f == 'cos' and g == 'cos':
                self.cc[k, l] = -x
            elif f == 'cos' and g == 'sin':
                self.cs[k, l] = -x
            elif f == 'sin' and g == 'cos':
                self.sc[k, l] = -x
            elif f == 'sin' and g == 'sin':
                self.ss[k, l] = -x

    def load(self, aa, f_name):
        """
        Loads set of expansion coefficents from isd+

        @param aa: Amino acid type
        @param f_name: File containing ramachandran definition
        """
        import os
        from numpy import reshape, array
        from csb.io import load

        f_name = os.path.expanduser(f_name)
        params, _energies = load(f_name)
        params = params[self.n]

        a, b, c, d = params[aa]
        a, b, c, d = reshape(array(a), (self.n, self.n)).astype('d'), \
                  reshape(array(b), (self.n, self.n)).astype('d'), \
                  reshape(array(c), (self.n, self.n)).astype('d'), \
                  reshape(array(d), (self.n, self.n)).astype('d')
        # Not a typo, I accidently swichted cos*sin and sin*cos
        self.cc, self.cs, self.sc, self.ss = -a, -c, -b, -d


    def periodicities(self):
        from numpy import arange

        return arange(self.n)

    def __call__(self, x, y):
        """
        returns the energy at positions (x,y) 

        @param x: x-coordinates for evaluation
        @type x: array-like

        @param y: y-coordinates for evaluation
        @type y: array-like
        """
        return self.energy(x,y)
        

    def set(self, coef):
        """
        Sets the fourier expansion coefficents and calculations the 
        new partation function

        @param coef: expansion coefficents
        @type coef: array like, with shape (4,n,n)
        """
        from numpy import reshape

        self.cc[:, :], self.ss[:, :], self.cs[:, :], self.sc[:, :] = \
                      reshape(coef, (4, self.n, self.n))
        self.normalize()

    def get(self):
        """
        returns current expansion coefficents
        """
        from numpy import array

        return array([self.cc, self.ss, self.cs, self.sc])

    
    def energy(self, x, y=None):
        """
        returns the energy at positions (x,y) 

        @param x: x-coordinates for evaluation
        @type x: array-like

        @param y: y-coordinates for evaluation
        @type y: array-like
        """
        from numpy import sin, cos, dot, multiply

        k = self.periodicities()
        cx, sx = cos(multiply.outer(k, x)), sin(multiply.outer(k, x))
        if y is not None:
            cy, sy = cos(multiply.outer(k, y)), sin(multiply.outer(k, y))
        else:
            cy, sy = cx, sx

        return dot(dot(cx.T, self.cc), cy) + \
               dot(dot(cx.T, self.cs), sy) + \
               dot(dot(sx.T, self.sc), cy) + \
               dot(dot(sx.T, self.ss), sy)

    def sample_weights(self):
        """
        creates a random set of expansion coefficents
        """
        from numpy import add, sqrt
        from numpy.random import standard_normal

        k = self.periodicities()
        k = add.outer(k ** 2, k ** 2)
        self.set([standard_normal(k.shape) / sqrt(k) for i in range(4)])    #@UnusedVariable
        self.normalize(True)

    def prob(self, x, y):
        """
        Returns the probability of the configurations x cross y
        """
        from csb.math import exp
        return exp(-self.beta * self(x, y))


    def z(self):
        """
        Calculate the partion function 
        """
        from scipy.integrate import dblquad
        from numpy import pi

        return dblquad(self.prob, 0., 2 * pi, lambda x: 0., lambda x: 2 * pi)


    def log_z(self, n=500, integration = 'simpson'):
        """
        Calculate the log partion function 
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
                'Choise from trapezoidal and simpson-rule Integration')
        

    def entropy(self, n=500):
        """
        Calculate the entropy of the model

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

        return - average_energy + log_z

    def calculate_statistics(self, data):
        """
        Calculate the sufficient statistics for the data
        """
        from numpy import cos, sin, dot, multiply

        k = self.periodicities()
        cx = cos(multiply.outer(k, data[:, 0]))
        sx = sin(multiply.outer(k, data[:, 0]))
        cy = cos(multiply.outer(k, data[:, 1]))
        sy = sin(multiply.outer(k, data[:, 1]))

        return dot(cx, cy.T), dot(sx, sy.T), dot(cx, sy.T), dot(sx, cy.T)

    def normalize(self, normalize_full=False):
        """
        Remove parameter, which do not have any influence on the model
        and compute the partition function 

        @param normalize_full: compute partition function
        @type normalize_full: boolean
        """
        self.cc[0, 0] = 0.
        self.ss[:, 0] = 0.
        self.ss[0, :] = 0.
        self.cs[:, 0] = 0.
        self.sc[0, :] = 0.

        if normalize_full:
            self.cc[0, 0] = self.log_z()

class Posterior:
    """
    Object to hold and calcuate the posterior (log)probability
    given an exponential family model and correspondingdata
    """

    def __init__(self, model, data):

        from numpy import array

        self.model = model
        self.data = array(data)
        self.stats = self.model.calculate_statistics(self.data)

        self.likelihoods = []

    def __call__(self, weights=None, n=100):
        """
        Returns the log posterior

        @param weights: optional expansion coeffients of the model,
                         if none are specified those of the model are used
        @param n: number of integration point for calculating the partition function
        """
        from numpy import sum

        if weights is not None:
            self.model.set(weights)

        a = sum(self.stats[0] * self.model.cc)
        b = sum(self.stats[1] * self.model.ss)
        c = sum(self.stats[2] * self.model.cs)
        d = sum(self.stats[3] * self.model.sc)

        log_z = self.data.shape[0] * self.model.log_z(n=n)

        likelihood = -self.model.beta * (a + b + c + d) - log_z

        self.likelihoods.append(likelihood)
        return likelihood


if __name__ == '__main__':
    from scipy.optimize import fmin_powell

    aa, ss = 'ALA', 'all'
    k = 2
    data = load_data(aa, ss)[:10000]
    n = len(data)
    # Setup model
    model = MaxentModel(k)
    model.sample_weights()
    posterior = Posterior(model, data)
    x = model.get() * 1.
    x0 = posterior.model.get().flatten()
    target = lambda w: -posterior(w)
    # Learn parameters which minimze posterior
    x = fmin_powell(target, x0)
    posterior.model.set(x)
    posterior.model.normalize(True)

