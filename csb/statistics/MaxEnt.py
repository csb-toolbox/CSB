
def load_data(aa, ss='all', path='~/data/dunbrack', isd=True):
    """
    Load raw rotamer data

    @param aa: Three-letter aminoacid code
    @param ss: Secondary structure element (H,E,all)
    @param path: search path for data files
    @paran isd map
    """
    from numpy import array, compress, transpose, where, pi, fmod
    import sys, os
    from csb.math import degree2radian

    sys.path.insert(0, os.path.expanduser(path))
    from pystartup import Load

    data = array(Load(path + '/%s_%s' % (aa, ss)).tolist())
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


def radian2torsionDeg(radianAngle):
    from numpy import pi
    ang = radianAngle
    if ang > 2. * pi or ang < 0.:
        raise ValueError, 'radian angle must be 0<ang<2pi, value: %f' % ang
    if ang > pi:
        ang = -((2. * pi) - ang)
    return ang * 180. / pi

def torsionDeg2radian(torsionDegAngle):
    from numpy import pi
    ang = torsionDegAngle
    if ang > 181. or ang < -181.:
        raise ValueError, 'degree angle must be -180<ang<180, value: %f' % ang

    if ang < 0:
        ang = 360. + ang
    return ang * pi / 180.


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

    def show(self):
        """
        Displays the probability distribution of a model in 3D
        
        """

        from numpy import arange, meshgrid, pi, array
        from csb.math import exp
        try:
            import matplotlib.cm as cm
            import matplotlib.pyplot as plt
        except:
            raise ImportError('Matplotlib not installed')

        delta = 0.025
        x = y = arange(0, 2 * pi, delta)
        X, Y = meshgrid(x, y)
        Z = array([self(xx, yy) for xx, yy in zip(X.ravel(), Y.ravel())])
        Z = Z.reshape((len(x), len(x)))
        im = plt.imshow(exp(-Z), interpolation='bilinear', cmap=cm.jet,
                        origin='lower', extent=[0, 2 * pi, 0, 2 * pi])
        plt.show()
        return

    def show3D(self):
        from numpy import arange, meshgrid, pi, array
        try:
            import matplotlib.cm as cm
            import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d import Axes3D
        except:
            raise ImportError('Matplotlib not installed')

        delta = 0.025
        x = y = arange(0, 2 * pi, delta)
        X, Y = meshgrid(x, y)
        z = array([self(xx, yy) for xx, yy in zip(X.ravel(), Y.ravel())])
        z = z.reshape((len(x), len(x)))
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.contourf(X, Y, clip(-z, -6.26, 20), 100, extend3d=True)
        fig.show()
        return


    def load_old(self,
                 aa, path='~/mnt/tay/mechelke/projects/isd/toppar/maxent_potentials'):
        import os
        from numpy import arange, zeros
        from csb.io import Load


        params, energies = eval(open(os.path.expanduser(path)).read())
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

    def load(self, aa,
             path='~/mnt/tay/mechelke/projects/isd/toppar/maxent_potentials_v3'):
        import os
        from numpy import arange, zeros, reshape, array
        from csb.io import Load

        path = os.path.expanduser(path)
        params, energies = Load(path)
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

        from numpy import sin, cos, concatenate, dot

        k = self.periodicities()
        cx, sx = cos(k * x), sin(k * x)
        cy, sy = cos(k * y), sin(k * y)
        return dot(cx, dot(self.cc, cy)) + dot(sx, dot(self.ss, sy)) + \
               dot(sx, dot(self.sc, cy)) + dot(cx, dot(self.cs, sy))

    def set(self, params):

        from numpy import reshape

        self.cc[:, :], self.ss[:, :], self.cs[:, :], self.sc[:, :] = \
                      reshape(params, (4, self.n, self.n))
        self.normalize()

    def get(self):

        from numpy import array

        return array([self.cc, self.ss, self.cs, self.sc])

    def energy(self, x, y=None):

        from numpy import sin, cos, concatenate, dot, multiply

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

        from numpy import add, sqrt, concatenate
        from numpy.random import standard_normal

        k = self.periodicities()
        k = add.outer(k ** 2, k ** 2)

        self.set([standard_normal(k.shape) / sqrt(k) for i in range(4)])
        self.normalize(True)

    def prob(self, x, y):
        from csb.math import exp
        return exp(-self.beta * self(x, y))

    def Z(self):

        from scipy.integrate import dblquad
        from numpy import pi

        return dblquad(self.prob, 0., 2 * pi, lambda x: 0., lambda x: 2 * pi)

    def logZ2(self, n=500):

        from numpy import pi, linspace, sum, array
        from csb.math import log_sum_exp, log, exp

        x = linspace(0., 2 * pi, n)
        dx = x[1] - x[0]
        E = self.energy(x).flatten()
        return log_sum_exp(-self.beta * E) + 2 * log(dx)

    def logZ3(self, n=500):
        """
        using trapezoidal rule
        """
        from numpy import pi, linspace, sum, array
        from csb.numeric import trapezoidal2D
        from csb.math import log_sum_exp, log, exp

        x = linspace(0., 2 * pi, n)
        dx = x[1] - x[0]

        f = -self.beta * self.energy(x)
        f_max = f.max()
        f -= f_max

        I = trapezoidal2D(exp(f))

        return log(I) + f_max + 2 * log(dx)

    def logZ(self, n=500):
        """
        using Simpson rule
        """
        from numpy import pi, linspace, sum, array
        from csb.numeric import simpson2D
        from csb.math import log_sum_exp, log, exp

        x = linspace(0., 2 * pi, 2 * n + 1)
        dx = x[1] - x[0]

        f = -self.beta * self.energy(x)
        f_max = f.max()
        f -= f_max

        I = simpson2D(exp(f))

        return log(I) + f_max + 2 * log(dx)

    def entropy(self, n=500):
        from csb.numeric import trapezoidal2D
        from numpy import pi, linspace, sum, array
        from csb.math import log_sum_exp, log, exp

        x = linspace(0., 2 * pi, n)
        dx = x[1] - x[0]

        f = -self.beta * self.energy(x)
        f_max = f.max()

        logZ = log(trapezoidal2D(exp(f - f_max))) + f_max + 2 * log(dx)
        E_av = trapezoidal2D(f * exp(f - f_max)) * exp(f_max + 2 * log(dx) - logZ)

        return - E_av + logZ

    def calculate_statistics(self, data):

        from numpy import cos, sin, dot, multiply

        k = self.periodicities()

        cx = cos(multiply.outer(k, data[:, 0]))
        sx = sin(multiply.outer(k, data[:, 0]))
        cy = cos(multiply.outer(k, data[:, 1]))
        sy = sin(multiply.outer(k, data[:, 1]))

        return dot(cx, cy.T), dot(sx, sy.T), dot(cx, sy.T), dot(sx, cy.T)

    def normalize(self, normalize_full=False):

        self.cc[0, 0] = 0.
        self.ss[:, 0] = 0.
        self.ss[0, :] = 0.
        self.cs[:, 0] = 0.
        self.sc[0, :] = 0.

        if normalize_full:
            self.cc[0, 0] = self.logZ()

class Posterior:

    def __init__(self, model, data):

        from numpy import array

        self.model = model
        self.data = array(data)
        self.stats = self.model.calculate_statistics(self.data)

        self.L = []

    def __call__(self, weights=None, n=500):

        from numpy import sum

        if weights is not None: self.model.set(weights)

        a = sum(self.stats[0] * self.model.cc)
        b = sum(self.stats[1] * self.model.ss)
        c = sum(self.stats[2] * self.model.cs)
        d = sum(self.stats[3] * self.model.sc)

        logZ = self.data.shape[0] * self.model.logZ(n=n)

        L = -self.model.beta * (a + b + c + d) - logZ

        self.L.append(L)

        return L

    def grad(self, weights=None):

        from numpy import sum

        if weights is not None: self.model.set(weights)

class MetropolisSampler:

    def __init__(self, posterior):

        self.posterior = posterior
        self.eps = 1e-1
        self.E = None
        self.accept = []
        self.beta = 1.0

    def propose(self, x):

        from numpy.random import random

        return x + 2 * (random(x.shape) - 0.5) * self.eps

    def sample(self, x):

        from numpy.random import random
        from numpy import exp

        if self.E is None:
            E = -self.posterior(x)
        else:
            E = self.E * 1.

        y = self.propose(x)
        E2 = -self.posterior(y)

        dE = E2 - E

        accept = random() < exp(-self.beta * dE)
        self.accept.append(1. * int(accept))

        if accept:
            self.E = E2
            return y

        else:
            self.E = E
            return x

    def generate_samples(self, x, n_samples=100):

        self.samples = []
        self.accept = []
        self.energies = []
        self.E = -self.posterior(x)

        for i in range(n_samples):

            y = self.sample(x)
            print i, self.accept[-1]
            self.samples.append(y)
            self.energies.append(self.E)


if __name__ == '__main__':
    from numpy import *
    from csb.statistics import Posterior, MetropolisSampler, MaxentModel
    from csb.statistics.MaxEnt import load_data
    from scipy.optimize import fmin_powell
    aa, ss = 'ALA', 'H'
    k = 5
    data = load_data(aa, ss)
    n = len(data)
    # Setup model
    model = MaxentModel(k)
    model.sample_weights()
    posterior = Posterior(model, data)
    x = model.get() * 1.
    x0 = posterior.model.get().flatten()
    target = lambda w:-posterior(w)
    # Learn parameters which minimze posterior
    x = fmin_powell(target, x0)
    posterior.model.set(x)
    posterior.model.normalize(True)
    z = posterior.model.show()
