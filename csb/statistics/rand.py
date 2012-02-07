# TODO test sampling non-visually

def probability_transform(shape, inv_cum, cum_min=0., cum_max=1.):
    """
    Generic sampler based on the probability transform.

    @param shape: shape of the random sample
    @param inv_cum: inversion of the cumulative density function from which one seeks to sample
    @param cum_min: lower value of the cumulative distribution
    @param cum_max: upper value of the cumulative distribution
    @return: random variates of the PDF implied by the inverse cumulative distribution
    """
    from numpy.random import random
    
    return inv_cum(cum_min + random(shape) * (cum_max - cum_min))

def truncated_gamma(shape=None, alpha=1., beta=1., x_min=None, x_max=None):
    """
    Generates random variates from a lower-and upper-bounded gamma distribution

    @param shape: shape of the random sample
    @param alpha: shape parameter (alpha > 0.)
    @param beta:  scale parameter (beta >= 0.)
    @param x_min: lower bound of variate
    @param x_max: upper bound of variate    
    @return: random variates of lower-bounded gamma distribution
    """
    from scipy.special import gammainc, gammaincinv
    from numpy.random import gamma
    from numpy import inf

    if x_min is None and x_max is None:
        return gamma(alpha, 1/beta, shape)
    elif x_min is None:
        x_min = 0.
    elif x_max is None:
        x_max = inf
        
    x_min = max(0.,x_min)
    x_max = min(1e300,x_max)

    a = gammainc(alpha, beta*x_min)
    b = gammainc(alpha, beta*x_max)

    return probability_transform(shape,
                                 lambda x,alpha=alpha: gammaincinv(alpha, x),
                                 a, b) / beta


    
def truncated_normal(shape=None, mu=0., sigma=1., x_min=None, x_max=None):
    """
    Generates random variates from a lower-and upper-bounded normal distribution

    @param shape: shape of the random sample
    @param mu:    location parameter 
    @param sigma: width of the distribution (sigma >= 0.)
    @param x_min: lower bound of variate
    @param x_max: upper bound of variate    
    @return: random variates of lower-bounded normal distribution
    """
    from scipy.special import erf, erfinv
    from numpy.random import standard_normal
    from numpy import inf, sqrt

    if x_min is None and x_max is None:
        return standard_normal(shape) * sigma + mu
    elif x_min is None:
        x_min = -inf
    elif x_max is None:
        x_max = inf
        
    x_min = max(-1e300, x_min)
    x_max = min(+1e300,x_max)
    var   = sigma**2 + 1e-300
    sigma = sqrt(2 * var)
    
    a = erf((x_min-mu)/sigma)
    b = erf((x_max-mu)/sigma)

    return probability_transform(shape, erfinv, a, b) * sigma + mu


def sample_dirichlet(alpha, n_samples = 1):
    """
    sample points from a dirichlet distribution with parameter alpha
    """
    from numpy import array, sum, transpose, ones
    from numpy.random import gamma

    alpha = array(alpha,ndmin=1)
    X = gamma(alpha,
              ones(len(alpha)),
              [n_samples,len(alpha)])
     
    return transpose(transpose(X) / sum(X,-1))


def sample_sphere3d(radius = 1. ,n_samples = 1 ):
    """
    sample points from 3D sphere

    @param radius: radius of the sphere
    @type radius: float

    @param n_samples: number of samples to return
    @type n_samples: int

    @return: n_samples times random cartesian coordinates inside the sphere
    @rtype: numpy array
    """
    from numpy.random  import random
    from numpy import arccos, transpose, cos, sin, pi, power

    r = radius * power(random(n_samples), 1/3.)
    theta = arccos(2. * (random(n_samples) - 0.5))
    phi = 2 * pi * random(n_samples)

    x = cos(phi) * sin(theta) * r
    y = sin(phi) * sin(theta) * r
    z = cos(theta) * r

    return transpose([x, y, z])



def sample_from_histogram(p, n_samples = 1):
    """
    returns the indice of bin according to the histogram p

    @param p: histogram
    @type p: numpy.array
    @param n_samples: number of samples to generate
    @type n_samples: integer
    """
    
    from numpy import add, less, argsort, take, arange
    from numpy.random import random

    indices = argsort(p)
    indices = take(indices, arange(len(p)-1,-1,-1))

    c = add.accumulate(take(p,indices)) / add.reduce(p)

    return indices[add.reduce(less.outer(c, random(n_samples)),0)]


def test_truncated_gamma(alpha = 2., beta  = 1., x_min = 0.1, x_max = 5.):

    from numpy import histogram, linspace

    from csb.math import exp, log_sum_exp, log
    
    x = truncated_gamma(10000, alpha, beta, x_min, x_max)

    hy, hx = histogram(x, 100)
    hx = 0.5*(hx[1:]+hx[:-1])
    hy = hy.astype('d')
    hy /= (hx[1]-hx[0]) * hy.sum()
    
    x = linspace(x_min, x_max, 1000)
    p = (alpha - 1) * log(x) - beta * x
    p -= log_sum_exp(p)
    p = exp(p) / (x[1]-x[0])
    


def test_truncated_normal(mu = 2., sigma  = 1., x_min = 0.1, x_max = 5.):

    from numpy import histogram, linspace

    from csb.math import exp, log_sum_exp
    
    x = truncated_normal(10000, mu, sigma, x_min, x_max)

    hy, hx = histogram(x, 100)
    hx = 0.5*(hx[1:]+hx[:-1])
    hy = hy.astype('d')
    hy /= (hx[1]-hx[0]) * hy.sum()
    
    x = linspace(mu - 5 * sigma, mu + 5 * sigma, 1000)
    p = - 0.5 * (x-mu)**2 / sigma**2
    p -= log_sum_exp(p)
    p = exp(p) / (x[1]-x[0])
    



def _sample_beta(kappa, n=1):
    from numpy import arccos,clip
    from csb.math import log, exp
    from numpy.random import random
    
    u = random(n)

    if kappa != 0.:
        x = clip(1 + 2 * log(u + (1-u) * exp(-kappa)) / kappa,-1.,1.)
    else:
        x = 2*u - 1

    if n == 1:
        return arccos(x)[0]
    else:
        return arccos(x)


def random_rotation(A, n_iter=10, initial_values = None):
    """
    generate random rotation R from

    exp(trace(dot(transpose(A),R)))
    
    """
    from numpy import cos, sin, dot, pi, clip
    from numpy.linalg import svd, det    
    from random import vonmisesvariate, randint
    from csb.math import euler

    U, L, V = svd(A)

    if det(U) < 0:
        L[2] *= -1
        U[:,2] *= -1
    if det(V) < 0:
        L[2] *= -1
        V[2] *= -1

    if initial_values is None:
        beta = 0.
    else:
        alpha, beta, gamma = initial_values

    for _i in range(n_iter):

        ## sample alpha and gamma
        phi = vonmisesvariate(0.,clip(cos(beta/2)**2 * (L[0]+L[1]),1e-308,1e10))
        psi = vonmisesvariate(pi, sin(beta/2)**2 * (L[0]-L[1]))
        u   = randint(0,1)
        
        alpha = 0.5 * (phi + psi) + pi * u
        gamma = 0.5 * (phi - psi) + pi * u

        ## sample beta
        kappa = cos(phi) * (L[0]+L[1]) + cos(psi) * (L[0]-L[1]) + 2 * L[2]
        beta  = _sample_beta(kappa)

    return dot(U, dot(euler(alpha, beta, gamma), V))
