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
    from csb.math import log
    from numpy.random import random, gamma
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

    return probability_transform(shape, lambda x,alpha=alpha: gammaincinv(alpha, x), a, b) / beta

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
    from csb.math import log
    from numpy.random import random, standard_normal
    from numpy import inf, sqrt

    if x_min is None and x_max is None:
        return standard_normal(shape) * sigma + mu
    elif x_min is None:
        x_min =-inf
    elif x_max is None:
        x_max = inf
        
    x_min = max(-1e300, x_min)
    x_max = min(+1e300,x_max)
    var   = sigma**2 + 1e-300
    sigma = sqrt(2 * var)
    
    a = erf((x_min-mu)/sigma)
    b = erf((x_max-mu)/sigma)

    return probability_transform(shape, erfinv, a, b) * sigma + mu

def test_truncated_gamma(alpha = 2., beta  = 1., x_min = 0.1, x_max = 5.):

    from numpy import histogram, linspace
    from gnuplot import plot
    from scipy.special import gammaln
    from csb.math import exp, log_sum_exp, log
    
    x = truncated_gamma(10000, alpha, beta, x_min, x_max)

    hy, hx = histogram(x, 100)
    hx = 0.5*(hx[1:]+hx[:-1])
    hy = hy.astype('d')
    hy/= (hx[1]-hx[0]) * hy.sum()
    
    x = linspace(x_min, x_max, 1000)
    p = (alpha - 1) * log(x) - beta * x
    p-= log_sum_exp(p)
    p = exp(p) / (x[1]-x[0])
    
    plot(zip(hx,hy),zip(x,p))

def test_truncated_normal(mu = 2., sigma  = 1., x_min = 0.1, x_max = 5.):

    from numpy import histogram, linspace
    from gnuplot import plot
    from scipy.special import gammaln
    from csb.math import exp, log_sum_exp, log
    
    x = truncated_normal(10000, mu, sigma, x_min, x_max)

    hy, hx = histogram(x, 100)
    hx = 0.5*(hx[1:]+hx[:-1])
    hy = hy.astype('d')
    hy/= (hx[1]-hx[0]) * hy.sum()
    
    x = linspace(mu - 5 * sigma, mu + 5 * sigma, 1000)
    p = - 0.5 * (x-mu)**2 / sigma**2
    p-= log_sum_exp(p)
    p = exp(p) / (x[1]-x[0])
    
    plot(zip(hx,hy),zip(x,p))
