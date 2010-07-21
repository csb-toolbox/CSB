"""
@note: doc in progress

"""

def geometric_mean(x, axis=None):
    """
    @param x:
    @param axis: compute the geometric mean along this axis 

    @return: geometric mean of x
    """
    from numpy import exp, log, clip, mean

    return exp(mean(log(clip(x, 1e-300, 1e300)), axis))

def harmonic_mean(x, axis=None):
    """
    @param x:
    @param axis: compute the harmonic mean along this axis 

    @return: harmonic mean of x
    """
    from numpy import mean

    return 1 / mean(1 / x, axis)

def kurtosis(x, axis=None):
    """
    @param x: random variables
    @param axis: compute the kurtosis along this axis

    @return: Sample kurtosis of x
    """
    from numpy import mean, std

    m = x.mean(axis)
    a = mean((x-m)**4,axis)
    s = std(x,axis)

    return a / s**4 - 3

def skewness(x, axis=None):
    """
    @param x: random variables
    @param axis: compute the skewness along this axis

    @return: Sample skewness of x
    """
    from numpy import mean, std

    s = std(x)
    return mean((x-x.mean())**3,axis) / s**3
    
def autocorrelation(x,n):
    """
    auto-correlation of a times series

    @param x: time series
    @type x: numpy.array
    @param n: Maximal lag for which to compute the auto-correlation
    @type n: int
    """
    from numpy import array, mean
    x = x - x.mean()
    return array([mean(x[i:]*x[:len(x)-i]) for i in range(n)])

def probabilistic_and(p, axis=0):
    """
    probabilistic version of AND
    """
    from numpy import array, multiply
    return multiply.reduce(array(p), axis=axis)

def probabilistic_or(p, axis=0):
    """
    probabilistic version of OR
    """
    from numpy import array
    return 1 - probabilistic_and(1-array(p), axis)

def probabilistic_xor(p, axis=0):
    """
    probabilistic version of XOR
    works only for axis=0
    """
    from numpy import array

    p = array(p)
    p_not = 1 - p
    
    P = []

    for i in range(p.shape[axis]):
        x = p_not * 1
        x[i] = p[i]
        P.append(probabilistic_and(x,0))

    return probabilistic_or(P,0)

def principal_coordinates(D, nd = None):
    """
    Reconstruction of a multidimensional configuration that
    optimally reproduces the input distance matrix.

    See: Gower, J (1966)
    """
    from numpy import clip, sqrt, take, argsort, sort
    from mulch.math import reverse
    from scipy.linalg import eigh
    
    ## calculate centered similarity matrix
    
    B = - clip(D2, 1e-150, 1e150)**2 / 2.

    b = B.mean(0)

    B = B - b
    B = (B.T - b).T
    B+= b.mean()

    ## calculate spectral decomposition

    v, U = eigh(B)
    v = v.real
    U = U.real

    U = take(U,argsort(v),1)
    v = sort(v)

    U = reverse(U,1)
    v = reverse(v)
    
    if nd is None: nd = len(v)

    X = U[:,:nd] * sqrt(clip(v[:nd], 0., 1e300))

    return X

def entropy(p):
    """
    Calculates the entropy of p 
    @return: entropy of p
    """
    from csb.math import log
    from numpy import sum
    
    return -sum(p * log(p))

def histogram2D(x, nbins=100, axes=None, nbatch=1000, normalize=True):
    """
    Non-greedy histogram two-dimemsional histogram

    @param x: input array of rank two
    @type x: numpy array
    @param nbins: number of bins
    @type nbins: integer
    @param axes: x- and y-axes used for binning the data (if provided this will be used instead of <nbins>)
    @type axes: tuple of two one-dimensional numpy arrays
    @param nbatch: size of batch that is used to sort the data into the 2D grid
    @type nbatch: integer
    @param normalize: specifies whether histogram should be normalized
    @type normalize: boolean

    @return: 2-rank array storing histogram, tuple of x- and y-axis
    """
    from numpy import linspace, zeros, argmin, fabs, subtract, transpose
    
    if axes is None:
        
        lower, upper = x.min(0), x.max(0)
        axes = [linspace(lower[i], upper[i], nbins) for i in range(lower.shape[0])]

    H = zeros((len(axes[0]),len(axes[1])))

    while len(x):

        y = x[:nbatch]
        x = x[nbatch:]

        I = transpose([argmin(fabs(subtract.outer(y[:,i],axes[i])),1) for i in range(2)])

        for i, j in I: H[i,j] += 1

    if normalize:
        H = H / H.sum() / (axes[0][1]-axes[0][0]) / (axes[1][1]-axes[1][0])

    return H, axes

def histogram_nd(x, nbins=100, axes=None, nbatch=1000, normalize=True):
    """
    Non-greedy histogram n-dimemsional histogram

    @param x: input array of rank (-1,n)
    @type x: numpy array
    @param nbins: number of bins
    @type nbins: integer
    @param axes: axes used for binning the data (if provided this will be used instead of <nbins>)
    @type axes: tuple of two one-dimensional numpy arrays
    @param nbatch: size of batch that is used to sort the data into the nD grid
    @type nbatch: integer
    @param normalize: specifies whether histogram should be normalized
    @type normalize: boolean

    @return: n-rank array storing histogram, tuple of axes
    """
    import numpy as np

    d = x.shape[1]
    
    if axes is None:
        
        lower, upper = x.min(0), x.max(0)
        axes = [np.linspace(lower[i], upper[i], nbins) for i in range(d)]

    shape = tuple(map(len, axes))

    H = np.zeros(shape)
    s = np.multiply.accumulate(np.array((1,) + H.shape[:-1]))[::-1]
    H = H.flatten()
    
    while len(x):

        y = x[:nbatch]
        x = x[nbatch:]

        I = np.transpose([np.argmin(np.fabs(np.subtract.outer(y[:,i],axes[i])),1)
                          for i in range(d)])
        I = np.dot(I,s)
        I = np.sort(I)

        i = list(set(I.tolist()))
        n = np.equal.outer(I,i).sum(0)
        
        H[i] += n
        
    if normalize:
        H = H / H.sum() / np.multiply.reduce([axes[i][1]-axes[i][0] for i in range(d)])

    H = np.reshape(H, shape)

    return H, axes

def density(x, nbins, normalize=True):
    """
    Histogram of univariate input data: basically calls numpy's histogram method and
    does a proper normalization.

    @param x: input numpy array
    @param nbins: number of bins
    @type nbins: integer
    @param normalize: if true, histogram will be normalized
    """
    from numpy import histogram
    
    hy, hx = histogram(x, nbins)
    hx = 0.5 * (hx[1:] + hx[:-1])
    hy = hy.astype('d')
    if normalize:
        hy /= (hx[1]-hx[0]) * hy.sum()

    return hx, hy

if __name__ == '__main__':

    from numpy import array

    raise

    for f in [probabilistic_and, probabilistic_or, probabilistic_xor]:

        print f

        for i in range(2):
            for j in range(2):
                p = array([i,j])
                print i, j, f(p)

