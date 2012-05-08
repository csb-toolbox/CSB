"""
Computational utility functions. 
"""

import numpy
import numpy.random


def fit(X, Y):
    """
    returns the translation vector and the rotation matrix
    minimizing the RMSD between two sets of vectors, i.e.
    if 
    >>> R,t = fit(X,Y)
    then
    >>> Y = dot(Y, transpose(R)) + t
    will be the fitted configuration.

    @param X: 3 x n input vector
    @type X: numpy array

    @param Y: 3 x n  input vector
    @type Y: numpy array

    @return: 3 x 3 rotation matrix and 3 x 1 translation vector
    @rtype: tuple
    """

    from numpy.linalg import svd, det
    from numpy import dot

    ## center configurations

    x = X.mean(0)
    y = Y.mean(0)

    ## SVD of correlation matrix

    V, _L, U = svd(dot((X-x).T, Y-y))

    ## calculate rotation and translation

    R = dot(V, U)

    if det(R) < 0.:
        U[-1] *= -1
        R = dot(V, U)

    t = x - dot(R, y)

    return R, t

def wfit(X, Y, w):
    """
    returns the translation vector and the rotation matrix
    minimizing the weighted RMSD between two sets of vectors, i.e.
    if 
    >>> R,t = wfit(X,Y,w)
    then
    >>> Y = dot(Y, transpose(R)) + t
    will be the fitted configuration.

    @param X: 3 x n input vector
    @type X: numpy array

    @param Y: 3 x n  input vector
    @type Y: numpy array

    @param w: input weights
    @type w: numpy array

    @return: 3 x 3 rotation matrix and 3 x 1 translation vector
    @rtype: tuple
    """

    from numpy.linalg import svd, det
    from numpy import dot, transpose, average

    ## center configurations

    norm = sum(w)
    x = dot(w, X) / norm
    y = dot(w, Y) / norm

    ## SVD of correlation matrix

    V, _L, U = svd(dot(transpose(X - x) * w, Y - y))

    ## calculate rotation and translation

    R = dot(V, U)

    if det(R) < 0.:
        U[2] *= -1
        R = dot(V, U)

    t = x - dot(R, y)

    return R, t

def probabilistic_fit(X,Y, w = None, niter = 10):
    from csb.statistics.rand import random_rotation
    from numpy import dot, transpose, average
    
    """
    Generates a superposition of X,Y where
    R ~ exp(trace(dot(transpose(dot(transpose(X-t),Y)),R)))
    t ~ N(t_opt, 1/sqrt(N))
    """
    if w is None:
        R, t = fit(X,Y)
    else:
        R, t = wfit(X,Y,w)
            
    N = len(X)

    for i in range(niter):
        ## sample rotation
        if w is None:
            A = dot(transpose(X-t),Y)
        else:
            A = dot(transpose(X-t)*w,Y)

        R = random_rotation(A)

        ## sample translation (without prior so far)
        if w is None:
            mu = average(X - dot(Y, transpose(R)),0)
            t = numpy.random.standard_normal(len(mu)) / numpy.sqrt(N) + mu
        else:
            mu = dot(w, X - dot(Y, transpose(R))) / numpy.sum(w)
            t = numpy.random.standard_normal(len(mu)) / numpy.sqrt(numpy.sum(w)) + mu

    return R,t
    
    

def fit_wellordered(X, Y, n_iter=None, n_stdv=2, tol_rmsd=.5,
                    tol_stdv=0.05, full_output=False):
    """
    Matches two arrays onto each other by iteratively throwing out
    highly deviating entries

    (Reference: Nilges et al.: Delineating well-ordered regions in
    protein structure ensembles).


    @param X: 3 x n input vector
    @type X: numpy array

    @param Y: 3 x n  input vector
    @type Y: numpy array

    @param n_stdv: number of standard deviations above which points are considered to be outliers

    @param tol_rmsd: tolerance in rmsd

    @param tol_stdv: tolerance in standard deviations

    @param full_output: also return full history of values calculated by the algorithm
    """

    from numpy import ones, compress, dot, transpose, sqrt, sum, nonzero, std, average


    rmsd_list = []

    rmsd_old = 0.
    stdv_old = 0.

    n = 0
    converged = False

    mask = ones(X.shape[0])

    while not converged:

        ## find transformation for best match

        R, t = fit(compress(mask, X, 0), compress(mask, Y, 0))

        ## calculate RMSD profile

        d = sqrt(sum((X - dot(Y, transpose(R)) - t) ** 2, 1))

        ## calculate rmsd and stdv

        rmsd = sqrt(average(compress(mask, d) ** 2, 0))
        stdv = std(compress(mask, d))

        ## check conditions for convergence

        if stdv < 1e-10: break

        d_rmsd = abs(rmsd - rmsd_old)
        d_stdv = abs(1 - stdv_old / stdv)

        if d_rmsd < tol_rmsd:
            if d_stdv < tol_stdv:
                converged = 1
            else:
                stdv_old = stdv
        else:
            rmsd_old = rmsd
            stdv_old = stdv

        ## store result

        perc = average(1. * mask)

        ## throw out non-matching rows

        new_mask = mask * (d < rmsd + n_stdv * stdv)
        outliers = nonzero(mask - new_mask)
        rmsd_list.append([perc, rmsd, outliers])

        mask = new_mask

        if n_iter and n >= n_iter: break

        n += 1

    if full_output:
        return (R, t), rmsd_list
    else:
        return (R, t)


def rmsd(X, Y):
    """
    Calculate the root mean squared deviation (RMSD) using Kabsch' formula 

    @param X: 3 x n input vector
    @type X: numpy array
    @param Y: 3 x n  input vector
    @type Y: numpy array

    @return: rmsd value between the input vectors
    @rtype: float
    """

    from numpy import sum, dot, transpose, sqrt, clip, average
    from numpy.linalg import svd

    X = X - average(X, 0)
    Y = Y - average(Y, 0)

    R_x = sum(X ** 2)
    R_y = sum(Y ** 2)

    L = svd(dot(transpose(Y), X))[1]

    return sqrt(clip(R_x + R_y - 2 * sum(L), 0., 1e300) / len(X))

def wrmsd(X, Y, w):
    """
    Calculate the weighted root mean squared deviation (wRMSD) using Kabsch' formula 

    @param X: 3 x n input vector
    @type X: numpy array
    @param Y: 3 x n  input vector
    @type Y: numpy array
    @param w: input weights
    @type w: numpy array

    @return: rmsd value between the input vectors
    @rtype: float
    """

    from numpy import sum, dot, transpose, sqrt, clip, average
    from numpy.linalg import svd

    ## normalize weights

    w = w / w.sum()

    X = X - dot(w,X)
    Y = Y - dot(w,Y)

    R_x = sum(X.T**2 * w)
    R_y = sum(Y.T**2 * w)

    L = svd(dot(transpose(Y)*w, X))[1]

    return sqrt(clip(R_x + R_y - 2 * sum(L), 0., 1e300))

def torsion_rmsd(x, y):
    """
    Compute the circular RMSD of two phi/psi angle sets.
    
    @param x: query phi/psi angles (Nx2 array, in radians)
    @type x: array
    @param y: subject phi/psi angles (Nx2 array, in radians)
    @type y: array
    
    @rtype: float
    """
    from numpy import array, sin, cos, sqrt
    
    phi, psi = (x - y).T
    assert len(phi) == len(psi)
    
    r = sin(phi).sum()**2 + cos(phi).sum()**2 + sin(psi).sum()**2 + cos(psi).sum()**2
    return 1 - (1.0 / len(phi)) * sqrt(r / 2.0)

def _tm_d0(Lmin):
    
    from numpy import power
     
    if Lmin > 15:
        d0 = 1.24 * power(Lmin - 15.0, 1.0 / 3.0) - 1.8
    else:
        d0 = 0.5

    return d0

def tm_score(x, y):
    """
    Evaluate the TM-score of two conformations as they are (no fitting is done)
    
    @param x: 3 x n input array
    @param x: 3 x n input array
    
    @return: computed TM-score
    @rtype: float
    """
    from numpy import sum

    d = sum((x - y) ** 2, 1) ** 0.5
    d0 = _tm_d0(len(x))
    
    return sum(1 / (1 + (d / d0) ** 2)) / len(x)

def tm_superimpose(x, y, fit_method=fit):
    """
    Compute the TM-score of two protein coordinate vector sets. 
    Reference:  hang.bioinformatics.ku.edu/TM-score
    
    @param x: 3 x n input vector
    @type x: numpy.array
    @param y: 3 x n  input vector
    @type y: numpy.array
    @param fit_method: a reference to a proper fitting function, e.g. fit
                       or fit_wellordered    
    @type fit_method: function
    
    @return: rotation matrix, translation vector, TM-score
    @rtype: tuple
    """
    from numpy import array, sum, dot, compress, power, ones, argmax

    mask = ones(len(x)).astype('i')
    d0 = _tm_d0(len(x))
    
    scores = []
    transformations = []

    for i in range(3):

        R, t = fit_method(compress(mask, x, 0), compress(mask, y, 0))

        scores.append(tm_score(x, dot(y, R.T) + t))
        transformations.append((R, t))

        cutoff = min(max(4.5, d0), 8.0) + i
        d = sum((x - dot(y, R.T) - t) ** 2, 1) ** 0.5

        while True:
            
            mask = (d < cutoff).astype('i')
            if sum(mask) >= 3 or 3 >= len(mask):
                break
            cutoff += 0.5

    R, t = transformations[argmax(scores)]
    score = max(scores)
    
    return R, t, score

def center_of_mass(x, m=None):
    """
    Compute the mean of a set of (optionally weighted) points

    @param x: array of rank (n,d) where n is the number of points
              and d the dimension
    @type x: numpy.array
    @param m: rank (n,) array of masses / weights
    @type m: numpy.array

    @return: center of mass
    @rtype: (d,) numpy.array
    """
    if m is None:
        return x.mean(0)
    else:
        from numpy import dot

        return dot(m, x) / m.sum()
    
def radius_of_gyration(x, m=None):
    """
    Compute the radius of gyration of a set of (optionally weighted) points

    @param x: array of rank (n,d) where n is the number of points
              and d the dimension
    @type x: numpy.array
    @param m: rank (n,) array of masses / weights
    @type m: numpy.array

    @return: center of mass
    @rtype: (d,) numpy.array
    """
    from numpy import sqrt, dot
    
    x = x - center_of_mass(x, m)
    if m is None:
        return sqrt((x**2).sum() / len(x))
    else:
        return sqrt(dot(m,(x**2).sum(1)) / m.sum())

def second_moments(x, m=None):
    """
    Compute the tensor second moments of a set of (optionally weighted) points

    @param x: array of rank (n,d) where n is the number of points
              and d the dimension
    @type x: numpy.array
    @param m: rank (n,) array of masses / weights
    @type m: numpy.array

    @return: second moments
    @rtype: (d,d) numpy.array
    """
    from numpy import dot

    x = (x - center_of_mass(x, m)).T
    
    if m is not None:
        return dot(x*m, x.T)
    else:
        return dot(x, x.T)

def inertia_tensor(x, m=None):
    """
    Compute the inertia tensor of a set of (optionally weighted) points

    @param x: array of rank (n,d) where n is the number of points
              and d the dimension
    @type x: numpy.array
    @param m: rank (n,) array of masses / weights
    @type m: numpy.array

    @return: inertia tensor
    @rtype: (d,d) numpy.array
    """
    from numpy import dot, eye

    x = (x - center_of_mass(x, m)).T
    r2= (x**2).sum(0)
    
    if m is not None:
        return eye(x.shape[0]) * dot(m,r2) - dot(x*m,x.T)
    else:
        return eye(x.shape[0]) * r2.sum() - dot(x, x.T)

def distance_matrix(X, Y=None):
    """
    Calculates a matrix of pairwise distances

    @param X: m x n input vector
    @type X: numpy array

    @param Y: k x n input vector or None, which defaults to Y=X
    @type Y: numpy array

    @return: m x k distance matrix
    @rtype: numpy array
    """
    from numpy import add, clip, sqrt, dot, transpose, sum

    if Y is None: Y = X

    if X.ndim < 2: X = X.reshape((1, -1))
    if Y.ndim < 2: Y = Y.reshape((1, -1))

    C = dot(X, transpose(Y))
    S = add.outer(sum(X**2, 1), sum(Y**2, 1))

    return sqrt(clip(S - 2 * C, 0., 1e300))

# vi:expandtab:smarttab:sw=4
