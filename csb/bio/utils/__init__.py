"""
Computational utility functions. 
"""

import numpy
import numpy.random

import csb.numeric


def fit(X, Y):
    """
    Return the translation vector and the rotation matrix
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

    V, _L, U = svd(dot((X - x).T, Y - y))

    ## calculate rotation and translation

    R = dot(V, U)

    if det(R) < 0.:
        U[-1] *= -1
        R = dot(V, U)

    t = x - dot(R, y)

    return R, t

def wfit(X, Y, w):
    """
    Return the translation vector and the rotation matrix
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

def probabilistic_fit(X, Y, w=None, niter=10):
    """
    Generate a superposition of X,Y where::
    
        R ~ exp(trace(dot(transpose(dot(transpose(X-t), Y)), R)))
        t ~ N(t_opt, 1 / sqrt(N))
    """
        
    from csb.statistics.rand import random_rotation
    from numpy import dot, transpose, average
    
    if w is None:
        R, t = fit(X, Y)
    else:
        R, t = wfit(X, Y, w)
            
    N = len(X)

    for i in range(niter):
        ## sample rotation
        if w is None:
            A = dot(transpose(X - t), Y)
        else:
            A = dot(transpose(X - t) * w, Y)

        R = random_rotation(A)

        ## sample translation (without prior so far)
        if w is None:
            mu = average(X - dot(Y, transpose(R)), 0)
            t = numpy.random.standard_normal(len(mu)) / numpy.sqrt(N) + mu
        else:
            mu = dot(w, X - dot(Y, transpose(R))) / numpy.sum(w)
            t = numpy.random.standard_normal(len(mu)) / numpy.sqrt(numpy.sum(w)) + mu

    return R, t
    
def fit_wellordered(X, Y, n_iter=None, n_stdv=2, tol_rmsd=.5,
                    tol_stdv=0.05, full_output=False):
    """
    Matche two arrays onto each other by iteratively throwing out
    highly deviating entries.
    
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
    Calculate the root mean squared deviation (RMSD) using Kabsch' formula.

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

def rmsd_cur(X, Y):
    """
    Calculate the RMSD of two conformations as they are (no fitting is done).
    For details, see L{rmsd}.

    @return: rmsd value between the input vectors
    @rtype: float
    """
    return distance_sq(X, Y).mean() ** 0.5

def wrmsd(X, Y, w):
    """
    Calculate the weighted root mean squared deviation (wRMSD) using Kabsch'
    formula.

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

    X = X - dot(w, X)
    Y = Y - dot(w, Y)

    R_x = sum(X.T ** 2 * w)
    R_y = sum(Y.T ** 2 * w)

    L = svd(dot(transpose(Y) * w, X))[1]

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
    
    r = sin(phi).sum() ** 2 + cos(phi).sum() ** 2 + sin(psi).sum() ** 2 + cos(psi).sum() ** 2
    return 1 - (1.0 / len(phi)) * sqrt(r / 2.0)

def _tm_d0(Lmin):
    
    from numpy import power
     
    if Lmin > 15:
        d0 = 1.24 * power(Lmin - 15.0, 1.0 / 3.0) - 1.8
    else:
        d0 = 0.5

    return max(0.5, d0)

def tm_score(x, y, L=None, d0=None):
    """
    Evaluate the TM-score of two conformations as they are (no fitting is done).
    
    @param x: 3 x n input array
    @type x: numpy array
    @param y: 3 x n input array
    @type y: numpy array
    @param L: length for normalization (default: C{len(x)})
    @type L: int
    @param d0: d0 in Angstroms (default: calculate from C{L})
    @type d0: float
    
    @return: computed TM-score
    @rtype: float
    """
    from numpy import sum

    if not L:
        L = len(x)
    if not d0:
        d0 = _tm_d0(L)
    d = distance(x, y)

    return sum(1 / (1 + (d / d0) ** 2)) / L

def tm_superimpose(x, y, fit_method=fit, L=None, d0=None, L_ini_min=4, iL_step=1):
    """
    Compute the TM-score of two protein coordinate vector sets. 
    Reference:  http://zhanglab.ccmb.med.umich.edu/TM-score
    
    @param x: 3 x n input vector
    @type x: numpy.array
    @param y: 3 x n  input vector
    @type y: numpy.array
    @param fit_method: a reference to a proper fitting function, e.g. fit
                       or fit_wellordered    
    @type fit_method: function
    @param L: length for normalization (default: C{len(x)})
    @type L: int
    @param d0: d0 in Angstroms (default: calculate from C{L})
    @type d0: float
    @param L_ini_min: minimum length of initial alignment window (increase
    to speed up but loose precision, a value of 0 disables local alignment
    initialization)
    @type L_ini_min: int
    @param iL_step: initial alignment window shift (increase to speed up
    but loose precision)
    @type iL_step: int
    
    @return: rotation matrix, translation vector, TM-score
    @rtype: tuple
    """
    from numpy import asarray, sum, dot, zeros, clip

    x, y = asarray(x), asarray(y)
    if not L:
        L = len(x)
    if not d0:
        d0 = _tm_d0(L)
    d0_search = clip(d0, 4.5, 8.0)
    best = None, None, 0.0

    L_ini_min = min(L, L_ini_min) if L_ini_min else L
    L_ini = [L_ini_min] + filter(lambda x: x > L_ini_min,
            [L // (2 ** n_init) for n_init in range(6)])

    # the outer two loops define a sliding window of different sizes for the
    # initial local alignment (disabled with L_ini_min=0)
    for L_init in L_ini:
        for iL in range(0, L - L_init + 1, min(L_init, iL_step)):
            mask = zeros(L, bool)
            mask[iL:iL + L_init] = True

            # refine mask until convergence, similar to fit_wellordered
            for i in range(20):
                R, t = fit_method(x[mask], y[mask])

                d = distance(x, dot(y, R.T) + t)
                score = sum(1 / (1 + (d / d0) ** 2)) / L

                if score > best[-1]:
                    best = R, t, score

                cutoff = d0_search + (-1 if i == 0 else 1)
                while True:
                    mask_prev = mask
                    mask = d < cutoff
                    if sum(mask) >= 3 or 3 >= len(mask):
                        break
                    cutoff += 0.5

                if (mask == mask_prev).all():
                    break

    return best

def center_of_mass(x, m=None):
    """
    Compute the mean of a set of (optionally weighted) points.

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
    Compute the radius of gyration of a set of (optionally weighted) points.

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
        return sqrt((x ** 2).sum() / len(x))
    else:
        return sqrt(dot(m, (x ** 2).sum(1)) / m.sum())

def second_moments(x, m=None):
    """
    Compute the tensor second moments of a set of (optionally weighted) points.

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
        return dot(x * m, x.T)
    else:
        return dot(x, x.T)

def inertia_tensor(x, m=None):
    """
    Compute the inertia tensor of a set of (optionally weighted) points.

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
    r2 = (x ** 2).sum(0)
    
    if m is not None:
        return eye(x.shape[0]) * dot(m, r2) - dot(x * m, x.T)
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
    S = add.outer(sum(X ** 2, 1), sum(Y ** 2, 1))

    return sqrt(clip(S - 2 * C, 0., 1e300))

def distance_sq(X, Y):
    """
    Squared distance between C{X} and C{Y} along the last axis. For details, see L{distance}.

    @return: scalar or array of length m
    @rtype: (m,) numpy.array
    """
    return ((numpy.asarray(X) - Y) ** 2).sum(-1)

def distance(X, Y):
    """
    Distance between C{X} and C{Y} along the last axis.

    @param X: m x n input vector
    @type X: numpy array
    @param Y: m x n input vector
    @type Y: numpy array
    @return: scalar or array of length m
    @rtype: (m,) numpy.array
    """
    return distance_sq(X, Y) ** 0.5

def average_structure(X):
    """
    Calculate an average structure from an ensemble of structures
    (i.e. X is a rank-3 tensor: X[i] is a (N,3) configuration matrix).

    @param X: m x n x 3 input vector
    @type X: numpy array
    
    @return: average structure
    @rtype: (n,3) numpy.array
    """
    from numpy.linalg import eigh

    B = csb.numeric.gower_matrix(X)
    v, U = eigh(B)
    if numpy.iscomplex(v).any():
        v = v.real
    if numpy.iscomplex(U).any():
        U = U.real

    indices = numpy.argsort(v)[-3:]
    v = numpy.take(v, indices, 0)
    U = numpy.take(U, indices, 1)
        
    x = U * numpy.sqrt(v)
    i = 0
    while is_mirror_image(x, X[0]) and i < 2:
        x[:, i] *= -1
        i += 1
    return x


def is_mirror_image(X, Y):
    """
    Check if two configurations X and Y are mirror images
    (i.e. their optimal superposition involves a reflection).

    @param X: n x 3 input vector
    @type X: numpy array
    @param Y: n x 3 input vector
    @type Y: numpy array
 
    @rtype: bool
    """
    from numpy.linalg import det, svd
    
    ## center configurations

    X = X - numpy.mean(X, 0)
    Y = Y - numpy.mean(Y, 0)

    ## SVD of correlation matrix

    V, L, U = svd(numpy.dot(numpy.transpose(X), Y))             #@UnusedVariable

    R = numpy.dot(V, U)

    return det(R) < 0

# vi:expandtab:smarttab:sw=4
