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
    from numpy import dot, transpose, average

    ## center configurations

    x = average(X, 0)
    y = average(Y, 0)

    ## SVD of correlation matrix

    V, _L, U = svd(dot(transpose(X - x), Y - y))

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
    """

    from numpy import sum, dot, transpose, sqrt, clip, average
    from numpy.linalg import svd

    X = X - average(X, 0)
    Y = Y - average(Y, 0)

    R_x = sum(X ** 2)
    R_y = sum(Y ** 2)

    L = svd(dot(transpose(Y), X))[1]

    return sqrt(clip(R_x + R_y - 2 * sum(L), 0., 1e300) / len(X))

def tm_d0(Lmin):
    
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
    @rtype: float
    """
    from numpy import sum

    d = sum((x - y) ** 2, 1) ** 0.5
    d0= tm_d0(len(x))
    
    return sum(1 / (1 + (d / d0) ** 2))

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
    
    @return: computed TM-Score between the vectors
    @rtype: float
    """
    from numpy import array, sum, dot, compress, power, ones, argmax

    x = array(x)
    y = array(y)

    mask = ones(len(x)).astype('i')

    scores = []
    transformations = []

    for i in range(3):

        R, t = fit_method(compress(mask, x, 0), compress(mask, y, 0))

        scores.append(tm_score(x, dot(y, R.T) + t))
        transformations.append((R,t))

        cutoff = min(max(4.5, tm_d0(len(x))), 8.0) + i
        d = sum((x - dot(y, R.T) - t) ** 2, 1) ** 0.5

        while True:
            
            mask = (d < cutoff).astype('i')
            if sum(mask) > 3: break
            cutoff += 0.5

    return max(scores), transformations[argmax(scores)]
            

