from numpy import vstack, arange, array, reshape
from mulch.math import log_sum_exp, log

def log_trapezoidal(log_y, x = None):
    """
    computes the logarithm of the 1D integral of x, using trepezoidal approximation

    assumes x is monotonically increasing
    """
    if x is not None:
        log_x_diff = log(x[1:]-x[:-1])
    else:
        log_x_diff = 0.

    log_y_add = log_sum_exp(vstack((log_y[:-1],log_y[1:])),0)    

    return log_sum_exp(log_y_add + log_x_diff) - log(2)


def log_midpoint_rule_2d(log_f, x, y):
    x_delta = x[:-1]-x[1:]
    y_delta = y[:-1]-y[1:]

    z = array([log_f[:,1:] ,log_f[:,:-1]])

    y_hat = log_sum_exp(z.reshape((2,-1)),0 )
    y_hat = reshape(y_hat, (len(x),len(y)-1))
    y_hat += log(y_delta) - log(2)

    return log_sum_exp(y_hat + log(x_delta)) - log(2.)
    

def log_trapezoidal_2d(log_f,x = None, y = None):
    """
    computes the logarithm of the 1D integral of x, using trepezoidal approximation

    assumes x and y is monotonically increasing
    """
    int_y = array([log_trapezoidal(log_f[i,:],y) for i in range(len(y))])
    
    return log_trapezoidal(int_y, x)
    
    
def trapezoidal(x, y):
    from numpy import dot

    return 0.5 * dot(x[1:]-x[:-1],y[1:]+y[:-1])
    

def trapezoidal_2d(f):
    """
    Approximate the integral of f from a to b in two dimensions, using trepezoidal approximation

    @param f: 2D numpy array of function values at equally spaces positions
    @return: approximation of the definit integral
    """

    I = f[0, 0] + f[-1, -1] + f[0, -1] + f[-1, 0]
    I += 2 * (f[1:-1, (0, -1)].sum() + f[(0, -1), 1:-1].sum())
    I += 4 * f[1:-1, 1:-1].sum()

    return I / 4.

def simpson_2d(f):
    """
    Approximate the integral of f from a to b in two dimensions, using Composite Simpson's rule

    @param f: 2D numpy array of function values
    @return: approximation of the definit integral
    """

    n = (f.shape[0] - 1) / 2
    i = 2 * arange(1, n + 1) - 1
    j = 2 * arange(1, n)

    I = f[0, 0] + f[-1, -1] + f[0, -1] + f[-1, 0]
    I += 4 * (f[0, i].sum() + f[-1, i].sum() + f[0, j].sum() + f[-1, j].sum())
    I += 4 * (f[i, 0].sum() + f[i, -1].sum() + f[j, 0].sum() + f[j, -1].sum())
    I += 16 * f[i][:, i].sum() + 8 * (f[i][:, j].sum() + f[j][:, i].sum())
    I += 4 * f[j][:, j].sum()

    return I / 9.
