def trapezoidal(x, y):

    from numpy import dot

    return 0.5 * dot(x[1:]-x[:-1],y[1:]+y[:-1])
    
def trapezoidal_2d(f):
    """
    Approximate the integral of f from a to b in two dimensions, using trepezoidal approximation

    @param f: 2D numpy array of function values
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
    from numpy import arange

    n = (f.shape[0] - 1) / 2
    i = 2 * arange(1, n + 1) - 1
    j = 2 * arange(1, n)

    I = f[0, 0] + f[-1, -1] + f[0, -1] + f[-1, 0]
    I += 4 * (f[0, i].sum() + f[-1, i].sum() + f[0, j].sum() + f[-1, j].sum())
    I += 4 * (f[i, 0].sum() + f[i, -1].sum() + f[j, 0].sum() + f[j, -1].sum())
    I += 16 * f[i][:, i].sum() + 8 * (f[i][:, j].sum() + f[j][:, i].sum())
    I += 4 * f[j][:, j].sum()

    return I / 9.
