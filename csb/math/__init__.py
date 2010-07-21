EXP_MIN = -1e308
EXP_MAX = +709
LOG_MIN = 1e-300
LOG_MAX = 1e+300

def log(x, x_min=LOG_MIN, x_max=LOG_MAX):
    """
    Save version of  log, clips argument such that overflow does not occur.

    @param x: input
    @type x: numpy array or float or int

    @param x_min: lower value for clipping
    @type x_min: float

    @param x_max: upper value for clipping
    @type x_max: float
    """    
    from numpy import log, clip

    x_min = max(x_min, LOG_MIN)
    x_max = min(x_max, LOG_MAX)

    return log(clip(x, x_min, x_max))

def exp(x, x_min=EXP_MIN, x_max=EXP_MAX):
    """
    Save version of exp, clips argument such that overflow does not occur.

    @param x: input
    @type x: numpy array or float or int

    @param x_min: lower value for clipping
    @type x_min: float

    @param x_max: upper value for clipping
    @type x_max: float
    """
    from numpy import exp, clip

    x_min = max(x_min, EXP_MIN)
    x_max = min(x_max, EXP_MAX)

    return exp(clip(x, x_min, x_max))

def sign(x):
    """
    Return the sign of the input.
    """
    return 2 * int(x>=0.) - 1

def isreal(x, tol=1e-14):
    """
    Checks if input array has no imaginary part.

    @param x: input array
    @type x: numpy array

    @param tol: tolerance to check for equality zero
    @type tol: float
    """
    return not hasattr(x,'real') or abs(x.imag) < tol

def log_sum_exp(x, axis=0):
    """
    Returns the logarithm of the sum of exponentials.

    @type x: Numpy array
    """
    xmax = x.max(axis)

    return log(exp(x-xmax).sum(axis)) + xmax

def radian2degree(x):
    """
    Converts radians angles to torsion angles

    @param x: radian angle
    @return: torsion angle of x
    """
    from numpy import pi, putmask

    x = x % (2*pi)
    putmask(x, x>pi, x-2*pi)
    return x * 180. / pi

def degree2radian(x):
    """
    Converts randian angles to torsion angles

    @param x: torsion angle
    @return: radian angle of x
    """
    from numpy import pi, putmask

    putmask(x, x<0., x+360.)
    return x * pi / 180.

