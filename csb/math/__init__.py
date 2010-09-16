EXP_MIN = -1e308
EXP_MAX = +709
LOG_MIN = 1e-300
LOG_MAX = 1e+300

from numpy import array, log, clip, exp, pi, putmask, cos, sin, pi, arccos, arctan2, cross, linalg



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

    x_min = max(x_min, EXP_MIN)
    x_max = min(x_max, EXP_MAX)

    return exp(clip(x, x_min, x_max))

def sign(x):
    """
    Return the sign of the input.
    """
    return 2 * int(x >= 0.) - 1

def isreal(x, tol=1e-14):
    """
    Checks if input array has no imaginary part.

    @param x: input array
    @type x: numpy array

    @param tol: tolerance to check for equality zero
    @type tol: float
    """
    return not hasattr(x, 'real') or abs(x.imag) < tol

def log_sum_exp(x, axis=0):
    """
    Returns the logarithm of the sum of exponentials.

    @type x: Numpy array
    """
    xmax = x.max(axis)

    return log(exp(x - xmax).sum(axis)) + xmax

def radian2degree(x):
    """
    Converts radians angles to torsion angles

    @param x: radian angle
    @return: torsion angle of x
    """

    x = x % (2 * pi)
    putmask(x, x > pi, x - 2 * pi)
    return x * 180. / pi

def degree2radian(x):
    """
    Converts randian angles to torsion angles

    @param x: torsion angle
    @return: radian angle of x
    """

    putmask(x, x < 0., x + 360.)
    return x * pi / 180.


def euler_angles(r):
    """
    Calculate the euler angles from a three dimensional rotation matrix

    @param r: 3x3 Rotation matrix, 
    """

    a = arctan2(r[2,1],r[2,0]) % (2*pi)
    b = arctan2((r[2,0]+r[2,1]) / (cos(a) + sin(a)), r[2,2]) % (2*pi)
    c = arctan2(r[1,2],-r[0,2]) % (2*pi)

    return a, b, c


def euler(a, b, c):
    """
    Calculate  a three dimensional rotation matrix from the euler angles 

    @param a: alpha, angle between the x-axis and the line of nodes
    @param b: beta, angle between the z axis of the different coordinate systems
    @param c: gamma,  angle between the line of nodes and the X-axis
    """
    from numpy import cos, sin, array

    ca, cb, cc = cos(a), cos(b), cos(c)
    sa, sb, sc = sin(a), sin(b), sin(c)

    return array([[ cc*cb*ca-sc*sa,  cc*cb*sa+sc*ca, -cc*sb],
                  [-sc*cb*ca-cc*sa, -sc*cb*sa+cc*ca,  sc*sb],
                  [     sb*ca,            sb*sa,        cb ]])


def norm(x):
    """
    Calculates the Eucledian norm of a d-dimensional vector.

    @param x: vector (i.e. rank one array)
    @return: length of vector
    """
    return (x**2).sum()**0.5


def polar(x):
    """
    Polar coordinate representation of a d-dimensional vector.

    @param x: vector (i.e. rank one array)
    @return: polar coordinates (radius and polar angles)
    """
    from numpy import zeros, arctan2, pi, cos, array

    (d,) = x.shape
    phi = zeros(d)
    for i in range(1,d)[::-1]:
        phi[i-1] = arctan2(x[i]/cos(phi[i]),x[i-1])
        
    return array([norm(x)] + phi[:-1].tolist())

def from_polar(x):
    """
    Reconstruct d-dimensional vector from polar coordinates

    @param x: vector (i.e. rank one array)
    @return: position in d-dimensional space
    """
    from numpy import cos, sin, zeros
    
    (d,) = x.shape

    c = cos(x[1:])
    s = sin(x[1:])
    r = x[0]

    x = zeros(d)
    x[0] = r
    
    for i in range(d-1):

        x[i+1] = x[i] * s[i]
        x[i] *= c[i]

    return x

def polar3d(x):
    """
    Polar coordinate representation of a three-dimensional vector.

    @param x: vector (i.e. rank one array)
    @return: polar coordinates (radius and polar angles)
    """

    assert x.shape == (3,)

    
    r = norm(x)
    theta = arccos(x[2]/r)
    phi = arctan2(x[1], x[0])

    return array([r, theta, phi])

def from_polar3d(x):
    """
    Reconstruct 3-dimensional vector from polar coordinates

    @param x: vector (i.e. rank one array)
    @return: position in 3-dimensional space
    """
    assert x.shape == (3,)


    r, theta, phi = x[:]
    s = sin(theta)
    c = cos(theta)
    S = sin(phi)
    C = cos(phi)
    
    return r * array([s*C, s*S, c])

def dihedral_angle(a,b,c,d):
    """
    Calculate the dihedral angle between 4 vectors
    representing 4 connected points. The angle is in  ]-180, 180].

    @param a, b, c, d: the four points that define the dihedral angle
    @type a, b, c, d: Numpy Array

    @return: angle in  ]-180, 180]
    """
    from numpy import fmod
    
    v = b - c
    m = cross((a - b),v)
    m /=  linalg.norm(m)
    n = cross((d - c),v)
    n /=  linalg.norm(m)

    c = dot(m,n)
    s = dot(cross(n,m),v) / linalg.norm(v)
    
    angle = math.degrees(math.atan2(s, c))        

    if angle > 0:
        return fmod(angles + 180, 360) - 180)
    else:
        return fmod(angles - 180, 360) + 180)

    
    
    
if __name__ == '__main__':

    from numpy.random import random
    import numpy as np
    
    d = 3

    x = random(d)

    print polar(x)
    print polar3d(x)

    tol = 1e-10

    while 1:

        x = random(d)
        
        if np.fabs(x-from_polar3d(polar3d(x))).sum() > tol: break
        if np.fabs(x-from_polar(polar(x))).sum() > tol: break
    
