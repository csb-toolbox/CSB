"""
Low level numeric / math utility functions.
"""

import sys
import math
import numpy


EXP_MIN = -308
EXP_MAX = +709
    
LOG_MIN = 1e-308
LOG_MAX = 1e+308

## Euler-Mascheroni constant
EULER_MASCHERONI = 0.57721566490153286060651209008240243104215933593992


def log(x, x_min=LOG_MIN, x_max=LOG_MAX):
    """
    Safe version of log, clips argument such that overflow does not occur.

    @param x: input
    @type x: numpy array or float or int

    @param x_min: lower value for clipping
    @type x_min: float

    @param x_max: upper value for clipping
    @type x_max: float
    """

    x_min = max(x_min, LOG_MIN)
    x_max = min(x_max, LOG_MAX)

    return numpy.log(numpy.clip(x, x_min, x_max))

def exp(x, x_min=EXP_MIN, x_max=EXP_MAX):
    """
    Safe version of exp, clips argument such that overflow does not occur.

    @param x: input
    @type x: numpy array or float or int

    @param x_min: lower value for clipping
    @type x_min: float

    @param x_max: upper value for clipping
    @type x_max: float
    """

    x_min = max(x_min, EXP_MIN)
    x_max = min(x_max, EXP_MAX)

    return numpy.exp(numpy.clip(x, x_min, x_max))

def sign(x):
    """
    Return the sign of the input.
    """
    return numpy.sign(x)

def isreal(x, tol=1e-14):
    """
    Check if input array has no imaginary part.

    @param x: input array
    @type x: numpy array

    @param tol: tolerance to check for equality zero
    @type tol: float
    """
    return not hasattr(x, 'real') or abs(x.imag) < tol

def log_sum_exp(x, axis=0):
    """
    Return the logarithm of the sum of exponentials.

    @type x: Numpy array
    """
    xmax = x.max(axis)
    return log(exp(x - xmax).sum(axis)) + xmax

def log_sum_exp_accumulate(x, axis=0):
    """
    Return the logarithm of the accumulated sums of exponentials.

    @type x: Numpy array
    """
    xmax = x.max(axis)
    return log(numpy.add.accumulate(exp(x - xmax), axis)) + xmax

def radian2degree(x):
    """
    Convert radians angles to torsion angles.

    @param x: radian angle
    @return: torsion angle of x
    """
    x = x % (2 * numpy.pi)
    numpy.putmask(x, x > numpy.pi, x - 2 * numpy.pi)
    return x * 180. / numpy.pi

def degree2radian(x):
    """
    Convert randian angles to torsion angles.

    @param x: torsion angle
    @return: radian angle of x
    """
    numpy.putmask(x, x < 0., x + 360.)
    return x * numpy.pi / 180.

def euler_angles(r):
    """
    Calculate the euler angles from a three dimensional rotation matrix.

    @param r: 3x3 Rotation matrix
    """
    a = numpy.arctan2(r[2, 1], r[2, 0]) % (2 * numpy.pi)
    b = numpy.arctan2((r[2, 0] + r[2, 1]) / (numpy.cos(a) + numpy.sin(a)), r[2, 2]) % (2 * numpy.pi)
    c = numpy.arctan2(r[1, 2], -r[0, 2]) % (2 * numpy.pi)

    return a, b, c

def euler(a, b, c):
    """
    Calculate a three dimensional rotation matrix from the euler angles.

    @param a: alpha, angle between the x-axis and the line of nodes
    @param b: beta, angle between the z axis of the different coordinate systems
    @param c: gamma, angle between the line of nodes and the X-axis
    """
    from numpy import cos, sin, array

    ca, cb, cc = cos(a), cos(b), cos(c)
    sa, sb, sc = sin(a), sin(b), sin(c)

    return array([[ cc * cb * ca - sc * sa, cc * cb * sa + sc * ca, -cc * sb],
                  [-sc * cb * ca - cc * sa, -sc * cb * sa + cc * ca, sc * sb],
                  [     sb * ca, sb * sa, cb ]])

def rotation_matrix(axis, angle):
    """
    Calculate a three dimensional rotation matrix for a rotation around
    the given angle and axis.

    @type axis: (3,) numpy array
    @param angle: angle in radians
    @type angle: float

    @rtype: (3,3) numpy.array
    """
    axis = numpy.asfarray(axis) / norm(axis)
    assert axis.shape == (3,)

    c = math.cos(angle)
    s = math.sin(angle)

    r = (1.0 - c) * numpy.outer(axis, axis)
    r.flat[[0,4,8]] += c
    r.flat[[5,6,1]] -= s * axis
    r.flat[[7,2,3]] += s * axis

    return r

def axis_and_angle(r):
    """
    Calculate axis and angle of rotation from a three dimensional
    rotation matrix.

    @param r: 3x3 Rotation matrix

    @return: axis unit vector as (3,) numpy.array and angle in radians as float
    @rtype: tuple
    """
    from scipy.linalg import logm

    B = logm(r).real
    a = numpy.array([-B[1,2], B[0,2], -B[0,1]])
    n = norm(a)

    return a / n, n

def norm(x):
    """
    Calculate the Eucledian norm of a d-dimensional vector.

    @param x: vector (i.e. rank one array)
    @return: length of vector
    """
    return numpy.linalg.norm(x)

def reverse(array, axis=0):
    """
    Reverse the order of elements in an array.
    """
    from numpy import take, arange
    return take(array, arange(array.shape[axis] - 1, -1, -1), axis)

def polar(x):
    """
    Polar coordinate representation of a d-dimensional vector.

    @param x: vector (i.e. rank one array)
    @return: polar coordinates (radius and polar angles)
    """
    
    (d,) = x.shape
    phi = numpy.zeros(d)
    for i in reversed(range(1, d)):
        phi[i - 1] = numpy.arctan2(x[i] / numpy.cos(phi[i]), x[i - 1])
        
    return numpy.array([norm(x)] + phi[:-1].tolist())

def from_polar(x):
    """
    Reconstruct d-dimensional vector from polar coordinates.

    @param x: vector (i.e. rank one array)
    @return: position in d-dimensional space
    """
    
    (d,) = x.shape

    c = numpy.cos(x[1:])
    s = numpy.sin(x[1:])
    r = x[0]

    x = numpy.zeros(d)
    x[0] = r
    
    for i in range(d - 1):

        x[i + 1] = x[i] * s[i]
        x[i] *= c[i]

    return x

def polar3d(x):
    """
    Polar coordinate representation of a three-dimensional vector.

    @param x: vector (i.e. rank one array)
    @return: polar coordinates (radius and polar angles)
    """

    if x.shape != (3,):
        raise ValueError(x)
    
    r = norm(x)
    theta = numpy.arccos(x[2] / r)
    phi = numpy.arctan2(x[1], x[0])

    return numpy.array([r, theta, phi])

def from_polar3d(x):
    """
    Reconstruct 3-dimensional vector from polar coordinates.

    @param x: vector (i.e. rank one array)
    @return: position in 3-dimensional space
    """
    assert x.shape == (3,)

    r, theta, phi = x[:]
    s = numpy.sin(theta)
    c = numpy.cos(theta)
    S = numpy.sin(phi)
    C = numpy.cos(phi)
    
    return r * numpy.array([s * C, s * S, c])

def dihedral_angle(a, b, c, d):
    """
    Calculate the dihedral angle between 4 vectors,
    representing 4 connected points. The angle is in range [-180, 180].

    @param a: the four points that define the dihedral angle
    @type a: array

    @return: angle in [-180, 180]
    """
    
    v = b - c
    m = numpy.cross((a - b), v)
    m /= norm(m)
    n = numpy.cross((d - c), v)
    n /= norm(m)

    c = numpy.dot(m, n)
    s = numpy.dot(numpy.cross(n, m), v) / norm(v)
    
    angle = math.degrees(math.atan2(s, c))        

    if angle > 0:
        return numpy.fmod(angle + 180, 360) - 180
    else:
        return numpy.fmod(angle - 180, 360) + 180

def psi(x):
    """
    Digamma function
    """
    from numpy import inf, log, sum, exp

    coef = [-1. / 12., 1. / 120., -1. / 252., 1. / 240., -1. / 132.,
            691. / 32760., -1. / 12.]

    if x == 0.:
        return -inf
    elif x < 0.:
        raise ValueError('not defined for negative values')
    elif x < 6.:
        return psi(x + 1) - 1. / x
    else:
        logx = log(x)
        res = logx - 0.5 / x
        res += sum([coef[i] * exp(-2 * (i + 1) * logx) for i in range(7)])
        return res

def approx_psi(x):
    from numpy import log, clip, where

    if type(x) == numpy.ndarray:
        y = 0. * x
        y[where(x < 0.6)] = -EULER_MASCHERONI - 1. / clip(x[where(x < 0.6)], 1e-154, 1e308)
        y[where(x >= 0.6)] = log(x[where(x >= 0.6)] - 0.5)
        return y
    else:
        if x < 0.6:
            return -EULER_MASCHERONI - 1 / clip(x, 1e-154, 1e308)
        else:
            return log(x - 0.5)

def d_approx_psi(x):
    from numpy import clip, where
    
    if type(x) == numpy.ndarray:
        y = 0. * x
        y[where(x < 0.6)] = 1. / clip(x[where(x < 0.6)], 1e-154, 1e308) ** 2
        y[where(x >= 0.6)] = 1. / (x[where(x >= 0.6)] - 0.5)
        return y
    else:

        if x < 0.6:
            return 1 / clip(x, 1e-154, 1e308) ** 2
        else:
            return 1 / (x - 0.5)
    
def inv_psi(y, tol=1e-10, n_iter=100, psi=psi):
    """
    Inverse digamma function
    """
    from numpy import exp
    from scipy.special import digamma
    ## initial value

    if y < -5 / 3. - EULER_MASCHERONI:
        x = -1 / (EULER_MASCHERONI + y)
    else:
        x = 0.5 + exp(y)

    ## Newton root finding

    for dummy in range(n_iter):

        z = digamma(x) - y

        if abs(z) < tol:
            break
        
        x -= z / d_approx_psi(x)

    return x

def log_trapezoidal(log_y, x=None):
    """
    Compute the logarithm of the 1D integral of x, using trepezoidal approximation.
    Assumes x is monotonically increasing.
    """
    if x is not None:
        log_x_diff = log(x[1:] - x[:-1])
    else:
        log_x_diff = 0.

    log_y_add = log_sum_exp(numpy.vstack((log_y[:-1], log_y[1:])), 0)    

    return log_sum_exp(log_y_add + log_x_diff) - log(2)

def log_midpoint_rule_2d(log_f, x, y):
    x_delta = x[:-1] - x[1:]
    y_delta = y[:-1] - y[1:]

    z = numpy.array([log_f[:, 1:] , log_f[:, :-1]])

    y_hat = log_sum_exp(z.reshape((2, -1)), 0)
    y_hat = numpy.reshape(y_hat, (len(x), len(y) - 1))
    y_hat += log(y_delta) - log(2)

    return log_sum_exp(y_hat + log(x_delta)) - log(2.)
    
def log_trapezoidal_2d(log_f, x=None, y=None):
    """
    Compute the logarithm of the 1D integral of x, using trepezoidal approximation.
    Assumes x and y is monotonically increasing.
    """
    int_y = numpy.array([log_trapezoidal(log_f[i, :], y) for i in range(len(y))])
    
    return log_trapezoidal(int_y, x)  
    
def trapezoidal(x, y):
    return 0.5 * numpy.dot(x[1:] - x[:-1], y[1:] + y[:-1])

def trapezoidal_2d(f):
    """
    Approximate the integral of f from a to b in two dimensions,
    using trepezoidal approximation.

    @param f: 2D numpy array of function values at equally spaces positions
    @return: approximation of the definit integral
    """

    I = f[0, 0] + f[-1, -1] + f[0, -1] + f[-1, 0]
    I += 2 * (f[1:-1, (0, -1)].sum() + f[(0, -1), 1:-1].sum())
    I += 4 * f[1:-1, 1:-1].sum()

    return I / 4.

def simpson_2d(f):
    """
    Approximate the integral of f from a to b in two dimensions,
    using Composite Simpson's rule.

    @param f: 2D numpy array of function values
    @return: approximation of the definit integral
    """

    n = int((f.shape[0] - 1) / 2)
    i = 2 * numpy.arange(1, n + 1) - 1
    j = 2 * numpy.arange(1, n)

    I = f[0, 0] + f[-1, -1] + f[0, -1] + f[-1, 0]
    I += 4 * (f[0, i].sum() + f[-1, i].sum() + f[0, j].sum() + f[-1, j].sum())
    I += 4 * (f[i, 0].sum() + f[i, -1].sum() + f[j, 0].sum() + f[j, -1].sum())
    I += 16 * f[i][:, i].sum() + 8 * (f[i][:, j].sum() + f[j][:, i].sum())
    I += 4 * f[j][:, j].sum()

    return I / 9.

def pad(x, s):
    """
    Add layers of zeros around grid.
    """
    s = numpy.array(s) - 1
    y = numpy.zeros(numpy.array(x.shape) + s)
    s /= 2
    slices = [slice(s[i], -s[i]) for i in range(len(s))]
    y[slices] = x

    return y

def trim(x, s):
    """
    Remove additional layers.
    """
    s = numpy.array(s) - 1
    s /= 2

    slices = [slice(s[i], -s[i]) for i in range(len(s))]

    return x[slices]

def zerofill(x, s):
    
    y = numpy.zeros(s)
    slices = [slice(-x.shape[i], None) for i in range(len(s))]
    y[slices] = x

    return y
              
def convolve(x, f):

    from numpy import fft, all

    sx = numpy.array(x.shape)
    sf = numpy.array(f.shape)

    if not all(sx >= sf): return convolve(f, x)
        
    y = fft.ifftn(fft.fftn(x) * fft.fftn(f, sx)).real
    slices = [slice(sf[i] - 1, sx[i]) for i in range(len(sf))]

    return y[slices]

def correlate(x, y):

    from numpy import fft

    sx = numpy.array(x.shape)
    sy = numpy.array(y.shape)

    if (sx >= sy).sum():

        slices = [slice(None, sx[i] - sy[i] + 1) for i in range(len(sx))]
        
        X = fft.fftn(x)
        Y = fft.fftn(zerofill(y, sx))

    else:

        sf = sx + sy - 1
        slices = [slice(None, sf[i]) for i in range(len(sf))]
        
        X = fft.fftn(x, sf)
        Y = fft.fftn(zerofill(y, sf), sf)

    return fft.ifftn(X.conjugate() * Y)[slices].real


def gower_matrix(X):
    """
    Gower, J.C. (1966). Some distance properties of latent root
    and vector methods used in multivariate analysis.
    Biometrika 53: 325-338

    @param X: ensemble coordinates
    @type X: (m,n,k) numpy.array

    @return: symmetric dissimilarity matrix
    @rtype: (n,n) numpy.array
    """
    X = numpy.asarray(X)

    B = sum(numpy.dot(x, x.T) for x in X) / float(len(X))
    b = B.mean(1)
    bb = b.mean()

    return (B - numpy.add.outer(b, b)) + bb

class MatrixInitError(Exception):
    pass

class InvertibleMatrix(object):
    """
    Matrix object which is intended to save time in MCMC sampling algorithms involving
    repeated integration of Hamiltonian equations of motion and frequent draws from
    multivariate normal distributions involving mass matrices as covariance matrices.
    It can be initialized either with the matrix one wants to use or its inverse.
    The main feature compared to standard numpy matrices / arrays is that it has also
    a property "inverse", which gives the inverse matrix. If the matrix (its inverse) is
    changed, the inverse (regular matrix) is calculated only when needed. This avoids costly
    matrix inversions.

    @param matrix: matrix-like object with whose values the Matrix object is supposed to
                   hold
    @type matrix:  invertible (n,n)-shaped numpy.ndarray

    @param inverse_matrix: matrix-like object with whose inverse the Matrix object is supposed
                           to hold
    @type inverse_matrix:  invertible (n,n)-shaped numpy.ndarray
    """
    
    def __init__(self, matrix=None, inverse_matrix=None):

        if (matrix is None and inverse_matrix is None):
            raise MatrixInitError("At least one matrix argument has to be specified")

        self._matrix = None
        self._inverse_matrix = None

        self._matrix_updated = False
        self._inverse_matrix_updated = False
        
        self._is_unity_multiple = False

        if matrix is not None and inverse_matrix is not None:
            if type(matrix) != numpy.ndarray or type(inverse_matrix) != numpy.ndarray:
                raise TypeError("Arguments have to be of type numpy.ndarray!")
            matrix = matrix.copy()
            inverse_matrix = inverse_matrix.copy()
            self._check_equal_shape(matrix, inverse_matrix)
            self._matrix = matrix
            self._inverse_matrix = inverse_matrix
            self._is_unity_multiple = self._check_unity_multiple(self._matrix)
            self._matrix_updated = True
            self._inverse_matrix_updated = True
        else:
            if matrix is not None:
                if type(matrix) != numpy.ndarray:
                    raise TypeError("Arguments have to be of type numpy.ndarray!")
                matrix = matrix.copy()
                self._check_square(matrix)
                self._matrix = matrix
                self._matrix_updated = True
                self._inverse_matrix_updated = False
                self._is_unity_multiple = self._check_unity_multiple(self._matrix)
            else:
                if type(inverse_matrix) != numpy.ndarray:
                    raise TypeError("Arguments have to be of type numpy.ndarray!")
                inverse_matrix = inverse_matrix.copy()
                self._check_square(inverse_matrix)
                self._inverse_matrix = inverse_matrix
                self._matrix_updated = False
                self._inverse_matrix_updated = True
                self._is_unity_multiple = self._check_unity_multiple(self._inverse_matrix)
            
    def _check_diagonal(self, matrix):

        return (matrix.T == matrix).all()

    def _check_unity_multiple(self, matrix):

        diagonal_elements_equal = numpy.all([matrix[i][i] == matrix[i+1][i+1] 
                                             for i in range(len(matrix) - 1)])

        return self._check_diagonal(matrix) and diagonal_elements_equal

    def _check_square(self, matrix):

        if matrix.shape[0] != matrix.shape[1]:
            raise ValueError("Matrix " + matrix.__name__ + " must be a square matrix!")

    def _check_equal_shape(self, matrix1, matrix2):

        if not numpy.all(matrix1.shape == matrix2.shape):
            raise ValueError("Matrices " + matrix1.__name__ + " and " + matrix2.__name__ + 
                             " must have equal shape!")
                    
    def __getitem__(self, loc):
        if self._matrix_updated == False:
            if self.is_unity_multiple:
                self._matrix = numpy.diag(1. / self._inverse_matrix.diagonal())
            else:
                self._matrix = numpy.linalg.inv(self._inverse_matrix)
            self._matrix_updated = True
            
        return self._matrix[loc]

    def __setitem__(self, i, value):
        if type(value) != numpy.ndarray:
            raise TypeError("Arguments have to be of type numpy.ndarray!")
        self._matrix[i] = numpy.array(value)
        self._matrix_updated = True
        self._inverse_matrix_updated = False
        self._is_unity_multiple = self._check_unity_multiple(self._matrix)

    def __array__(self, dtype=None):
        return self._matrix

    def __mul__(self, value):
        if type(value) in (float, int):
            v = float(value)
            if self._matrix_updated:
                return InvertibleMatrix(v * self._matrix)
            else:
                return InvertibleMatrix(inverse_matrix = self._inverse_matrix / v)

        else:
            raise ValueError("Only float and int are supported for multiplication!")
    
    __rmul__ = __mul__

    def __truediv__(self, value):
        if type(value) in (float, int):
            v = float(value)
            if self._matrix_updated:
                return InvertibleMatrix(self._matrix / v)
            else:
                return InvertibleMatrix(inverse_matrix=self._inverse_matrix * v)

        else:
            raise ValueError("Only float and int are supported for division!")
        
    __div__ = __truediv__

    def __imul__(self, value):
        if type(value) in (float, int):
            if self._matrix_updated:
                self._matrix *= float(value)
                self._deprecate_inverse_matrix()
            else:
                self._inverse_matrix /= float(value)
                self._inverse_matrix_updated = True
                self._deprecate_matrix()

            return self
        else:
            raise ValueError("Only float and int are supported for multiplication!")

    def __itruediv__(self, value):
        if type(value) in (float, int):
            if self._matrix_updated:
                self._matrix /= float(value)
                self._matrix_updated = True
                self._deprecate_inverse_matrix()
            else:
                self._inverse_matrix *= float(value)
                self._inverse_matrix_updated = True
                self._deprecate_matrix()

            return self
        else:
            raise ValueError("Only float and int are supported for division!")

    __idiv__ = __itruediv__

    def __str__(self):
        if self._matrix is not None and self._inverse_matrix is not None:
            return "csb.numeric.InvertibleMatrix object holding the following numpy matrix:\n"\
              +self._matrix.__str__() +"\n and its inverse:\n"+self._inverse_matrix.__str__()
        else:
            if self._matrix is None:
                return "csb.numeric.InvertibleMatrix object holding only the inverse matrix:\n"\
                  +self._inverse_matrix.__str__()
            else:
                return "csb.numeric.InvertibleMatrix object holding the matrix:\n"\
                  +self._matrix.__str__()

    def __len__(self):
        if self._matrix is not None:
            return len(self._matrix)
        else:
            return len(self._inverse_matrix)
    
    def _deprecate_matrix(self):
        self._matrix_updated = False

    def _deprecate_inverse_matrix(self):
        self._inverse_matrix_updated = False

    def _update_matrix(self, matrix=None):
        """
        Updates the _matrix field given a new matrix or by setting
        it to the inverse of the _inverse_matrix field.

        @param matrix: matrix-like object which the Matrix object
                       is supposed to represent
        @type matrix:  (n,n)-shaped numpy.ndarray or list
        """
        
        if matrix is None:
            self._matrix = numpy.linalg.inv(self._inverse_matrix)
        else:
            self._matrix = matrix
            
        self._matrix_updated = True
        self._deprecate_inverse_matrix()
            

    def _update_inverse_matrix(self, inverse=None):
        """
        Updates the __inverse_matrix field given a new matrix or by setting
        it to the inverse of the _matrix field.

        @param inverse: matrix-like object which the Matrix object
                        is supposed to represent
        @type inverse:  (n,n)-shaped numpy.ndarray or list
        """

        if inverse is None:
            if self.is_unity_multiple:
                self._inverse_matrix = numpy.diag(1. / self._matrix.diagonal())
            else:
                self._inverse_matrix = numpy.linalg.inv(self._matrix)
        else:
            self._inverse_matrix = inverse
            
        self._inverse_matrix_updated = True
        
    @property
    def inverse(self):
        if self._inverse_matrix_updated == False:
            self._update_inverse_matrix()
            
        return self._inverse_matrix.copy()
    @inverse.setter
    def inverse(self, value):
        if type(value) != numpy.ndarray:
            raise TypeError("Arguments have to be of type numpy.ndarray!")
        self._check_equal_shape(value, self._matrix)
        self._update_inverse_matrix(numpy.array(value))
        self._deprecate_matrix()
        self._is_unity_multiple = self._check_unity_multiple(self._inverse_matrix)

    @property
    def is_unity_multiple(self):
        """
        This property can be used to save computation time when drawing
        from multivariate normal distributions with the covariance matrix
        given by an instance of this object. By probing this property,
        one can instead draw from normal distributions and rescale the samples
        afterwards to avoid costly matrix inversions
        """
        return self._is_unity_multiple
