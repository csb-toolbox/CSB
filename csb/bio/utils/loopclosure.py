
def loop_closure(moving, fixed, min_rmsd=0.1, maxit=10000):
    """
    Implementation of the Full cyclic coordinate descent algorithm

    by Boomsma and Hamelryck (2005)

    This algorithm trys to close a loop given by superimposing the last three
    coordinates of the moving set onto the fixed set of size three.
    The first coordinat of the moving set does stay constant

    @param moving: Set of the coordinates, which should form the loop.The first coordinate wont be moved and therefore be the anchor point. The last three coordinates will be matched onto the three fixed  points
    @type moving: array_like

    @param fixed: Anchor points for the end of the loop
    @type fixed array like

    @param min_rmsd: goal rmsd between fixed and last three moving coordinates, if this is reached the algorithm terminates
    @type min_rmsd: float

    @param maxit: Maximal number of iterations to perform
    @type maxit: it
    
    """

    from numpy import zeros, dot, transpose, average, array, sqrt
    from numpy.linalg import svd, det

    if len(moving) < 6:
        raise ValueError("Moving too short")

    if len(fixed) != 3:
        raise ValueError("Fixed should have length 3")

    n = len(moving)
    # Copy stuff
    moving = array(moving)
    fixed = array(fixed)
    result = zeros(moving.shape)
    result[:, :] = moving[:, :]
    converged = False
    it = 0
    while not converged:
        for i in xrange(1, n - 3):
            # new pivot position around which the moving residues are rotated
            pivot = result[i, :]

            # Solve for rotation by svd
            V, _L, U = svd(dot(transpose(fixed - pivot),
                              result[-3:, :] - pivot))

            # Calc rotation 
            R = dot(V, U)
            # Check for rotoinversion
            if det(R) < 0.:
                U[-1] *= -1
                R = dot(V, U)

            # Apply rotation
            result[(i + 1):] = dot(result[(i + 1):] - pivot, transpose(R))\
                               + pivot
            # Calculate rmsd of fixed and moved residues
            rmsd = sqrt(average((result[-3:] - fixed[-3:]) ** 2))

        it += 1
        # Check for convergence 
        if rmsd < min_rmsd or it > maxit:
            converged = True

    return result, rmsd, it


def make_random_chain(n=12):
    """
    Return a list of random vectors, each with distance
    3.8 from the previous vector.

    @param n: length of chain
    @type n: int
    """

    from numpy import zeros
    from numpy.random import random
    from numpy.linalg import norm

    v = zeros(3)
    l = [v]
    for i in range(0, n - 1):               #@UnusedVariable
        nv = random(3)
        nv /= norm(nv)
        nv = v + nv ** 3.8
        l.append(nv)
        v = nv
    return l


def rotate_last_three(chain):
    """
    Take the last three vectors of chain, copy them,
    and apply a random rotation

    @param chain: coordinates of chain
    @type chain:  array_like
    """

    from numpy.random import random
    from numpy import pi, arccos, dot, transpose
    from csb.math import euler

    alpha = random() * 2 * pi
    gamma = random() * 2 * pi
    beta = arccos(2 * (random() - 0.5))

    R = euler(alpha, beta, gamma)
    return dot(chain[-3:], transpose(R))

if __name__ == "__main__":

    for i in range(100):
        # Moving segment
        moving = make_random_chain()
        # Fixed segment 
        # Last three residues of the moving segment
        # after applying a random rotation/translation
        fixed = rotate_last_three(moving)

        # Close the gap
        results, rmsd, it = loop_closure(moving, fixed, 0.1, 100)

        # Print result
        print "Final RMSD ", rmsd
        print "Number of iterations ", it
