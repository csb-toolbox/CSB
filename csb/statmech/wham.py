"""
This module provides algorithms to estimate free energies and density of states
from tempered ensembles using histogram reweighting
"""

import numpy


from csb.math import log, log_sum_exp
from csb.statistics import histogram_nd

class AbstractWHAM(object):
    """
    Abstract base class
    """

    def __init__(self, ensembles, raw_energies, n):
        self._f = numpy.zeros(len(ensembles))
        self._e = raw_energies
        self._n = n
        self._L = []
        self._log_g = None
        self._ensembles = ensembles


    def log_g(self, normalize = True):
        """
        Returns the Density of states (DOS)

        @param normalize: Ensure that the density of states sums to one

        @rtype float
        """
        if normalize:
            return self._log_g - log_sum_exp(self._log_g)
        else:
            return self._log_g 
    
    @property
    def free_energies(self):
        """
        Free energies
        """
        return self._f

    
    
    def _stop_criterium(self, tol = 1e-10):
        """
        general stop criterium; if the relative difference between
        sequential negative log likelihoods is less than a predefined
        tolerance
        
        @param tol: tolerance
        @type tol: float

        @rtype: boolean
        """
        L = self._L
        return  tol is not None and len(L) > 1 and \
                   abs((L[-2] - L[-1]) /(L[-2] + L[-1])) < tol
        
    
    def estimate(self, *params):
        """
        Estimate the density of states
        """
        pass

    def log_z(self, beta = 1., ensembles = None):
        """
        Compute the partition function for an ensemble at inverse temperature
        beta or for a defined ensemble

        @param beta: Inverse Temperature
        @type beta: float or list

        @param beta: List of ensembles for which the partition function should be evaluated
        @type beta: List of ensembles

        @rtype: float or array
        """
        
        pass
 
class WHAM(AbstractWHAM):
    """
    Implementation of the original WHAM methods based on histograms
    """

    def __init__(self, ensembles, raw_energies, n):
        super(WHAM, self).__init__(ensembles, raw_energies, n)

        self._ex = None
        self._h = None
        
    def estimate(self, n_bins = 100,  n_iter = 10000, tol = 1e-10):

        self._L = []        
        h, e = histogram_nd(self._e, nbins = n_bins, normalize = False)
        self._ex = e = numpy.array(e)
        self._h = h
        f = self._f
        
        log_h = log(h)
        log_g = h * 0.0
        log_g -= log_sum_exp(log_g)
        log_n = log(self._n)

        e_ij = -numpy.squeeze(numpy.array([ensemble.energy(e)
                                           for ensemble in self._ensembles] )).T

        for _i in range(n_iter):

            ## update density of states
            y = log_sum_exp(numpy.reshape((e_ij - f + log_n).T,
                                          (len(f),-1)),0)
            log_g = log_h - numpy.reshape(y, log_g.shape)
            log_g -= log_sum_exp(log_g)

            ## update free energies
            f = log_sum_exp(numpy.reshape(e_ij.T + log_g.flatten(),
                                          (len(f), -1)).T, 0)
            self._L.append((self._n * f).sum() - (h * log_g).sum())

            self._log_g = log_g
            self._f = f

            if self._stop_criterium(tol):
                break

        return f, log_g


    def log_z(self, beta=1., ensembles = None):
        """
        use trapezoidal rule to evaluate the partition function
        TODO does so far only work for the one 1D
        """
        from numpy import array, multiply, reshape

        is_float = False

        if type(beta) == float:
            beta = reshape(array(beta),(-1,))
            is_float = True

        for i in range(self._ex.shape[0]):
            x = self._ex[i,1:] - self._ex[i,:-1]
            y = - multiply.outer(beta, self._ex[i]) + self._log_g
            y = reshape(array([y.T[1:], y.T[:-1]]), (2, len(beta)*len(x)))
            y = log_sum_exp(y, 0) - log(2)
            y = reshape(y, (len(x), len(beta))).T + log(x)

        log_z = log_sum_exp(y.T, 0)

        if is_float:
            return float(log_z)
        else:
            return log_z



    
class ImplicitWHAM(AbstractWHAM):
    """
    Implementation of the implicit WHAM outlined in Habeck 2012, in which histograms
    are reduced to delta peaks, this allows to use energies samples at different orders 
    of magnitude, improving the accuracy of the DOS estimates
    """

    def estimate(self, n_iter = 10000, tol = 1e-10):

        e_ij =  numpy.array([ensemble.energy(self._e)
                             for ensemble in self._ensembles]).T

        f = self._f
        log_n = log(self._n)
        self._L = []
        for _i in range(n_iter):

            ## update density of states
            log_g = - log_sum_exp((-e_ij - f + log_n).T, 0)
            log_g -= log_sum_exp(log_g)

            ## update free energies            
            f = log_sum_exp((- e_ij.T + log_g).T, 0)
            self._L.append((self._n * f).sum() - log_g.sum())

            self._f = f
            self._log_g = log_g

            if self._stop_criterium(tol):
                break

        return f, log_g

 
    def log_g(self, normalize=True):

        e_ij = - numpy.array([ensemble.energy(self._e)
                              for ensemble in self._ensembles]).T
        
        log_g = - e_ij - self._f + log(self._n)
        log_g = - log_sum_exp(log_g.T, 0)

        if normalize:
            log_g -= log_sum_exp(log_g)

        return log_g


    def log_z(self, beta=1., ensembles = None):

        from numpy import multiply

        e_ij = numpy.array([ensemble.energy(self._e)
                              for ensemble in self._ensembles]).T

        x = log_sum_exp((- e_ij - self._f + log(self._n)).T, 0)

        if ensembles is not None:
            e_ij_prime = numpy.array([ensemble.energy(self._e)
                                      for ensemble in ensembles])
        else:
            e_ij_prime =  multiply.outer(beta, self._e)
        
        
        log_z = log_sum_exp((- e_ij_prime - x).T, 0)

        return log_z
