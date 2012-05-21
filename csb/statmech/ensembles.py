"""
Statistical Ensembles
"""

from csb.numeric import log, exp
from abc import ABCMeta, abstractmethod


class StatisticalEnsemble(object):
    
    __metaclass__ = ABCMeta
  
    def __call__(self, raw_energies):
        return exp(-self.energy(raw_energies))

    def log_prob(self, raw_energies):
        return -self.energy(raw_energies)

    @abstractmethod
    def energy(self, raw_energies):
        """
        Transforms the raw energies as if they were observed
        in this statistical ensemble
        """
        pass

    def gradient(self, raw_energies):
        raise NotImplementedError()


class BoltzmannEnsemble(StatisticalEnsemble):

    def __init__(self, beta=1.):
        
        self._beta = float(beta)

    @property
    def beta(self):
        """
        Inverse temperature
        """
        return self._beta
    @beta.setter
    def beta(self, value):
        value = float(value)
        if value <= 0.:
            raise ValueError("Inverse temperature {0} < 0".formate(value))
        self._beta = value

    def energy(self, raw_energies):
        return raw_energies * self._beta
        

class TsallisEnsemble(StatisticalEnsemble):

    def __init__(self, q=1., e_min=0.):

        self._q = q
        self._e_min = e_min
    
    @property
    def q(self):
        """
        q-analoge of the temperature
        """
        return self._q
    @q.setter
    def q(self, value):
        if value <= 0.:
            raise ValueError("Inverse temperature {0} < 0".formate(value))
        self._q = value

    @property
    def e_min(self):
        """
        lower bound of the energy
        """
        return self._e_min
    @e_min.setter
    def e_min(self, value):
        self._e_min = value

    def energy(self, raw_energies):
        q = self.q
        e_min = self.e_min
        
        if (q < 1 + 1e-10):
            return raw_energies * q
        else:
            return log(1 + (raw_energies - e_min) * (q - 1)) * q / (q - 1) + e_min


class CompositeEnsemble(StatisticalEnsemble):

    def __init__(self, ensembles=[]):

        self._ensembles = ensembles

    @property
    def ensembles(self):
        """
        Collection of statistical ensembles
        """
        return self._ensembles
    @ensembles.setter
    def ensembles(self, value):
        if not isinstance(value, list):
            if len(value) > 0:
                if not isinstance(value[0], StatisticalEnsemble):
                    raise  ValueError("Not a list of statistical ensembles")
                else:
                    self._enesmbles = value
            else:
                self._enesmbles = value

    def energy(self, raw_energies):
        return sum([self._ensembles[i].energy(raw_energies[i])
                    for i in range(len(self.ensembles))], 0)
    
    def gradient(self, raw_energies):
        return sum([self._ensembles[i].gradient(raw_energies[i])
                    for i in range(len(self.ensembles))], 0)
