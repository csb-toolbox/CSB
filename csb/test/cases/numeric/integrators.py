import numpy as np
import csb.test as test

from math import cos

from csb.numeric.integrators import LeapFrog, VelocityVerlet, AbstractGradient
from csb.statistics.samplers import State


@test.functional
class TestIntegrators(test.Case):

    def setUp(self):
        
        super(TestIntegrators, self).setUp()
        
        self.dt = 0.1
        self.grad = self._createGradient(1.)
        self.nsteps = 100
        self.state = State(np.array([1.]), np.array([0.]))
        
    def _createGradient(self, sigma):
        
        class Grad(AbstractGradient):
            def evaluate(self, q, t):
                return q / (sigma ** 2)
            
        return Grad()        

    def _run(self, algorithm):
        
        result = algorithm.integrate(self.state, self.nsteps).last.position
        self.assertWithinDelta(result, cos(self.nsteps * self.dt), 0.1)
        
    def testLeapFrog(self):
        
        algorithm = LeapFrog(self.dt, self.grad)
        self._run(algorithm)

    def testVelocityVerlet(self):
        
        algorithm = VelocityVerlet(self.dt, self.grad)
        self._run(algorithm)


if __name__ == '__main__':

    test.Console()
