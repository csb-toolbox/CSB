import numpy as np
import csb.test as test

from csb.statistics.pdf import Normal

from csb.numeric.integrators import AbstractGradient, VelocityVerlet, LeapFrog, FastLeapFrog

from csb.statistics.samplers import State
from csb.statistics.samplers.mc import MDRENSSwapParameterInfo, ThermostattedMDRENSSwapParameterInfo
from csb.statistics.samplers.mc import AlternatingAdjacentSwapScheme
from csb.statistics.samplers.mc.singlechain import HMCSampler, RWMCSampler
from csb.statistics.samplers.mc.multichain import ReplicaExchangeMC, MDRENS, ThermostattedMDRENS


class SamplePDF(Normal):
    
    def log_prob(self, x):
        return sum(map(super(SamplePDF, self).log_prob, x))
    
    
@test.functional
class TestSinglechain(test.Case):

    def setUp(self):
        
        super(TestSinglechain, self).setUp()

        self.pdf = SamplePDF()
        self.grad = self._createGradient(1.)
        self.dt = 1.2
        self.nsteps = 25
        self.nits = 15000
        self.state = State(np.random.normal(size=1))
        
    def _createGradient(self, sigma):
        
        class Grad(AbstractGradient):
            def evaluate(self, q, t):
                return q / (sigma ** 2)
            
        return Grad()        

    def _run(self, sampler):
        states = []

        for i in range(self.nits):
            sampler.sample()
            states.append(sampler.state.position[0])

        self.assertWithinDelta(np.array(states).mean(), 0., 1e-1)
        self.assertWithinDelta(np.array(states).var(), 1., 1e-1)

    def testHMCSampler(self):
        
        sampler = HMCSampler(self.pdf, self.state, self.grad,
                                  self.dt, self.nsteps, VelocityVerlet)
        self._run(sampler)

    def testRWMCSampler(self):

        sampler = RWMCSampler(self.pdf, self.state, 1.2)
        self._run(sampler)
        
        
@test.functional
class TestMultichain(test.Case):

    def setUp(self):
        
        super(TestMultichain, self).setUp()

        self.sigma1, self.sigma2, self.sigma3 = 1., 1. / np.sqrt(3), 1. / np.sqrt(5)
        init_state = State(np.random.normal(size=1))

        self.samplers = [RWMCSampler(SamplePDF(sigma=self.sigma1), init_state, 2),
                         RWMCSampler(SamplePDF(sigma=self.sigma1), init_state, 1.5),
                         RWMCSampler(SamplePDF(sigma=self.sigma1), init_state, 1.0)]
        
        self.timesteps = [0.05, 0.05]
        self.nsteps = [100, 100]
        self.nits = 5000
        self.algorithm = None
        self.temperatures = [lambda l: 1., lambda l: 1.]
        
        self.gradients = [self._createGradientL(self.sigma1, self.sigma2),
                          self._createGradientL(self.sigma2, self.sigma3)]
        
        self.params = [ThermostattedMDRENSSwapParameterInfo(self.samplers[0], self.samplers[1],
                                                            self.timesteps[0], self.nsteps[0],
                                                            self.gradients[0]),
                       ThermostattedMDRENSSwapParameterInfo(self.samplers[1], self.samplers[2],
                                                            self.timesteps[1], self.nsteps[1],
                                                            self.gradients[1])]

    def _createGradientL(self, sigma1, sigma2):
        
        class Grad(AbstractGradient):
            def evaluate(self, q, l):
                return q * (l / (sigma2 ** 2) + (1. - l) / (sigma1 ** 2))

        return Grad()

    def _run(self, algorithm):
        
        samples = []
        swapper = AlternatingAdjacentSwapScheme(algorithm)

        for i in range(self.nits):
            if i % 5 == 0 and i > 0:
                swapper.swap_all()
                samples.append(algorithm.state)

            samples.append(algorithm.sample())

        states = [x[0].position[0] for x in samples]

        self.assertWithinDelta(np.array(states).mean(), 0., 1.5e-1)
        self.assertWithinDelta(np.array(states).var(), 1., 1.5e-1)

    def testReplicaExchangeMC(self):
        
        algorithm = ReplicaExchangeMC(self.samplers, self.params)
        self._run(algorithm)

    def testMDRENS(self):
        
        algorithm = MDRENS(self.samplers, self.params, FastLeapFrog)
        self._run(algorithm)

    def testThermostattedMDRens(self):
        
        algorithm = ThermostattedMDRENS(self.samplers, self.params, VelocityVerlet)
        self._run(algorithm)
        
        
if __name__ == '__main__':

    test.Console()
