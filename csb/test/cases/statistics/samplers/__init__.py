import numpy as np
import csb.test as test

from csb.statistics.pdf import Normal, AbstractDensity

from csb.numeric.integrators import AbstractGradient, VelocityVerlet, LeapFrog, FastLeapFrog

from csb.statistics.samplers import State
from csb.statistics.samplers.mc import MDRENSSwapParameterInfo, ThermostattedMDRENSSwapParameterInfo
from csb.statistics.samplers.mc import RESwapParameterInfo
from csb.statistics.samplers.mc import AlternatingAdjacentSwapScheme
from csb.statistics.samplers.mc.singlechain import HMCSampler, RWMCSampler
from csb.statistics.samplers.mc.multichain import ReplicaExchangeMC, MDRENS, ThermostattedMDRENS
from csb.statistics.samplers.mc.propagators import RWMCPropagator, HMCPropagator


class SamplePDF(Normal):
    
    def log_prob(self, x):
        return sum(map(super(SamplePDF, self).log_prob, x))

class MultimodalPDF(AbstractDensity):

    def log_prob(self, x):
        return sum(-2.5 * np.cos(2.5 * x) - 0.04 * x ** 2)

    def grad(self, x, t):
        return -6.25 * np.sin(2.5 * x) + 0.08 * x

    
@test.functional
class TestMCPropagators(test.Case):

    def setUp(self):
        
        super(TestMCPropagators, self).setUp()

        self.pdf = SamplePDF()
        self.gradient = self._createGradient(1.)
        self.timestep = 1.2
        self.stepsize = 1.2
        self.nsteps = 25
        self.nits = 10000
        self.state = State(np.random.normal(size=1))

    def _createGradient(self, sigma):
        
        class Grad(AbstractGradient):
            def evaluate(self, q, t):
                return q / (sigma ** 2)
        return Grad()       

    def check_result(self, trajectory):

        dim = len(trajectory[0].position)
        for i in range(dim):
            states = [state.position[i] for state in trajectory]
            self.assertAlmostEqual(np.array(states).mean(), 0., delta=1e-1)
            self.assertAlmostEqual(np.array(states).var(), 1., delta=1e-1)

    def testRWMCPropagator(self):

        gen = RWMCPropagator(self.pdf, self.stepsize)

        self.check_result(gen.generate(self.state, self.nits))

    def testHMCPropagator(self):

        gen = HMCPropagator(self.pdf, self.gradient, self.timestep, self.nsteps)

        self.check_result(gen.generate(self.state, self.nits))

    def testHMCPropagator_mm(self):

        mm = np.array([[1., 0.], [0., 2.]])
        init_state = State(np.random.normal(size=2))
        gen = HMCPropagator(self.pdf, self.gradient, self.timestep, self.nsteps)
        self.check_result(gen.generate(init_state, self.nits))
        
        
@test.functional
class TestMultichain(test.Case):

    def setUp(self):
        
        super(TestMultichain, self).setUp()

        self.temperatures = [0.4, 2.0]

        init_state = State(np.random.uniform(low=-3.0, high=3.0, size=1))
        
        self.samplers = [RWMCSampler(MultimodalPDF(), init_state, 0.5, 
                                     temperature=self.temperatures[0]),
                         RWMCSampler(MultimodalPDF(), init_state, 5.5, 
                                     temperature=self.temperatures[1])]
        
        self.nits = 10000
        self.algorithm = None
        self.Ts = [lambda l: l * self.temperatures[1] + (1. - l) * self.temperatures[0]]

        self.grad = self.samplers[0]._pdf.grad

    def set_1p_samplerparams(self):
        init_state = State(np.random.uniform(low=-3.0, high=3.0, size=1))
        self.samplers[0].stepsize = 0.5
        self.samplers[1].stepsize = 5.5
        self.samplers[0].state = init_state
        self.samplers[1].state = init_state

    def set_2p_samplerparams(self):
        init_state = State(np.random.uniform(low=-3.0, high=3.0, size=2))
        self.samplers[0].stepsize = 0.3
        self.samplers[1].stepsize = 2.5
        self.samplers[0].state = init_state
        self.samplers[1].state = init_state
        
    def _run(self, algorithm):

        xmin1 = -2.5
        xmax1 = 0.0
        xmin2 = 0.0
        xmax2 = 2.5
        p_occ = 0.382

        swapper = AlternatingAdjacentSwapScheme(algorithm)

        n_occ1 = 0
        n_occ2 = 0
        
        for i in range(self.nits):
            if i % 5 == 0:
                swapper.swap_all()

            else:
                algorithm.sample()

            x = self.samplers[0].state.position[0]
            
            if x > xmin1 and x < xmax1:
                n_occ1 += 1
            if x > xmin2 and x < xmax2:
                n_occ2 += 1

        p_occ_sampled1 = float(n_occ1) / float(self.nits)
        p_occ_sampled2 = float(n_occ2) / float(self.nits)

        # Assert by comparison with real occupation probabilities and a tolerance of
        # four standard deviations of a run with n=15000 samples and 100 iterations
        
        self.assertAlmostEqual(p_occ_sampled1, p_occ, delta=4.0 * 0.035)
        self.assertAlmostEqual(p_occ_sampled2, p_occ, delta=4.0 * 0.035)
        
    def testReplicaExchangeMC(self):

        params = [RESwapParameterInfo(self.samplers[0], self.samplers[1])]
        algorithm = ReplicaExchangeMC(self.samplers, params)
        self._run(algorithm)

    def testMDRENS(self):

        self.set_1p_samplerparams()
        params = [MDRENSSwapParameterInfo(self.samplers[0], self.samplers[1],
                                          0.025, 15, self.grad)]
        algorithm = MDRENS(self.samplers, params)
        self._run(algorithm)

    def testThermostattedMDRens(self):
        
        self.set_1p_samplerparams()
        params = [ThermostattedMDRENSSwapParameterInfo(self.samplers[0], self.samplers[1],
                                                       0.05, 15, self.grad,
                                                       temperature=self.Ts[0])]
        algorithm = ThermostattedMDRENS(self.samplers, params)
        self._run(algorithm)

    def testThermostattedMDRens_mm(self):
        mm = np.array([[1., 0.], [0., 3.]])
        self.set_2p_samplerparams()
        params = [ThermostattedMDRENSSwapParameterInfo(self.samplers[0], self.samplers[1],
                                                       0.025, 75, self.grad,
                                                       temperature=self.Ts[0],
                                                       mass_matrix=mm)]
        algorithm = ThermostattedMDRENS(self.samplers, params)
        self._run(algorithm)
                
if __name__ == '__main__':

    test.Console()
