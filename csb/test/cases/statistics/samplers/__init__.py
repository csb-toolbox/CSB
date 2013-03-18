import numpy as np
import csb.test as test

from csb.statistics.pdf import Normal, AbstractDensity

from csb.numeric.integrators import AbstractGradient, VelocityVerlet, LeapFrog, FastLeapFrog
from csb.numeric import InvertibleMatrix

from csb.statistics.samplers import State
from csb.statistics.samplers.mc.multichain import MDRENSSwapParameterInfo, MDRENS
from csb.statistics.samplers.mc.multichain import ThermostattedMDRENSSwapParameterInfo
from csb.statistics.samplers.mc.multichain import RESwapParameterInfo, AlternatingAdjacentSwapScheme
from csb.statistics.samplers.mc.multichain import ReplicaExchangeMC, ThermostattedMDRENS
from csb.statistics.samplers.mc.multichain import HMCStepRENS, HMCStepRENSSwapParameterInfo
from csb.statistics.samplers.mc.singlechain import HMCSampler, RWMCSampler, AbstractNCMCSampler
from csb.statistics.samplers.mc.propagators import RWMCPropagator, HMCPropagator
from csb.statistics.samplers.mc.propagators import AbstractNCMCPropagator
from csb.statistics.samplers.mc.neqsteppropagator import ReducedHamiltonian, HamiltonianSysInfo
from csb.statistics.samplers.mc.neqsteppropagator import MDPropagation, MDPropagationParam
from csb.statistics.samplers.mc.neqsteppropagator import Protocol, Step
from csb.statistics.samplers.mc.neqsteppropagator import ReducedHamiltonianPerturbation

class SamplePDF(Normal):
    
    def log_prob(self, x):
        return sum(map(super(SamplePDF, self).log_prob, x))

    def grad(self, x, t):
        return x / (self.sigma ** 2)

class MultimodalPDF(AbstractDensity):

    def log_prob(self, x):
        return sum(-2.5 * np.cos(2.5 * x) - 0.04 * x ** 2)

    def grad(self, x, t):
        return -6.25 * np.sin(2.5 * x) + 0.08 * x

class Multimodal2DPDF(AbstractDensity):

    k = 0.5

    def _E1(self, x):
        return 2.5 * np.cos(2.5 * x[0]) + 0.04 * x[0] ** 2

    def _E2(self, x):
        return self.k * x[1] ** 2

    def log_prob(self, x):
        return -self._E1(x) - self._E2(x)

    def grad(self, x, t):
        return np.array([(-6.25 * np.sin(2.5 * x[0]) + 0.08 * x[0]) * self._E2(x),
                         self._E1(x) * self.k * x[1]])

    
@test.functional
class TestMCPropagators(test.Case):

    def setUp(self):
        
        super(TestMCPropagators, self).setUp()

        self.pdf = SamplePDF()
        self.gradient = self._createGradient(1.)
        self.timestep = 1.2
        self.stepsize = 1.2
        self.nsteps = 15
        self.nits = 10000
        self.state = State(np.random.normal(size=1))

    def _createGradient(self, sigma):
        
        class Grad(AbstractGradient):
            def evaluate(self, q, t):
                return q / (sigma ** 2)
        return Grad()       

    def checkResult(self, trajectory):

        dim = len(trajectory[0].position)
        for i in range(dim):
            states = [state.position[i] for state in trajectory]
            self.assertAlmostEqual(np.array(states).mean(), 0., delta=0.15)
            self.assertAlmostEqual(np.array(states).var(), 1., delta=0.15)

    def testRWMCPropagator(self):

        gen = RWMCPropagator(self.pdf, self.stepsize)

        self.checkResult(gen.generate(self.state, self.nits))

    def testHMCPropagator(self):

        gen = HMCPropagator(self.pdf, self.gradient, self.timestep, self.nsteps)

        self.checkResult(gen.generate(self.state, self.nits))

    def testHMCPropagatorMM(self):

        mm = InvertibleMatrix(np.array([[1., 0.], [0.,  2.]]))
        init_state = State(np.random.normal(size=2))
        gen = HMCPropagator(self.pdf, self.gradient, self.timestep * 1.5, self.nsteps, mass_matrix=mm)

        self.checkResult(gen.generate(init_state, self.nits))

    def testNCMCPropagator(self):

        Nhalf = 30
        dt = 0.3
        md_tl = 1
        
        ks = np.linspace(1.0, 0.2, Nhalf).tolist()
        sigmas = [1/np.sqrt(k) for k in ks]
        sigmas += sigmas[::-1][1:]
        N = len(sigmas)
        pdfs = [SamplePDF(sigma=s) for s in sigmas]
        hamiltonians = [ReducedHamiltonian(pdfs[i].log_prob, pdfs[i].grad) for i in range(N)]
        sys_infos = [HamiltonianSysInfo(hamiltonians[i]) for i in range(N)]
        
        steps = [Step(ReducedHamiltonianPerturbation(sys_infos[i], sys_infos[i+1],
                                                     evaluate_work=False),
                      MDPropagation(sys_infos[i+1], MDPropagationParam(dt, md_tl, pdfs[i+1].grad),
                                    evaluate_heat=False))
                      for i in range(N - 1)]
        rv_steps = [Step(ReducedHamiltonianPerturbation(sys_infos[i], sys_infos[i+1],
                                                        evaluate_work=False),
                         MDPropagation(sys_infos[i], MDPropagationParam(dt, md_tl, pdfs[i].grad),
                                       evaluate_heat=False))
                   for i in range(N - 1)]
        
        for s in rv_steps:
            s.set_propagation_first()
        protocol = Protocol(steps)
        rv_protocol = Protocol(rv_steps)

        class MDProbStepNCMCSampler(AbstractNCMCSampler):
            def _calc_pacc(self, proposal_communicator):
                return np.exp(-proposal_communicator.traj.deltaH)

        class MDPropStepNCMCPropagator(AbstractNCMCPropagator):
            def _init_sampler(self, init_state):
                self._sampler = MDProbStepNCMCSampler(init_state, self.protocol,
                                                      self.reverse_protocol)

        gen = MDPropStepNCMCPropagator(protocol, rv_protocol)

        init_state = State(np.array([1.0]), np.array([-0.5]))
        
        traj = gen.generate(init_state, self.nits, return_trajectory=True)

        self.checkResult(traj)

        self.assertAlmostEqual(np.mean([x.momentum[0] for x in traj]), 0.0, delta=0.15)
        self.assertAlmostEqual(np.std([x.momentum[0] for x in traj]), 1.0, delta=0.15)




        
@test.functional
class TestMultichain(test.Case):

    def setUp(self):
        super(TestMultichain, self).setUp()

        self.samplers = None

    def set1pParams(self):
        init_state = State(np.random.uniform(low=-3.0, high=3.0, size=1))
        self.temperatures = [0.4, 2.0]
        self.samplers = [RWMCSampler(MultimodalPDF(), init_state, 0.5,
                                     temperature=self.temperatures[0]),
                         RWMCSampler(MultimodalPDF(), init_state, 5.5,
                                     temperature=self.temperatures[1])]
        self.grad = self.samplers[0]._pdf.grad
        self.nits = 10000
        self.Ts = [lambda l: l * self.temperatures[i+1] + (1. - l) * self.temperatures[i]
                   for i in range(len(self.samplers) - 1)]
        
    def set2pParams(self):
        init_state = State(np.random.uniform(low=-3.0, high=3.0, size=2))
        pdf = Multimodal2DPDF()
        self.temperatures = [0.4, 1.0, 2.0]
        self.samplers = [RWMCSampler(pdf, init_state, 0.2,
                                     temperature=self.temperatures[0]),
                         RWMCSampler(pdf, init_state, .8,
                                     temperature=self.temperatures[1]),
                         RWMCSampler(pdf, init_state, 2.,
                                     temperature=self.temperatures[2])]
        self.grad = self.samplers[0]._pdf.grad
        self.nits = 20000
        self.Ts = [lambda l: l * self.temperatures[i+1] + (1. - l) * self.temperatures[i]
                   for i in range(len(self.samplers) - 1)]
        
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
        self.set1pParams()
        params = [RESwapParameterInfo(self.samplers[0], self.samplers[1])]
        algorithm = ReplicaExchangeMC(self.samplers, params)
        self._run(algorithm)

    def testMDRENS(self):

        self.set1pParams()
        params = [MDRENSSwapParameterInfo(self.samplers[0], self.samplers[1],
                                          0.025, 15, self.grad)]
        algorithm = MDRENS(self.samplers, params, integrator=VelocityVerlet)
        self._run(algorithm)

    def testThermostattedMDRens(self):
        
        self.set1pParams()
        params = [ThermostattedMDRENSSwapParameterInfo(self.samplers[0], self.samplers[1],
                                                       0.05, 15, self.grad,
                                                       temperature=self.Ts[0])]
        algorithm = ThermostattedMDRENS(self.samplers, params)
        self._run(algorithm)

    def testThermostattedMDRensMM(self):

        self.set2pParams()
        mm1 = InvertibleMatrix(np.array([[1.0, 0.0], [0.0, 5.0]]))
        mm2 = InvertibleMatrix(np.array([[.5, 0.0], [0.0, 10.0]]))
        pdf = Multimodal2DPDF()
        params = [ThermostattedMDRENSSwapParameterInfo(self.samplers[0], self.samplers[1],
                                                       0.01, 15, pdf.grad,
                                                       temperature=self.Ts[0],
                                                       mass_matrix=mm1),
                  ThermostattedMDRENSSwapParameterInfo(self.samplers[1], self.samplers[2],
                                                       0.1, 15, pdf.grad,
                                                       temperature=self.Ts[1],
                                                       mass_matrix=mm2)]
        algorithm = ThermostattedMDRENS(self.samplers, params)

        self._run(algorithm)

    def testHMCStepRENS(self):

        self.set1pParams()
        params = [HMCStepRENSSwapParameterInfo(self.samplers[0], self.samplers[1], 0.05, 3, 1,
                                               self.grad, 5)]

        algorithm = HMCStepRENS(self.samplers, params)

        self._run(algorithm)


if __name__ == '__main__':

    test.Console()
