import numpy as np

import csb.test as test
import csb.numeric

from csb.statistics.pdf import Normal, BaseDensity

from csb.numeric.integrators import AbstractGradient, VelocityVerlet, LeapFrog, FastLeapFrog
from csb.numeric import InvertibleMatrix

from csb.statistics.samplers import State
from csb.statistics.samplers.mc import Trajectory
from csb.statistics.samplers.mc.multichain import MDRENSSwapParameterInfo, MDRENS
from csb.statistics.samplers.mc.multichain import ThermostattedMDRENSSwapParameterInfo
from csb.statistics.samplers.mc.multichain import RESwapParameterInfo, AlternatingAdjacentSwapScheme
from csb.statistics.samplers.mc.multichain import ReplicaExchangeMC, ThermostattedMDRENS
from csb.statistics.samplers.mc.multichain import HMCStepRENS, HMCStepRENSSwapParameterInfo
from csb.statistics.samplers.mc.multichain import AbstractSwapCommunicator, AbstractExchangeMC
from csb.statistics.samplers.mc.multichain import AbstractSwapParameterInfo, ReplicaHistory
from csb.statistics.samplers.mc.singlechain import HMCSampler, RWMCSampler, AbstractNCMCSampler
from csb.statistics.samplers.mc.singlechain import AbstractSingleChainMC
from csb.statistics.samplers.mc.propagators import RWMCPropagator, HMCPropagator, MDPropagator
from csb.statistics.samplers.mc.propagators import AbstractNCMCPropagator, AbstractPropagator
from csb.statistics.samplers.mc.neqsteppropagator import ReducedHamiltonian, HamiltonianSysInfo
from csb.statistics.samplers.mc.neqsteppropagator import PlainMDPropagation, PlainMDPropagationParam
from csb.statistics.samplers.mc.neqsteppropagator import AbstractMDPropagation, HMCPropagation
from csb.statistics.samplers.mc.neqsteppropagator import Protocol, Step, AbstractPerturbation
from csb.statistics.samplers.mc.neqsteppropagator import ReducedHamiltonianPerturbation, AbstractPropagation
from csb.statistics.samplers.mc.neqsteppropagator import NonequilibriumStepPropagator
from csb.statistics.samplers.mc.neqsteppropagator import NonequilibriumTrajectory
from csb.statistics.samplers.mc.neqsteppropagator import HMCPropagationParam

class SamplePDF(Normal):
    
    def log_prob(self, x):
        return sum(map(super(SamplePDF, self).log_prob, x))

    def grad(self, x, t):
        return x / (self.sigma ** 2)

class MultimodalPDF(BaseDensity):

    def log_prob(self, x):
        return sum(-2.5 * np.cos(2.5 * x) - 0.04 * x ** 2)

    def grad(self, x, t):
        return -6.25 * np.sin(2.5 * x) + 0.08 * x

class Multimodal2DPDF(BaseDensity):

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

    @test.skip("Takes quite a long time to run.")
    def testNCMCPropagator(self):

        Nhalf = 5
        dt = 0.1
        md_tl = 5
        
        ks = np.linspace(1.0, 0.2, Nhalf).tolist()
        sigmas = [1/np.sqrt(k) for k in ks]
        sigmas += sigmas[::-1][1:]
        N = len(sigmas)
        pdfs = [SamplePDF(sigma=s) for s in sigmas]
        hamiltonians = [ReducedHamiltonian(pdfs[i].log_prob, pdfs[i].grad) for i in range(N)]
        sys_infos = [HamiltonianSysInfo(hamiltonians[i]) for i in range(N)]
        
        steps = [Step(ReducedHamiltonianPerturbation(sys_infos[i], sys_infos[i+1],
                                                     evaluate_work=False),
                      PlainMDPropagation(sys_infos[i+1], 
                                         PlainMDPropagationParam(dt, md_tl, pdfs[i+1].grad),
                                         evaluate_heat=False))
                 for i in range(N - 1)]
        rv_steps = [Step(ReducedHamiltonianPerturbation(sys_infos[i], sys_infos[i+1],
                                                        evaluate_work=False),
                         PlainMDPropagation(sys_infos[i],
                                            PlainMDPropagationParam(dt, md_tl, pdfs[i].grad),
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

        init_state = State(np.array([1.0]))
        traj = gen.generate(init_state, self.nits, return_trajectory=True)
        self.checkResult(traj)

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
    
    @test.skip("Takes some time, rendered optional by a unit test.")
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

class MockSwapCommunicator(AbstractSwapCommunicator):
    
    pass

class MockSwapParameterInfo(AbstractSwapParameterInfo):

    pass

class MockSampler(AbstractSingleChainMC):

    def __init__(self, pdf, state, temperature=1.0):
        
        self._state = state
        self._pdf = pdf
        self._temperature = temperature

    def _propose(self):
        
        pass

    def _calc_pacc(self):

        pass

class MockedAbstractExchangeMC(AbstractExchangeMC):
    
    def _propose_swap(self, param_info):

        return MockSwapCommunicator(param_info, Trajectory([State(np.array([1.0])),
                                                           State(np.array([2.0]))]),
                                                Trajectory([State(np.array([2.0])),
                                                           State(np.array([1.0]))]))

    def _calc_pacc_swap(self, swapcom):
        
        swapcom.acceptance_probability = 0.75

        return swapcom

@test.unit
class TestAbstractExchangeMC(test.Case):

    def setUp(self):

        self.samplers = [MockSampler(None, State(np.array([3.0]))),
                         MockSampler(None, State(np.array([5.0])))]

        self.param_info = MockSwapParameterInfo(self.samplers[0], self.samplers[1])

        self.algo = MockedAbstractExchangeMC(self.samplers, [self.param_info])

    
    def testAcceptSwap(self):

        swapcom = MockSwapCommunicator(self.param_info,
                                       Trajectory([State(np.array([1.0])),
                                                   State(np.array([2.0]))]),
                                       Trajectory([State(np.array([2.0])),
                                                   State(np.array([1.0]))]))

        np.random.seed(5)

        swapcom.acceptance_probability = 0.75
        res = self.algo._accept_swap(swapcom)
        assert(res)

        swapcom.acceptance_probability = 0.15
        res = self.algo._accept_swap(swapcom)
        assert(not res)

    def testSwap(self):

        np.random.seed(5)

        res = self.algo.swap(0)

        assert(res)
        self.assertEqual(self.samplers[0].state.position[0], 1.0)
        self.assertEqual(self.samplers[1].state.position[0], 2.0)
        self.assertEqual(self.algo.statistics.stats[0].total_swaps, 1)
        self.assertEqual(self.algo.statistics.stats[0].accepted_swaps, 1)

        np.random.seed(4)

        res = self.algo.swap(0)

        assert(not res)
        self.assertEqual(self.samplers[0].state.position[0], 1.0)
        self.assertEqual(self.samplers[1].state.position[0], 2.0)
        self.assertEqual(self.algo.statistics.stats[0].total_swaps, 2)
        self.assertEqual(self.algo.statistics.stats[0].accepted_swaps, 1)

@test.unit
class TestReplicaExchangeMC(test.Case):

    def setUp(self):

        pdf1 = HO()
        pdf2 = HO(k1=2.0, k2=2.0)

        self.samplers = [MockSampler(pdf1, State(np.array([3.0]))),
                         MockSampler(pdf2, State(np.array([5.0])))]

        self.param_info = RESwapParameterInfo(self.samplers[0], self.samplers[1])

        self.algo = ReplicaExchangeMC(self.samplers, [self.param_info])

    def testProposeSwap(self):
        
        res = self.algo._propose_swap(self.param_info)
        self.assertEqual(res.traj12.initial.position[0], 3.0)
        self.assertEqual(res.traj12.final.position[0], 3.0)
        self.assertEqual(res.traj21.initial.position[0], 5.0)
        self.assertEqual(res.traj21.final.position[0], 5.0)

    def testCalcPaccSwap(self):

        swapcom = self.algo._propose_swap(self.param_info)
        res = self.algo._calc_pacc_swap(swapcom)

        self.assertEqual(res.acceptance_probability, csb.numeric.exp(-12.5 + 4.5 - 9.0 + 25.0))

class HO(object):

    def __init__(self, k1=1.0, k2=1.0, x1=0.0, x2=0.0, tau=1.0):
        self.k1 = k1
        self.k2 = k2
        self.x1 = x1
        self.x2 = x2
        self.tau = tau
        self.kt = lambda t: self.k2 * t / self.tau + (1 - t / self.tau) * self.k1
        self.xt = lambda t: self.x2 * t / self.tau + (1 - t / self.tau) * self.x1

    def log_prob(self, x, t=0.0):
        return -0.5 * self.kt(t) * sum((x - self.xt(t)) ** 2)

    def gradient(self, x, t):
        return self.kt(t) * (x - self.xt(t))

class MockPropagator(AbstractPropagator):

    def __init__(self):
        pass

    def generate(self, init_state, length, return_trajectory=False):

        final_state = State(init_state.position * 2, init_state.momentum * 2)
        
        return Trajectory([init_state, final_state])
    
class PlainMDPropagationMocked(PlainMDPropagation):

    def _propagator_factory(self):

        return MockPropagator()

class HMCPropagationMocked(HMCPropagation):

    def _propagator_factory(self):

        return MockPropagator()

class MockPerturbation(AbstractPerturbation):

    @property
    def sys_before(self):

        pdf = HO()

        return HamiltonianSysInfo(ReducedHamiltonian(pdf.log_prob, pdf.gradient))

    @property
    def sys_after(self):

        pdf = HO()

        return HamiltonianSysInfo(ReducedHamiltonian(pdf.log_prob, pdf.gradient))

    def __init__(self):
        pass
    
    def _run_perturbator(self, state):

        final = State(state.position * 2, state.momentum * 2)
        
        return Trajectory([state, final])

    def _calculate_work(self, traj):

        return 42.0

    def _calculate_jacobian(self, traj):

        return 1.1

    
class MockPropagation(AbstractPropagation):

    def __init__(self):
        pass

    @property
    def sys(self):

        pdf = HO()

        return HamiltonianSysInfo(pdf.log_prob, pdf.gradient)

    def _run_propagator(self, state):

        final = State(state.position * 2, state.momentum * 2)

        return Trajectory([state, final])

    def _calculate_heat(self, traj):

        return -42.0

    def _propagator_factory(self):

        return None


class MockStep(Step):

    def __init__(self, return_momentum=True):

        self._return_momentum = return_momentum
        self._perform = None
        self.perform = self._perform_pert_prop

    @property
    def perturbation(self):

        return MockPerturbation()

    def _perform_pert_prop(self, state, extra_info=None):

        if self._return_momentum == True:
            final = State(state.position * 2, state.momentum * 2)
        else:
            final = State(state.position * 2)

        res = NonequilibriumTrajectory([state, final], heat=-42.0, work=42.0, jacobian=1.1)

        return res, None, None

    def _perform_prop_pert(self, state, extra_info=None):

        if self._return_momentum == True:
            final = State(state.position * 2, state.momentum * 2)
        else:
            final = State(state.position * 2)

        res = NonequilibriumTrajectory([state, final], heat=42.0, work=-42.0, jacobian=1.1)
        
        return res, None, None


class MockProtocol(Protocol):

    def __init__(self, momentum=True):

        self._momentum = momentum
        self.steps = [MockStep(self._momentum), MockStep(self._momentum)]


@test.unit
class TestNeqsteppropagator(test.Case):

    def testReducedHamiltonian(self):
        pdf = HO(k1=2.0, k2=2.0)
        init = State(np.array([2.0]), np.array([-2.0]))
        ham = ReducedHamiltonian(lambda x: pdf.log_prob(x, 0.0), pdf.gradient, temperature=4.0)

        self.assertEqual(4.0, ham.E(init.position))
        self.assertEqual(2.0, ham.kinetic_energy(init.momentum))
        self.assertEqual(0.0, ham.kinetic_energy(None))
        self.assertEqual(-1.0, ham.rlog_prob(init.position))
        self.assertEqual(0.5, ham.rkinetic_energy(init.momentum))
        self.assertEqual(1.5, ham(init))

    def testHMCPropagation(self):

        pdf = HO()
        sys = HamiltonianSysInfo(ReducedHamiltonian(pdf.log_prob, pdf.gradient))
        param = HMCPropagationParam(None, None, None)
        hmcprop = HMCPropagationMocked(sys, param)

        init = State(np.array([2.0]), np.array([2.0]))

        ## Test _set_mass_matrix
        d = len(init.position)
        param = HMCPropagationParam(None, None, None, mass_matrix=InvertibleMatrix(np.eye(d)))
        hmcprop = HMCPropagationMocked(sys, param)
        hmcprop._set_mass_matrix(init)
        self.assertEqual(hmcprop.param.mass_matrix,
                         InvertibleMatrix(np.eye(len(init.position))))

        
        param = HMCPropagationParam(None, None, None)
        hmcprop = HMCPropagationMocked(sys, param)
        hmcprop._set_mass_matrix(init)
        self.assertEqual(hmcprop.param.mass_matrix,
                         InvertibleMatrix(np.eye(len(init.position))))
        
        ## Test _calculate_heat
        final = State(init.position * 2, init.momentum * 2)
        traj = Trajectory([init, final])
        self.assertEqual(hmcprop._calculate_heat(traj), 6.0)

        ## Test __call__
        result = hmcprop(init)
        self.assertEqual(init.position, result.initial.position)
        self.assertEqual(init.momentum, result.initial.momentum)
        self.assertEqual(result.final.position, init.position * 2)
        self.assertEqual(result.final.momentum, init.momentum * 2)
        self.assertEqual(result.heat, 6.0)

    def testPlainMDPropagation(self):

        pdf = HO()
        sys = HamiltonianSysInfo(ReducedHamiltonian(pdf.log_prob, pdf.gradient))

        init = State(np.array([2.0]), np.array([2.0]))

        ## Test _set_mass_matrix
        d = len(init.position)
        param = PlainMDPropagationParam(None, None, None, 
                                        mass_matrix=InvertibleMatrix(np.eye(d)))
        mdprop = PlainMDPropagationMocked(sys, param)
        mdprop._set_mass_matrix(init)
        self.assertEqual(mdprop.param.mass_matrix,
                         InvertibleMatrix(np.eye(d)))

        param = PlainMDPropagationParam(None, None, None)
        mdprop = PlainMDPropagationMocked(sys, param)
        mdprop._set_mass_matrix(init)
        self.assertEqual(mdprop.param.mass_matrix,
                         InvertibleMatrix(np.eye(d)))

        ## Test _calculate_heat
        final = State(init.position * 2, init.momentum * 2)
        traj = Trajectory([init, final])
        self.assertEqual(mdprop._calculate_heat(traj), 12.0)

        ## Test __call__
        result = mdprop(init)
        self.assertEqual(init.position, result.initial.position)
        self.assertEqual(init.momentum, result.initial.momentum)
        self.assertEqual(result.final.position, init.position * 2)
        self.assertEqual(result.final.momentum, init.momentum * 2)
        self.assertEqual(result.heat, 12.0)
        
        
    def testReducedHamiltonianPerturbation(self):

        pdf = HO(k1=1.0, k2=2.0)
        redham1 = ReducedHamiltonian(lambda x: pdf.log_prob(x, 0.0))
        redham2 = ReducedHamiltonian(lambda x: pdf.log_prob(x, 1.0))
        sys1 = HamiltonianSysInfo(redham1)
        sys2 = HamiltonianSysInfo(redham2)
        init = State(np.array([2.0]), np.array([2.0]))
        traj = Trajectory([init, init])
        
        hampert = ReducedHamiltonianPerturbation(sys1, sys2)

        ## Test _calculate_work
        self.assertEqual(hampert._calculate_work(traj), 2.0)

        ## Test __call__
        result = hampert(init)
        self.assertEqual(result.initial.position[0], init.position[0])
        self.assertEqual(result.initial.momentum[0], init.momentum[0])
        self.assertEqual(result.initial.position[0], result.final.position[0])
        self.assertEqual(result.initial.momentum[0], result.final.momentum[0])
        self.assertEqual(result.work, 2.0)
        self.assertEqual(result.jacobian, 1.0)
        

    def testStep(self):

        step = Step(MockPerturbation(), MockPropagation())
        init = State(np.array([2.0]), np.array([2.0]))

        ## Test step with first perturbation, then propagation
        res = step.perform(init)[0]

        self.assertEqual(res.final.position, init.position * 4)
        self.assertEqual(res.final.momentum, init.momentum * 4)
        self.assertEqual(res.heat, -42.0)
        self.assertEqual(res.work, 42.0)
        self.assertEqual(res.jacobian, 1.1)

        ## Test step with first perturbation, then propagation
        step.set_propagation_first()
        res = step.perform(init)[0]

        self.assertEqual(step.perform, step._perform_prop_pert)
        self.assertEqual(res.final.position, init.position * 4)
        self.assertEqual(res.final.momentum, init.momentum * 4)
        self.assertEqual(res.heat, -42.0)
        self.assertEqual(res.work, 42.0)
        self.assertEqual(res.jacobian, 1.1)
        
    def testNonequilibriumStepPropagator(self):

        protocol = Protocol([MockStep(True) for i in range(10)])

        gen = NonequilibriumStepPropagator(protocol)

        ## Test generate()
        init = State(np.array([2.0]), np.array([2.0]))

        res = gen.generate(init)

        self.assertEqual(res.final.position, init.position * (2 ** 10))
        self.assertEqual(res.final.momentum, init.momentum * (2 ** 10))    
        self.assertEqual(res.work, 10 * 42)
        self.assertEqual(res.heat, -10 * 42)
        self.assertEqual(res.jacobian, 1.1 ** 10)


class MockedNCMCSampler(AbstractNCMCSampler):

    def _calc_pacc(self, proposal_communicator):
        
        return proposal_communicator.traj.final.position[0]


@test.unit
class TestNCMCSampler(test.Case):

    def testProposeWithMomentum(self):

        self.protocol = MockProtocol(True)
        self.reverse_protocol = MockProtocol(True)
        for s in self.reverse_protocol.steps:
            s.set_propagation_first()
        
        ## Test with momentum
        init = State(np.array([2.0]), np.array([2.0]))

        sampler = MockedNCMCSampler(init, self.protocol, self.reverse_protocol)

        ## Test _propose
        # Make sure the first random number is < 0.5
        np.random.seed(5)
        result = sampler._propose()

        self.assertEqual(result.traj.heat, - 2 * 42)
        self.assertEqual(result.traj.work, 2 * 42)
        self.assertEqual(result.traj.initial.position[0], init.position[0])
        self.assertEqual(result.traj.final.position[0],  init.position[0] * 4)
        self.assertEqual(result.traj.initial.momentum[0], init.momentum[0])
        self.assertEqual(result.traj.final.momentum[0],  init.momentum[0] * 4)

        # Make sure the first random number is > 0.5
        np.random.seed(4)
        result = sampler._propose()

        self.assertEqual(result.traj.heat, 2 * 42)
        self.assertEqual(result.traj.work, - 2 * 42)
        self.assertEqual(result.traj.initial.position[0], init.position[0])
        self.assertEqual(result.traj.final.position[0],  init.position[0] * 4)
        self.assertEqual(result.traj.initial.momentum[0], init.momentum[0])
        self.assertEqual(result.traj.final.momentum[0],  init.momentum[0] * 4)


    def testProposeWithoutMomentum(self):

        self.protocol = MockProtocol(False)
        self.reverse_protocol = MockProtocol(False)
        for s in self.reverse_protocol.steps:
            s.set_propagation_first()
        
        ## Test without momentum
        init = State(np.array([2.0]))

        sampler = MockedNCMCSampler(init, self.protocol, self.reverse_protocol)

        ## Test _propose
        # Make sure the first random number is < 0.5
        np.random.seed(5)
        result = sampler._propose()

        self.assertEqual(result.traj.heat, - 2 * 42)
        self.assertEqual(result.traj.work, 2 * 42)
        self.assertEqual(result.traj.initial.position[0], init.position[0])
        self.assertEqual(result.traj.final.position[0],  init.position[0] * 4)
        self.assertEqual(result.traj.initial.momentum, None)
        self.assertEqual(result.traj.final.momentum, None)

        # Make sure the first random number is > 0.5
        np.random.seed(4)
        result = sampler._propose()

        self.assertEqual(result.traj.heat, 2 * 42)
        self.assertEqual(result.traj.work, - 2 * 42)
        self.assertEqual(result.traj.initial.position[0], init.position[0])
        self.assertEqual(result.traj.final.position[0],  init.position[0] * 4)
        self.assertEqual(result.traj.initial.momentum, None)
        self.assertEqual(result.traj.final.momentum, None)


class HState(State):

    def clone(self):
        s = super(HState, self).clone()
        s.history = self.history
        
        return s


@test.functional
class TestReplicaHistory(test.Case):

    def setUp(self):

        pass
    
    def _runSimulation(self, n_replicas, swap_interval, first_swap):

        temperatures = np.linspace(1.0, 5.0, n_replicas)

        init_states = [HState(np.array([1.0])) for T in temperatures]

        for i, x in enumerate(init_states):
            x.history = []

        samplers = [RWMCSampler(SamplePDF(), init_states[i], stepsize=1.0, temperature=T) 
                    for i, T in enumerate(temperatures)]

        params = [RESwapParameterInfo(samplers[i], samplers[i+1]) for i in range(len(samplers) - 1)]

        algo = ReplicaExchangeMC(samplers, params)
        swapper = AlternatingAdjacentSwapScheme(algo)
        samples = []

        for i in range(500):
            if (i - first_swap) % swap_interval == 0 and i > 0 and i >= first_swap:
                swapper.swap_all()
            else:
                algo.sample()

            for j, s in enumerate(algo.state):
                s.history.append(j)

            samples.append(algo.state)

        return samples


    def _assertIdenticalHistories(self, samples, interval, first_swap=None):

        rh = ReplicaHistory(samples, interval, first_swap)

        for i in range(len(samples[0])):
            h = rh.calculate_history(i)
            for j, x in enumerate(samples[-1]):
                if x.history == h:
                    return True

        return False

    def _assertIdenticalProjTrajs(self, samples, interval, first_swap=None):

        rh = ReplicaHistory(samples, interval, first_swap)

        ## Calculate projected trajectories directly from test data history
        trajs1 = [Trajectory([x for x in [y[j] for y in samples] if x.history[0] == j]) 
                  for j in range(len(samples[0]))]
        
        ok = []
        for i in range(len(samples[0])):
            trajs2 = rh.calculate_projected_trajectories(i)
            ok.append(True in [np.all(np.array(t1) == np.array(t2)) for t1 in trajs1
                               for t2 in trajs2])
            
        return np.all(ok)

    def testTwoReplicas(self):

        swap_interval = 5
        first_swap = 5

        samples = self._runSimulation(2, swap_interval, first_swap)

        assert(self._assertIdenticalHistories(samples, swap_interval))
        assert(self._assertIdenticalProjTrajs(samples, swap_interval))

    def testFourReplicas(self):

        swap_interval = 5
        first_swap = 5

        samples = self._runSimulation(4, swap_interval, first_swap)

        assert(self._assertIdenticalHistories(samples, swap_interval))
        assert(self._assertIdenticalProjTrajs(samples, swap_interval))

    def testFiveReplicas(self):

        swap_interval = 5
        first_swap = 5
        
        samples = self._runSimulation(5, swap_interval, first_swap)

        assert(self._assertIdenticalHistories(samples, swap_interval))
        assert(self._assertIdenticalProjTrajs(samples, swap_interval))

    def testFiveReplicasOffset(self):

        swap_interval = 6
        first_swap = 7

        samples = self._runSimulation(5, swap_interval, first_swap)

        assert(self._assertIdenticalHistories(samples, swap_interval, first_swap))
        assert(self._assertIdenticalProjTrajs(samples, swap_interval, first_swap))


if __name__ == '__main__':

    test.Console()
