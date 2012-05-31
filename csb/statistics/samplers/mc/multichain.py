"""
Implements several extended-ensemble Monte Carlo sampling algorithms.
"""

import csb.numeric

from csb.statistics.samplers.mc import AbstractExchangeMC, AbstractRENS, RESwapCommunicator
from csb.statistics.samplers.mc.propagators import MDPropagator, ThermostattedMDPropagator
from csb.numeric.integrators import AbstractGradient


class InterpolationFactory(object):
    """
    Produces interpolations for functions, changed during non-equilibrium
    trajectories.
    
    @param protocol: protocol to be used to generate non-equilibrium trajectories
    @type protocol: function mapping t to [0...1] for fixed tau
    @param tau: switching time
    @type tau: float    
    """
    
    def __init__(self, protocol, tau):
        
        self._protocol = None
        self._tau = None
        
        self.protocol = protocol
        self.tau = tau
        
    @property
    def protocol(self):
        return self._protocol
    @protocol.setter
    def protocol(self, value):
        if not hasattr(value, '__call__'):
            raise TypeError(value)
        self._protocol = value
                    
    @property
    def tau(self):
        return self._tau
    @tau.setter
    def tau(self, value):
        self._tau = float(value)
        
    def build_gradient(self, gradient):
        """
        Create a gradient instance with according to given protocol
        and switching time.
        
        @param gradient: gradient with G(0) = G_1 and G(1) = G_2
        @type gradient: callable    
        """
        return Gradient(gradient, self._protocol, self._tau)
    
    def build_temperature(self, temperature):
        """
        Create a temperature function according to given protocol and
        switching time.

        @param temperature: temperature with T(0) = T_1 and T(1) = T_2
        @type temperature: callable        
        """
        return lambda t: temperature(self.protocol(t, self.tau))
        
class Gradient(AbstractGradient):
    
    def __init__(self, gradient, protocol, tau):
        
        self._protocol = protocol
        self._gradient = gradient
        self._tau = tau
    
    def evaluate(self, q, t):
        return self._gradient(q, self._protocol(t, self._tau))

class ReplicaExchangeMC(AbstractExchangeMC):
    """
    Replica Exchange (RE, Swendsen & Yang 1986) implementation.
    """
        
    def _propose_swap(self, param_info):
        
        return RESwapCommunicator(param_info, param_info.sampler2.state, param_info.sampler1.state)
    
    def _calc_pacc_swap(self, swapcom):
        
        E1 = lambda x:-swapcom.sampler1._pdf.log_prob(x)
        E2 = lambda x:-swapcom.sampler2._pdf.log_prob(x)
        
        state1 = swapcom.sampler1.state
        state2 = swapcom.sampler2.state
        
        proposal1 = swapcom.proposal1
        proposal2 = swapcom.proposal2
        
        swapcom.acceptance_probability = csb.numeric.exp(-E1(proposal1.position) + E1(state1.position) - \
                                                 E2(proposal2.position) + E2(state2.position))

class MDRENS(AbstractRENS):
    """
    Replica Exchange with Nonequilibrium Switches (RENS, Ballard & Jarzynski 2009)
    with Molecular Dynamics (MD) trajectories.

    @param samplers: Samplers which sample their
                         respective equilibrium distributions
    @type samplers: list of L{AbstractSingleChainMC}

    @param param_infos: ParameterInfo instance holding
                        information required to perform a MDRENS swap
    @type param_infos: list of L{MDRENSSwapParameterInfo}

    @param integrator: Subclass of L{AbstractIntegrator} to be used to
                       calculate the non-equilibrium trajectories
    @type integrator: type
    """

    def __init__(self, samplers, param_infos, integrator):
        
        super(MDRENS, self).__init__(samplers, param_infos)
        
        self._integrator = integrator
        
    def _run_traj_generator(self, traj_info):
        
        tau = traj_info.param_info.traj_length * traj_info.param_info.timestep
        factory = InterpolationFactory(traj_info.protocol, tau)
                
        gen = MDPropagator(factory.build_gradient(traj_info.param_info.gradient),
                           traj_info.param_info.timestep, self._integrator)
        
        traj = gen.generate(traj_info.init_state, int(traj_info.param_info.traj_length))
        return traj

class ThermostattedMDRENS(MDRENS):
    """
    Replica Exchange with Nonequilibrium Switches (RENS, Ballard & Jarzynski, 2009)
    with Andersen-thermostatted Molecular Dynamics (MD) trajectories.

    @param samplers: Samplers which sample their
                         respective equilibrium distributions
    @type samplers: list of L{AbstractSingleChainMC}

    @param param_infos: ParameterInfo instance holding
                        information required to perform a MDRENS swap
    @type param_infos: list of L{ThermostattedMDRENSSwapParameterInfo}

    @param integrator: Subclass of L{AbstractIntegrator} to be used to
                       calculate the non-equilibrium trajectories
    @type integrator: type
    """

    def __init__(self, samplers, param_infos, integrator):
        
        super(ThermostattedMDRENS, self).__init__(samplers, param_infos, integrator)

    def _run_traj_generator(self, traj_info):
        
        tau = traj_info.param_info.traj_length * traj_info.param_info.timestep
        factory = InterpolationFactory(traj_info.protocol, tau)
        
        grad = factory.build_gradient(traj_info.param_info.gradient)
        temp = factory.build_temperature(traj_info.param_info.temperature)
        
        gen = ThermostattedMDPropagator(grad,
                                        traj_info.param_info.timestep, self._integrator,
                                        temp,
                                        traj_info.param_info.collision_probability,
                                        traj_info.param_info.collision_interval)
        
        traj = gen.generate(traj_info.init_state, traj_info.param_info.traj_length)

        return traj
