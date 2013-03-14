"""
Propagator class employing stepwise trajectories as used in the NCMC
algorithm (Nilmeier et al., "Nonequilibrium candidate Monte Carlo is 
an efficient tool for equilibrium simulation", PNAS 2011)
"""

import csb

import numpy

from abc import ABCMeta, abstractmethod
from csb.statistics.samplers.mc import TrajectoryBuilder, Trajectory
from csb.statistics.samplers.mc.propagators import AbstractPropagator, MDPropagator, HMCPropagator
from csb.numeric import InvertibleMatrix
from csb.numeric.integrators import FastLeapFrog

class NonequilibriumTrajectory(Trajectory):
    """
    Trajectory holding additional information about energy difference
    the Jacobian.

    @param items: sequence of trajectory states
    @type items: list of L{State}s

    @param heat: heat produced during the trajectory
    @type heat: float

    @param work: work expended during the trajectory
    @type work: float

    @param deltaH: energy difference between initial and final states
    @type deltaH: float

    @param jacobian: product of Jacobians of perturbations applied  in the
                     calculation of the trajectory
    @type jacobian: float
    """

    def __init__(self, items, heat=0.0, work=0.0, deltaH=0.0, jacobian=1.0):
        
        super(NonequilibriumTrajectory, self).__init__(items, heat=heat, work=work)

        self._deltaH = None
        self.deltaH = deltaH
        self._jacobian = None
        self.jacobian = jacobian

    @property
    def jacobian(self):
        return self._jacobian
    @jacobian.setter
    def jacobian(self, value):
        self._jacobian = value

    @property
    def deltaH(self):
        return self._deltaH
    @deltaH.setter
    def deltaH(self, value):
        self._deltaH = value


class AbstractSystemInfo(object):
    """
    Subclasses hold all information describing a current system setup
    (Hamiltonian, boundaries, ...)
    """
    
    pass
        
class PerturbationResult(object):
    """
    Instances hold the result of a perturbation.

    @param state: state resulting from perturbation
    @type state: L{State}
    
    @param perturbed_sys: L{AbstractSystemInfo} instance 
                          describing the perturbed system
    @type perturbed_sys: L{AbstractSystemInfo}
    
    @param work: work performed on the system during perturbation
    @type work: float
    
    @param jacobian: jacobian of the perturbation
    @type jacobian: float
    """
    
    def __init__(self, state, perturbed_sys, work=0., jacobian=1.):
        self._state = None
        self.state = state
        self._work = None
        self.work = work
        self._jacobian = None
        self.jacobian = jacobian
        self._perturbed_sys = None
        self.perturbed_sys = perturbed_sys

class Protocol(object):
    """
    Describes a stepwise protocol as in Nilmeier et al. (2011).

    @param steps: the steps making up the protocol
    @type steps: list of L{Step}s
    """

    def __init__(self, steps):

        self._steps = None
        self.steps = steps

    def reverse(self):
        """
        Reverses the protocol, that is, reverses the order of
        propagation and perturbation in each step.
        """

        raise NotImplementedError()

    @property
    def steps(self):
        """
        The steps making up the protocol
        """
        return self._steps
    @steps.setter
    def steps(self, value):
        self._steps = value
    
class Step(object):
    '''
    Defines a step in an NCMC-like stepwise protocol.

    @param perturbation: The perturbation of the system
    @type perturbation: L{AbstractPerturbation}

    @param propagation: The propagation of the perturbed system
    @type propagation: L{AbstractPropagation}
    '''
    
    def __init__(self, perturbation, propagation):

        self._perturbation = None
        self.perturbation = perturbation
        self._propagation = None
        self.propagation = propagation
        self._perform = None
        self.perform = self._perform_pert_prop

    def _perform_pert_prop(self, state):
        perturbation_result = self.perturbation(state)
        propagation_result = self.propagation(perturbation_result.state)

        result_state = propagation_result.final

        return NonequilibriumTrajectory([state, result_state],
                                        heat=propagation_result.heat,
                                        work=perturbation_result.work,
                                        jacobian=perturbation_result.jacobian)

    def _perform_prop_pert(self, state):
        propagation_result = self.propagation(state)
        perturbation_result = self.perturbation(propagation_result.final)

        result_state = perturbation_result.state

        return NonequilibriumTrajectory([state, result_state],
                                        heat=propagation_result.heat,
                                        work=perturbation_result.work,
                                        jacobian=perturbation_result.jacobian)

    def set_perturbation_first(self):
        """
        Perform first perturbation, then propagation
        """
        
        self.perform = self._perform_pert_prop

    def set_propagation_first(self):
        """
        Perform first propagation, then perturbation
        """
        
        self.perform = self._perform_prop_pert
        
    @property
    def perturbation(self):
        return self._perturbation
    @perturbation.setter
    def perturbation(self, value):
        self._perturbation = value

    @property
    def propagation(self):
        return self._propagation
    @propagation.setter
    def propagation(self, value):
        self._propagation = value
        

class ReducedHamiltonian(object):
    """
    Describes a reduced Hamiltonian (Hamiltonian, its position gradient
    and the system temperature)

    @param log_prob: log probability of the PDF under consideration, that is,
                     the negative potential energy of the system
    @type log_prob: callable

    @param gradient: gradient of the negative log probability of the PDF under
                     consideration, that is, the gradient of the potential energy;
                     function of position array and time
    @type gradient: callable

    @param temperature: system temperature
    @type temperature: float
    
    @param mass_matrix: system mass matrix
    @type mass_matrix: L{InvertibleMatrix}
    """                 

    def __init__(self, log_prob, gradient=None, temperature=1.0, mass_matrix=1.0):
        self._log_prob = None
        self.log_prob = log_prob
        self._gradient = None
        self.gradient = gradient
        self._temperature = None
        self.temperature = temperature
        self._mass_matrix = None
        self.mass_matrix = mass_matrix

    def E(self, x):
        """
        Potential energy of the system, aka negative log probability

        @param x: position vector
        @type x: 1D numpy array
        """
        
        return -self.log_prob(x)

    def kinetic_energy(self, p):
        """
        Kinetic energy of the system

        @param p: system momentum vector
        @type p: 1D numpy array
        """
        
        if p is not None:
            if self.mass_matrix == 1.0:
                return 0.5 * sum(p ** 2)
            else:
                return 0.5 * numpy.dot(p, numpy.dot(self.mass_matrix.inverse, p))
        else:
            return 0.0

    def rlog_prob(self, x):
        """
        Reduced log probability

        @param x: position vector
        @type x: 1D numpy array
        """
        
        return self.log_prob(x) / self.temperature

    def rkinetic_energy(self, p):
        """
        Reduced kinetic energy

        @param p: system momentum vector
        @type p: 1D numpy array
        """
        
        return self.kinetic_energy(p) / self.temperature

    def __call__(self, x):
        """
        Evaluates the reduced Hamiltionian at the state x

        @param x: system state
        @type x: L{State}
        """
        
        return -self.rlog_prob(x.position) + self.rkinetic_energy(x.momentum)

    @property
    def log_prob(self):
        return self._log_prob
    @log_prob.setter
    def log_prob(self, value):
        self._log_prob = value
    
    @property
    def gradient(self):
        return self._gradient
    @gradient.setter
    def gradient(self, value):
        self._gradient = value

    @property
    def temperature(self):
        return self._temperature
    @temperature.setter
    def temperature(self, value):
        self._temperature = value

    @property
    def mass_matrix(self):
        return self._mass_matrix
    @mass_matrix.setter
    def mass_matrix(self, value):
        self._mass_matrix = value

        
class AbstractPerturbation(object):
    """
    Describes an abstract system perturbation

    @param sys_before: information about the system before the perturbation
    @type sys_before: L{AbstractSystemInfo}

    @param sys_after: information about the system after the perturbation
    @type sys_after: L{AbstractSystemInfo}

    @param param: parameters neccessary for system perturbation
    @type param: L{AbstractPerturbationParam}
    """

    __metaclass__ = ABCMeta

    def __init__(self, sys_before, sys_after, param=None):
        self._sys_before = None
        self.sys_before = sys_before
        self._sys_after = None
        self.sys_after = sys_after
        self.param = param

    @abstractmethod
    def _evaluate(self, state):
        """
        Performs the perturbation of the system and / or the state

        @param state: system state
        @type state: L{State}
        """

        pass

    def __call__(self, state):
        """
        Performs the perturbation of the system and / or the state

        @param state: system state
        @type state: L{State}
        """
        
        return self._evaluate(state)

    @abstractmethod
    def reverse(self):
        """
        Reverses the perturbation of the system
        """

        pass

    @property
    def sys_before(self):
        return self._sys_before
    @sys_before.setter
    def sys_before(self, value):
        self._sys_before = value

    @property
    def sys_after(self):
        return self._sys_after
    @sys_after.setter
    def sys_after(self, value):
        self._sys_after = value
        
    @property
    def param(self):
        return self._param
    @param.setter
    def param(self, value):
        self._param = value
        
        
class AbstractPropagation(object):
    """
    Describes an abstract system propagation

    @param sys: information about the current system setup
    @type sys: L{AbstractSystemInfo}

    @param param: parameters neccessary for propagating the system
    @type param: L{AbstractPropagationParam}
    """

    def __init__(self, sys, param):

        self._sys = None
        self.sys = sys
        self._param = None
        self.param = param

    def _evaluate(self, state):
        """
        Performs the propagation of a state in the specified system

        @param state: system state
        @type state: L{State}
        """
        
        pass

    def __call__(self, state):
        """
        Performs the propagation of a state in the specified system

        @param state: system state
        @type state: L{State}
        """
        
        return self._evaluate(state)

    @property
    def sys(self):
        return self._sys
    @sys.setter
    def sys(self, value):
        self._sys = value
    
    @property
    def param(self):
        return self._param
    @param.setter
    def param(self, value):
        self._param = value


class ReducedHamiltonianPerturbation(AbstractPerturbation):
    """
    System perturbation by changing the reduced Hamiltonian

    @param sys_before: information about the system before the perturbation
    @type sys_before: L{AbstractSystemInfo}

    @param sys_after: information about the system after the perturbation
    @type sys_after: L{AbstractSystemInfo}
    """

    def __init__(self, sys_before, sys_after):

        super(ReducedHamiltonianPerturbation, self).__init__(sys_before, sys_after)

    def _evaluate(self, state):

        return PerturbationResult(state, self.sys_after,
                                  work=self.sys_after.hamiltonian(state) - \
                                       self.sys_before.hamiltonian(state))

    def reverse(self):

        self.sys_before, self.sys_after = self.sys_after, self.sys_before


class MDPerturbation(AbstractPerturbation):
    """
    System perturbation by letting the system evolve deterministically
    under the action of a parametrically switched Hamiltonian at constant
    temperature.
        
    @param sys_before: information about the system before the perturbation
    @type sys_before: L{AbstractSystemInfo}

    @param sys_after: information about the system after the perturbation
    @type sys_after: L{AbstractSystemInfo}

    @param param: parameters neccessary for system perturbation
    @type param: L{MDPerturbationParam}
    """
    
    def __init__(self, sys_before, sys_after, param):

        super(MDPerturbation, self).__init__(sys_before, sys_after)

        self._param = None
        self.param = param
        
    def _evaluate(self, state):

        init_state = state.clone()
        
        if self.param.mass_matrix is None:
            d = len(init_state.position)
            self.param.mass_matrix = InvertibleMatrix(numpy.eye(d), numpy.eye(d))
        
        if init_state.momentum is None:
            momentum = None
            if not self.param.mass_matrix.is_unity_multiple:
                covariance_matrix = self.sys_before.hamiltonian.temperature * self.param.mass_matrix
                momentum = numpy.random.multivariate_normal(mean=numpy.zeros(len(state.position)),
                                                            cov=covariance_matrix)
            else:
                variance = self.sys_before.hamiltonian.temperature * self.param.mass_matrix[0][0]
                momentum = numpy.random.normal(size=len(state.position),
                                               scale=numpy.sqrt(variance))

            init_state.momentum = momentum


        gen = MDPropagator(self.param.gradient, self.param.timestep,
                           mass_matrix=self.param.mass_matrix,
                           integrator=self.param.integrator)
        
        traj = gen.generate(init_state, self.param.traj_length)

        final = traj.final

        return PerturbationResult(final, self.sys_after,
                                  work=self.sys_after.hamiltonian(final) - \
                                       self.sys_before.hamiltonian(state))

    def reverse(self):

        self.sys_before, self.sys_after = self.sys_after, self.sys_before
        self.param.gradient = lambda x, t: self.param.gradient(x,
                                                               self.traj_length * self.timestep - t)

    @property
    def param(self):
        return self._param
    @param.setter
    def param(self, value):
        self._param = value

    
class HMCPropagation(AbstractPropagation):
    """
    System propagation by HMC

    @param sys: information about the current system setup
    @type sys: L{AbstractSystemInfo}

    @param param: parameters neccessary for propagating the system
    @type param: L{HMCPropagationParam}
    """
    
    def __init__(self, sys, param):

        self._param = None
        self.param = param
        self._sys = None
        self.sys = sys
        if self.param.gradient is None:
            self.param.gradient = self.sys.hamiltonian.gradient

    def _evaluate(self, state):

        from csb.statistics.samplers.mc import PropagationResult
        
        if self.param.mass_matrix is None:
            d = len(state.position)
            self.param.mass_matrix = InvertibleMatrix(numpy.eye(d), numpy.eye(d))
            
        gen = HMCPropagator(self.sys.hamiltonian, self.param.gradient,
                            self.param.timestep, self.param.traj_length,
                            temperature=self.sys.hamiltonian.temperature,
                            integrator=self.param.integrator, mass_matrix=self.param.mass_matrix)
        
        final = gen.generate(state, self.param.iterations).final

        return PropagationResult(state, final, heat=self.sys.hamiltonian.E(final.position) - \
                                                    self.sys.hamiltonian.E(state.position))

    @property
    def param(self):
        return self._param
    @param.setter
    def param(self, value):
        self._param = value

class AbstractPerturbationParam(object):
    """
    Subclasses hold informations required for different kinds
    of system perturbation
    """

    pass


class AbstractPropagationParam(object):
    """
    Subclasses hold informations required for different kinds
    of system propagation
    """

    pass


class MDParam(object):
    """
    Holds all required information for calculating a MD trajectory.

    @param timestep: integration timestep
    @type timestep: float

    @param traj_length: MD trajectory length
    @type traj_length: int

    @param gradient: gradient governing the equations of motion, function of
                     position array and time
    @type gradient: callable

    @param temperature: System temperature for drawing momenta from the
                        Maxwell distribution
    @type temperature: float
    
    @param integrator: Integration scheme to be utilized
    @type integrator: L{AbstractIntegrator}

    @param mass_matrix: mass matrix for kinetic energy definition
    @type mass_matrix: L{InvertibleMatrix}
    """
    
    def __init__(self, timestep, traj_length, gradient, temperature=1.0,
                 integrator=FastLeapFrog, mass_matrix=None):

        self._timestep = None
        self.timestep = timestep
        self._traj_length = None
        self.traj_length = traj_length
        self._gradient = None
        self.gradient = gradient
        self._temperature = None
        self.temperature = temperature
        self._integrator = None
        self.integrator = integrator
        self._mass_matrix = None
        self.mass_matrix = mass_matrix

    @property
    def timestep(self):
        return self._timestep
    @timestep.setter
    def timestep(self, value):
        self._timestep = value

    @property
    def traj_length(self):
        return self._traj_length
    @traj_length.setter
    def traj_length(self, value):
        self._traj_length = value

    @property
    def gradient(self):
        return self._gradient
    @gradient.setter
    def gradient(self, value):
        self._gradient = value

    @property
    def temperature(self):
        return self._temperature
    @temperature.setter
    def temperature(self, value):
        self._temperature = value
 
    @property
    def integrator(self):
        return self._integrator
    @integrator.setter
    def integrator(self, value):
        self._integrator = value
    
    @property
    def mass_matrix(self):
        return self._mass_matrix
    @mass_matrix.setter
    def mass_matrix(self, value):
        self._mass_matrix = value


class MDPerturbationParam(MDParam, AbstractPerturbationParam):
    """
    Holds all required information for perturbing a system by
    calculating a MD trajectory with time-dependent Hamiltonian, but
    constant temperature defined in 
    MDPerturbation.sys_before.hamiltonian.temperature.

    @param timestep: integration timestep
    @type timestep: float

    @param traj_length: MD trajectory length
    @type traj_length: int

    @param gradient: gradient governing the equations of motion, function of
                     position array and time
    @type gradient: callable
    
    @param integrator: Integration scheme to be utilized
    @type integrator: L{AbstractIntegrator}

    @param mass_matrix: mass matrix for kinetic energy definition
    @type mass_matrix: L{InvertibleMatrix}
    """

    def __init__(self, timestep, traj_length, gradient, integrator=FastLeapFrog,
                 mass_matrix=None):

        super(MDPerturbationParam, self).__init__(self, timestep, traj_length,
                                                  gradient, integrator=integrator,
                                                  mass_matrix=mass_matrix)

    pass
        
                                             
class HMCPropagationParam(MDParam, AbstractPropagationParam):
    """
    Holds all required information for propagating a system by HMC.
    The system temperature is taken from the 
    HMCPropagation.sys.hamiltonian.temperature property.
    
    @param timestep: integration timestep
    @type timestep: float

    @param traj_length: MD trajectory length
    @type traj_length: int

    @param gradient: gradient governing the equations of motion, function of
                     position array and time
    @type gradient: callable
    
    @param iterations: number of HMC iterations to be performed
    @type iterations: int

    @param integrator: Integration scheme to be utilized
    @type integrator: l{AbstractIntegrator}

    @param mass_matrix: mass matrix for kinetic energy definition
    @type mass_matrix: L{InvertibleMatrix}
    """
    
    def __init__(self, timestep, traj_length, gradient, iterations=1,
                 integrator=FastLeapFrog, mass_matrix=None):

        super(HMCPropagationParam, self).__init__(timestep, traj_length, gradient,
                                                  integrator=integrator, mass_matrix=mass_matrix)

        self._iterations = None
        self.iterations = iterations

    @property
    def iterations(self):
        return self._iterations
    @iterations.setter
    def iterations(self, value):
        self._iterations = value


class HamiltonianSysInfo(AbstractSystemInfo):
    """
    Holds information describing a system by its Hamiltonian only.
    
    @param hamiltonian: the Hamiltonian of the system to be described
    @type hamiltonian: L{ReducedHamiltonian}
    """
    
    def __init__(self, hamiltonian):
        
        self._hamiltonian = None
        self.hamiltonian = hamiltonian
    
    @property
    def hamiltonian(self):
        return self._hamiltonian
    @hamiltonian.setter
    def hamiltonian(self, value):
        self._hamiltonian = value


class NonequilibriumStepPropagator(AbstractPropagator):
    """
    Propagator class which propagates a system using NCMC-like
    stepwise trajectories
    
    @param protocol: stepwise protocol to be followed
    @type protocol: L{Protocol}
    """
    
    def __init__(self, protocol):

        self._protocol = None
        self.protocol = protocol

    def generate(self, init_state, return_trajectory=False):

        builder = TrajectoryBuilder.create(full=return_trajectory)
        builder.add_initial_state(init_state)
            
        work = 0.
        heat = 0.
        total_jacobian = 1.

        estate = init_state.clone()

        for i in range(len(self.protocol.steps)):
            shorttraj = self.protocol.steps[i].perform(estate)

            heat += shorttraj.heat
            work += shorttraj.work
            total_jacobian *= shorttraj.jacobian

            estate = shorttraj.final
            if i != len(self.protocol.steps) - 1:
                builder.add_intermediate_state(estate)
            else:
                builder.add_final_state(estate)

        traj = builder.product
        result = None

        red_heat = heat
        red_work = work
        red_deltaH = self.protocol.steps[-1].perturbation.sys_after.hamiltonian(traj.final) -\
                     self.protocol.steps[0].perturbation.sys_before.hamiltonian(traj.initial)

        if return_trajectory:
            result = NonequilibriumTrajectory([x for x in traj], heat=red_heat, work=red_work,
                                              deltaH=red_deltaH, jacobian=total_jacobian)
        else:
            result = NonequilibriumTrajectory([traj.initial, traj.final], heat=red_heat,
                                              work=red_work, deltaH=red_deltaH, 
                                              jacobian=total_jacobian)

        return result

    @property
    def protocol(self):
        return self._protocol
    @protocol.setter
    def protocol(self, value):
        self._protocol = value

    @property
    def verbose(self):
        return self._verbose
    @verbose.setter
    def verbose(self, value):
        self._verbose = value
