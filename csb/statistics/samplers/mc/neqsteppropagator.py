"""
Propagator class employing stepwise trajectories as used in the NCMC
algorithm (Nilmeier et al., "Nonequilibrium candidate Monte Carlo is 
an efficient tool for equilibrium simulation", PNAS 2011)
"""

import csb

import numpy

from abc import ABCMeta, abstractmethod
from csb.statistics.samplers.mc import TrajectoryBuilder, Trajectory, augment_state, PropagationResult
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
        
class PerturbationResult(Trajectory):
    """
    Instances hold the result of a perturbation.

    @param items: list of states defining a phase-space trajectory
    @type items: list of L{AbstractState}s
    
    @param work: work performed on the system during perturbation
    @type work: float
    
    @param jacobian: jacobian of the perturbation
    @type jacobian: float
    
    @param perturbed_sys: L{AbstractSystemInfo} instance 
                          describing the perturbed system
    @type perturbed_sys: L{AbstractSystemInfo}
    """
    
    def __init__(self, items, perturbed_sys, work, heat=0.0, jacobian=1.0):

        super(PerturbationResult, self).__init__(items, heat, work)

        self._jacobian = None
        self.jacobian = jacobian
        self._perturbed_sys = None
        self.perturbed_sys = perturbed_sys

    @property
    def jacobian(self):
        return self._jacobian
    @jacobian.setter
    def jacobian(self, value):
        self._jacobian = value

    @property
    def perturbed_sys(self):
        return self._perturbed_sys
    @perturbed_sys.setter
    def perturbed_sys(self, value):
        self._perturbed_sys = value

        
class Protocol(object):
    """
    Describes a stepwise protocol as in Nilmeier et al. (2011).

    @param steps: the steps making up the protocol
    @type steps: list of L{Step}s
    """

    def __init__(self, steps):

        self._steps = None
        self.steps = steps

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
    """
    Defines a step in an NCMC-like stepwise protocol.

    @param perturbation: The perturbation of the system
    @type perturbation: L{AbstractPerturbation}

    @param propagation: The propagation of the perturbed system
    @type propagation: L{AbstractPropagation}
    """
        
    def __init__(self, perturbation, propagation):

        self._perturbation = None
        self.perturbation = perturbation
        self._propagation = None
        self.propagation = propagation
        self._perform = None
        self.perform = self._perform_pert_prop

    def _perform_pert_prop(self, state):

        perturbation_result = self.perturbation(state)
        propagation_result = self.propagation(perturbation_result.final)

        result_state = propagation_result.final

        return NonequilibriumTrajectory([state, result_state],
                                        heat=propagation_result.heat,
                                        work=perturbation_result.work,
                                        jacobian=perturbation_result.jacobian)

    def _perform_prop_pert(self, state):

        propagation_result = self.propagation(state)
        perturbation_result = self.perturbation(propagation_result.final)
        result_state = perturbation_result.final

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

    def __init__(self, log_prob, gradient=None, temperature=1.0, mass_matrix=None):
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
            if self.mass_matrix is None:
                return 0.5 * sum(p ** 2)
            else:
                if self.mass_matrix.is_unity_multiple:
                    return 0.5 * sum(p ** 2) / self.mass_matrix[0][0]
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

    @param evaluate_work: Allows to switch off the work evaluation,
                          which might not always be needed, in order to
                          save computation time.
    @type evaluate_work: boolean
    """

    __metaclass__ = ABCMeta

    def __init__(self, sys_before, sys_after, param=None, evaluate_work=True):
        self._sys_before = None
        self.sys_before = sys_before
        self._sys_after = None
        self.sys_after = sys_after
        self.param = param
        self._evaluate_work = None
        self.evaluate_work = evaluate_work

    @abstractmethod
    def _run_perturbator(self, state):
        """
        Calculates the trajectory of the system while it is being perturbed.

        @param state: The initial system state
        @type state: L{State}

        @return: The trajectory of the system while it is being perturbed
        @rtype: L{Trajectory}
        """

        pass

    @abstractmethod
    def _calculate_work(self, traj):
        """
        Calculates the work expended during perturbation of the system.

        @param traj: The trajectory of the system while being perturbed
        @type traj: L{Trajectory}

        @return: The work expended during perturbation
        @rtype: float
        """

        pass

    @abstractmethod
    def _calculate_jacobian(self, traj):
        """
        Calculates the Jacobian determinant which reflects phase
        space compression during perturbation.

        @param traj: The trajectory of the system while being perturbed
        @type traj: L{Trajectory}

        @return: The Jacobian determinant
        @rtype: float
        """

        pass

    def _evaluate(self, state):
        """
        Performs the perturbation of the system and / or the state

        @param state: system state
        @type state: L{State}
        """

        traj = self._run_perturbator(state)
        work = self._calculate_work(traj)
        jacobian = self._calculate_jacobian(traj)

        return PerturbationResult([traj.initial, traj.final], self.sys_after, 
                                  work, jacobian=jacobian)
        
    def __call__(self, state):
        """
        Performs the perturbation of the system and / or the state

        @param state: system state
        @type state: L{State}
        """
        
        return self._evaluate(state)

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

    @property
    def evaluate_work(self):
        return self._evaluate_work
    @evaluate_work.setter
    def evaluate_work(self, value):
        self._evaluate_work = value
        
        
class AbstractPropagation(object):
    """
    Describes an abstract system propagation

    @param sys: information about the current system setup
    @type sys: L{AbstractSystemInfo}

    @param param: parameters neccessary for propagating the system
    @type param: L{AbstractPropagationParam}

    @param evaluate_heat: Allows to switch off the heat evaluation,
                          which might not always be needed, in order to
                          save computation time.
    @type evaluate_heat: boolean
    """

    __metaclass__ = ABCMeta

    def __init__(self, sys, param, evaluate_heat=True):

        self._sys = None
        self.sys = sys
        self._param = None
        self.param = param
        self._evaluate_heat = None
        self.evaluate_heat = evaluate_heat

    @abstractmethod
    def _propagator_factory(self):
        """
        Factory method which returns the propagator to be used for
        propagating the system.

        @return: Some propagator object
        @rtype: L{AbstractPropagator}
        """

    @abstractmethod
    def _run_propagator(self, state):
        """
        Propagates the system using the propagator instance returned
        by _propagator_factory().

        @param state: Initial state
        @type state: L{State}

        @return: The result of the propagation
        @rtype: L{PropagationResult}
        """
        
        pass

    @abstractmethod
    def _calculate_heat(self, traj):
        """
        Calculates the heat resulting from system propagation.

        @param traj: The trajectory of the system during propagation
        @type traj: L{Trajectory}

        @return: The heat resulting from system propagation.
        @rtype: float
        """
        
        pass

    def _evaluate(self, state):
        """
        Performs the propagation of a state in the specified system

        @param state: system state
        @type state: L{State}
        """
        
        traj = self._run_propagator(state)
        heat = self._calculate_heat(traj)
        
        return PropagationResult(traj.initial, traj.final, heat=heat)

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
    
    @property
    def evaluate_heat(self):
        return self._evaluate_heat
    @evaluate_heat.setter
    def evaluate_heat(self, value):
        self._evaluate_heat = value

    
class ReducedHamiltonianPerturbation(AbstractPerturbation):
    """
    System perturbation by changing the reduced Hamiltonian

    @param sys_before: information about the system before the perturbation
    @type sys_before: L{AbstractSystemInfo}

    @param sys_after: information about the system after the perturbation
    @type sys_after: L{AbstractSystemInfo}

    @param evaluate_work: Allows to switch off the work evaluation,
                          which might not always be needed, in order to
                          save computation time.
    @type evaluate_work: boolean
    """

    def __init__(self, sys_before, sys_after, evaluate_work=True):

        super(ReducedHamiltonianPerturbation, self).__init__(sys_before, sys_after,
                                                             evaluate_work=evaluate_work)

    def _calculate_work(self, traj):

        work = 0.0
        if self.evaluate_work == True:
            work = self.sys_after.hamiltonian(traj.final) - \
                   self.sys_before.hamiltonian(traj.initial)

        return work

    def _calculate_jacobian(self, traj):

        return 1.0

    def _run_perturbator(self, state):

        return Trajectory([state, state])


class AbstractMCPropagation(AbstractPropagation):
    """
    System propagation by some MC algorithm.

    @param sys: information about the current system setup
    @type sys: L{AbstractSystemInfo}

    @param param: parameters neccessary for propagating the system
    @type param: L{AbstractPropagationParam}

    @param evaluate_heat: Allows to switch off the heat evaluation,
                          which might not always be needed, in order to
                          save computation time.
    @type evaluate_heat: boolean
    """

    ## Not neccessary, but otherwise epydoc complains
    def __init__(self, sys, param, evaluate_heat=True):

        super(AbstractMCPropagation, self).__init__(sys, param, evaluate_heat=True)

    def _calculate_heat(self, traj):
        
        heat = 0.0
        if self.evaluate_heat == True:
            heat = self.sys.hamiltonian.E(traj.final.position) - \
                   self.sys.hamiltonian.E(traj.initial.position)

        return heat

    def _run_propagator(self, state):
            
        gen = self._propagator_factory()
        
        return gen.generate(state, self.param.iterations)
    
class HMCPropagation(AbstractMCPropagation):
    """
    System propagation by HMC

    @param sys: information about the current system setup
    @type sys: L{AbstractSystemInfo}

    @param param: parameters neccessary for propagating the system
    @type param: L{HMCPropagationParam}

    @param evaluate_heat: Allows to switch off the heat evaluation,
                          which might not always be needed, in order to
                          save computation time.
    @type evaluate_heat: boolean
    """
    
    def __init__(self, sys, param, evaluate_heat=True):

        super(HMCPropagation, self).__init__(sys, param, evaluate_heat)
        
        if self.param.gradient is None:
            self.param.gradient = self.sys.hamiltonian.gradient

    def _set_mass_matrix(self, state):
        """
        Sets the mass matrix in the param object.

        @param state: The initial state which is used to determine the dimension
                      of the mass matrix
        @type state: L{State}
        """

        if self.param.mass_matrix is None:
            d = len(state.position)
            self.param.mass_matrix = InvertibleMatrix(numpy.eye(d))

    def _propagator_factory(self):

        return HMCPropagator(self.sys.hamiltonian, self.param.gradient,
                             self.param.timestep, self.param.traj_length,
                             temperature=self.sys.hamiltonian.temperature,
                             integrator=self.param.integrator, mass_matrix=self.param.mass_matrix)

    def _evaluate(self, state):
        
        self._set_mass_matrix(state)

        return super(HMCPropagation, self)._evaluate(state)

    @property
    def param(self):
        return self._param
    @param.setter
    def param(self, value):
        self._param = value


class AbstractMDPropagation(AbstractPropagation):
    """
    System propagation by some MD algorithm

    @param sys: information about the current system setup
    @type sys: L{AbstractSystemInfo}

    @param param: parameters neccessary for propagating the system
    @type param: L{MDPropagationParam}

    @param evaluate_heat: Allows to switch off the heat evaluation,
                          which might not always be needed, in order to
                          save computation time.
    @type evaluate_heat: boolean
    """

    __metaclass__ = ABCMeta

    ## Not neccessary, but otherwise epydoc complains
    def __init__(self, sys, param, evaluate_heat=True):

        super(AbstractMDPropagation, self).__init__(sys, param, evaluate_heat=True)
    
    def _set_mass_matrix(self, state):
        """
        Sets the mass matrix in the param object.

        @param state: The initial state which is used to determine the dimension
                      of the mass matrix
        @type state: L{State}
        """

        if self.param.mass_matrix is None:
            d = len(state.position)
            self.param.mass_matrix = InvertibleMatrix(numpy.eye(d))

    def _augment_state(self, state):
        """
        Augments the initial state by a momentum if none is defined.

        @param state: Initial state
        @type state: L{State}
        """

        if state.momentum == None:
            state = augment_state(state, self.sys.hamiltonian.temperature,
                                  self.param.mass_matrix)

        return state

    def _run_propagator(self, state):
        
        gen = self._propagator_factory()
        state = self._augment_state(state)
        traj = gen.generate(state, self.param.traj_length)

        return traj

        
class PlainMDPropagation(AbstractMDPropagation):
    """
    System propagation by plain, microcanonical MD

    @param sys: information about the current system setup
    @type sys: L{AbstractSystemInfo}

    @param param: parameters neccessary for propagating the system
    @type param: L{PlainMDPropagationParam}

    @param evaluate_heat: Allows to switch off the heat evaluation,
                          which might not always be needed, in order to
                          save computation time.
    @type evaluate_heat: boolean
    """

    ## Not neccessary, but otherwise epydoc is complaining
    def __init__(self, sys, param, evaluate_heat=True):

        super(PlainMDPropagation, self).__init__(sys, param, evaluate_heat=evaluate_heat)

    def _propagator_factory(self):

        return MDPropagator(self.param.gradient, self.param.timestep,
                            mass_matrix=self.param.mass_matrix,
                            integrator=self.param.integrator)

    def _calculate_heat(self, traj):
        
        heat = 0.0
        if self.evaluate_heat == True:
            heat = self.sys.hamiltonian(traj.final) - self.sys.hamiltonian(traj.initial)

        return heat

    def _evaluate(self, state):

        self._set_mass_matrix(state)

        return super(PlainMDPropagation, self)._evaluate(state)
        
        
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


class MDPropagationParam(MDParam, AbstractPropagationParam):

    pass


class PlainMDPropagationParam(MDParam, AbstractPropagationParam):
    """
    Holds all required information for propagating a system by MD.
    The system temperature is taken from the 
    MDPropagation.sys.hamiltonian.temperature property.
    
    @param timestep: integration timestep
    @type timestep: float

    @param traj_length: MD trajectory length
    @type traj_length: int

    @param gradient: gradient governing the equations of motion, function of
                     position array and time
    @type gradient: callable

    @param integrator: Integration scheme to be utilized
    @type integrator: l{AbstractIntegrator}

    @param mass_matrix: mass matrix for kinetic energy definition
    @type mass_matrix: L{InvertibleMatrix}
    """
    
    def __init__(self, timestep, traj_length, gradient,
                 integrator=FastLeapFrog, mass_matrix=None):

        super(PlainMDPropagationParam, self).__init__(timestep, traj_length, gradient,
                                                 integrator=integrator, mass_matrix=mass_matrix)


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

    def _calculate_deltaH(self, traj):
        """
        Calculate the difference of the Hamiltonian between the initial and
        the final state of a NCMC trajectory.

        @param traj: The NCMC trajectory between whose initial and final states
                     the Hamiltonian difference should be calculated
        @type traj: L{NonequilibriumTrajectory}
        """

        return self.protocol.steps[-1].perturbation.sys_after.hamiltonian(traj.final) - \
               self.protocol.steps[0].perturbation.sys_before.hamiltonian(traj.initial)
        

    def generate(self, init_state, return_trajectory=False):

        estate = init_state.clone()
        
        reduced_work = 0.
        reduced_heat = 0.

        builder = TrajectoryBuilder.create(full=return_trajectory)
            
        total_jacobian = 1.

        for i in range(len(self.protocol.steps)):
            shorttraj = self.protocol.steps[i].perform(estate)
            if i == 0:
                builder.add_initial_state(shorttraj.initial)

            reduced_heat += shorttraj.heat
            reduced_work += shorttraj.work
            total_jacobian *= shorttraj.jacobian

            estate = shorttraj.final
            if i != len(self.protocol.steps) - 1:
                builder.add_intermediate_state(estate)
            else:
                builder.add_final_state(estate)

        traj = builder.product

        reduced_deltaH = self._calculate_deltaH(traj)
        
        if init_state.momentum is None:
            for s in traj:
                s.momentum = None
        
        result = NonequilibriumTrajectory([x for x in traj],
                                           heat=reduced_heat,
                                           work=reduced_work,
                                           deltaH=reduced_deltaH,
                                           jacobian=total_jacobian)
    
        return result

    @property
    def protocol(self):
        return self._protocol
    @protocol.setter
    def protocol(self, value):
        self._protocol = value
