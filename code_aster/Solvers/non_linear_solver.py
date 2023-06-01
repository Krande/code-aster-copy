# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
# This file is part of code_aster.
#
# code_aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# code_aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------

from ..Messages import MessageLog
from ..Objects import NonLinearResult
from ..Supervis import ConvergenceError, IntegrationError
from ..Utilities import DEBUG, logger, no_new_attributes, profile
from .solver_features import SolverFeature
from .solver_features import SolverOptions as SOP

"""
Notations to use for mechanics.

We assume that we compute the solution at N time steps (t^1, ..., t^N).
Moreover, we assume that the physical state is given at t^0 (and is entierely known).

The solution is known at (t^0, ..., t^{n-1}) and we search the solution at t^{n}.
At t^{n-1}, this is the previous solution (the variables end with '_prev').
At t^{n}, this is the current solution (the variables end with '_curr').

The [0, I] incremental solutions beetwen t^{n} and t^{n} are the incremental solution
(the variables end with '_incr'). It is typically the solution of Newton or Picard
solver at a given iteration.
The total different between t^{n-1} and t^{n} is the step solution (the variables end
with '_step'). This the sum of incremental solution.

For example for the displacement field and a Newton solver.

- displ^{n-1}: previous displacement (displ_prev)
- displ^{n}: current displacement (displ_curr)
- delta_u = -K^{-1}*R(u): incremental displacement (displ_incr)
- Delta_u = sum_{i=0}^I delta_u^i ( = u^{n} - u^{n-1} at convergence): step
  displacement (displ_step)
"""


class NonLinearSolver(SolverFeature):
    """Main object to solve a non linear problem."""

    provide = SOP.ProblemSolver
    required_features = [
        SOP.PhysicalProblem,
        SOP.PhysicalState,
        SOP.Storage,
        SOP.StepSolver,
        SOP.TimeStepper,
        SOP.Keywords,
    ]
    optional_features = [SOP.PostStepHook]

    step_rank = current_matrix = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self):
        super().__init__()

    # convenient shortcuts properties
    @property
    def stepper(self):
        """:py:class:`~.time_stepper.TimeStepper`: object to be used."""
        return self.get_feature(SOP.TimeStepper)

    @property
    def param(self):
        """dict: object to be used."""
        return self.get_feature(SOP.Keywords)

    def getResult(self):
        """Get the Result object.

        Returns:
            NonLinearResult: the Result object.
        """
        return self.get_feature(SOP.Storage).getResult()

    def isFinished(self):
        """Tell if there are steps to be computed.

        Returns:
            bool: *True* if there is no step to be computed, *False* otherwise.
        """
        return self.stepper.isFinished()

    def _storeRank(self, time):
        """Store the current physical state.

        Arguments:
            time (float): current (pseudo)-time.
        """
        storage_manager = self.get_feature(SOP.Storage)
        storage_manager.storeState(time, self.phys_pb, self.phys_state)
        storage_manager.completed()

    @profile
    def initialize(self):
        """Initialize run"""
        self.check_features()
        args = self.get_feature(SOP.Keywords)
        # essential to be called enough soon (may change the size of VARI field)
        phys_pb = self.phys_pb
        phys_pb.computeBehaviourProperty(args["COMPORTEMENT"], "NON", 2)
        phys_pb.computeListOfLoads()
        phys_pb.computeDOFNumbering()
        if phys_pb.getMaterialField().hasExternalStateVariableForLoad():
            phys_pb.computeReferenceExternalStateVariables()
        self.setInitialState()
        self.step_rank = 0
        self._storeRank(self.phys_state.time)
        # register observers
        for source in self.get_childs(SOP.IncrementalSolver | SOP.EventSource):
            source.add_observer(self.stepper)
        for source in self.get_childs(SOP.ConvergenceCriteria | SOP.EventSource):
            source.add_observer(self.stepper)

    @profile
    def setInitialState(self):
        """Initialize the physical state."""
        phys_state = self.phys_state
        phys_state.zeroInitialState(self.phys_pb)
        init_time = self.stepper.getInitial()
        init_state = self._get("ETAT_INIT")
        if init_state:
            if "INST_ETAT_INIT" in init_state:
                init_time = init_state.get("INST_ETAT_INIT")
            if "EVOL_NOLI" in init_state:
                resu = init_state.get("EVOL_NOLI")
                assert isinstance(resu, NonLinearResult), resu
                extract_time = init_state.get("INST")
                if extract_time is None:
                    extract_time = resu.getLastTime()
                if init_time is None:
                    init_time = extract_time
                phys_state.primal = resu.getField("DEPL", para="INST", value=extract_time)
                phys_state.stress = resu.getField("SIEF_ELGA", para="INST", value=extract_time)
                phys_state.internVar = resu.getField("VARI_ELGA", para="INST", value=extract_time)
            if "DEPL" in init_state:
                phys_state.primal = init_state.get("DEPL")
            if "SIGM" in init_state:
                phys_state.stress = init_state.get("SIGM")
            if "VARI" in init_state:
                phys_state.internVar = init_state.get("VARI")

            if init_time is not None:
                self.stepper.setInitial(init_time)

        phys_state.time = init_time

    @profile
    def run(self):
        """Solve the problem."""
        self.initialize()

        solv = self.get_feature(SOP.StepSolver)

        matr_update_step = self._get("NEWTON", "REAC_ITER", 1)

        # Create field for external state variables
        if self.phys_pb.getMaterialField().hasExternalStateVariable():
            self.phys_state.externVar = self.phys_pb.getExternalStateVariables(self.phys_state.time)

        # Solve nonlinear problem
        while not self.isFinished():
            timeEndStep = self.stepper.getCurrent()
            self.phys_state.time_step = timeEndStep - self.phys_state.time
            if self.stepper.splitting_level <= 0:
                logger.info(MessageLog.GetText("I", "MECANONLINE6_7", valr=timeEndStep))
            else:
                args = dict(valr=timeEndStep, vali=self.stepper.splitting_level)
                logger.info(MessageLog.GetText("I", "MECANONLINE6_5", **args))

            if self.phys_pb.getMaterialField().hasExternalStateVariable():
                self.phys_state.externVar_next = self.phys_pb.getExternalStateVariables(timeEndStep)

            solv.initialize()
            if (self.step_rank + 1) % matr_update_step:
                solv.current_matrix = self.current_matrix

            if logger.getEffectiveLevel() <= DEBUG:
                self.phys_state.debugPrint("<t-> ")
            self.phys_state.stash()
            try:
                solv.solve()
            except (ConvergenceError, IntegrationError) as exc:
                logger.warning(exc.message)
                self.stepper.failed(exc)
            else:
                # DELTA_GRANDEUR, COLLISION, INTERPENETRATION, INSTABILITE (post_hook ?)
                if not self.stepper.check_event(self.phys_state):
                    # + reset current_matrix to None (REAC_INCR)
                    self.phys_state.revert()
                    continue
                self.phys_state.commit()
                self._storeRank(timeEndStep)
                self.stepper.completed()
                self.current_matrix = solv.current_matrix
                self.post_hooks()
                self.step_rank += 1

    def post_hooks(self):
        """Call post hooks"""
        for hook in self.get_features(SOP.PostStepHook):
            hook(self.phys_state)

    def _get(self, keyword, parameter=None, default=None):
        """ "Return a keyword value"""
        args = self.param
        if parameter is not None:
            if args.get(keyword) is None:
                return default
            return _F(args[keyword])[0].get(parameter, default)

        return args.get(keyword, default)
