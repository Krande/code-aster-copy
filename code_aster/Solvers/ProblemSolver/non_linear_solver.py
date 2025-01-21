# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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

from ...Messages import MessageLog, UTMESS
from ...Objects import NonLinearResult, ThermalResult, HHO
from ...Supervis import ConvergenceError, IntegrationError, SolverError
from ...Utilities import DEBUG, logger, no_new_attributes, profile
from ..Basics import ProblemType as PBT
from ..Basics import SolverFeature
from ..Basics import SolverOptions as SOP


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
        self.step_rank = None
        self.current_matrix = None

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

    def _storeState(self, state, ignore_policy=False):
        """Store the physical state.

        Arguments:
            time (float): current (pseudo)-time.
            ignore_policy (bool): ignore storing-policy.

        Returns:
            bool: *True* if it was actually stored, else *False*.
        """
        storage_manager = self.get_feature(SOP.Storage)
        return storage_manager.storeState(
            self.step_rank,
            state.time_curr,
            self.phys_pb,
            state,
            is_final_time=self.isFinished(),
            ignore_policy=ignore_policy,
        )

    @profile
    def initialize(self):
        """Initialize run"""
        self.check_features()
        args = self.get_feature(SOP.Keywords)
        # essential to be called enough soon (may change the size of VARI field)
        phys_pb = self.phys_pb
        if self._get("ETAT_INIT"):
            phys_pb.computeBehaviourProperty(args["COMPORTEMENT"], "OUI", 2)
        else:
            phys_pb.computeBehaviourProperty(args["COMPORTEMENT"], "NON", 2)
        phys_pb.computeListOfLoads()
        phys_pb.computeDOFNumbering()
        if phys_pb.getMaterialField().hasExternalStateVariableForLoad():
            phys_pb.computeReferenceExternalStateVariables()
        self.step_rank = 0
        self.setInitialState()
        self._storeState(self.phys_state)
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
        init_state = self._get("ETAT_INIT")
        nume_equa = self.phys_pb.getDOFNumbering().getEquationNumbering()

        if init_state:

            def _msginit(field, result=None):
                if result:
                    logger.info(MessageLog.GetText("I", "ETATINIT_32", valk=(field, result)))
                else:
                    logger.info(MessageLog.GetText("I", "ETATINIT_33", valk=(field,)))

            def _extract_param(init_state, resu):
                """Extract parameters for getField()."""
                extract_time = init_state.get("INST")
                if extract_time is None:
                    extract_time = resu.getLastTime()
                if init_state.get("NUME_ORDRE"):
                    para, value = "NUME_ORDRE", init_state["NUME_ORDRE"]
                else:
                    para, value = "INST", extract_time
                return para, value

            def _raise_elga_error():
                """Raise an error to make sure that initial_state model and
                field model are the same."""
                UTMESS("F", "MECANONLINE_9")

            def _extract_resu_field_and_check_model(resu, para, val, name_field, model):
                """Extract the field from the result, then check that the model
                of the field is the same as the model given in parameter.

                Returns the extracted field if the model is the same"""
                field = resu.getField(name_field, para=para, value=val)
                if model is field.getModel():
                    return field
                else:
                    _raise_elga_error()

            def _get_field_and_check_model(state, name_field, model):
                """Extract the field from the initial state, then check that
                the model of the field is the same as the model given in parameter.

                Returns the field if the model is the same"""
                field = state.get(name_field)
                if model is field.getModel():
                    return field
                else:
                    _raise_elga_error()

            model = self.phys_pb.getModel()

            if "EVOL_NOLI" in init_state:
                resu = init_state.get("EVOL_NOLI")
                assert isinstance(resu, NonLinearResult), resu
                para, value = _extract_param(init_state, resu)

                phys_state.primal_curr = resu.getField(
                    "DEPL", para=para, value=value
                ).copyUsingDescription(nume_equa, False)
                _msginit("DEPL", resu.userName)

                if phys_state.pb_type == PBT.MecaDyna:
                    phys_state.current.dU = resu.getField(
                        "VITE", para=para, value=value
                    ).copyUsingDescription(nume_equa)
                    _msginit("VITE", resu.userName)

                    phys_state.current.d2U = resu.getField(
                        "ACCE", para=para, value=value
                    ).copyUsingDescription(nume_equa)
                    _msginit("ACCE", resu.userName)

                phys_state.stress = _extract_resu_field_and_check_model(
                    resu=resu, para=para, val=value, name_field="SIEF_ELGA", model=model
                )
                _msginit("SIEF_ELGA", resu.userName)

                phys_state.internVar = _extract_resu_field_and_check_model(
                    resu=resu, para=para, val=value, name_field="VARI_ELGA", model=model
                )
                _msginit("VARI_ELGA", resu.userName)

                list_of_loads = self.phys_pb.getListOfLoads()
                if list_of_loads.hasDifferentialLoads():
                    nume_didi = init_state.get("NUME_DIDI")
                    if nume_didi:
                        displ = resu.getField("DEPL", nume_didi).copyUsingDescription(nume_equa)
                    else:
                        displ = phys_state.primal_curr
                    list_of_loads.setDifferentialDisplacement(displ)

            if "EVOL_THER" in init_state:
                resu = init_state.get("EVOL_THER")
                assert isinstance(resu, ThermalResult), resu
                para, value = _extract_param(init_state, resu)

                phys_state.primal_curr = resu.getField(
                    "TEMP", para=para, value=value
                ).copyUsingDescription(nume_equa)

            if "CHAM_NO" in init_state:
                phys_state.primal_curr = init_state.get("CHAM_NO").copyUsingDescription(nume_equa)

            if "DEPL" in init_state:
                phys_state.primal_curr = init_state.get("DEPL").copyUsingDescription(
                    nume_equa, False
                )
                _msginit("DEPL")

            if "SIGM" in init_state:
                phys_state.stress = _get_field_and_check_model(
                    state=init_state, name_field="SIGM", model=model
                )
                _msginit("SIEF_ELGA")

            if "VARI" in init_state:
                phys_state.internVar = _get_field_and_check_model(
                    state=init_state, name_field="VARI", model=model
                )
                _msginit("VARI_ELGA")

            if "VITE" in init_state:
                phys_state.current.dU = init_state.get("VITE").copyUsingDescription(nume_equa)
                _msginit("VITE")

            if "ACCE" in init_state:
                phys_state.current.d2U = init_state.get("ACCE").copyUsingDescription(nume_equa)
                _msginit("ACCE")

            if "VALE" in init_state:
                if model.existsHHO():
                    phys_state.primal_curr = HHO(self.phys_pb).projectOnHHOSpace(init_state["VALE"])
                else:
                    phys_state.primal_curr = phys_state.createPrimal(
                        self.phys_pb, value={"TEMP": init_state.get("VALE")}
                    )

        init_time = self.stepper.getInitial()
        self.computeExternalStateVariables(init_time)
        phys_state.time_curr = init_time

        if init_state:
            if init_state.get("STAT") == "OUI":
                solv = self.get_feature(SOP.StepSolver)
                solv.initialize()
                args = {"valr": phys_state.time_curr, "vali": self.stepper.splitting_level}
                logger.info(MessageLog.GetText("I", "MECANONLINE6_5", **args))
                solv.solve()
                if (
                    self.stepper.size() == 1
                    and self.stepper.getCurrent() == self.stepper.getPrevious()
                ):
                    self.stepper.completed()

            self.post_hooks()

        phys_state.commit()

    @profile
    def run(self):
        """Solve the problem."""
        self.initialize()
        matr_update_step = self._get("NEWTON", "REAC_INCR", 1)

        # Solve nonlinear problem
        solv = self.get_feature(SOP.StepSolver)
        last_stored = False
        while not self.isFinished():
            self.phys_state.time_curr = self.stepper.getCurrent()
            self.phys_state.time_step = self.phys_state.time_curr - self.phys_state.time_prev
            if self.stepper.splitting_level <= 0:
                logger.info(
                    MessageLog.GetText("I", "MECANONLINE6_7", valr=self.phys_state.time_curr)
                )
            else:
                args = dict(valr=self.phys_state.time_curr, vali=self.stepper.splitting_level)
                logger.info(MessageLog.GetText("I", "MECANONLINE6_5", **args))

            self.computeExternalStateVariables(self.phys_state.time_curr)
            solv.initialize()

            if matr_update_step == 0 or (self.step_rank + 1) % matr_update_step:
                solv.current_matrix = self.current_matrix
            else:
                solv.current_matrix = None

            if logger.getEffectiveLevel() <= DEBUG:
                self.phys_state.debugPrint("<t-> ")
            self.phys_state.stash()
            try:
                solv.solve()
            except (ConvergenceError, IntegrationError, SolverError) as exc:
                logger.warning(exc.format("I"))
                try:
                    self.stepper.failed(exc)
                except (ConvergenceError, IntegrationError, SolverError):
                    # an error occurred, ensure that the previous step was stored
                    logger.warning(
                        "An error occurred, ensure that the last converged step is saved"
                    )
                    self._storeState(self.phys_state.getState(-1), ignore_policy=True)
                    raise
            else:
                if not self.stepper.check_event(self.phys_state):
                    # + reset current_matrix to None (REAC_INCR)
                    self.phys_state.revert()
                    continue
                self.post_hooks()
                self.phys_state.commit()
                self.stepper.completed()
                self.current_matrix = solv.current_matrix
                self.step_rank += 1
                last_stored = self._storeState(self.phys_state)
        # ensure that last step was stored
        if not last_stored:
            self._storeState(self.phys_state, ignore_policy=True)

    def post_hooks(self):
        """Call post hooks"""
        for hook in self.get_features(SOP.PostStepHook):
            hook(self)

    def computeExternalStateVariables(self, current_time):
        """Compute and set external variables in the physical state.

        Arguments:
            current_time (float): Current time value.
        """
        if self.phys_pb.getMaterialField().hasExternalStateVariable():
            self.phys_state.externVar = self.phys_pb.getExternalStateVariables(current_time)

    def _get(self, keyword, parameter=None, default=None):
        """ "Return a keyword value"""
        args = self.param
        if parameter is not None:
            if args.get(keyword) is None:
                return default
            return _F(args[keyword][0]).get(parameter, default)

        return args.get(keyword, default)
