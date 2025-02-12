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

from libaster import deleteTemporaryObjects, resetFortranLoggingLevel, setFortranLoggingLevel

from ...Messages import UTMESS, MessageLog
from ...Objects import HHO, NonLinearResult, ThermalResult
from ...Supervis import ConvergenceError, IntegrationError, SolverError
from ...Utilities import (
    DEBUG,
    INFO,
    WARNING,
    ExecutionParameter,
    Options,
    logger,
    no_new_attributes,
    profile,
)
from ..Basics import ContextMixin
from ..Basics import ProblemType as PBT
from ..StepSolvers import BaseStepSolver
from .storage_manager import StorageManager
from .time_stepper import TimeStepper


class NonLinearOperator(ContextMixin):
    """Solver for linear and non linear problem.

    Arguments:
        main (*NonLinearFeature*): Main object.
        result (*misc*): The result object.
    """

    __needs__ = (
        "contact",
        "keywords",
        "linear_solver",
        "oper",
        "problem",
        "problem_type",
        "result",
        "state",
    )

    _stepper = _store = _step_solver = None
    _verb = None
    # FIXME: add _ prefix?
    _step_idx = current_matrix = None
    __setattr__ = no_new_attributes(object.__setattr__)

    @classmethod
    def builder(cls, context):
        """Builder of a NonLinearOperator object.

        Args:
            context (Context): Context of the problem.

        Returns:
            instance: New object.
        """
        # same as constructor
        return cls(context)

    def __init__(self, context) -> None:
        super().__init__()
        self.context = context
        self._step_idx = None
        self.current_matrix = None
        self._verb = logger.getEffectiveLevel(), ExecutionParameter().option & Options.ShowSyntax

    # convenient shortcuts properties to init and access subobjects
    @property
    def stepper(self):
        """:py:class:`~.time_stepper.TimeStepper`: object to be used."""
        if not self._stepper:
            logger.debug("+++ init Stepper")
            self._stepper = TimeStepper.from_keywords(**self.keywords["INCREMENT"])
        return self._stepper

    @property
    def store(self):
        if not self._store:
            logger.debug("+++ init StorageManager")
            kwds = self.keywords
            reuse = kwds.get("REUSE")
            self._store = StorageManager(
                self.result, kwds.get("ARCHIVAGE"), reused=reuse is self.result
            )
            if reuse:
                init_state = kwds.get("ETAT_INIT")
                assert init_state
                if "EVOL_NOLI" in init_state:
                    # Pour l'instant, on se restreint au cas où la sd passée
                    # par reuse est la même que celle passée dans ETAT_INIT
                    assert init_state["EVOL_NOLI"] is reuse
                # FIXME: access to _eps!
                init_index = reuse.getIndexFromParameter(
                    "INST", self.stepper.getInitial(), "RELATIF", self.stepper._eps
                )
                # if stepper is None ?! should be happen
                # init_index = reuse.getLastIndex()
                self._store.setFirstStorageIndex(init_index + 1)
        return self._store

    @property
    def step_solver(self):
        if self._step_solver:
            return self._step_solver
        logger.debug("+++ init StepSolver")
        self._step_solver = solv = BaseStepSolver.factory(self.context)
        return solv

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
        return self.store.storeState(
            self._step_idx,
            state.time_curr,
            self.problem,
            state,
            is_final_time=self.isFinished(),
            ignore_policy=ignore_policy,
        )

    @profile
    def initialize(self):
        """Initialize run"""
        phys_pb = self.problem
        kwds = self.keywords
        # essential to be called enough soon (may change the size of VARI field)
        if self.get_keyword("ETAT_INIT"):
            phys_pb.computeBehaviourProperty(kwds["COMPORTEMENT"], "OUI", 2)
        else:
            phys_pb.computeBehaviourProperty(kwds["COMPORTEMENT"], "NON", 2)
        phys_pb.computeListOfLoads()
        phys_pb.computeDOFNumbering()
        if phys_pb.getMaterialField().hasExternalStateVariableForLoad():
            phys_pb.computeReferenceExternalStateVariables()
        self._step_idx = 0
        self.setInitialState()
        self._storeState(self.state)
        # FIXME: register observers
        # for source in self.get_childs(SOP.IncrementalSolver | SOP.EventSource):
        #     source.add_observer(self.stepper)
        # for source in self.get_childs(SOP.ConvergenceCriteria | SOP.EventSource):
        #     source.add_observer(self.stepper)

    # FIXME: mixin by problem_type / factory
    @profile
    def setInitialState(self):
        """Initialize the physical state."""
        self.state.zeroInitialState(self.problem)
        init_state = self.get_keyword("ETAT_INIT")
        nume_equa = self.problem.getDOFNumbering().getEquationNumbering()
        if init_state:
            model = self.problem.getModel()

            if "EVOL_NOLI" in init_state:
                resu = init_state.get("EVOL_NOLI")
                assert isinstance(resu, NonLinearResult), resu
                para, value = _extract_param(init_state, resu)

                self.state.primal_curr = resu.getField(
                    "DEPL", para=para, value=value
                ).copyUsingDescription(nume_equa, False)
                _msginit("DEPL", resu.userName)

                if self.state.pb_type == PBT.MecaDyna:
                    self.state.current.dU = resu.getField(
                        "VITE", para=para, value=value
                    ).copyUsingDescription(nume_equa)
                    _msginit("VITE", resu.userName)

                    self.state.current.d2U = resu.getField(
                        "ACCE", para=para, value=value
                    ).copyUsingDescription(nume_equa)
                    _msginit("ACCE", resu.userName)

                self.state.stress = _extract_resu_field_and_check_model(
                    resu=resu, para=para, val=value, name_field="SIEF_ELGA", model=model
                )
                _msginit("SIEF_ELGA", resu.userName)

                self.state.internVar = _extract_resu_field_and_check_model(
                    resu=resu, para=para, val=value, name_field="VARI_ELGA", model=model
                )
                _msginit("VARI_ELGA", resu.userName)

            if "EVOL_THER" in init_state:
                resu = init_state.get("EVOL_THER")
                assert isinstance(resu, ThermalResult), resu
                para, value = _extract_param(init_state, resu)

                self.state.primal_curr = resu.getField(
                    "TEMP", para=para, value=value
                ).copyUsingDescription(nume_equa)

            if "CHAM_NO" in init_state:
                self.state.primal_curr = init_state.get("CHAM_NO").copyUsingDescription(nume_equa)

            if "DEPL" in init_state:
                self.state.primal_curr = init_state.get("DEPL").copyUsingDescription(
                    nume_equa, False
                )
                _msginit("DEPL")

            if "SIGM" in init_state:
                self.state.stress = _get_field_and_check_model(
                    state=init_state, name_field="SIGM", model=model
                )
                _msginit("SIEF_ELGA")

            if "VARI" in init_state:
                self.state.internVar = _get_field_and_check_model(
                    state=init_state, name_field="VARI", model=model
                )
                _msginit("VARI_ELGA")

            if "VITE" in init_state:
                self.state.current.dU = init_state.get("VITE").copyUsingDescription(nume_equa)
                _msginit("VITE")

            if "ACCE" in init_state:
                self.state.current.d2U = init_state.get("ACCE").copyUsingDescription(nume_equa)
                _msginit("ACCE")

            if "VALE" in init_state:
                if model.existsHHO():
                    self.state.primal_curr = HHO(self.problem).projectOnHHOSpace(init_state["VALE"])
                else:
                    self.state.primal_curr = self.state.createPrimal(
                        self.problem, value={"TEMP": init_state.get("VALE")}
                    )

        init_time = self.stepper.getInitial()
        self.computeExternalStateVariables(init_time)
        self.state.time_curr = init_time

        if init_state:
            if init_state.get("STAT") == "OUI":
                self.step_solver.initialize()
                args = {"valr": self.state.time_curr, "vali": self.stepper.splitting_level}
                logger.info(MessageLog.GetText("I", "MECANONLINE6_5", **args))
                self.step_solver.solve()
                if (
                    self.stepper.size() == 1
                    and self.stepper.getCurrent() == self.stepper.getPrevious()
                ):
                    self.stepper.completed()

            self.post_hooks()

        self.state.commit()

    def run(self):
        """Solve the problem.

        Returns:
            *misc*: result object.
        """
        try:
            self._setLoggingLevel(self.get_keyword("INFO", default=1))
            self.run_()
        finally:
            self._resetLoggingLevel()
            deleteTemporaryObjects()
        return self.result

    @profile
    def run_(self):
        """Solve the problem."""
        self.initialize()
        matr_update_step = self.get_keyword("NEWTON", "REAC_INCR", 1)

        # Solve nonlinear problem
        solv = self.step_solver
        state = self.state
        last_stored = False
        while not self.isFinished():
            state.time_curr = self.stepper.getCurrent()
            state.time_step = state.time_curr - state.time_prev
            if self.stepper.splitting_level <= 0:
                logger.info(MessageLog.GetText("I", "MECANONLINE6_7", valr=state.time_curr))
            else:
                args = dict(valr=state.time_curr, vali=self.stepper.splitting_level)
                logger.info(MessageLog.GetText("I", "MECANONLINE6_5", **args))

            self.computeExternalStateVariables(state.time_curr)
            solv.initialize()

            if matr_update_step == 0 or (self._step_idx + 1) % matr_update_step:
                solv.current_matrix = self.current_matrix
            else:
                solv.current_matrix = None

            if logger.getEffectiveLevel() <= DEBUG:
                state.debugPrint("<t-> ")
            state.stash()
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
                    self._storeState(state.getState(-1), ignore_policy=True)
                    raise
            else:
                if not self.stepper.check_event(state):
                    # + reset current_matrix to None (REAC_INCR)
                    state.revert()
                    continue
                self.post_hooks()
                state.commit()
                self.stepper.completed()
                self.current_matrix = solv.current_matrix
                self._step_idx += 1
                last_stored = self._storeState(state)
        # ensure that last step was stored
        if not last_stored:
            self._storeState(state, ignore_policy=True)

    def post_hooks(self):
        """Call post hooks"""
        # FIXME: todo
        # for hook in self.get_features(SOP.PostStepHook):
        #     hook(self)

    def computeExternalStateVariables(self, current_time):
        """Compute and set external variables in the physical state.

        Arguments:
            current_time (float): Current time value.
        """
        if self.problem.getMaterialField().hasExternalStateVariable():
            self.state.externVar = self.problem.getExternalStateVariables(current_time)

    def _setLoggingLevel(self, level):
        """Set logging level.

        Arguments:
            level (int): verbosity level (meaning INFO keyword).
        """
        info = {0: WARNING, 1: INFO, 2: DEBUG, 3: DEBUG, 4: DEBUG}
        if level is None:
            level = 1
        setFortranLoggingLevel(level)
        logger.setLevel(info[level])
        # Disable printing of python command
        if level < 3:
            ExecutionParameter().disable(Options.ShowSyntax)

    def _resetLoggingLevel(self):
        """Reset logging level."""
        level, show = self._verb
        resetFortranLoggingLevel()
        logger.setLevel(level)
        if show:
            ExecutionParameter().enable(Options.ShowSyntax)


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
