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

from ...Cata.Language.SyntaxObjects import _F
from ...Messages import UTMESS, MessageLog
from ...Objects import HHO, LinearSolver, NonLinearResult, ThermalResult
from ...Supervis import AsterError, ConvergenceError, IntegrationError, SolverError
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
from ..Basics import ContextMixin, PhysicalState, SolverFeature
from ..Basics import ProblemType as PBT
from ..Basics import SolverOptions as SOP
from ..OperatorsManager import BaseOperatorsManager
from ..StepSolvers import BaseStepSolver
from .convergence_manager import ConvergenceManager
from .incremental_solver import IncrementalSolver
from .line_search import LineSearch
from .newton_solver import NewtonSolver
from .raspen_solver import RASPENSolver
from .snes_solver import SNESSolver
from .storage_manager import StorageManager
from .time_stepper import TimeStepper


class ProblemSolver(SolverFeature, ContextMixin):
    """Solver for linear and non linear problem.

    Arguments:
        main (*NonLinearFeature*): Main object.
        result (*misc*): The result object.
    """

    required_features = [
        SOP.Keywords,
        SOP.LinearSolver,
        SOP.PhysicalProblem,
        SOP.PhysicalState,
        SOP.StepSolver,
        SOP.Storage,
        SOP.TimeStepper,
    ]
    optional_features = [
        SOP.Contact,
        SOP.ConvergenceCriteria,
        SOP.ConvergenceManager,
        SOP.IncrementalSolver,
        SOP.LineSearch,
        SOP.PostStepHook,
    ]

    _stepper = _store = _step_solver = None
    _verb = None
    # FIXME: add _ prefix?
    step_rank = current_matrix = None
    __setattr__ = no_new_attributes(object.__setattr__)

    @classmethod
    def builder(cls, _):
        """Disabled for this object."""
        raise NotImplementedError("ProblemSolver builds the Context itself!")

    def __init__(self, phys_pb, problem_type, result, keywords) -> None:
        super().__init__()
        self.problem = phys_pb
        self.problem_type = problem_type
        self.result = result
        self.keywords = keywords
        self.state = PhysicalState(problem_type, size=1)
        self.oper = BaseOperatorsManager.factory(self.context)
        self.step_rank = None
        self.current_matrix = None
        self._verb = logger.getEffectiveLevel(), ExecutionParameter().option & Options.ShowSyntax

    # convenient shortcuts properties to init and access subobjects
    @property
    def stepper(self):
        """:py:class:`~.time_stepper.TimeStepper`: object to be used."""
        if not self._stepper:
            logger.debug("+++ init Stepper")
            self._stepper = TimeStepper.from_keywords(**self.keywords["INCREMENT"])
            self.use(self._stepper)
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
            self.use(self._store)
        return self._store

    @property
    def step_solver(self):
        if self._step_solver:
            return self._step_solver
        logger.debug("+++ init StepSolver")
        self._step_solver = solv = BaseStepSolver.factory(self.context)
        solv.setParameters(self.keywords)
        # # FIXME: todo
        # for feat, required in solv.undefined():
        #     solv.use(self._getF(feat, required))
        # self.use(solv)
        solv.setup()
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
            self.step_rank,
            state.time_curr,
            self.problem,
            state,
            is_final_time=self.isFinished(),
            ignore_policy=ignore_policy,
        )

    # @profile
    def initialize(self):
        """Initialize run"""
        phys_pb = self.problem
        kwds = self.keywords
        # essential to be called enough soon (may change the size of VARI field)
        if self._get("ETAT_INIT"):
            phys_pb.computeBehaviourProperty(kwds["COMPORTEMENT"], "OUI", 2)
        else:
            phys_pb.computeBehaviourProperty(kwds["COMPORTEMENT"], "NON", 2)
        phys_pb.computeListOfLoads()
        phys_pb.computeDOFNumbering()
        if phys_pb.getMaterialField().hasExternalStateVariableForLoad():
            phys_pb.computeReferenceExternalStateVariables()
        self.step_rank = 0
        self.setInitialState()
        self._storeState(self.state)
        # register observers
        for source in self.get_childs(SOP.IncrementalSolver | SOP.EventSource):
            source.add_observer(self.stepper)
        for source in self.get_childs(SOP.ConvergenceCriteria | SOP.EventSource):
            source.add_observer(self.stepper)

    # FIXME: mixin by problem_type / factory
    # @profile
    def setInitialState(self):
        """Initialize the physical state."""
        self.state.zeroInitialState(self.problem)
        init_state = self._get("ETAT_INIT")
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

    # @profile
    def run(self):
        """Solve the problem."""
        self.initialize()
        matr_update_step = self._get("NEWTON", "REAC_INCR", 1)

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

            if matr_update_step == 0 or (self.step_rank + 1) % matr_update_step:
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
                self.step_rank += 1
                last_stored = self._storeState(state)
        # ensure that last step was stored
        if not last_stored:
            self._storeState(state, ignore_policy=True)

    def post_hooks(self):
        """Call post hooks"""
        # FIXME: todo
        for hook in self.get_features(SOP.PostStepHook):
            hook(self)

    def computeExternalStateVariables(self, current_time):
        """Compute and set external variables in the physical state.

        Arguments:
            current_time (float): Current time value.
        """
        if self.problem.getMaterialField().hasExternalStateVariable():
            self.state.externVar = self.problem.getExternalStateVariables(current_time)

    # FIXME: self.keywords.get(...) or self.get_keyword(...) / make 'keywords' an object
    def _get(self, keyword, parameter=None, default=None):
        """ "Return a keyword value"""
        args = self.keywords
        if parameter is not None:
            if args.get(keyword) is None:
                return default
            # FIXME: to be checked: use _F() or not, take [0] or not...
            return _F(args[keyword])[0].get(parameter, default)

        return args.get(keyword, default)

    # ---- OLD ----
    def setKeywords(self, **args):
        """Set parameters from user keywords.

        Arguments:
            args (dict) : user keywords.
        """
        args = _F(args)
        self.use(args, SOP.Keywords)

    def _initialize(self):
        """Initialize default objects when required."""
        # args = self.get_feature(SOP.Keywords, optional=True)
        # self._setLoggingLevel(args.get("INFO", 1))
        # # required to build other default objects
        # self.use(self.problem)
        # self.use(self.state)
        # self.use(args, SOP.Keywords)
        # self.use(self._get_stepper())
        # self.use(self._get_storage())
        # self.use(self._get_step_solver())
        # self.use(self.get_features(SOP.PostStepHook), provide=SOP.PostStepHook)
        # self.check_features()

    def run_ops(self):
        """Solve the problem.

        Returns:
            *misc*: result object.
        """
        try:
            self._initialize()
            self.run()
        finally:
            self._resetLoggingLevel()
            deleteTemporaryObjects()
        return self.result

    def _setLoggingLevel(self, level):
        """Set logging level.

        Arguments:
            level (int): verbosity level (meaning INFO keyword).
        """
        info = {0: WARNING, 1: INFO, 2: DEBUG, 3: DEBUG, 4: DEBUG}
        if level is None:
            level = 1
        level = 2
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

    # methods that build and return objects (allowing dependencies)
    def _get_storage(self):
        logger.debug("+++ get StorageManager")
        store = self.get_feature(SOP.Storage, optional=True)
        if not store:
            args = self.get_feature(SOP.Keywords)
            reuse = args.get("REUSE")
            store = StorageManager(self.result, args.get("ARCHIVAGE"), reused=reuse is self.result)
            if reuse:
                init_state = args.get("ETAT_INIT")
                assert init_state
                if "EVOL_NOLI" in init_state:
                    # Pour l'instant, on se restreint au cas où la sd passée
                    # par reuse est la même que celle passée dans ETAT_INIT
                    assert init_state["EVOL_NOLI"] is reuse
                init_index = None
                stepper = self.get_feature(SOP.TimeStepper, optional=True)
                if stepper:
                    init_index = reuse.getIndexFromParameter(
                        "INST", stepper.getInitial(), "RELATIF", stepper._eps
                    )
                else:
                    init_index = reuse.getLastIndex()
                store.setFirstStorageIndex(init_index + 1)
        self.use(store)
        return store

    def _get_stepper(self):
        logger.debug("+++ get Stepper")
        stepper = self.get_feature(SOP.TimeStepper, optional=True)
        if not stepper:
            args = self.get_feature(SOP.Keywords)
            stepper = TimeStepper.from_keywords(**args["INCREMENT"])
        self.use(stepper)
        return stepper

    def _get_linear_solver(self):
        logger.debug("+++ get LinearSolver")
        if not self.has_feature(SOP.LinearSolver):
            args = self.get_feature(SOP.Keywords)
            self.use(LinearSolver.factory("STAT_NON_LINE", args["SOLVEUR"]))
        return self.get_feature(SOP.LinearSolver)

    def _get_line_search(self):
        logger.debug("+++ get LineSearch")
        if not self.has_feature(SOP.LineSearch):
            args = self.get_feature(SOP.Keywords)
            line = LineSearch(args.get("RECH_LINEAIRE"))
        for feat, required in line.undefined():
            line.use(self._getF(feat, required))
        line.setup()
        self.use(line)
        return self.get_feature(SOP.LineSearch)

    def _get_conv_manager(self):
        logger.debug("+++ get ConvergenceManager")
        converg = self.get_feature(SOP.ConvergenceManager, optional=True)
        if not converg:
            args = self.get_feature(SOP.Keywords)
            converg = ConvergenceManager()
            for crit in ("RESI_GLOB_RELA", "RESI_GLOB_MAXI", "ITER_GLOB_MAXI"):
                value = args["CONVERGENCE"].get(crit)
                if value is not None:
                    converg.setdefault(crit, value)
            if args.get("CONTACT"):
                converg.setdefault("RESI_GEOM", args["CONTACT"].get("RESI_GEOM"))
        for feat, required in converg.undefined():
            converg.use(self._getF(feat, required))
        self.use(converg)
        return converg

    def _get_incremental_solver(self):
        logger.debug("+++ get IncrementalSolver")
        incr_solver = self.get_feature(SOP.IncrementalSolver, optional=True)
        if not incr_solver:
            incr_solver = IncrementalSolver()
        for feat, required in incr_solver.undefined():
            incr_solver.use(self._getF(feat, required))
        self.use(incr_solver)
        return incr_solver

    def _get_step_conv_solver(self):
        logger.debug("+++ get ConvergenceCriteria")
        step_conv_solv = self.get_feature(SOP.ConvergenceCriteria, optional=True)
        args = self.get_feature(SOP.Keywords)
        use_local_solver = False
        if not step_conv_solv:
            method = args.get("METHODE", "NEWTON")
            if method == "NEWTON":
                step_conv_solv = NewtonSolver()
            elif method == "SNES":
                step_conv_solv = SNESSolver()
            elif method == "RASPEN":
                step_conv_solv = RASPENSolver()
                local_solver = SNESSolver(local=True)
                step_conv_solv.local_solver = local_solver
                use_local_solver = True
            else:
                raise AsterError(f"Unkwown method {method}")
        # FIXME replace by 'use(*, Keywords)'
        if not step_conv_solv.param:
            step_conv_solv.setParameters(args)
            if use_local_solver:
                local_solver.setParameters(args)
        for feat, required in step_conv_solv.undefined():
            step_conv_solv.use(self._getF(feat, required))
            if use_local_solver:
                local_solver.use(self._getF(feat, required))
        self.use(step_conv_solv)
        if use_local_solver:
            self.use(local_solver)
        return step_conv_solv

    def _get_step_solver(self):
        logger.debug("+++ get StepSolver")
        step_solver = self.get_feature(SOP.StepSolver, optional=True)
        if not step_solver:
            step_solver = BaseStepSolver.create(self.get_feature(SOP.Keywords))
        # FIXME replace by 'use(*, Keywords)'
        if not step_solver.param:
            step_solver.setParameters(self.get_feature(SOP.Keywords))
        for feat, required in step_solver.undefined():
            step_solver.use(self._getF(feat, required))
        self.use(step_solver)
        step_solver.setup()
        return step_solver

    def _getF(self, option, required):
        if option & SOP.PhysicalProblem:
            return self.problem
        if option & SOP.PhysicalState:
            return self.state
        if option & SOP.Storage:
            return self.store
        if option & SOP.LinearSolver:
            return self._get_linear_solver()
        if option & SOP.LineSearch:
            return self._get_line_search()
        if option & SOP.Contact:
            return self.contact
        if option & SOP.ConvergenceManager:
            return self._get_conv_manager()
        if option & SOP.IncrementalSolver:
            return self._get_incremental_solver()
        if option & SOP.ConvergenceCriteria:
            return self._get_step_conv_solver()
        if option & SOP.StepSolver:
            return self.step_solver
        if option & SOP.TimeStepper:
            return self.stepper
        if required:
            raise NotImplementedError(f"unsupported feature id: {option}")


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
