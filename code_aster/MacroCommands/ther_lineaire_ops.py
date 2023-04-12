# coding: utf-8

# Copyright (C) 1991 - 2023  EDF R&D                www.code-aster.org
#
# This file is part of Code_Aster.
#
# Code_Aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Code_Aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.

from libaster import deleteCachedObjects, resetFortranLoggingLevel, setFortranLoggingLevel

from ..Commands import CALC_CHAMP
from ..Objects import (
    HHO,
    AssemblyMatrixTemperatureReal,
    DiscreteComputation,
    FieldOnNodesReal,
    LinearSolver,
    ParallelThermalLoadFunction,
    ParallelThermalLoadReal,
    PhysicalProblem,
    ThermalDirichletBC,
    ThermalLoadFunction,
    ThermalLoadReal,
    ThermalResult,
)
from ..Solvers import PhysicalState, StorageManager, TimeStepper
from ..Utilities import SearchList, logger, print_stats, profile


def _checkArgs(args):

    if args.get("RESULTAT") is not None and args.get("reuse") is not None:
        assert args.get("RESULTAT") is args.get("reuse")

    if args.get("RESULTAT") is not None or args.get("reuse") is not None:
        assert args.get("ETAT_INIT") is not None
        assert args.get("ETAT_INIT").get("STAT") is None


def _hasExchangeFields(args):
    has_fields = False
    if args.get("EXCIT") is not None:
        for loadkws in args["EXCIT"]:
            load = loadkws["CHARGE"]
            if isinstance(load, (ThermalLoadFunction, ThermalLoadReal)):
                has_fields = (
                    has_fields
                    or load.hasLoadResult()
                    or load.hasLoadField("COEFH")
                    or load.hasLoadField("HECHP")
                )
    return has_fields


@profile
def _addLoads(phys_pb, args):
    """Add loads from cata inplace

    Arguments:
        phys_pb (PhysicalProblem): Physical problem to add load
        **args (dict): User's keywords.

    Returns:
        (PhysicalProblem): Modifel physical problem inplace
    """

    if args.get("EXCIT") is not None:
        for load in args["EXCIT"]:
            if isinstance(
                load["CHARGE"],
                (
                    ThermalLoadFunction,
                    ThermalLoadReal,
                    ThermalDirichletBC,
                    ParallelThermalLoadFunction,
                    ParallelThermalLoadReal,
                ),
            ):
                phys_pb.addLoadFromDict(load)
            else:
                raise RuntimeError("Unknown load")

    phys_pb.computeListOfLoads()

    return phys_pb


@profile
def _setupInitialField(phys_pb, args):

    logger.debug("<THER_LINEAIRE><ETAT_INIT>: Start")

    initial_field = None
    is_stat_init = False
    initial_state = args.get("ETAT_INIT")

    if args["TYPE_CALCUL"] == "STAT" or "STAT" in initial_state:
        logger.debug(
            "<THER_LINEAIRE><ETAT_INIT>: Stationnary Computation initialized with a null field"
        )
        initial_field = FieldOnNodesReal(phys_pb.getDOFNumbering())
        initial_field.setValues(0.0)
        is_stat_init = True
    elif "CHAM_NO" in initial_state:
        logger.debug(
            "<THER_LINEAIRE><ETAT_INIT>: Initialized with given field '%s'"
            % initial_state.get("CHAM_NO").getName()
        )
        initial_field = initial_state.get("CHAM_NO")
    elif "VALE" in initial_state:
        logger.debug(
            "<THER_LINEAIRE><ETAT_INIT>: Initialized with constant field with value %s"
            % initial_state["VALE"]
        )
        # For HHO, there is a projection to do.
        if phys_pb.getModel().existsHHO():
            initial_field = HHO(phys_pb).projectOnHHOSpace(initial_state["VALE"])
        else:
            initial_field = FieldOnNodesReal(phys_pb.getDOFNumbering())
            initial_field.setValues(initial_state["VALE"])
    elif "EVOL_THER" in initial_state:
        resu_ther = initial_state.get("EVOL_THER")
        para = resu_ther.getAccessParameters()

        index = initial_state.get("NUME_ORDRE") or para["NUME_ORDRE"][-1]

        if initial_state.get("INST") is not None:
            timelist = SearchList(
                para["INST"], initial_state["PRECISION"], initial_state["CRITERE"]
            )

            index = timelist.index(initial_state["INST"])

        initial_field = resu_ther.getField("TEMP", index).duplicate()
        logger.debug(
            "<THER_LINEAIRE><ETAT_INIT>: Initialized with field from '%s' at index '%s'"
            % (resu_ther.getName(), index)
        )
    else:
        assert False

    assert initial_field is not None
    logger.debug("<THER_LINEAIRE><ETAT_INIT>: Finish")
    return initial_field, is_stat_init


@profile
def _createTimeStepper(args):
    """Create time stepper from catalogues

    Arguments:
        **args (dict): User's keywords.

    Returns:
        (TimeStepper): a time stepper.
    """

    logger.debug("<THER_LINEAIRE><TIMESTEPPER>: Start")

    # Create result
    result = args.get("RESULTAT") or args.get("reuse")
    if result is None:
        last_prev_inst = None
        if args["TYPE_CALCUL"] == "STAT":
            first_index = 1
        else:
            first_index = 0
    else:
        para = result.getAccessParameters()
        inst_prev = para["INST"]
        last_prev_inst = inst_prev[-1]
        index_prev = para["NUME_ORDRE"]

    increment = args.get("INCREMENT")

    time_values = [0.0]
    if increment is not None:
        listInst = increment.get("LIST_INST")
        prec = increment.get("PRECISION")

        list_values = listInst.getValues()

        timelist_current = SearchList(list_values, prec)

        logger.debug("<THER_LINEAIRE><TIMESTEPPER>: list_values = %s" % list_values)

        nume_inst_init = increment.get("NUME_INST_INIT") or 0
        nume_inst_fin = increment.get("NUME_INST_FIN") or len(list_values)

        if increment.get("INST_INIT") is None:
            inst_init = last_prev_inst or list_values[0]
        else:
            inst_init = increment.get("INST_INIT")

        inst_fin = increment.get("INST_FIN") or list_values[-1]

        logger.debug("<THER_LINEAIRE><TIMESTEPPER>: nume_inst_init = %s" % nume_inst_init)
        logger.debug("<THER_LINEAIRE><TIMESTEPPER>: nume_inst_fin = %s" % nume_inst_fin)
        logger.debug("<THER_LINEAIRE><TIMESTEPPER>: inst_init = %s" % inst_init)
        logger.debug("<THER_LINEAIRE><TIMESTEPPER>: inst_fin = %s" % inst_fin)

        assert nume_inst_init >= 0
        assert nume_inst_fin <= len(list_values)
        assert nume_inst_init < nume_inst_fin

        assert inst_init >= list_values[0]
        assert inst_fin <= list_values[-1]
        assert inst_init <= inst_fin

        assert timelist_current.unique(inst_init)
        assert timelist_current.unique(inst_fin)

        time_values = [
            time
            for time in list_values[nume_inst_init : nume_inst_fin + 1]
            if ((inst_init - prec) <= time <= (inst_fin + prec))
        ]

        assert len(time_values) > 0

        # first index to use
        if last_prev_inst is not None:
            assert result is not None
            timelist_prev = SearchList(inst_prev, prec)
            idx = timelist_prev.index(last_prev_inst)
            first_index = index_prev[idx]

    logger.debug("<THER_LINEAIRE><TIMESTEPPER>: first_index = %s" % first_index)
    logger.debug("<THER_LINEAIRE><TIMESTEPPER>: time_values = %s" % time_values)
    return first_index, TimeStepper(time_values)


@profile
def _computeMatrix(disr_comp, matrix, is_evol, time_value, time_delta, time_theta):
    """Compute and assemble the thermal matrix

    Arguments:
        disr_comp (DiscreteComputation): to compute discrete quantities
        matrix (AssemblyMatrixTemperatureReal): matrix to compute and assemble inplace
        time_value (float): Current time
        time_delta (float): Time increment
        time_theta (float): Theta parameter for time scheme

    Returns:
        AssemblyMatrixTemperatureReal: matrix computed and assembled
    """
    logger.debug("<THER_LINEAIRE><MATRIX>: Start")

    phys_pb = disr_comp.getPhysicalProblem()

    logger.debug("<THER_LINEAIRE><MATRIX>: Linear Conductivity")
    matr_elem_rigi = disr_comp.getLinearStiffnessMatrix(time_value, with_dual=False)
    matrix.addElementaryMatrix(matr_elem_rigi, time_theta)

    matr_elem_exch = disr_comp.getExchangeThermalMatrix(time_value)
    matrix.addElementaryMatrix(matr_elem_exch, time_theta)

    if phys_pb.getDOFNumbering().useLagrangeMultipliers():
        logger.debug("<THER_LINEAIRE><MATRIX>: Dual Conductivity")
        matr_elem_dual = disr_comp.getDualLinearConductivityMatrix()
        matrix.addElementaryMatrix(matr_elem_dual)

    if is_evol:
        logger.debug("<THER_LINEAIRE><MATRIX>: Linear Capacity")

        matr_elem_capa = disr_comp.getLinearCapacityMatrix(time_value)
        matrix.addElementaryMatrix(matr_elem_capa, 1.0 / time_delta)

    matrix.assemble(True)

    logger.debug("<THER_LINEAIRE><MATRIX>: Finish")

    return matrix


@profile
def _computeRhs(disr_comp, is_evol, time_value, time_delta, time_theta, previousPrimalField):
    """Compute and assemble the right hand side

    Arguments:
        disr_comp (DiscreteComputation): to compute discrete quantities
        time_value (float): Current time
        time_delta (float): Time increment
        time_theta (float): Theta parameter for integration
        previousPrimalField (fieldOnNodesReal): solution field at previous time

     Returns:
         FieldOnNodesReal: vector of load
    """

    logger.debug("<THER_LINEAIRE><RHS>: Start")

    # compute imposed temperature with Lagrange
    rhs = disr_comp.getImposedDualBC(time_value, time_delta, time_theta)
    logger.debug("<THER_LINEAIRE><RHS>: Nodal BC")
    # compute neumann forces
    rhs += disr_comp.getNeumannForces(time_value, time_delta, time_theta, previousPrimalField)
    logger.debug("<THER_LINEAIRE><RHS>: Neumann BC")

    if is_evol:
        rhs += disr_comp.getTransientThermalForces(
            time_value, time_delta, time_theta, previousPrimalField
        )
        logger.debug("<THER_LINEAIRE><RHS>: Transient Load BC")

    logger.debug("<THER_LINEAIRE><RHS>: Finish")
    return rhs


def ther_lineaire_ops(self, **args):
    """Execute the command THER_LINEAIRE.

    Arguments:
        **args (dict): User's keywords.

    Returns:
        ThermalResult: result for linear thermal problem
    """

    logger.debug("<THER_LINEAIRE>: Initialization")
    logger.debug("<THER_LINEAIRE>: Args : %s" % args)

    _checkArgs(args)

    verbosity = args.get("INFO") or 1
    setFortranLoggingLevel(verbosity)

    is_evol = args["TYPE_CALCUL"] == "TRAN"

    # Create result
    result = args.get("RESULTAT") or args.get("reuse")
    if result is None:
        result = ThermalResult()
        result.setTitle(args.get("TITRE") or "RESU_THER_LINEAIRE")

    # Create physical problem
    model = args["MODELE"]
    phys_pb = PhysicalProblem(model, args["CHAM_MATER"], args.get("CARA_ELEM"))
    logger.debug("<THER_LINEAIRE>: Physical Problem created")

    # for HHO model
    hho = HHO(phys_pb)

    # Add loads
    phys_pb = _addLoads(phys_pb, args)
    has_exchange_fields = _hasExchangeFields(args)
    logger.debug("<THER_LINEAIRE>: Loads added")

    # Compute numbering
    phys_pb.computeDOFNumbering()
    logger.debug("<THER_LINEAIRE>: DOFNumbering computed")

    # Setup ETAT_INIT
    initial_field, is_stat_init = _setupInitialField(phys_pb, args)

    # Create time stepper
    first_index, timeStepper = _createTimeStepper(args)

    # do not erase initial state if there are equals
    save_initial_state = True
    if (
        is_evol
        and "reuse" in args
        and "EVOL_THER" in args["ETAT_INIT"]
        and args["reuse"] == args["ETAT_INIT"]["EVOL_THER"]
    ):
        first_index += 1
        save_initial_state = False

    # Create linear solver
    linear_solver = LinearSolver.factory("THER_LINEAIRE", args["SOLVEUR"])
    if (model.getMesh().isParallel()) and (not linear_solver.supportParallelMesh()):
        raise RuntimeError("ParallelMesh not allowed with this linear solver")

    linear_solver.build()

    # Create storage manager
    storage_manager = StorageManager(result, args["ARCHIVAGE"])
    storage_manager.setInitialIndex(first_index)

    # Define main objects
    phys_state = PhysicalState()
    disc_comp = DiscreteComputation(phys_pb)

    # we define the matrix before to have an unique name
    # because of a bug with LDLT_SP
    matrix = AssemblyMatrixTemperatureReal(phys_pb)
    # the matrix depends on times or external variables
    is_const = phys_pb.getCodedMaterial().constant()

    # Detect external state variables
    hasExternalStateVariable = phys_pb.getMaterialField().hasExternalStateVariable()

    # Compute reference value vector for external state variables
    if phys_pb.getMaterialField().hasExternalStateVariableWithReference():
        phys_pb.computeReferenceExternalStateVariables()

    # Run computation
    logger.debug("<THER_LINEAIRE>: Start computation")

    phys_state.primal = initial_field
    time_delta_prev = timeStepper.null_increment

    # Compute initial state
    if is_evol:
        phys_state.time = timeStepper.getNext()
        time_theta = 1.0
        time_delta = timeStepper.null_increment
        if is_stat_init:
            matrix = _computeMatrix(
                disc_comp, matrix, False, phys_state.time, time_delta, time_theta
            )
            profile(linear_solver.factorize)(matrix)

            rhs = _computeRhs(
                disc_comp, False, phys_state.time, time_delta, time_theta, phys_state.primal
            )

            # solve linear system
            diriBCs = profile(disc_comp.getDirichletBC)(phys_state.time)
            phys_state.primal = profile(linear_solver.solve)(rhs, diriBCs)

        if save_initial_state:
            storage_manager.storeState(
                phys_state.time, phys_pb, phys_state, param={"PARM_THETA": time_theta}
            )
            if model.existsHHO():
                hho_field = hho.projectOnLagrangeSpace(phys_state.primal)
                storage_manager.storeField(hho_field, "HHO_TEMP", phys_state.time)

            storage_manager.completed()
        timeStepper.completed()

    # Loop on time step
    while not timeStepper.hasFinished():
        phys_state.time = timeStepper.getNext()

        if is_evol:
            time_theta = args.get("PARM_THETA")
            time_delta = timeStepper.getIncrement()
        else:
            time_theta = 1.0
            time_delta = timeStepper.null_increment

        logger.debug("<THER_LINEAIRE>:     IS_EVOL %s" % is_evol)
        logger.debug("<THER_LINEAIRE>:     IS_CONST = %s" % is_const)
        logger.debug("<THER_LINEAIRE>:     HAS_EXT_STATE_VAR = %s" % hasExternalStateVariable)
        logger.debug("<THER_LINEAIRE>:     CURRENT TIME %s" % phys_state.time)
        logger.debug("<THER_LINEAIRE>:     TIME_VALUE %s" % phys_state.time)
        logger.debug("<THER_LINEAIRE>:     TIME_DELTA %s" % time_delta)
        logger.debug("<THER_LINEAIRE>:     TIME_THETA %s" % time_theta)

        if (
            not is_const
            or (is_const and time_delta is timeStepper.null_increment)
            or (is_const and has_exchange_fields)
            or (is_const and abs(time_delta - time_delta_prev) > 1.0e-12)
            or (is_const and hasExternalStateVariable)
        ):

            matrix = _computeMatrix(
                disc_comp, matrix, is_evol, phys_state.time, time_delta, time_theta
            )
            profile(linear_solver.factorize)(matrix)

        rhs = _computeRhs(
            disc_comp, is_evol, phys_state.time, time_delta, time_theta, phys_state.primal
        )

        # solve linear system
        diriBCs = profile(disc_comp.getDirichletBC)(phys_state.time)
        phys_state.primal = profile(linear_solver.solve)(rhs, diriBCs)

        if storage_manager.hasToBeStored(phys_state.time):
            storage_manager.storeState(
                phys_state.time, phys_pb, phys_state, param={"PARM_THETA": time_theta}
            )
            if model.existsHHO():
                hho_field = hho.projectOnLagrangeSpace(phys_state.primal)
                storage_manager.storeField(hho_field, "HHO_TEMP", phys_state.time)

            storage_manager.completed()

        timeStepper.completed()
        time_delta_prev = time_delta

    logger.debug("<THER_LINEAIRE>: Finish computation")
    # delete factorized matrix - free memory
    linear_solver.deleteFactorizedMatrix()

    # store field
    result = storage_manager.getResult()

    if model.isXfem():
        result = CALC_CHAMP(RESULTAT=result, reuse=result, THERMIQUE="TEMP_ELGA")

    if verbosity > 1:
        print_stats()
    resetFortranLoggingLevel()

    # cleaning
    deleteCachedObjects()

    return result
