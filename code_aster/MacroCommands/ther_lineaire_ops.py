# coding: utf-8

# Copyright (C) 1991 - 2022  EDF R&D                www.code-aster.org
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

from libaster import deleteCachedObjects, setFortranLoggingLevel, resetFortranLoggingLevel

from ..Commands import CALC_CHAMP
from ..Objects import (
    AssemblyMatrixTemperatureReal,
    DiscreteComputation,
    ThermalResult,
    LinearSolver,
    ThermalDirichletBC,
    ThermalLoadFunction,
    ThermalLoadReal,
    ParallelThermalLoadFunction,
    ParallelThermalLoadReal,
    PhysicalProblem,
    FieldOnNodesReal,
    PostProcessing
)
from ..Utilities import logger, print_stats, profile
from .NonLinearSolver import PhysicalState, StorageManager, TimeStepper


def get_index(inst, linst, prec=1.E-6, crit="RELATIF"):

    assert crit in ("RELATIF", "ABSOLU")

    if crit == "RELATIF":
        min_v = inst*(1-prec)
        max_v = inst*(1+prec)
    else:
        min_v = inst-prec
        max_v = inst+prec

    min_v, max_v = sorted((min_v, max_v))
    test = tuple((i >= min_v and i <= max_v) for i in linst)
    return tuple(i for i, v in enumerate(test) if v is True)


def exists_unique(*args, **kwargs):
    return len(get_index(*args, **kwargs)) == 1


def get_unique_index(*args, **kwargs):
    assert(exists_unique(*args, **kwargs))
    return get_index(*args, **kwargs)[0]


def _checkArgs(args):

    if (args.get("RESULTAT") is not None and args.get("reuse") is not None):
        assert (args.get("RESULTAT") is args.get("reuse"))

    if args.get("RESULTAT") is not None:
        assert(args.get("ETAT_INIT") is not None)
        assert(args.get("ETAT_INIT").get("STATIONNAIRE") is None)


def _checkIfStat(args):

    if args.get("ETAT_INIT") is None:
        return True
    else:
        if "STATIONNAIRE" in args.get("ETAT_INIT"):
            return True
    return False


def _hasExchangeFields(args):
    has_fields = False
    if args.get("EXCIT") is not None:
        for loadkws in args["EXCIT"]:
            load = loadkws["CHARGE"]
            if isinstance(load, (ThermalLoadFunction, ThermalLoadReal)):
                has_fields = has_fields or load.hasLoadResult() or load.hasLoadField(
                    "COEFH") or load.hasLoadField("HECHP")
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
def _setupInitialField(is_stat, phys_pb, args):

    logger.debug("<THER_LINEAIRE><ETAT_INIT>: Start")
    initial_field = None

    if is_stat:
        logger.debug(
            "<THER_LINEAIRE><ETAT_INIT>: Stationnary Computation initialized with a null field")
        initial_field = FieldOnNodesReal(phys_pb.getDOFNumbering())
        initial_field.setValues(0.0)
    else:

        ls_kws = args.get("ETAT_INIT")
        assert(ls_kws is not None)
        kws = ls_kws[0]

        if kws.get("CHAM_NO") is not None:
            logger.debug("<THER_LINEAIRE><ETAT_INIT>: Initialized with given field '%s'" % kws.get(
                "CHAM_NO").getName())
            initial_field = kws.get("CHAM_NO")

        elif kws.get("EVOL_THER") is not None:
            resu_ther = kws.get("EVOL_THER")
            para = resu_ther.getAccessParameters()

            irank = kws.get("NUME_ORDRE") or para["NUME_ORDRE"][-1]

            if kws.get("INST") is not None:
                irank = get_unique_index(kws.get("INST"), para["INST"],
                                         kws.get("PRECISION"), kws.get("CRITERE"))

            initial_field = resu_ther.getFieldOnNodesReal("TEMP", irank)
            logger.debug("<THER_LINEAIRE><ETAT_INIT>: Initialized with field from '%s' at rank '%s'" % (resu_ther.getName(),
                                                                                                        irank))

        elif kws.get("VALE") is not None:
            # For HHO, there is a projection to do.
            if phys_pb.getModel().existsHHO():
                assert abs(kws.get("VALE")) == 0.0

            initial_field = FieldOnNodesReal(phys_pb.getDOFNumbering())
            initial_field.setValues(kws.get("VALE"))
            logger.debug(
                "<THER_LINEAIRE><ETAT_INIT>: Initialized with constant field with value %s" % kws.get("VALE"))
        else:
            assert(False)

    assert(initial_field is not None)
    logger.debug("<THER_LINEAIRE><ETAT_INIT>: Finish")
    return initial_field


@profile
def _setupArchivage(timestepper_values, args):
    logger.debug("<THER_LINEAIRE><ARCHIVAGE>: Start")

    arch_args = args["ARCHIVAGE"]
    arch_prec = arch_args.get("PRECISION")
    arch_crit = arch_args.get("CRITERE")
    if "LIST_INST" in arch_args:
        arch_times = arch_args["LIST_INST"].getValues()
    else:
        arch_times = timestepper_values

    for t in arch_times:
        assert(exists_unique(t, timestepper_values,
               prec=arch_prec, crit=arch_crit))

    logger.debug("<THER_LINEAIRE><ARCHIVAGE>: arch_times %s " % arch_times)
    logger.debug("<THER_LINEAIRE><ARCHIVAGE>: Finish")
    return arch_times


@profile
def _createTimeStepper(is_stat, args):
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
        first_rank = 0
    else:
        para = result.getAccessParameters()
        inst_prev = para["INST"]
        last_prev_inst = inst_prev[-1]
        rank_prev = para["NUME_ORDRE"]

    timeValues = None
    increment = args.get("INCREMENT")

    if increment is None:
        time_values = [0.0]
    else:
        listInst = increment.get("LIST_INST")
        prec = increment.get("PRECISION")

        if listInst is not None:
            list_values = listInst.getValues()
            logger.debug(
                "<THER_LINEAIRE><TIMESTEPPER>: list_values = %s" % list_values)

            nume_inst_init = increment.get("NUME_INST_INIT") or 0
            nume_inst_fin = increment.get("NUME_INST_FIN") or len(list_values)

            if increment.get("INST_INIT") is None:
                inst_init = last_prev_inst or list_values[0]
            else:
                inst_init = increment.get("INST_INIT")

            inst_fin = increment.get("INST_FIN") or list_values[-1]

            logger.debug(
                "<THER_LINEAIRE><TIMESTEPPER>: nume_inst_init = %s" % nume_inst_init)
            logger.debug(
                "<THER_LINEAIRE><TIMESTEPPER>: nume_inst_fin = %s" % nume_inst_fin)
            logger.debug(
                "<THER_LINEAIRE><TIMESTEPPER>: inst_init = %s" % inst_init)
            logger.debug(
                "<THER_LINEAIRE><TIMESTEPPER>: inst_fin = %s" % inst_fin)

            assert nume_inst_init >= 0
            assert nume_inst_fin <= len(list_values)
            assert nume_inst_init < nume_inst_fin

            assert inst_init >= list_values[0]
            assert inst_fin <= list_values[-1]
            assert inst_init <= inst_fin

            assert exists_unique(inst_init, list_values, prec)
            assert exists_unique(inst_fin, list_values, prec)

            time_values = [time for time in list_values[nume_inst_init:nume_inst_fin+1]
                           if ((inst_init - prec) <= time <= (inst_fin + prec))]

    assert len(time_values) > 0
    is_evol = len(time_values) > 1
    if not is_evol:
        assert(is_stat)

    # first rank to use
    if last_prev_inst is not None:
        assert(result is not None)
        idx = get_unique_index(last_prev_inst, inst_prev, prec)
        first_rank = rank_prev[idx] + 1

    logger.debug("<THER_LINEAIRE><TIMESTEPPER>: first_rank = %s" % first_rank)
    logger.debug("<THER_LINEAIRE><TIMESTEPPER>: time_values = %s" %
                 time_values)
    return is_evol, first_rank, TimeStepper(time_values)


@profile
def _computeMatrix(disr_comp, matrix,
                   is_stat, time_value, time_delta, time_theta):
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

    if not is_stat:
        logger.debug("<THER_LINEAIRE><MATRIX>: Linear Capacity")

        matr_elem_capa = disr_comp.getLinearCapacityMatrix(time_value)
        matrix.addElementaryMatrix(matr_elem_capa, 1.0/time_delta)

    matrix.assemble(True)

    logger.debug("<THER_LINEAIRE><MATRIX>: Finish")

    return matrix


@profile
def _computeRhs(disr_comp,
                is_evol, time_value, time_delta, time_theta,
                previousPrimalField):
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
    rhs += disr_comp.getNeumannForces(time_value, time_delta,
                             time_theta, previousPrimalField)
    logger.debug("<THER_LINEAIRE><RHS>: Neumann BC")

    if is_evol:
        rhs += disr_comp.getTransientThermalForces(time_value, time_delta, time_theta,
                                              previousPrimalField)
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

    is_stat = _checkIfStat(args)

    # Create result
    result = args.get("RESULTAT") or args.get("reuse")
    if result is None:
        result = ThermalResult()
        result.setTitle(args.get("TITRE") or "RESU_THER_LINEAIRE")

    # Create physical problem
    model = args["MODELE"]
    phys_pb = PhysicalProblem(model, args["CHAM_MATER"], args.get("CARA_ELEM"))
    logger.debug("<THER_LINEAIRE>: Physical Problem created")

    postpro = PostProcessing(phys_pb)

    # Add loads
    phys_pb = _addLoads(phys_pb, args)
    has_exchange_fields = _hasExchangeFields(args)
    logger.debug("<THER_LINEAIRE>: Loads added")

    # Compute numbering
    phys_pb.computeDOFNumbering()
    logger.debug("<THER_LINEAIRE>: DOFNumbering computed")

    # Setup ETAT_INIT
    initial_field = _setupInitialField(is_stat, phys_pb, args)

    # Create time stepper
    is_evol, rank, timeStepper = _createTimeStepper(is_stat, args)

    # Archivage
    arch_args = args["ARCHIVAGE"]
    arch_prec = arch_args.get("PRECISION")
    arch_crit = arch_args.get("CRITERE")
    arch_times = _setupArchivage(timeStepper.times, args)

    # Create linear solver
    linear_solver = LinearSolver.factory("THER_LINEAIRE", args["SOLVEUR"])
    if (model.getMesh().isParallel()) and (not linear_solver.supportParallelMesh()):
        raise RuntimeError("ParallelMesh not allowed with this linear solver")

    linear_solver.build()

    # Create storage manager
    storage_manager = StorageManager(result)

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

    is_first = True
    phys_state.primal = initial_field
    time_delta_prev = timeStepper.null_increment

    while not timeStepper.hasFinished():
        phys_state.time = timeStepper.getNext()

        if is_stat:
            time_theta = 1.0
            time_delta = timeStepper.null_increment
        else:
            time_theta = args.get("PARM_THETA")
            time_delta = timeStepper.getIncrement()

        logger.debug("<THER_LINEAIRE>:     IS_EVOL %s" % is_evol)
        logger.debug("<THER_LINEAIRE>:     IS_STAT %s" % is_stat)
        logger.debug("<THER_LINEAIRE>:     IS_CONST = %s" % is_const)
        logger.debug("<THER_LINEAIRE>:     HAS_EXT_STATE_VAR = %s" %
                     hasExternalStateVariable)
        logger.debug("<THER_LINEAIRE>:     CURRENT TIME %s" % phys_state.time)
        logger.debug("<THER_LINEAIRE>:     TIME_VALUE %s" % phys_state.time)
        logger.debug("<THER_LINEAIRE>:     TIME_DELTA %s" % time_delta)
        logger.debug("<THER_LINEAIRE>:     TIME_THETA %s" % time_theta)

        if not (is_first and not is_stat):
            if (not is_const
                or (is_const and time_delta is timeStepper.null_increment)
                or (is_const and has_exchange_fields)
                or (is_const and abs(time_delta - time_delta_prev) > 1.e-12)
                    or (is_const and hasExternalStateVariable)):
                matrix = _computeMatrix(disc_comp, matrix,
                                        is_stat, phys_state.time, time_delta, time_theta)
                profile(linear_solver.factorize)(matrix)

            rhs = _computeRhs(disc_comp,
                              is_evol, phys_state.time, time_delta, time_theta,
                              phys_state.primal)

            # solve linear system
            diriBCs = profile(disc_comp.getDirichletBC)(phys_state.time)
            phys_state.primal = profile(linear_solver.solve)(rhs, diriBCs)

        if (rank == 0) or not is_first:
            if exists_unique(phys_state.time, arch_times, arch_prec, arch_crit):
                storage_manager.storeState(rank, phys_state.time, phys_pb, phys_state,
                                           theta=time_theta)
                if model.existsHHO():
                    proj_hho = postpro.projectHHO(phys_state.primal)
                    storage_manager.storeField(proj_hho, "HHO_TEMP", rank)
                rank += 1

        timeStepper.completed()
        is_first = False
        is_stat = False
        time_delta_prev = time_delta

    logger.debug("<THER_LINEAIRE>: Finish computation")
    # delete factorized matrix - free memory
    linear_solver.deleteFactorizedMatrix()

    # store field
    result = storage_manager.getResult()

    if model.isXfem():
        result = CALC_CHAMP(RESULTAT=result, reuse=result,
                            THERMIQUE="TEMP_ELGA")

    if verbosity > 1:
        print_stats()
    resetFortranLoggingLevel()

    # cleaning
    deleteCachedObjects()

    return result
