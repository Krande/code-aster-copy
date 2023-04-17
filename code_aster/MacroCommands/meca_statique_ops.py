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

from libaster import deleteTemporaryObjects, resetFortranLoggingLevel, setFortranLoggingLevel

from ..Commands import CALC_CHAMP
from ..Objects import (
    AssemblyMatrixDisplacementReal,
    DiscreteComputation,
    ElasticResult,
    LinearSolver,
    MechanicalDirichletBC,
    MechanicalLoadFunction,
    MechanicalLoadReal,
    ParallelMechanicalLoadFunction,
    ParallelMechanicalLoadReal,
    PhysicalProblem,
)
from ..Solvers import PhysicalState, StorageManager, TimeStepper
from ..Utilities import logger, print_stats, profile


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
                    MechanicalLoadFunction,
                    MechanicalLoadReal,
                    ParallelMechanicalLoadFunction,
                    ParallelMechanicalLoadReal,
                    MechanicalDirichletBC,
                ),
            ):
                phys_pb.addLoadFromDict(load)
            else:
                raise RuntimeError("Unknown load")

    phys_pb.computeListOfLoads()

    return phys_pb


@profile
def _createTimeStepper(args):
    """Create time stepper from catalogues

    Arguments:
        **args (dict): User's keywords.

    Returns:
        (TimeStepper): a time stepper.
    """

    timeValues = [0.0]
    inst = args.get("INST")
    if inst is not None:
        timeValues = [inst]
    else:
        listInst = args.get("LIST_INST")
        if listInst is not None:
            timeValues = listInst.getValues()
            inst_fin = args.get("INST_FIN")
            if inst_fin is not None:
                timeValues = [time for time in timeValues if time <= (inst_fin + 1.0e-6)]

            if "RESULTAT" in args:
                nbIndex = args["RESULTAT"].getNumberOfIndexes()
                inst_deb = args["RESULTAT"].getTimeValue(nbIndex)
                timeValues = [time for time in timeValues if time >= (inst_deb + 1.0e-6)]

    return TimeStepper(timeValues)


@profile
def _computeMatrix(disr_comp, matrix, time):
    """Compute and assemble the elastic matrix

    Arguments:
        disr_comp (DiscreteComputation): to compute discrete quantities
        matrix (AssemblyMatrixDisplacementReal): matrix to compute and assemble inplace
        time (float): current time

    Returns:
        AssemblyMatrixDisplacementReal: matrix computed and assembled
    """

    matr_elem = disr_comp.getLinearStiffnessMatrix(time=time, with_dual=True)
    matrix.addElementaryMatrix(matr_elem)

    profile(matrix.assemble)(True)

    return matrix


@profile
def _computeRhs(phys_pb, disr_comp, time):
    """Compute and assemble the right hand side

    Arguments:
         phys_pb (PhysicalProblem): physical problem
         disr_comp (DiscreteComputation): to compute discrete quantities
         time (float): current time

     Returns:
         FieldOnNodesReal: vector of load
    """

    # compute imposed displacement with Lagrange
    rhs = disr_comp.getImposedDualBC(time)

    # compute neumann forces
    rhs += disr_comp.getNeumannForces(time)

    if phys_pb.getMaterialField().hasExternalStateVariableForLoad():
        rhs += disr_comp.getExternalStateVariablesForces(time)

    return rhs


@profile
def _computeStress(phys_pb, result):
    """Compute SIEF_ELGA en STRX_ELGA

    Arguments:
        phys_pb (PhysicalProblem): phisical problem
        result (ElasticResult): result to fill in (in place)

    Returns:
        ElasticResult: result with stress fields
    """

    option = ["SIEF_ELGA"]

    if phys_pb.getModel().existsMultiFiberBeam():
        option.append("STRX_ELGA")

    result = CALC_CHAMP(reuse=result, RESULTAT=result, CONTRAINTE=option)

    return result


def meca_statique_ops(self, **args):
    """Execute the command MECA_STATIQUE.

    Arguments:
        **args (dict): User's keywords.

    Returns:
        ElasticResult: result for linear elasticity problem
    """

    logger.debug("<MECA_STATIQUE>: Initialization")

    verbosity = args["INFO"]
    setFortranLoggingLevel(verbosity)

    # Create result
    result = args.get("RESULTAT")
    if result is None:
        result = ElasticResult()
        title = args.get("TITRE")
        if title is not None:
            result.setTitle(title)

    # Create physical problem
    model = args["MODELE"]
    phys_pb = PhysicalProblem(model, args["CHAM_MATER"], args.get("CARA_ELEM"))

    # Add loads
    phys_pb = _addLoads(phys_pb, args)

    # Compute numbering
    phys_pb.computeDOFNumbering()

    # Create linear solver
    linear_solver = LinearSolver.factory("MECA_STATIQUE", args["SOLVEUR"])
    if (model.getMesh().isParallel()) and (not linear_solver.supportParallelMesh()):
        raise RuntimeError("ParallelMesh not allowed with this linear solver")

    if model.xfemPreconditioningEnable():
        linear_solver.enableXfem()
    linear_solver.build()

    # Create time stepper
    timeStepper = _createTimeStepper(args)

    # Create storage manager
    storage_manager = StorageManager(result)

    # Define main objects
    phys_state = PhysicalState()
    disc_comp = DiscreteComputation(phys_pb)

    # we define the matrix before to have an unique name
    # because of a bug with LDLT_SP
    matrix = AssemblyMatrixDisplacementReal(phys_pb)
    # the matrix depends on times or external variables
    isConst = phys_pb.getCodedMaterial().constant()
    isFirst = True

    # Compute reference value vector for external state variables
    if phys_pb.getMaterialField().hasExternalStateVariableWithReference():
        phys_pb.computeReferenceExternalStateVariables()

    # first index to use
    storage_manager.setInitialIndex(result.getNumberOfIndexes() + 1)

    # Run computation
    logger.debug("<MECA_STATIQUE>: Run computation")
    while not timeStepper.hasFinished():
        phys_state.time = timeStepper.getCurrent()

        # compute matrix and factorize it
        if not isConst or isFirst:
            matrix = _computeMatrix(disc_comp, matrix, phys_state.time)
            profile(linear_solver.factorize)(matrix)

        # compute rhs
        rhs = _computeRhs(phys_pb, disc_comp, phys_state.time)

        # solve linear system
        diriBCs = profile(disc_comp.getDirichletBC)(phys_state.time)
        phys_state.primal = profile(linear_solver.solve)(rhs, diriBCs)

        # store field
        storage_manager.storeState(phys_state.time, phys_pb, phys_state)

        timeStepper.completed()
        storage_manager.completed()
        isFirst = False

    # delete factorized matrix - free memory
    linear_solver.deleteFactorizedMatrix()

    # store field
    result = storage_manager.getResult()

    # cleaning because some objects are still on VOLATILE
    deleteTemporaryObjects()

    # compute stress if requested
    logger.debug("<MECA_STATIQUE>: Compute stress")
    if args["OPTION"] == "SIEF_ELGA":
        result = _computeStress(phys_pb, result)

    if verbosity > 1:
        print_stats()
    resetFortranLoggingLevel()

    return result
