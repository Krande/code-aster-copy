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

from libaster import deleteTemporaryObjects, setFortranLoggingLevel, resetFortranLoggingLevel

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
from ..Utilities import logger, print_stats, profile
from .NonLinearSolver import PhysicalState, StorageManager, TimeStepper


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

    return TimeStepper(timeValues)


@profile
def _computeMatrix(disr_comp, matrix, time, externVar):
    """Compute and assemble the elastic matrix

    Arguments:
        disr_comp (DiscreteComputation): to compute discrete quantities
        matrix (AssemblyMatrixDisplacementReal): matrix to compute and assemble inplace
        time (float): current time
        externVar (field): current external state variables

    Returns:
        AssemblyMatrixDisplacementReal: matrix computed and assembled
    """

    matrix.clearElementaryMatrix()
    matr_elem = disr_comp.elasticStiffnessMatrix(time, externVarField=externVar)
    matrix.addElementaryMatrix(matr_elem)
    matrix.assemble()

    return matrix


@profile
def _computeRhs(phys_pb, disr_comp, time, timeField, externVarField):
    """Compute and assemble the right hand side

    Arguments:
         phys_pb (PhysicalProblem): physical problem
         disr_comp (DiscreteComputation): to compute discrete quantities
         time (float): current time
         timeField (ConstantFieldOnCell): field with value of current time
         externVarField (fieldOnCellsReal): external state variable at current time

     Returns:
         FieldOnNodesReal: vector of load
    """

    # compute imposed displacement with Lagrange
    rhs = disr_comp.imposedDisplacement(time)

    # compute neumann forces
    rhs += disr_comp.neumann([time, 0.0, 0.0], externVarField=externVarField)

    if phys_pb.getMaterialField().hasExternalStateVariableForLoad():
        rhs += disr_comp.computeExternalStateVariablesLoad(time, timeField, externVarField)

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

    # Les noms des champs sont "en dur" à cause de nmvcex/vectme => à supprimer après 31903
    externVarName = phys_pb.getMaterialField().getName() + "      .TOUT"
    timeFieldName = phys_pb.getMaterialField().getName() + "      .INST"
    externVarRefeName = phys_pb.getModel().getName()[0:8] + ".CHVCREF   "

    # Detect external state variables
    hasExternalStateVariable = phys_pb.getMaterialField().hasExternalStateVariable()

    # Compute reference value vector for external state variables
    externVarRefe = None
    if phys_pb.getMaterialField().hasExternalStateVariableWithReference():
        externVarRefe = disc_comp.computeExternalStateVariablesReference(externVarRefeName)
        phys_pb.setExternalStateVariablesReference(externVarRefe)

    # first rank to use
    rank = result.getNumberOfRanks() + 1

    # Run computation
    logger.debug("<MECA_STATIQUE>: Run computation")
    while not timeStepper.hasFinished():
        phys_state.time = timeStepper.getNext()

        # Update external state variable if required
        if hasExternalStateVariable:
            phys_state.externVar = disc_comp.createExternalStateVariablesField(
                externVarName, phys_state.time
            )

        # Update time field
        timeField = disc_comp.createTimeField(timeFieldName, phys_state.time)

        # compute matrix and factorize it
        if not isConst or isFirst:
            matrix = _computeMatrix(disc_comp, matrix, phys_state.time, phys_state.externVar)
            profile(linear_solver.factorize)(matrix)

        # compute rhs
        rhs = _computeRhs(phys_pb, disc_comp, phys_state.time, timeField, phys_state.externVar)

        # solve linear system
        diriBCs = profile(disc_comp.dirichletBC)(phys_state.time)
        phys_state.displ = profile(linear_solver.solve)(rhs, diriBCs)

        # store rank
        storage_manager.storeState(rank, phys_state.time, phys_pb, phys_state)

        timeStepper.completed()
        rank += 1
        isFirst = False

        # Suppression obligatoire car nom en "dur" => à supprimer après 31903
        del timeField
        phys_state.externVar = None

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
