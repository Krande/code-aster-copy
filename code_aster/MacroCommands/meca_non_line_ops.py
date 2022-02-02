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


from ..Cata.Syntax import _F
from ..Objects import (
    MechanicalDirichletBC,
    MechanicalLoadFunction,
    MechanicalLoadReal,
    ParallelMechanicalLoadFunction,
    ParallelMechanicalLoadReal,
    NonLinearResult,
)
from ..Messages import UTMESS
from ..Utilities import print_stats
from .NonLinearSolver import NonLinearSolver, TimeStepper


def meca_non_line_ops(self, **args):
    """Execute the command.

    Arguments:
        **args (dict): User's keywords.
    """

    UTMESS('A', 'QUALITY1_2')

    args = _F(args)

    snl = NonLinearSolver()
    snl.setLoggingLevel(args["INFO"])
    snl.setPhysicalProblem(args["MODELE"], args["CHAM_MATER"], args["CARA_ELEM"])

    # Add parameters
    snl.setKeywords(
        CONVERGENCE=args["CONVERGENCE"],
        NEWTON=args["NEWTON"],
        ETAT_INIT=args["ETAT_INIT"],
        INCREMENT=args["INCREMENT"],
        INFO=args["INFO"],
    )

    # Add behaviour
    snl.setBehaviourProperty(args["COMPORTEMENT"])

    # Add loads
    if args["EXCIT"] is not None:
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
                snl.phys_pb.addLoadFromDict(load)
            else:
                raise RuntimeError("Unknown load")

    # Add linear solver
    snl.setLinearSolver(keywords=args["SOLVEUR"])

    # Add stepper
    tini = None
    if args["ETAT_INIT"] is not None:
        if "EVOL_NOLI" in args["ETAT_INIT"]:
            resu = args["ETAT_INIT"].get("EVOL_NOLI")
            assert isinstance(resu, NonLinearResult)
            tini = resu.getTimeValue(resu.getNumberOfRanks() - 1)
            if "INST_ETAT_INIT" in args["ETAT_INIT"]:
                tini = args["ETAT_INIT"].get("INST_ETAT_INIT")
    timeStepper = TimeStepper(args["INCREMENT"]["LIST_INST"].getValues()[1::])
    if tini is not None:
        timeStepper.updateTimes(tini)
    snl.setStepper(timeStepper)

    # Run computation
    snl.run()
    print_stats()
    return snl.getResult()
