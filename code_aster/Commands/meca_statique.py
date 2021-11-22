# coding: utf-8

# Copyright (C) 1991 - 2021  EDF R&D                www.code-aster.org
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

# person_in_charge: nicolas.sellenet@edf.fr

from ..Objects import ( MechanicalLoadReal, DirichletBC, MechanicalLoadFunction,
                        ParallelMechanicalLoadReal, ParallelMechanicalLoadFunction,
                        LinearStaticAnalysis, ElasticResult)
from ..Supervis import ExecuteCommand
from ..Utilities import force_list
from .common_keywords import create_solver


class MechanicalSolver(ExecuteCommand):
    """Solver for static linear mechanical problems."""
    command_name = "MECA_STATIQUE"

    def create_result(self, keywords):
        """Create Result"""
        self._result = ElasticResult()

    @staticmethod
    def _addLoad(mechaSolv, fkw):
        load = fkw["CHARGE"]
        func = fkw.get("FONC_MULT")

        if func !=  None:
            if isinstance(load, DirichletBC):
                mechaSolv.addDirichletBC(load, func)
            elif isinstance(load, ( MechanicalLoadReal, MechanicalLoadFunction,
                                ParallelMechanicalLoadReal, ParallelMechanicalLoadFunction )):
                mechaSolv.addLoad(load, func)
            else:
                assert False
        else:
            if isinstance(load, DirichletBC):
                mechaSolv.addDirichletBC(load)
            elif isinstance(load, ( MechanicalLoadReal, MechanicalLoadFunction,
                                ParallelMechanicalLoadReal, ParallelMechanicalLoadFunction )):
                mechaSolv.addLoad(load)
            else:
                assert False

    def exec_(self, keywords):
        """Execute the command.

        Arguments:
            keywords (dict): User's keywords.
        """
        model = keywords["MODELE"]
        materialField = keywords["CHAM_MATER"]
        caraElem = keywords.get("CARA_ELEM")

        mechaSolv = None
        if caraElem != None:
            mechaSolv = LinearStaticAnalysis(model, materialField, caraElem)
        else:
            mechaSolv = LinearStaticAnalysis(model, materialField)

        inst = keywords.get("INST")
        if inst != None:
            mechaSolv.setTimeStepManager([inst])
        listInst = keywords.get("LIST_INST")
        if listInst != None:
            timeValues = listInst.getValues()
            inst_fin = keywords.get("INST_FIN")
            if inst_fin is not None:
                timeValues = [time for time in timeValues if time <= inst_fin]
            mechaSolv.setTimeStepManager(timeValues)

        fkw = force_list(keywords.get("EXCIT", []))
        for curDict in fkw:
            self._addLoad(mechaSolv, curDict)

        mechaSolv.setStressComputation( keywords["OPTION"] == "SIEF_ELGA" )

        solver = create_solver(keywords.get("SOLVEUR"))
        mechaSolv.setLinearSolver(solver)
        self._result = mechaSolv.execute( self._result )

    def post_exec(self, keywords):
        """Execute the command.

        Arguments:
            keywords (dict): User's keywords.
        """

        if self._result is not None and self._result.getNumberOfRanks() > 0:
            self._result.build()


    def add_dependencies(self, keywords):
        """Register input *DataStructure* objects as dependencies.

        Arguments:
            keywords (dict): User's keywords.
        """

        # Add no dependencies since everything is done in c++ directly

MECA_STATIQUE = MechanicalSolver.run
