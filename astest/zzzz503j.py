# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

import code_aster
from code_aster.Commands import *

code_aster.init("--test")

test = code_aster.TestCase()

###################################################################
#
#   Solve Helmholtz problem with HHO
#   Continuous:
#   (A * grad u, grad v) + (H * u, v) = (H * f, v)
#   with f given and A > 0, H > 0
#
#   HHO:
#   sum_{T \in Th} (A * GkT(huT), GkT(hvT))_T + (H * u_T, v_T) _T = (H * f, v_T)_T
#
###################################################################

mesh0 = code_aster.Mesh.buildSquare(refine=3)

mesh = CREA_MAILLAGE(MAILLAGE=mesh0,
                     MODI_HHO=_F(TOUT='OUI',),)

# define model
model = AFFE_MODELE(MAILLAGE=mesh,
                    AFFE=_F(TOUT='OUI',
                            MODELISATION='PLAN_HHO',
                            FORMULATION='LINEAIRE',
                            PHENOMENE='THERMIQUE')
                    )

# define material
H = 2.0
A = 4.0

coeff = DEFI_MATERIAU(THER=_F(LAMBDA=A,  RHO_CP=H))

mater = AFFE_MATERIAU(MAILLAGE=mesh,
                      AFFE=_F(TOUT='OUI',  MATER=coeff)
                      )

# define BC
f = 100.
load = AFFE_CHAR_THER(MODELE=model,
                      SOURCE=_F(GROUP_MA='SURFACE',  SOUR=H*f),
                      )

bc = AFFE_CHAR_CINE(MODELE=model,
                    THER_IMPO=_F(GROUP_MA=('RIGHT', 'LEFT',
                                 'TOP', 'BOTTOM',),  TEMP=0.0),
                    )

# define discrete object
phys_pb = code_aster.PhysicalProblem(model, mater)
phys_pb.addDirichletBC(bc)
phys_pb.addLoad(load)

disc_comp = code_aster.DiscreteComputation(phys_pb)

# compute DOF numbering
phys_pb.computeDOFNumbering()

# compute (A * GkT(huT), GkT(hvT))_T
matK = disc_comp.getLinearStiffnessMatrix()

# compute (H * u_T, v_T) _T
matM = disc_comp.getMassMatrix()

# compute (H * f, v_T)_T
rhs = disc_comp.getNeumannForces()

# compute BC
diriBCs = disc_comp.getDirichletBC()

# assemble matrix
matrix = code_aster.AssemblyMatrixTemperatureReal(phys_pb)

matrix.addElementaryMatrix(matK)
matrix.addElementaryMatrix(matM)
matrix.assemble()

# solve linear system
mySolver = code_aster.MumpsSolver()
mySolver.factorize(matrix)
solution = mySolver.solve(rhs, diriBCs)

test.assertAlmostEqual(solution.norm("NORM_2"), 28.660962846752938, delta=1e-6)

# project HHO solution
hho_field = code_aster.HHO(phys_pb).projectOnLagrangeSpace(solution)
test.assertAlmostEqual(hho_field.norm("NORM_2"), 32.18599383551358, delta=1e-6)

# save result
hho_field.printMedFile("hhoField.med")

test.printSummary()

FIN()
