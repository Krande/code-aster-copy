# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

mesh = CREA_MAILLAGE(MAILLAGE=mesh0, MODI_HHO=_F(TOUT="OUI"))

# define model
model = AFFE_MODELE(
    MAILLAGE=mesh,
    AFFE=_F(TOUT="OUI", MODELISATION="PLAN_HHO", FORMULATION="LINEAIRE", PHENOMENE="THERMIQUE"),
)

# define material
H = 2.0
A = 4.0

coeff = DEFI_MATERIAU(THER=_F(LAMBDA=1.0, RHO_CP=1.0))

mater = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=coeff))

# define BC
f = 100.0
bc = AFFE_CHAR_CINE(
    MODELE=model, THER_IMPO=_F(GROUP_MA=("RIGHT", "LEFT", "TOP", "BOTTOM"), TEMP=0.0)
)
load = AFFE_CHAR_THER(MODELE=model, SOURCE=_F(GROUP_MA="SURFACE", SOUR=H * f))

# define discrete object
phys_pb = code_aster.PhysicalProblem(model, mater)
phys_pb.addDirichletBC(bc)
phys_pb.addLoad(load)

disc_comp = code_aster.DiscreteComputation(phys_pb)

hho = code_aster.HHO(phys_pb)

# compute DOF numbering
phys_pb.computeDOFNumbering()

# compute (GkT(huT), GkT(hvT))_T
matEK = disc_comp.getLinearStiffnessMatrix()
matK = code_aster.AssemblyMatrixTemperatureReal(phys_pb)
matK.addElementaryMatrix(matEK)
matK.assemble()

# compute (u_T, v_T) _T
matEM = disc_comp.getMassMatrix()
matM = code_aster.AssemblyMatrixTemperatureReal(phys_pb)
matM.addElementaryMatrix(matEM)
matM.assemble()

# compute ( H * f, v_T)_T
form = FORMULE(VALE="X-X+Y-Y+100", NOM_PARA=["X", "Y"])
f_hho = hho.projectOnHHOSpace(form)
rhs2 = H * matM * f_hho
rhs = disc_comp.getVolumetricForces()

# test.assertAlmostEqual(rhs.norm("NORM_2"), rhs2.norm("NORM_2"), delta=1e-6)

# compute BC
diriBCs = disc_comp.getDirichletBC()

# lhs matrix
lhs = A * matK + H * matM

# solve linear system
mySolver = code_aster.MumpsSolver()
mySolver.factorize(lhs)
solution = mySolver.solve(rhs, diriBCs)

sol_ref = 28.544813405889784
test.assertAlmostEqual((solution.norm("NORM_2") - sol_ref) / sol_ref, 0.0, delta=1e-4)

# project HHO solution
h1_field = hho.projectOnLagrangeSpace(solution)
hho_field = hho.projectOnHHOSpace(h1_field)
h1_ref = 32.400164138793706
test.assertAlmostEqual((h1_field.norm("NORM_2") - h1_ref) / h1_ref, 0.0, delta=1e-6)

# save result
hho_field.printMedFile("hhoField.med")

# Non-linear process with A(u) = A * (1.1+max(u)_Omega), H(u) = H * (1 + max(u)_Omega)
u_hho = hho.projectOnHHOSpace(0.0)

print("Newton solver:")
for i in range(100):
    max_u = u_hho.norm("NORM_INFINITY")
    Au = A * (1.1 + max_u)
    Hu = H * (1 + max_u)
    Resi = Au * matK * u_hho + Hu * matM * (u_hho - f_hho)
    Jaco = Au * matK + Hu * matM

    print("*Iter %d: residual %f" % (i, Resi.norm("NORM_2")))
    if Resi.norm("NORM_2") < 10e-8:
        break

    mySolver.factorize(Jaco)
    du_hho = mySolver.solve(-Resi, diriBCs)
    u_hho += du_hho

u_hho_ref = 27.931912612339957
test.assertAlmostEqual((u_hho.norm("NORM_2") - u_hho_ref) / u_hho_ref, 0.0, delta=1e-4)

# test projection
const = FORMULE(VALE="X-X+Y-Y+100", NOM_PARA=["X", "Y"])
f_hho = hho.projectOnHHOSpace(const)
f_lagr = hho.projectOnLagrangeSpace(f_hho)
f2_hho = hho.projectOnHHOSpace(f_lagr)
diff = f_hho - f2_hho
test.assertAlmostEqual(diff.norm("NORM_2"), 0.0, delta=1e-4)

f0 = FORMULE(VALE="X+Y", NOM_PARA=["X", "Y"])
f_hho = hho.projectOnHHOSpace(f0)
f_lagr = hho.projectOnLagrangeSpace(f_hho)
f2_hho = hho.projectOnHHOSpace(f_lagr)
f2_lagr = hho.projectOnLagrangeSpace(f2_hho)
diff = f_lagr - f2_lagr
test.assertAlmostEqual(diff.norm("NORM_2"), 0.0, delta=1e-4)
diff = f_hho - f2_hho
test.assertAlmostEqual(diff.norm("NORM_2"), 0.0, delta=1e-4)

f1 = FORMULE(VALE="X*X+Y*Y", NOM_PARA=["X", "Y"])
f_hho = hho.projectOnHHOSpace(f1)
f_lagr = hho.projectOnLagrangeSpace(f_hho)
f2_hho = hho.projectOnHHOSpace(f_lagr)
f2_lagr = hho.projectOnLagrangeSpace(f2_hho)
diff = f_lagr - f2_lagr
test.assertAlmostEqual(diff.norm("NORM_2"), 0.0, delta=1e-4)
diff = f_hho - f2_hho
test.assertAlmostEqual(diff.norm("NORM_2"), 0.0, delta=5e-2)

hho_elga = hho.evaluateAtQuadraturePoints(f_hho)
test.assertAlmostEqual(max(hho_elga.getValues()), 1.8957, delta=1e-2)

f_proj = hho.projectOnHHOCellSpace(hho_elga)
test.assertAlmostEqual(f_proj.norm("NORM_2"), 6.376190569410008, delta=1e-2)
f3_lagr = hho.projectOnLagrangeSpace(f2_hho)
diff = f_lagr - f3_lagr
test.assertAlmostEqual(diff.norm("NORM_2"), 0.0, delta=1e-4)

# this forces deletion of MUMPS objects
del mySolver

test.printSummary()

FIN()
