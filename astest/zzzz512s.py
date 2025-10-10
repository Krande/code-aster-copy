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


from code_aster.Commands import *
from code_aster import CA

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

###################################################################
#
#   Solve Poisson problem with HHO
#
#   -div(A*grad(u))) = f
#   Continuous:
#   (A * grad u, grad v) = (f, v)
#   with f given and A > 0
#
#   Solution is a polynomial of order k
#   So method of order k should have a null error
#
#   The script to compute solution is given in zzzz512i.datg
#
#   HHO:
#   sum_{T \in Th} (A * GkT(huT), GkT(hvT))_T = (f, v_T)_T
#
###################################################################


# define material
A = 4.0

u = {
    "CONSTANTE": FORMULE(VALE="-1", NOM_PARA=("X", "Y", "Z")),
    "LINEAIRE": FORMULE(VALE="Z+X-1", NOM_PARA=("X", "Y", "Z")),
    "QUADRATIQUE": FORMULE(VALE="Y*Y-X*Y+Z+X-1", NOM_PARA=("X", "Y", "Z")),
    "CUBIQUE": FORMULE(VALE="Z*Z*X+Y*Y-X*Y+Z+X-1", NOM_PARA=("X", "Y", "Z")),
    "QUARTIQUE": FORMULE(VALE="-X*X*X*X+Z*Z*X+Y*Y-X*Y+Z+X-1", NOM_PARA=("X", "Y", "Z")),
}

f = {
    "CONSTANTE": FORMULE(VALE="0", NOM_PARA=("X", "Y", "Z")),
    "LINEAIRE": FORMULE(VALE="0", NOM_PARA=("X", "Y", "Z")),
    "QUADRATIQUE": FORMULE(VALE="-2*A", NOM_PARA=("X", "Y", "Z"), A=A),
    "CUBIQUE": FORMULE(VALE="-2*A*(X+1)", NOM_PARA=("X", "Y", "Z"), A=A),
    "QUARTIQUE": FORMULE(VALE="2*A*(6*X*X-X-1)", NOM_PARA=("X", "Y", "Z"), A=A),
}

mesh0 = LIRE_MAILLAGE(FORMAT="MED", UNITE=20)

mesh = CREA_MAILLAGE(MAILLAGE=mesh0, MODI_HHO=_F(TOUT="OUI"))

mesh = DEFI_GROUP(
    reuse=mesh, MAILLAGE=mesh, CREA_GROUP_MA=_F(NOM="3D", TOUT="OUI", TYPE_MAILLE="3D")
)

coeff = DEFI_MATERIAU(THER=_F(LAMBDA=1.0, RHO_CP=1.0))

mater = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=coeff))

formu = ["CONSTANTE", "LINEAIRE", "QUADRATIQUE", "CUBIQUE", "QUARTIQUE"]

for form in formu:
    # define model
    model = AFFE_MODELE(
        MAILLAGE=mesh,
        AFFE=_F(TOUT="OUI", MODELISATION="3D_HHO", FORMULATION=form, PHENOMENE="THERMIQUE"),
    )

    bc = AFFE_CHAR_CINE_F(MODELE=model, THER_IMPO=_F(TEMP=u[form], GROUP_MA="FACES"), INFO=1)
    load = AFFE_CHAR_THER_F(MODELE=model, SOURCE=_F(GROUP_MA="3D", SOUR=f[form]))

    # define discrete object
    phys_pb = CA.PhysicalProblem(model, mater)
    phys_pb.addDirichletBC(bc)
    phys_pb.addLoad(load)

    disc_comp = CA.DiscreteComputation(phys_pb)

    # compute DOF numbering
    phys_pb.computeDOFNumbering()

    # compute (GkT(huT), GkT(hvT))_T
    matEK = disc_comp.getLinearStiffnessMatrix()
    matK = CA.AssemblyMatrixTemperatureReal(phys_pb)
    matK.assemble(matEK, phys_pb.getListOfLoads())

    # compute (u_T, v_T) _T
    matEM = disc_comp.getMassMatrix()
    matM = CA.AssemblyMatrixTemperatureReal(phys_pb)
    matM.assemble(matEM, phys_pb.getListOfLoads())

    # compute ( f, v_T)_T
    rhs = disc_comp.getVolumetricForces()

    # lhs matrix
    lhs = A * matK

    # BC
    diriBC = disc_comp.getDirichletBC()

    # solve linear system
    mySolver = CA.MumpsSolver()
    mySolver.factorize(lhs)
    u_sol = mySolver.solve(rhs, diriBC)

    del mySolver

    hho = CA.HHO(phys_pb)

    # project function
    u_hho = hho.projectOnHHOSpace(u[form])

    u_diff = u_hho - u_sol

    test.assertAlmostEqual(u_diff.norm("NORM_2") / u_hho.norm("NORM_2"), 0.0, delta=1e-7)

    l2_diff = (matM * u_diff).dot(u_diff)
    l2_ref = (matM * u_hho).dot(u_hho)
    test.assertAlmostEqual(l2_diff / l2_ref, 0.0, delta=1e-10)

    # project HHO solution
    h1_field = hho.projectOnLagrangeSpace(u_sol)
    hho_field = hho.projectOnHHOSpace(h1_field)
    # not zero if k > 2 because quadratic cell
    # test.assertAlmostEqual((hho_field - u_sol).norm("NORM_2"), 0.0, delta=1e-6)


test.printSummary()

FIN()
