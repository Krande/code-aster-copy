# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2026 - EDF R&D - www.code-aster.org
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

# --------------------------------------------------------------------
# This file is part of code_aster.
#
# code_aster is free software: you can redistribute it and/or modifZ
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

CA.init("--test")

test = CA.TestCase()

###################################################################################
#
#   Patch test with analytical solution
#   Solution is a polynomial of order k
#   So method of order k should have a null error
#
#   The script to compute solution is given in zzzz512v.datg
#
####################################################################################

E = 200000.0
Nu = 0.3

lamb = E * Nu / (1 + Nu) / (1 - 2 * Nu)
mu = E / 2 / (1 + Nu)

uR = {
    "LINEAIRE": FORMULE(VALE="X", NOM_PARA=("X", "Y")),
    "QUADRATIQUE": FORMULE(VALE="X", NOM_PARA=("X", "Y")),
    "CUBIQUE": FORMULE(VALE="X", NOM_PARA=("X", "Y")),
    "QUARTIQUE": FORMULE(VALE="X", NOM_PARA=("X", "Y")),
}
uZ = {
    "LINEAIRE": FORMULE(VALE="Y+1", NOM_PARA=("X", "Y")),
    "QUADRATIQUE": FORMULE(VALE="Y+1", NOM_PARA=("X", "Y")),
    "CUBIQUE": FORMULE(VALE="Y+1", NOM_PARA=("X", "Y")),
    "QUARTIQUE": FORMULE(VALE="Y+1", NOM_PARA=("X", "Y")),
}

zero = FORMULE(VALE="0", NOM_PARA=("X", "Y"))

fR = {
    "LINEAIRE": zero,
    "QUADRATIQUE": zero,
    "CUBIQUE": FORMULE(VALE="0.0", NOM_PARA=("X", "Y"), lamb=lamb, mu=mu),
    "QUARTIQUE": FORMULE(VALE="0.0", NOM_PARA=("X", "Y"), lamb=lamb, mu=mu),
}
fZ = {
    "LINEAIRE": zero,
    "QUADRATIQUE": zero,
    "CUBIQUE": FORMULE(VALE="0.0", NOM_PARA=("X", "Y"), lamb=lamb, mu=mu),
    "QUARTIQUE": FORMULE(VALE="0.0", NOM_PARA=("X", "Y"), lamb=lamb, mu=mu),
}

mesh0 = LIRE_MAILLAGE(FORMAT="MED", UNITE=20)

mesh = CREA_MAILLAGE(MAILLAGE=mesh0, MODI_HHO=_F(TOUT="OUI"))

mesh = DEFI_GROUP(
    reuse=mesh, MAILLAGE=mesh, CREA_GROUP_MA=_F(NOM="2D", TOUT="OUI", TYPE_MAILLE="2D")
)

# define material
coeff = DEFI_MATERIAU(ELAS=_F(E=E, NU=Nu, RHO=1.0), HHO=_F(COEF_STAB=2 * mu))

mater = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=coeff))

for form in ["LINEAIRE", "QUADRATIQUE", "CUBIQUE", "QUARTIQUE"]:
    model = AFFE_MODELE(
        MAILLAGE=mesh,
        AFFE=_F(TOUT="OUI", MODELISATION="AXIS_HHO", FORMULATION=form, PHENOMENE="MECANIQUE"),
    )

    bc = AFFE_CHAR_CINE_F(
        MODELE=model, MECA_IMPO=_F(GROUP_MA="BOUNDARIES", DX=uR[form], DY=uZ[form])
    )

    load = AFFE_CHAR_MECA_F(MODELE=model, FORCE_INTERNE=_F(GROUP_MA="2D", FX=fR[form], FY=fZ[form]))

    # solve linear system
    LREEL = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=1, NOMBRE=1))

    resu = STAT_NON_LINE(
        MODELE=model,
        CHAM_MATER=mater,
        INCREMENT=_F(LIST_INST=LREEL),
        EXCIT=(_F(CHARGE=bc), _F(CHARGE=load)),
    )

    u_sol = resu.getField("DEPL", para="INST", value=1.0)

    # define discrete object
    phys_pb = CA.PhysicalProblem(model, mater)
    phys_pb.addDirichletBC(bc)
    phys_pb.computeDOFNumbering()

    hho = CA.HHO(phys_pb)

    # project function
    u_hho = hho.projectOnHHOSpace([uR[form], uZ[form]])

    u_diff = u_hho - u_sol

    test.assertAlmostEqual(u_diff.norm("NORM_2") / u_hho.norm("NORM_2"), 0.0, delta=5e-6)

    dc = CA.DiscreteComputation(phys_pb)
    mass = dc.getMassMatrix(assembly=True)
    l2_diff = (mass * u_diff).dot(u_diff)
    l2_ref = (mass * u_hho).dot(u_hho)
    test.assertAlmostEqual(l2_diff / l2_ref, 0.0, delta=1e-10)

    # compute (u_T, v_T) _T
    matEM = CALC_MATR_ELEM(MODELE=model, OPTION="MASS_MECA", CHAM_MATER=mater)

    matM = CA.AssemblyMatrixDisplacementReal(phys_pb)
    matM.assemble(matEM, phys_pb.getListOfLoads())

    l2_diff = (matM * u_diff).dot(u_diff)
    test.assertAlmostEqual(l2_diff / l2_ref, 0.0, delta=1e-10)

FIN()
