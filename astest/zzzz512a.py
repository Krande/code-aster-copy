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

CA.init("--test")

test = CA.TestCase()

###################################################################################
#
#   Patch test with analytical solution
#   Solution is a polynomial of order k
#   So method of order k should have a null error
#
#   The script to compute solution is given in zzzz512a.datg
#
####################################################################################

E = 200000.0
Nu = 0.3

lamb = E * Nu / (1 + Nu) / (1 - 2 * Nu)
mu = E / 2 / (1 + Nu)

uX = {
    "LINEAIRE": FORMULE(VALE="X+1", NOM_PARA=("X", "Y", "Z")),
    "QUADRATIQUE": FORMULE(VALE="Z*(X+1)", NOM_PARA=("X", "Y", "Z")),
    "CUBIQUE": FORMULE(VALE="Z*(X+1)*Y+X*X", NOM_PARA=("X", "Y", "Z")),
    "QUARTIQUE": FORMULE(VALE="Z*(X+1)*Y*X+X*X", NOM_PARA=("X", "Y", "Z")),
}
uY = {
    "LINEAIRE": FORMULE(VALE="Y+3", NOM_PARA=("X", "Y", "Z")),
    "QUADRATIQUE": FORMULE(VALE="X*(Y+3)", NOM_PARA=("X", "Y", "Z")),
    "CUBIQUE": FORMULE(VALE="X*(Y+3)*Z+Y*Y", NOM_PARA=("X", "Y", "Z")),
    "QUARTIQUE": FORMULE(VALE="X*(Y+3)*Z*Y+Y*Y", NOM_PARA=("X", "Y", "Z")),
}
uZ = {
    "LINEAIRE": FORMULE(VALE="Z-2", NOM_PARA=("X", "Y", "Z")),
    "QUADRATIQUE": FORMULE(VALE="Y*(Z-2)", NOM_PARA=("X", "Y", "Z")),
    "CUBIQUE": FORMULE(VALE="Y*(Z-2)*X+Z*Z", NOM_PARA=("X", "Y", "Z")),
    "QUARTIQUE": FORMULE(VALE="Y*(Z-2)*X*Z+Z*Z", NOM_PARA=("X", "Y", "Z")),
}

zero = FORMULE(VALE="0", NOM_PARA=("X", "Y", "Z"))
f0 = -(lamb + mu)
fc = FORMULE(VALE="f0", NOM_PARA=("X", "Y", "Z"), f0=f0)

fX = {
    "LINEAIRE": zero,
    "QUADRATIQUE": fc,
    "CUBIQUE": FORMULE(
        VALE="-lamb*(Y + Z + 2.0) - mu*Y - mu*Z - 4.0*mu",
        NOM_PARA=("X", "Y", "Z"),
        lamb=lamb,
        mu=mu,
    ),
    "QUARTIQUE": FORMULE(
        VALE="-6.0*lamb*Y*Z + 2.0*lamb*Y - 3.0*lamb*Z - 2.0*lamb - 8.0*mu*Y*Z + 2.0*mu*Y - 3.0*mu*Z - 4.0*mu",
        NOM_PARA=("X", "Y", "Z"),
        lamb=lamb,
        mu=mu,
    ),
}
fY = {
    "LINEAIRE": zero,
    "QUADRATIQUE": fc,
    "CUBIQUE": FORMULE(
        VALE="-lamb*(X + Z + 2.0) - mu*X - mu*Z - 4.0*mu",
        NOM_PARA=("X", "Y", "Z"),
        lamb=lamb,
        mu=mu,
    ),
    "QUARTIQUE": FORMULE(
        VALE="-6.0*lamb*X*Z + 2.0*lamb*X - lamb*Z - 2.0*lamb - 8.0*mu*X*Z + 2.0*mu*X - mu*Z - 4.0*mu",
        NOM_PARA=("X", "Y", "Z"),
        lamb=lamb,
        mu=mu,
    ),
}
fZ = {
    "LINEAIRE": zero,
    "QUADRATIQUE": fc,
    "CUBIQUE": FORMULE(
        VALE="-lamb*(X + Y + 2.0) - mu*X - mu*Y - 4.0*mu",
        NOM_PARA=("X", "Y", "Z"),
        lamb=lamb,
        mu=mu,
    ),
    "QUARTIQUE": FORMULE(
        VALE="-6.0*lamb*X*Y - 3.0*lamb*X - lamb*Y - 2.0*lamb - 8.0*mu*X*Y - 3.0*mu*X - mu*Y - 4.0*mu",
        NOM_PARA=("X", "Y", "Z"),
        lamb=lamb,
        mu=mu,
    ),
}

mesh0 = LIRE_MAILLAGE(FORMAT="MED", UNITE=20)

mesh = CREA_MAILLAGE(MAILLAGE=mesh0, MODI_HHO=_F(TOUT="OUI"))

mesh = DEFI_GROUP(
    reuse=mesh, MAILLAGE=mesh, CREA_GROUP_MA=_F(NOM="3D", TOUT="OUI", TYPE_MAILLE="3D")
)

# define material
coeff = DEFI_MATERIAU(ELAS=_F(E=E, NU=Nu, RHO=1.0), HHO=_F(COEF_STAB=2 * mu))

mater = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=coeff))

nbCells = mesh.getNumberOfCells()
if nbCells > 1500:
    formu = ["LINEAIRE", "QUADRATIQUE", "CUBIQUE"]
else:
    formu = ["LINEAIRE", "QUADRATIQUE", "CUBIQUE", "QUARTIQUE"]

for form in formu:
    model = AFFE_MODELE(
        MAILLAGE=mesh,
        AFFE=_F(TOUT="OUI", MODELISATION="3D_HHO", FORMULATION=form, PHENOMENE="MECANIQUE"),
    )

    bc = AFFE_CHAR_CINE_F(
        MODELE=model, MECA_IMPO=_F(GROUP_MA="FACES", DX=uX[form], DY=uY[form], DZ=uZ[form])
    )

    load = AFFE_CHAR_MECA_F(
        MODELE=model, FORCE_INTERNE=_F(GROUP_MA="3D", FX=fX[form], FY=fY[form], FZ=fZ[form])
    )

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
    u_hho = hho.projectOnHHOSpace([uX[form], uY[form], uZ[form]])

    u_diff = u_hho - u_sol

    test.assertAlmostEqual(u_diff.norm("NORM_2") / u_hho.norm("NORM_2"), 0.0, delta=5e-6)

    dc = CA.DiscreteComputation(phys_pb)
    mass = dc.getMassMatrix(assembly=True)
    l2_diff = (mass * u_diff).dot(u_diff)
    l2_ref = (mass * u_hho).dot(u_hho)
    test.assertAlmostEqual(l2_diff / l2_ref, 0.0, delta=1e-10)

    # print("u=", u_hho.getValues())
    # print("uh=", u_sol.getValues())
    # print("diif=", u_diff.getValues())

FIN()
