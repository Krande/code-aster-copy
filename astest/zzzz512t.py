# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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

###################################################################################
#
#   Analytical solution
#   u = (X, Y)
#   eps = Id
#   sig = (2*mu + d *lamda) * Id
#   f = (0, 0)
#
#
####################################################################################


mesh0 = CA.ParallelMesh.buildCube(refine=1)

mesh = CREA_MAILLAGE(MAILLAGE=mesh0, MODI_HHO=_F(TOUT="OUI"))

mesh = DEFI_GROUP(reuse=mesh, MAILLAGE=mesh, CREA_GROUP_NO=(_F(GROUP_MA="VOLUME", NOM="VOLUME"),))

# define material
coeff = DEFI_MATERIAU(ELAS=_F(E=100.0, NU=0.3))

mater = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=coeff))


for form in ("LINEAIRE", "QUADRATIQUE"):
    model = AFFE_MODELE(
        MAILLAGE=mesh,
        AFFE=_F(TOUT="OUI", MODELISATION="3D_HHO", FORMULATION=form, PHENOMENE="MECANIQUE"),
    )

    uX = FORMULE(VALE="X+1", NOM_PARA=("X", "Y", "Z"))
    uY = FORMULE(VALE="Y+3", NOM_PARA=("X", "Y", "Z"))
    uZ = FORMULE(VALE="Z-2", NOM_PARA=("X", "Y", "Z"))

    bc = AFFE_CHAR_CINE_F(
        MODELE=model,
        MECA_IMPO=_F(
            GROUP_MA=("RIGHT", "LEFT", "FRONT", "BACK", "TOP", "BOTTOM"), DX=uX, DY=uY, DZ=uZ
        ),
    )

    # define discrete object
    phys_pb = CA.PhysicalProblem(model, mater)
    phys_pb.addDirichletBC(bc)
    phys_pb.computeDOFNumbering()

    hho = CA.HHO(phys_pb)

    # project function
    u_hho = hho.projectOnHHOSpace([uX, uY, uZ])

    # solve linear system
    LREEL = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=1, NOMBRE=1))

    resu = STAT_NON_LINE(
        MODELE=model, CHAM_MATER=mater, INCREMENT=_F(LIST_INST=LREEL), EXCIT=(_F(CHARGE=bc),)
    )

    u_sol = resu.getField("DEPL", para="INST", value=1.0)

    u_diff = u_hho - u_sol

    test.assertAlmostEqual(u_diff.norm("NORM_2"), 0.0, delta=1e-8)

    # test getValues on FieldOnNodes
    u_depl = resu.getField("HHO_DEPL", para="INST", value=1.0)
    test.assertAlmostEqual(
        sum(u_depl.getValues(["DX", "DY", "DZ"], ["VOLUME"])),
        sum(u_depl.getValues(["DX", "DY", "DZ"])),
    )
    test.assertAlmostEqual(
        sum(u_depl.getValues(["DX", "DY", "DZ"], ["VOLUME"])), sum(u_depl.getValues())
    )
    test.assertAlmostEqual(sum(u_depl.getValues([], ["VOLUME"])), sum(u_depl.getValues()))
    test.assertEqual(len(u_depl.getValues(["HHO_F1"], ["VOLUME"])), 0)
    test.assertEqual(len(u_depl.getValues(["HHO_F1"], ["VOL"])), 0)

# print("f=", f_hho.getValues())
# print("u=", u_hho.getValues())
# print("diif=", (u_hho - u_sol).getValues())

FIN()
