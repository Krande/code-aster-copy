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

mesh0 = CA.Mesh.buildSquare(refine=0)

mesh = CREA_MAILLAGE(MAILLAGE=mesh0, MODI_HHO=_F(TOUT="OUI"))

# define model
model = AFFE_MODELE(
    MAILLAGE=mesh,
    AFFE=_F(TOUT="OUI", MODELISATION="D_PLAN_HHO", FORMULATION="LINEAIRE", PHENOMENE="MECANIQUE"),
)

# define material
coeff = DEFI_MATERIAU(ELAS=_F(E=1.0, NU=0.0))

mater = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=coeff))

# define discrete object
phys_pb = CA.PhysicalProblem(model, mater)
phys_pb.computeDOFNumbering()
hho = CA.HHO(phys_pb)

f_hho = hho.projectOnHHOSpace([1.0, 2.0])

ref = (1 + 2) * (1 + 4)
test.assertAlmostEqual(sum(f_hho.getValues()), ref, delta=5e-5)

fx = FORMULE(VALE="X", NOM_PARA=["X", "Y"])
fy = FORMULE(VALE="X", NOM_PARA=["X", "Y"])

f2_hho = hho.projectOnHHOSpace([fx, fy])
test.assertAlmostEqual(sum(f2_hho.getValues()), 4.422649730810375, delta=5e-5)

test.printSummary()

FIN()
