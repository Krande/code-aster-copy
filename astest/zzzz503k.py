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
from zzzz503k_func import *

code_aster.init("--test")

test = code_aster.TestCase()

###################################################################################
#
#   Solve coupling problem with HHO
#   u -> displacement, dx -> macro-damage, d -> micro-damage
#
#   Continuous:
#   R_u(u,dx,d) = (PK1(u,dx,d), grad v) = 0, \forall v
#   R_dx(u,dx,d) = (Ax * grad dx, grad px) + (Hx*(dx-d), px), \forall px
#
#   d(t+dt) = min(1, max(d(t), (2 phi(u) + Hx(dx+1/beta) ) / (2 phi(u) + Hx)))
#
#   HHO:
#   sum_{T \in Th} (PK(huT, d_T), GkT(hvT))_T + stab(huT, hvT) = 0
#   sum_{T \in Th} (Ax * GkT(hdxT), GkT(hpxT))_T + stab(hdxT, hpxT) +
#      (H * (dx_T - dT), px_T)_T = 0
#
####################################################################################


mesh0 = code_aster.Mesh.buildSquare(refine=3)

mesh = CREA_MAILLAGE(MAILLAGE=mesh0,
                     MODI_HHO=_F(TOUT='OUI',),)

# define model
model_ther = AFFE_MODELE(MAILLAGE=mesh,
                         AFFE=_F(TOUT='OUI',
                                 MODELISATION='PLAN_HHO',
                                 FORMULATION='LINEAIRE',
                                 PHENOMENE='THERMIQUE')
                         )

model_meca = AFFE_MODELE(MAILLAGE=mesh,
                         AFFE=_F(TOUT='OUI',
                                 MODELISATION='D_PLAN_HHO',
                                 FORMULATION='LINEAIRE',
                                 PHENOMENE='MECANIQUE')
                         )

THER = {"MODELE": model_ther}
MECA = {"MODELE": model_meca}


# define BC
bc_ther = AFFE_CHAR_CINE(MODELE=model_ther,
                         THER_IMPO=_F(GROUP_MA=('RIGHT', 'LEFT',
                                                'TOP', 'BOTTOM',),  TEMP=0.0),
                         )

THER["EXCIT"] = { "CHARGE" : bc_ther}


bc_meca = AFFE_CHAR_CINE(MODELE=model_meca,
                         MECA_IMPO=(_F(GROUP_MA=('RIGHT'),  DX=0.0, DY=0.0),
                                    _F(GROUP_MA=("LEFT"), DY=0.2),),
                         )
MECA["EXCIT"] = { "CHARGE" : bc_meca}

# load

THER["SOURCE"] = 100.

# material

MECA["MATER"] = DEFI_MATERIAU(ELAS=_F(E=200000.,
                                      NU=0.3,
                                      ALPHA=1.0,)
                              )

A = 10
H = 2
THER["MATER"] = {"LAMBDA": A, "RHO_CP": H}

cs = CoupledSolver(MECA, THER)

sol = cs.solve()

u_ref = 1028.577291151862
dx_ref = 1175.794430261763
d_ref = 800.0000000000003

test.assertAlmostEqual((sol.u.norm("NORM_2") - u_ref)/u_ref, 0, delta=1e-6)
test.assertAlmostEqual((sol.dx.norm("NORM_2") - dx_ref)/dx_ref, 0, delta=1e-6)
test.assertAlmostEqual((sol.d.norm("NORM_2") - d_ref)/d_ref, 0, delta=1e-6)

test.printSummary()

FIN()
