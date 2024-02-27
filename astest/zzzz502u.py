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

import os

from code_aster.Commands import *
from code_aster import CA
from code_aster.CA import MPI

rank = MPI.ASTER_COMM_WORLD.Get_rank()

current_dir = os.getcwd()

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

mesh = CA.ParallelMesh()

mesh.readMedFile("./zzzz502u.mmed")

# Definition du modele
model = AFFE_MODELE(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))

# Definition du materiau
mater = DEFI_MATERIAU(ELAS=_F(E=10e9, NU=0.3))

# Affectation du materiau sur le maillage
affectMat = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=mater))

# construction de l AFFE_CHAR_MECA
nb_noeuds = [3, 3, 3, 3, 5]
liaison = []
group_no = []
for i in range(len(nb_noeuds)):
    nb = nb_noeuds[i]
    nnode = nb
    coef = 1.0 / (nnode - 1)
    for ddl in ["X", "Y", "Z"]:
        group_no.append(("Liaison_" + str(i + 1)).strip())
        liaison.append(
            _F(
                GROUP_NO=("Liaison_" + str(i + 1)).strip(),
                DDL=["D" + ddl] * nnode,
                COEF_MULT=[coef] * (nnode - 1) + [-1],
                COEF_IMPO=0.0,
            )
        )

Liaison_ddl = AFFE_CHAR_MECA(MODELE=model, LIAISON_DDL=liaison, INFO=2)

# Definition des conditions aux limites
clim = AFFE_CHAR_CINE(
    MODELE=model,
    MECA_IMPO=(
        _F(GROUP_NO="Gauche", DX=0.0),
        _F(GROUP_NO="Gauche", DY=0.0),
        _F(GROUP_NO="Gauche", DZ=0.0),
        _F(GROUP_NO="Droite", DX=1.0),
    ),
)

COEF3 = DEFI_FONCTION(NOM_PARA="INST", PROL_DROITE="CONSTANT", VALE=(0.0, 1.0))

COEF0 = DEFI_FONCTION(NOM_PARA="INST", PROL_DROITE="CONSTANT", VALE=(1.0, 1.0))

# time increment
TEMPS = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=1.0, NOMBRE=1.0))

resu = STAT_NON_LINE(
    MODELE=model,
    CHAM_MATER=affectMat,
    EXCIT=(_F(CHARGE=clim, FONC_MULT=COEF3), _F(CHARGE=Liaison_ddl, FONC_MULT=COEF0)),
    COMPORTEMENT=_F(RELATION="ELAS"),
    NEWTON=_F(MATRICE="TANGENTE", PREDICTION="ELASTIQUE", REAC_ITER=1),
    CONVERGENCE=_F(RESI_GLOB_RELA=1.0e-6, ITER_GLOB_MAXI=30, ARRET="OUI"),
    SOLVEUR=_F(METHODE="MUMPS", NPREC=8),
    # SOLVEUR=_F(METHODE='PETSC', PRE_COND='JACOBI', ),
    INCREMENT=_F(LIST_INST=TEMPS),
)

depl = resu.getField("DEPL", 1)

if rank == 0:
    test.assertTrue(abs(depl[0] - 6.56591e-01) < 1e-6)
    test.assertTrue(abs(depl[15] - 1e0) < 1e-10)
else:
    test.assertTrue(abs(depl[0] - 1e0) < 1e-10)
    test.assertTrue(abs(depl[15] - 6.56591e-01) < 1e-6)

FIN()
