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

import os
import code_aster
from code_aster.Commands import *

current_dir = os.getcwd()

code_aster.init("--test")

test = code_aster.TestCase()

mesh = code_aster.Mesh()

mesh.readMedFile('./zzzz502u.mmed')

# Definition du modele
model=AFFE_MODELE(MAILLAGE=mesh,
                 AFFE=_F(TOUT='OUI',
                         PHENOMENE='MECANIQUE',
                         MODELISATION='3D',),);

# Definition du materiau
mater=DEFI_MATERIAU(ELAS=_F(E=10e9,
                            NU=0.3,),);

# Affectation du materiau sur le maillage
affectMat=AFFE_MATERIAU(MAILLAGE=mesh,
                    AFFE=_F(TOUT='OUI',
                            MATER=mater,),);

# construction de l AFFE_CHAR_MECA
nb_noeuds = [3, 3, 3, 3, 5]
liaison  = []
group_no = []
for i in range(len(nb_noeuds)):
    nb = nb_noeuds[i]
    nnode = nb
    coef  = 1./(nnode-1)
    for ddl in ['X', 'Y', 'Z'] :
        group_no.append(('Liaison_' + str(i+1)).strip())
        liaison.append(_F(
                GROUP_NO  = ('Liaison_' + str(i+1)).strip(),
                DDL       = ['D' + ddl]*nnode,
                COEF_MULT = [coef]*(nnode-1) + [-1],
                COEF_IMPO = 0.,
                )
            )

Liaison_ddl = AFFE_CHAR_MECA(
    MODELE     = model,
    LIAISON_DDL= liaison,
    INFO=2
    )

# Definition des conditions aux limites
clim = AFFE_CHAR_CINE(
    MODELE=model,
    MECA_IMPO=(
        _F(GROUP_NO='Gauche', DX=0., ),
        _F(GROUP_NO='Gauche', DY=0., ),
        _F(GROUP_NO='Gauche', DZ=0., ),
        _F(GROUP_NO='Droite', DX=1., ),
    ),
)

COEF3 = DEFI_FONCTION(NOM_PARA='INST',
                      PROL_DROITE='CONSTANT',
                      VALE=(0., 1.,), )

COEF0 = DEFI_FONCTION(NOM_PARA='INST',
                      PROL_DROITE='CONSTANT',
                      VALE=(1., 1.,), )

# time increment
TEMPS = DEFI_LIST_REEL(
    DEBUT      = 0.,
    INTERVALLE = _F(JUSQU_A=1., NOMBRE=1.),
    )

resu = STAT_NON_LINE(MODELE=model,
                   CHAM_MATER=affectMat,
                   EXCIT=(_F(CHARGE=clim, FONC_MULT=COEF3,),
                          _F(CHARGE=Liaison_ddl, FONC_MULT=COEF0,),
                          ),
                   COMPORTEMENT=_F(RELATION='ELAS', ),
                   NEWTON=_F(MATRICE='TANGENTE',
                             PREDICTION='ELASTIQUE',
                             REAC_ITER=1, ),
                   CONVERGENCE=_F(RESI_GLOB_RELA=1.E-6,
                                  ITER_GLOB_MAXI=30,
                                  ARRET='OUI', ),
                   SOLVEUR=_F(METHODE='MUMPS', NPREC=8, ),
                   #SOLVEUR=_F(METHODE='PETSC', PRE_COND='JACOBI', ),
                   INCREMENT=_F(LIST_INST=TEMPS, ), )



depl = resu.getFieldOnNodesReal('DEPL', 1)

mpi_value = depl.norm("NORM_1")

seq_value = 17.841069929231153
test.assertTrue( abs(seq_value - mpi_value)/abs(seq_value) < 1e-6)

FIN();
