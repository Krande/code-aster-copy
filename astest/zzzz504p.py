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
from code_aster.CA import MPI

DEBUT(CODE="OUI", ERREUR=_F(ALARME="EXCEPTION"))

nProc = MPI.ASTER_COMM_WORLD.Get_size()
parallel = nProc > 1
rank = MPI.ASTER_COMM_WORLD.Get_rank()

if parallel:
    MAIL = CA.ParallelMesh()
    MAIL.readMedFile("zzzz504p/%d.med" % rank, partitioned=True)

model = AFFE_MODELE(
    MAILLAGE=MAIL,
    AFFE=(
        _F(MODELISATION="3D", PHENOMENE="MECANIQUE", TOUT="OUI"),
        _F(MODELISATION="DKT", PHENOMENE="MECANIQUE", GROUP_MA="Surf"),
    ),
    VERI_PLAN="OUI",
)

CARAMECA = AFFE_CARA_ELEM(MODELE=model, COQUE=_F(GROUP_MA="Surf", EPAIS=0.1))

MATER1 = DEFI_MATERIAU(ELAS=_F(E=200000.0, NU=0.3, RHO=1000.0))

AFFMAT = AFFE_MATERIAU(MAILLAGE=MAIL, AFFE=_F(TOUT="OUI", MATER=MATER1))

load = AFFE_CHAR_CINE(
    MODELE=model,
    MECA_IMPO=(
        _F(GROUP_MA="Encast", DX=0.0),
        _F(GROUP_MA="Encast", DY=0.0),
        _F(GROUP_MA="Encast", DZ=0.0),
    ),
)

load2 = AFFE_CHAR_MECA(
    MODELE=model, PESANTEUR=_F(GRAVITE=1.0, DIRECTION=(0.0, -1.0, 0.0)), INFO=1, VERI_NORM="NON"
)

resu = MECA_STATIQUE(
    MODELE=model,
    CHAM_MATER=AFFMAT,
    CARA_ELEM=CARAMECA,
    EXCIT=(_F(CHARGE=load), _F(CHARGE=load2)),
    SOLVEUR=_F(METHODE="MUMPS"),
    INFO=2,
)

TEST_RESU(
    RESU=(
        _F(
            GROUP_MA=("Surf",),
            NOM_CHAM="SIEF_ELGA",
            NOM_CMP="SIYY",
            POINT=1,
            SOUS_POINT=2,
            NUME_ORDRE=1,
            RESULTAT=resu,
            CRITERE="RELATIF",
            VALE_CALC=-405.5059907903177,
        ),
    )
)


FIN()
