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

test = CA.TestCase()

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

nProc = CA.MPI.ASTER_COMM_WORLD.Get_size()
rank = CA.MPI.ASTER_COMM_WORLD.Get_rank()


pMesh = LIRE_MAILLAGE(UNITE=20, FORMAT="MED", PARTITIONNEUR="PTSCOTCH", INFO_MED=1)
pMesh = MODI_MAILLAGE(
    reuse=pMesh, MAILLAGE=pMesh, ORIE_PEAU=_F(GROUP_MA_PEAU=("HAUT", "BAS", "DROITE", "GAUCHE"))
)

cellNum = pMesh.getLocalToGlobalCellIds()
cellNumRef = [[0, 1, 2, 6, 7, 4, 9], [5, 7, 8, 2, 3, 9, 10]]
test.assertSequenceEqual(cellNum, cellNumRef[rank])


pmodel = AFFE_MODELE(
    MAILLAGE=pMesh,
    AFFE=(
        _F(MODELISATION="D_PLAN", PHENOMENE="MECANIQUE", TOUT="OUI"),
        _F(MODELISATION="D_PLAN_INCO_UPG", PHENOMENE="MECANIQUE", GROUP_MA="F3"),
    ),
)


MA = DEFI_MATERIAU(ELAS=_F(E=210000e6, NU=0.3))

pMATE = AFFE_MATERIAU(MAILLAGE=pMesh, AFFE=_F(TOUT="OUI", MATER=MA))

pload0 = AFFE_CHAR_MECA(MODELE=pmodel, DDL_IMPO=(_F(GROUP_MA="GAUCHE", DX=-1.0),))

pload1 = AFFE_CHAR_MECA(
    MODELE=pmodel,
    LIAISON_GROUP=(
        _F(
            GROUP_MA_1="GAUCHE",
            DDL_1="DX",
            GROUP_MA_2="DROITE",
            DDL_2="DX",
            COEF_MULT_1=1.0,
            COEF_MULT_2=1.0,
            COEF_IMPO=0.0,
        ),
    ),
)

pload2 = AFFE_CHAR_MECA(MODELE=pmodel, PRES_REP=(_F(GROUP_MA="HAUT", PRES=1.0),))


pload3 = AFFE_CHAR_CINE(MODELE=pmodel, MECA_IMPO=(_F(GROUP_MA="BAS", DY=0.0),))


pRESU = MECA_STATIQUE(
    MODELE=pmodel,
    CHAM_MATER=pMATE,
    EXCIT=(_F(CHARGE=pload0), _F(CHARGE=pload1), _F(CHARGE=pload2), _F(CHARGE=pload3)),
    INST=1.0,
    SOLVEUR=_F(METHODE="PETSC", RESI_RELA=1e-9),
)

TEST_RESU(
    RESU=(
        _F(
            GROUP_NO=("GN0",),
            NOM_CHAM="DEPL",
            NOM_CMP="DX",
            INST=1.0,
            PRECISION=0.00001,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=pRESU,
            CRITERE="ABSOLU",
            VALE_CALC=-0.33333333333333204,
            VALE_REFE=-0.33333333333333204,
        ),
        _F(
            GROUP_NO=("GN0",),
            NOM_CHAM="DEPL",
            NOM_CMP="DY",
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=pRESU,
            ORDRE_GRANDEUR=1e-6,
            CRITERE="ABSOLU",
            VALE_CALC=2.6988334396258498e-17,
            VALE_REFE=0.0,
        ),
    )
)


TEST_RESU(
    RESU=(
        _F(
            GROUP_NO=("GN1",),
            NOM_CHAM="DEPL",
            NOM_CMP="DX",
            INST=1.0,
            PRECISION=0.00001,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=pRESU,
            CRITERE="RELATIF",
            VALE_CALC=0.33333333333333515,
            VALE_REFE=0.33333333333333515,
        ),
        _F(
            CRITERE=("RELATIF",),
            GROUP_NO=("GN1",),
            NOM_CHAM="DEPL",
            NOM_CMP="PRES",
            INST=1.0,
            PRECISION=0.00001,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=pRESU,
            VALE_CALC=19999999999.380936,
            VALE_REFE=19999999999.380936,
        ),
    )
)


RAMPE = DEFI_FONCTION(NOM_PARA="INST", VALE=(0.0, 0.0, 1.0, 1.0))

LINSTC = DEFI_LIST_REEL(VALE=(0.0, 0.5, 0.75, 1.0))

nRESU = STAT_NON_LINE(
    MODELE=pmodel,
    CHAM_MATER=pMATE,
    EXCIT=(
        _F(CHARGE=pload0, FONC_MULT=RAMPE),
        _F(CHARGE=pload1, FONC_MULT=RAMPE),
        _F(CHARGE=pload2, FONC_MULT=RAMPE),
        _F(CHARGE=pload3, FONC_MULT=RAMPE),
    ),
    INCREMENT=_F(LIST_INST=LINSTC),
    SOLVEUR=_F(METHODE="PETSC", RESI_RELA=1e-9),
)

TEST_RESU(
    RESU=(
        _F(
            GROUP_NO=("GN0",),
            NOM_CHAM="DEPL",
            NOM_CMP="DX",
            INST=1.0,
            PRECISION=0.00001,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=nRESU,
            CRITERE="ABSOLU",
            VALE_CALC=-0.33333333333333204,
            VALE_REFE=-0.33333333333333204,
        ),
        _F(
            GROUP_NO=("GN0",),
            NOM_CHAM="DEPL",
            NOM_CMP="DY",
            INST=1.0,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=nRESU,
            ORDRE_GRANDEUR=1e-6,
            CRITERE="ABSOLU",
            VALE_CALC=2.6988334396258498e-17,
            VALE_REFE=0.0,
        ),
    )
)


TEST_RESU(
    RESU=(
        _F(
            GROUP_NO=("GN1",),
            NOM_CHAM="DEPL",
            NOM_CMP="DX",
            INST=1.0,
            PRECISION=0.00001,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=nRESU,
            CRITERE="RELATIF",
            VALE_CALC=0.33333333333333515,
            VALE_REFE=0.33333333333333515,
        ),
        _F(
            CRITERE=("RELATIF",),
            GROUP_NO=("GN1",),
            NOM_CHAM="DEPL",
            NOM_CMP="PRES",
            INST=1.0,
            PRECISION=0.00001,
            REFERENCE="AUTRE_ASTER",
            RESULTAT=nRESU,
            VALE_CALC=19999999999.380936,
            VALE_REFE=19999999999.380936,
        ),
    )
)


test.printSummary()

FIN()
