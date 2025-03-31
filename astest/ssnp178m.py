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

import os.path as osp

from code_aster.Commands import *
from code_aster.Utilities import SharedTmpdir

POURSUITE(CODE="OUI", ERREUR=_F(ALARME="EXCEPTION"))

tmpdir = SharedTmpdir("ssnp178m_")

medfile = osp.join(tmpdir.path, "resu.ssnp178m.med")
DEFI_FICHIER(UNITE=87, FICHIER=medfile, TYPE="BINARY")
IMPR_RESU(
    FICHIER_UNIQUE="OUI",
    FORMAT="MED",
    UNITE=87,
    INFO=1,
    RESU=_F(RESULTAT=resnonl),
    VERSION_MED="4.1.0",
)

# read std field and project on hpc field
mesh_std = LIRE_MAILLAGE(UNITE=87)

model_std = AFFE_MODELE(
    AFFE=_F(MODELISATION=("D_PLAN_INCO_UPG",), PHENOMENE="MECANIQUE", TOUT="OUI"), MAILLAGE=mesh_std
)

fieldmat_std = AFFE_MATERIAU(AFFE=_F(MATER=(mater,), TOUT="OUI"), MAILLAGE=mesh_std)

depl_std = AFFE_CHAR_CINE(
    MECA_IMPO=_F(DX=0.0, DY=0.0, GROUP_MA=("Encast",)), MODELE=model_std, INFO=1
)

load_std = AFFE_CHAR_MECA(FORCE_CONTOUR=_F(FY=0.1125, GROUP_MA=("load",)), MODELE=model_std)

sol_std = LIRE_RESU(
    MODELE=model_std,
    FORMAT="MED",
    UNITE=87,
    TYPE_RESU="EVOL_NOLI",
    COMPORTEMENT=_F(DEFORMATION="PETIT", RELATION="VMIS_ISOT_LINE", TOUT="OUI"),
    CHAM_MATER=fieldmat_std,
    EXCIT=(
        _F(CHARGE=depl_std, TYPE_CHARGE="FIXE_CSTE"),
        _F(CHARGE=load_std, FONC_MULT=rampe, TYPE_CHARGE="FIXE_CSTE"),
    ),
    FORMAT_MED=(
        _F(NOM_RESU="resnonl", NOM_CHAM="DEPL"),
        _F(NOM_RESU="resnonl", NOM_CHAM="SIEF_ELGA"),
        _F(NOM_RESU="resnonl", NOM_CHAM="VARI_ELGA"),
    ),
    TOUT_ORDRE="OUI",
)

DEFI_FICHIER(ACTION="LIBERER", UNITE=87)

del tmpdir

sol_proj = PROJ_CHAMP(
    RESULTAT=sol_std,
    INFO=2,
    NOM_CHAM="DEPL",
    VIS_A_VIS=(_F(TOUT_1="OUI", TOUT_2="OUI"),),
    MODELE_1=model_std,
    MODELE_2=model,
    DISTANCE_ALARME=0.01,
)

TEST_RESU(
    RESU=_F(
        INST=1.0,
        GROUP_NO=("A",),
        REFERENCE="AUTRE_ASTER",
        RESULTAT=sol_proj,
        NOM_CHAM="DEPL",
        NOM_CMP="DY",
        VALE_CALC=2.24084161927,
        VALE_REFE=2.13651,
        CRITERE="RELATIF",
        PRECISION=5.0e-2,
    )
)

TEST_RESU(
    RESU=_F(
        INST=1.0,
        GROUP_NO=("A",),
        REFERENCE="AUTRE_ASTER",
        RESULTAT=sol_proj,
        NOM_CHAM="DEPL",
        NOM_CMP="DX",
        VALE_CALC=-1.69640391995,
        VALE_REFE=-1.61818,
        CRITERE="RELATIF",
        PRECISION=5.0e-2,
    )
)


FIN()
