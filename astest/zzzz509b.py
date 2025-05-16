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

DEBUT(CODE="OUI", DEBUG=_F(SDVERI="OUI"), INFO=1)

test = CA.TestCase()

mesh = LIRE_MAILLAGE(UNITE=20, FORMAT="MED")


mo = AFFE_MODELE(
    MAILLAGE=mesh,
    AFFE=(
        _F(GROUP_MA=("cube", "xsup"), MODELISATION="3D", PHENOMENE="MECANIQUE"),
        _F(GROUP_MA="barre", MODELISATION="BARRE", PHENOMENE="MECANIQUE"),
        _F(GROUP_MA="grille", MODELISATION="GRILLE_MEMBRANE", PHENOMENE="MECANIQUE"),
    ),
)


cael = AFFE_CARA_ELEM(
    MODELE=mo,
    BARRE=_F(GROUP_MA="barre", SECTION="GENERALE", CARA="A", VALE=100.0),
    GRILLE=_F(GROUP_MA="grille", SECTION=2.0, VECT_1=(1.0, 0.0, 0.0)),
)


young_b = DEFI_FONCTION(NOM_PARA="TEMP", ABSCISSE=(10, 110), ORDONNEE=(30.0e3, 25.0e3))


young_a = DEFI_FONCTION(NOM_PARA="TEMP", ABSCISSE=(10, 110), ORDONNEE=(200.0e3, 200.0e3))


nu = DEFI_CONSTANTE(VALE=0.0)


alpha = DEFI_CONSTANTE(VALE=1.0e-5)


beton = DEFI_MATERIAU(
    ELAS_FO=_F(E=young_b, NU=nu, ALPHA=alpha, TEMP_DEF_ALPHA=10),
    ECRO_NL=_F(R0=1.0e30),
    VISC_ELAS=_F(K=50.0e3, TAU=1.0),
)


acier = DEFI_MATERIAU(
    ELAS_FO=_F(E=young_a, NU=nu, ALPHA=alpha, TEMP_DEF_ALPHA=10),
    ECRO_LINE=_F(SY=1.0e30, D_SIGM_EPSI=0),
)


bc = AFFE_CHAR_MECA(
    MODELE=mo,
    FACE_IMPO=(_F(GROUP_MA="xinf", DX=0), _F(GROUP_MA="yinf", DY=0), _F(GROUP_MA="zinf", DZ=0)),
)


rampe = DEFI_FONCTION(NOM_PARA="INST", ABSCISSE=(0, 1), ORDONNEE=(0, 1))


load = AFFE_CHAR_MECA(MODELE=mo, FORCE_FACE=_F(GROUP_MA="xsup", FX=4.6))


# --------------------------------------------------------------------------------------------------
# Temperature applique au liner
# --------------------------------------------------------------------------------------------------

temp_list = [
    CREA_CHAMP(
        OPERATION="AFFE",
        TYPE_CHAM="NOEU_TEMP_R",
        MODELE=mo,
        AFFE=(_F(GROUP_MA=("cube", "barre", "grille"), NOM_CMP="TEMP", VALE=temp),),
    )
    for temp in (10, 110)
]


temp = CREA_RESU(
    OPERATION="AFFE",
    TYPE_RESU="EVOL_THER",
    AFFE=[
        _F(NOM_CHAM="TEMP", INST=t, CHAM_GD=temp_field)
        for (t, temp_field) in zip((0, 1), temp_list)
    ],
)


mat = AFFE_MATERIAU(
    MODELE=mo,
    AFFE=(
        _F(GROUP_MA="cube", MATER=beton),
        _F(GROUP_MA="barre", MATER=acier),
        _F(GROUP_MA="grille", MATER=acier),
    ),
    AFFE_VARC=_F(
        GROUP_MA=("cube", "barre", "grille"),
        NOM_VARC="TEMP",
        EVOL=temp,
        NOM_CHAM="TEMP",
        VALE_REF=10,
    ),
)

stepping = DEFI_LIST_REEL(VALE=(0, 0.5, 1))

evolRefe = STAT_NON_LINE(
    MODELE=mo,
    CHAM_MATER=mat,
    CARA_ELEM=cael,
    EXCIT=(_F(CHARGE=bc), _F(CHARGE=load, FONC_MULT=rampe)),
    COMPORTEMENT=(
        _F(GROUP_MA="cube", RELATION="VMIS_ISOT_NL", REGU_VISC="OUI"),
        _F(GROUP_MA="grille", RELATION="GRILLE_ISOT_LINE"),
        _F(GROUP_MA="barre", RELATION="VMIS_ISOT_LINE"),
    ),
    INCREMENT=_F(LIST_INST=stepping),
    NEWTON=_F(MATRICE="TANGENTE", REAC_ITER=1),
    CONVERGENCE=_F(ITER_GLOB_MAXI=50, RESI_GLOB_MAXI=3.0e2 * 1.0e-3),
    MESURE=_F(TABLE="OUI"),
)


evol = MECA_NON_LINE(
    MODELE=mo,
    CHAM_MATER=mat,
    CARA_ELEM=cael,
    EXCIT=(_F(CHARGE=bc), _F(CHARGE=load, FONC_MULT=rampe)),
    COMPORTEMENT=(
        _F(GROUP_MA="cube", RELATION="VMIS_ISOT_NL", REGU_VISC="OUI"),
        _F(GROUP_MA="grille", RELATION="GRILLE_ISOT_LINE"),
        _F(GROUP_MA="barre", RELATION="VMIS_ISOT_LINE"),
    ),
    INCREMENT=_F(LIST_INST=stepping),
    NEWTON=_F(MATRICE="TANGENTE", REAC_ITER=1),
    CONVERGENCE=_F(ITER_GLOB_MAXI=50, RESI_GLOB_MAXI=3.0e2 * 1.0e-3),
    INFO=1,
)

fieldRefe = evol.getField("DEPL", 2)
field = evol.getField("DEPL", 2)

test.assertAlmostEqual(fieldRefe.norm("NORM_2"), field.norm("NORM_2"), 6)

FIN()
