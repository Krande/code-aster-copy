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

DEBUT(CODE="OUI")

test = CA.TestCase()

POUTRE0 = LIRE_MAILLAGE(FORMAT="MED")

DEFI_GROUP(reuse=POUTRE0, MAILLAGE=POUTRE0, CREA_GROUP_NO=_F(TOUT_GROUP_MA="OUI"))

MODELE = AFFE_MODELE(
    MAILLAGE=POUTRE0,
    AFFE=(
        _F(GROUP_MA=("Ligne1", "Ligne2"), PHENOMENE="MECANIQUE", MODELISATION="POU_D_T"),
        _F(GROUP_MA=("Point1", "Point2"), PHENOMENE="MECANIQUE", MODELISATION="DIS_TR"),
    ),
)

#
# ------------------------
# CARACTERISTIQUES POUTRE
# ------------------------
#

CARA = AFFE_CARA_ELEM(
    MODELE=MODELE,
    POUTRE=(
        _F(GROUP_MA="Ligne1", SECTION="CERCLE", CARA=("R", "EP"), VALE=(6.50e-03, 2.10e-03)),
        _F(GROUP_MA="Ligne2", SECTION="CERCLE", CARA=("R", "EP"), VALE=(6.50e-03, 2.20e-03)),
    ),
    DISCRET=(
        _F(GROUP_MA=("Point1", "Point2"), CARA="K_TR_D_N", VALE=(0.0, 0.0, 0.0, 6.29, 0.0, 6.29)),
        _F(
            GROUP_MA=("Point1", "Point2"),
            CARA="M_TR_D_N",
            VALE=(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        ),
    ),
)

#
# --------------------------
# CARACTERISTIQUES MATERIAU
# --------------------------
#

MAT = DEFI_MATERIAU(ELAS=_F(E=2.80e9, NU=0.3, RHO=1500.0))

#
# ---------------------
# AFFECTATION MATERIAU
# ---------------------
#

CHMAT = AFFE_MATERIAU(MAILLAGE=POUTRE0, AFFE=_F(TOUT="OUI", MATER=MAT))

#
# -------------------------
# CONDITIONS AUX LIMITES
# POUTRE ROTULEE - ROTULEE
# -------------------------
#

CHARG1 = AFFE_CHAR_MECA(
    MODELE=MODELE, DDL_IMPO=_F(GROUP_MA="Point1", DX=0.0, DY=0.0, DZ=0.0, DRX=0, DRY=0, DRZ=0)
)

CHARG2 = AFFE_CHAR_MECA(MODELE=MODELE, FORCE_POUTRE=_F(GROUP_MA=("Ligne1", "Ligne2"), FY=1.0))

RESU = MECA_STATIQUE(
    MODELE=MODELE, CHAM_MATER=CHMAT, CARA_ELEM=CARA, EXCIT=(_F(CHARGE=CHARG1), _F(CHARGE=CHARG2))
)

TEST_RESU(
    RESU=_F(
        CRITERE="ABSOLU",
        GROUP_NO="Point2",
        NOM_CHAM="DEPL",
        NOM_CMP="DY",
        NUME_ORDRE=1,
        PRECISION=1.0e-6,
        REFERENCE="AUTRE_ASTER",
        RESULTAT=RESU,
        VALE_CALC=2.227968e-02,
        VALE_REFE=2.227968e-02,
    )
)


test.assertTrue(True)
test.printSummary()

FIN()
