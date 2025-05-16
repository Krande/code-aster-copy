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

DEBUT(CODE="OUI", ERREUR=_F(ALARME="EXCEPTION"), DEBUG=_F(SDVERI="OUI"), INFO=1)

test = CA.TestCase()

mesh = LIRE_MAILLAGE(FORMAT="MED", UNITE=20)

model = AFFE_MODELE(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))
# Very high elasticity limit to simulate elasticity
acier = DEFI_MATERIAU(ELAS=_F(E=200000.0, NU=0.3), ECRO_LINE=_F(D_SIGM_EPSI=2000.0, SY=200000.0))

mater = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=acier))

# Build reference field with CREA_CHAMP

value = 2.0

refe1 = CREA_CHAMP(
    TYPE_CHAM="ELGA_SIEF_R",
    OPERATION="AFFE",
    MODELE=model,
    AFFE=_F(
        TOUT="OUI",
        NOM_CMP=("SIXX", "SIYY", "SIZZ", "SIXY", "SIYZ", "SIXZ"),
        VALE=(value, value, value, value, value, value),
    ),
)
refe2 = CREA_CHAMP(
    OPERATION="AFFE",
    TYPE_CHAM="ELGA_VARI_R",
    MODELE=model,
    PROL_ZERO="OUI",
    AFFE=_F(TOUT="OUI", NOM_CMP=("V1", "V2"), VALE=(value, value)),
)

# Using Python binding for behaviour
study = CA.PhysicalProblem(model, mater)

# With default values: no initial state, no implex and info=1
study.computeBehaviourProperty(COMPORTEMENT=(_F(RELATION="VMIS_ISOT_LINE", TOUT="OUI"),))
behav = study.getBehaviourProperty()
# Build testfield with model and compor
fieldSIEF = CA.FieldOnCellsReal(model, "ELGA", "SIEF_R")
fieldSIEF.setValues(value)

fieldVARI = CA.FieldOnCellsReal(model, "ELGA", "VARI_R", behav)
fieldVARI.setValues(value)


# Test
test.assertAlmostEqual(len(refe1.getValues()), len(fieldSIEF.getValues()))
test.assertAlmostEqual(len(refe2.getValues()), len(fieldVARI.getValues()))
test.assertAlmostEqual(refe1.getValues(), fieldSIEF.getValues())
test.assertAlmostEqual(refe2.getValues(), fieldVARI.getValues())
test.assertEqual(fieldSIEF.getLocalization(), "ELGA")
test.assertEqual(fieldSIEF.getNumberOfComponents(), 6)
test.assertEqual(fieldVARI.getNumberOfComponents(), 2)
test.assertSequenceEqual(
    ["SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ"], fieldSIEF.getComponents()
)
test.assertSequenceEqual(["V1", "V2"], fieldVARI.getComponents())

fno = fieldSIEF.toFieldOnNodes()
fsno = fieldSIEF.toSimpleFieldOnNodes()
felno = fieldSIEF.asLocalization("ELNO")

test.assertEqual(felno.getLocalization(), "ELNO")
test.assertSequenceEqual(["SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ"], felno.getComponents())


fs0 = fieldSIEF.toSimpleFieldOnCells()
fneuts = fs0.asPhysicalQuantity("DEPL_R", {"SIXX": "DZ", "SIXZ": "DX"})
test.assertSequenceEqual(["DZ", "DX"], fneuts.getComponents())
test.assertEqual(fneuts.getPhysicalQuantity(), "DEPL_R")
test.assertEqual(fneuts.getLocalization(), "ELGA")

fneut = fieldSIEF.asPhysicalQuantity("DEPL_R", {"SIXX": "DZ", "SIXZ": "DX"})
test.assertEqual(fneut.getPhysicalQuantity(), "DEPL_R")
test.assertEqual(fneut.getLocalization(), "ELGA")

sf2 = fieldVARI.toSimpleFieldOnCells()
test.assertSequenceEqual(["V1", "V2"], sf2.getComponents())
test.assertEqual(sf2.getPhysicalQuantity(), "VARI_R")
test.assertEqual(sf2.getLocalization(), "ELGA")
test.assertFalse(sf2.hasValue(200, 1, 4800))
test.assertTrue(sf2.hasValue(200, 1, 1))
test.assertTrue(len(sf2.getValuesOnCell(200)[0]) == sf2.getNumberOfPointsOfCell(200))

sfr2 = sf2.restrict(["V1"])
test.assertSequenceEqual(["V1"], sfr2.getComponents())
test.assertEqual(sfr2.getPhysicalQuantity(), "VARI_R")
test.assertEqual(sfr2.getLocalization(), "ELGA")
test.assertAlmostEqual(sf2.getValue(200, 0, 0), sfr2.getValue(200, 0, 0))

sfr3 = sf2.asPhysicalQuantity("NEUT_R", {"V1": "X2"})
test.assertSequenceEqual(["X2"], sfr3.getComponents())
test.assertEqual(sfr3.getPhysicalQuantity(), "NEUT_R")
test.assertEqual(sfr3.getLocalization(), "ELGA")
test.assertAlmostEqual(sfr2.getValue(200, 0, 0), sfr3.getValue(200, 0, 0))

fn3 = fno.asPhysicalQuantity("TEMP_R", {"SIYY": "TEMP"})
test.assertSequenceEqual(["TEMP"], fn3.getComponents())
test.assertEqual(fn3.getPhysicalQuantity(), "TEMP_R")
test.assertEqual(fn3.getLocalization(), "NOEU")

# Test constructeur avec le caraelem


MAIL = LIRE_MAILLAGE(UNITE=18, FORMAT="MED")

MODELE = AFFE_MODELE(
    MAILLAGE=MAIL, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="TUYAU_3M")
)

CAREL = AFFE_CARA_ELEM(
    MODELE=MODELE,
    POUTRE=(
        _F(
            SECTION="CERCLE",
            GROUP_MA="TOUT",
            CARA=("R", "EP"),
            VALE=(0.10, 0.05),
            TUYAU_NSEC=8,
            TUYAU_NCOU=3,
        ),
    ),
)


MAT = DEFI_MATERIAU(ELAS=_F(E=2.1e11, NU=0.3, RHO=7800.0))

CHMAT = AFFE_MATERIAU(MAILLAGE=MAIL, AFFE=_F(TOUT="OUI", MATER=MAT))

value = 10

refe1 = CREA_CHAMP(
    TYPE_CHAM="ELGA_SIEF_R",
    OPERATION="AFFE",
    AFFE_SP=_F(CARA_ELEM=CAREL),
    MODELE=MODELE,
    AFFE=_F(
        TOUT="OUI",
        NOM_CMP=("SIXX", "SIYY", "SIZZ", "SIXY", "SIYZ", "SIXZ"),
        VALE=(value, value, value, value, value, value),
    ),
)


study = CA.PhysicalProblem(MODELE, CHMAT, CAREL)

study.computeBehaviourProperty(COMPORTEMENT=(_F(RELATION="VMIS_ISOT_LINE", TOUT="OUI"),))
behav = study.getBehaviourProperty()
fieldSIEF = CA.FieldOnCellsReal(MODELE, "ELGA", "SIEF_R", behav, CAREL)

# TEST LENGTH EQUALITY
test.assertAlmostEqual(len(refe1.getValues()), len(fieldSIEF.getValues()))

test.printSummary()
CA.close()
