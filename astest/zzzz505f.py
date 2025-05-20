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

DEBUT(CODE="OUI", ERREUR=_F(ALARME="EXCEPTION"), DEBUG=_F(SDVERI="OUI"))

test = CA.TestCase()

fmt_raison = (
    "-" * 80
    + """

   Exception interceptee
   Raison : %s

"""
    + "-" * 80
    + "\n"
)

E = 200.0e9

rho = 8000.0

nu = 0.3


mesh = LIRE_MAILLAGE(UNITE=20, FORMAT="MED")

model = AFFE_MODELE(
    MAILLAGE=mesh, AFFE=_F(GROUP_MA="VOL", PHENOMENE="MECANIQUE", MODELISATION="3D")
)

steel = DEFI_MATERIAU(ELAS=_F(E=E, NU=nu, RHO=rho))

CHMAT = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(GROUP_MA="VOL", MATER=steel))

clamp = AFFE_CHAR_MECA(
    MODELE=model, DDL_IMPO=_F(GROUP_MA="ENCAS", BLOCAGE=("DEPLACEMENT", "ROTATION"))
)

gravity1 = AFFE_CHAR_MECA(MODELE=model, PESANTEUR=_F(GRAVITE=100.0, DIRECTION=(-1.0, 0, 1)))
gravity2 = AFFE_CHAR_MECA(MODELE=model, PESANTEUR=_F(GRAVITE=200.0, DIRECTION=(-1.0, 0, 1)))
gravity3 = AFFE_CHAR_MECA(MODELE=model, PESANTEUR=_F(GRAVITE=300.0, DIRECTION=(-1.0, 0, 1)))

elasMult = MACRO_ELAS_MULT(
    MODELE=model,
    CHAM_MATER=CHMAT,
    CHAR_MECA_GLOBAL=clamp,
    CAS_CHARGE=(
        _F(NOM_CAS="grav1", CHAR_MECA=gravity1),
        _F(NOM_CAS="grav2", CHAR_MECA=gravity2),
        _F(NOM_CAS="grav3", CHAR_MECA=gravity3),
    ),
    SOLVEUR=_F(RESI_RELA=1.0e-5),
)

grav1Vale = 0.06942126997803875
grav2Vale = 0.1388425399560775
grav3Vale = 0.20826380993411645

TEST_RESU(
    RESU=(
        _F(
            GROUP_NO="P",
            RESULTAT=elasMult,
            NOM_CHAM="DEPL",
            NOM_CMP="DZ",
            NOM_CAS="grav1",
            VALE_CALC=grav1Vale,
        ),
        _F(
            GROUP_NO="P",
            RESULTAT=elasMult,
            NOM_CHAM="DEPL",
            NOM_CMP="DZ",
            NOM_CAS="grav2",
            VALE_CALC=grav2Vale,
        ),
        _F(
            GROUP_NO="P",
            RESULTAT=elasMult,
            NOM_CHAM="DEPL",
            NOM_CMP="DZ",
            NOM_CAS="grav3",
            VALE_CALC=grav3Vale,
        ),
    )
)

# Test interface Python de la SD résultat
nbIndexes = elasMult.getNumberOfIndexes()
firstIndex = elasMult.getFirstIndex()
lastIndex = elasMult.getLastIndex()
test.assertEqual(nbIndexes, 3)
test.assertEqual(firstIndex, 1)
test.assertEqual(lastIndex, 3)

# Check error message for wrong access parameter
is_ok = 0
try:
    pouet = elasMult.getTime(1)
except CA.AsterError as err:
    if err.id_message == "RESULT2_9":
        is_ok = 1
test.assertEqual(is_ok, 1)


allFieldNames = elasMult.getFieldsNames()
test.assertSequenceEqual(allFieldNames, ["DEPL", "SIEF_ELGA"])

para = elasMult.getAccessParameters()
allParameters = list(para.keys())
test.assertSequenceEqual(allParameters, ["NUME_ORDRE", "NOM_CAS"])

ranks_index = elasMult.getAccessParameters()["NUME_ORDRE"]
ranks_caseName = elasMult.getAccessParameters()["NOM_CAS"]
test.assertSequenceEqual(ranks_index, [1, 2, 3])
test.assertSequenceEqual(ranks_caseName, ["grav1", "grav2", "grav3"])

# Get displacements
disp1 = elasMult.getField("DEPL", value="grav1", para="NOM_CAS")
disp2 = elasMult.getField("DEPL", value="grav2", para="NOM_CAS")
disp3 = elasMult.getField("DEPL", value="grav3", para="NOM_CAS")

# Create new result
syntElasMult = CA.MultipleElasticResult()
syntElasMult.allocate(nbIndexes)

syntElasMult.setField(disp1, "DEPL", value="grav1", para="NOM_CAS")
syntElasMult.setField(disp2, "DEPL", value="grav2", para="NOM_CAS")
syntElasMult.setField(disp3, "DEPL", value="grav3", para="NOM_CAS")
syntElasMult.setField(disp1, "DEPL", value="grav4", para="NOM_CAS")


TEST_RESU(
    RESU=(
        _F(
            GROUP_NO="P",
            RESULTAT=syntElasMult,
            NOM_CHAM="DEPL",
            NOM_CMP="DZ",
            NOM_CAS="grav1",
            VALE_CALC=grav1Vale,
        ),
        _F(
            GROUP_NO="P",
            RESULTAT=syntElasMult,
            NOM_CHAM="DEPL",
            NOM_CMP="DZ",
            NOM_CAS="grav2",
            VALE_CALC=grav2Vale,
        ),
        _F(
            GROUP_NO="P",
            RESULTAT=syntElasMult,
            NOM_CHAM="DEPL",
            NOM_CMP="DZ",
            NOM_CAS="grav3",
            VALE_CALC=grav3Vale,
        ),
        _F(
            GROUP_NO="P",
            RESULTAT=syntElasMult,
            NOM_CHAM="DEPL",
            NOM_CMP="DZ",
            NOM_CAS="grav4",
            VALE_CALC=grav1Vale,
        ),
    )
)

syntElasMult.setField(disp2, "DEPL", value="grav4", para="NOM_CAS")

TEST_RESU(
    RESU=(
        _F(
            GROUP_NO="P",
            RESULTAT=syntElasMult,
            NOM_CHAM="DEPL",
            NOM_CMP="DZ",
            NOM_CAS="grav4",
            VALE_CALC=grav2Vale,
        ),
    )
)

nbIndexes = syntElasMult.getNumberOfIndexes()
firstIndex = syntElasMult.getFirstIndex()
lastIndex = syntElasMult.getLastIndex()

test.assertEqual(nbIndexes, 4)
test.assertEqual(firstIndex, 1)
test.assertEqual(lastIndex, 4)

FIN(INFO_RESU="OUI")
