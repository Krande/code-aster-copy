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

DEBUT(CODE="OUI", DEBUG=_F(SDVERI="OUI"))

# From hsnv120a

mesh = LIRE_MAILLAGE(FORMAT="ASTER")

model = AFFE_MODELE(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))

COU_TRAC = DEFI_FONCTION(
    NOM_PARA="EPSI", PROL_DROITE="LINEAIRE", VALE=(0.005, 1000.0, 1.005, 3000.0)
)

ACIER_L = DEFI_MATERIAU(
    ELAS=_F(E=200000.0, NU=0.3, ALPHA=1.0e-4), ECRO_LINE=_F(D_SIGM_EPSI=2000.0, SY=1000.0)
)

ACIER_T = DEFI_MATERIAU(ELAS=_F(E=200000.0, NU=0.3, ALPHA=1.0e-4), TRACTION=_F(SIGM=COU_TRAC))

L_INST = DEFI_LIST_REEL(
    DEBUT=0.0,
    INTERVALLE=(
        _F(JUSQU_A=1.0, NOMBRE=1),
        _F(JUSQU_A=2.0, NOMBRE=1),
        _F(JUSQU_A=2.5, NOMBRE=5),
        _F(JUSQU_A=3.0, NOMBRE=1),
    ),
)

F_CHAR = DEFI_FONCTION(
    NOM_PARA="INST", PROL_DROITE="CONSTANT", VALE=(0.0, 0.0, 1.0, 0.0, 2.0, 1.0, 3.0, 0.0)
)

F_TEMP = DEFI_FONCTION(
    NOM_PARA="INST",
    NOM_RESU="TEMP",
    PROL_DROITE="CONSTANT",
    VALE=(0.0, 20.0, 1.0, 120.0, 2.0, 120.0, 3.0, 20.0),
)

CHP_TEMP = CREA_CHAMP(
    OPERATION="AFFE",
    TYPE_CHAM="NOEU_TEMP_F",
    MAILLAGE=mesh,
    AFFE=_F(TOUT="OUI", NOM_CMP="TEMP", VALE_F=F_TEMP),
)

TEMP = CREA_RESU(
    OPERATION="AFFE",
    TYPE_RESU="EVOL_THER",
    AFFE=_F(NOM_CHAM="TEMP", LIST_INST=L_INST, CHAM_GD=CHP_TEMP),
)


materialFieldTrac = AFFE_MATERIAU(
    MAILLAGE=mesh,
    AFFE=_F(TOUT="OUI", MATER=ACIER_T),
    AFFE_VARC=_F(TOUT="OUI", EVOL=TEMP, NOM_VARC="TEMP", NOM_CHAM="TEMP", VALE_REF=20.0),
)


CHR_LIAI = AFFE_CHAR_MECA(
    MODELE=model,
    DDL_IMPO=(
        _F(GROUP_NO="NO2", DX=0.0, DY=0.0, DZ=0.0),
        _F(GROUP_NO="NO6", DX=0.0, DY=0.0),
        _F(GROUP_NO="NO1", DX=0.0, DZ=0.0),
        _F(GROUP_NO=("NO9", "NO13", "NO14", "NO5", "NO17"), DX=0.0),
    ),
)

CHR_FORC = AFFE_CHAR_MECA(MODELE=model, FORCE_FACE=_F(GROUP_MA="MA2", FX=1298.0))

resultTrac = MECA_NON_LINE(
    MODELE=model,
    CHAM_MATER=materialFieldTrac,
    EXCIT=(_F(CHARGE=CHR_LIAI), _F(CHARGE=CHR_FORC, FONC_MULT=F_CHAR)),
    COMPORTEMENT=_F(RELATION="ELAS_VMIS_TRAC", DEFORMATION="GREEN_LAGRANGE"),
    INCREMENT=_F(LIST_INST=L_INST),
    CONVERGENCE=_F(RESI_GLOB_MAXI=1.0e-5, ITER_GLOB_MAXI=50),
)

resultTrac = CALC_CHAMP(reuse=resultTrac, INST=2.0, RESULTAT=resultTrac, FORCE="FORC_NODA")

# Value to test

timeToTest = 2.0
valeRefe = -1.08170000e8
valeCalc = -1.0816666666667e08

TEST_RESU(
    RESU=(
        _F(
            INST=timeToTest,
            RESULTAT=resultTrac,
            NOM_CHAM="FORC_NODA",
            GROUP_NO="NO8",
            NOM_CMP="DX",
            VALE_CALC=valeCalc,
            VALE_REFE=valeRefe,
            REFERENCE="ANALYTIQUE",
            PRECISION=1.0e-3,
        ),
    )
)

# ==================================================================================================
#
# Compute FORC_NODA
#
# ==================================================================================================


# Get fields from result
sief_elga = resultTrac.getField("SIEF_ELGA", para="INST", value=timeToTest)
disp = resultTrac.getField("DEPL", para="INST", value=timeToTest)


# Create the physical problem
phys_pb = CA.PhysicalProblem(model, materialFieldTrac)

# We have non-linear (GREEN-LAGRANGE) => we need behaviour Map !
phys_pb.computeBehaviourProperty(
    _F(RELATION="ELAS_VMIS_TRAC", DEFORMATION="GREEN_LAGRANGE"), "NON", 1
)

# Create object for discrete computation
disc_comp = CA.DiscreteComputation(phys_pb)

# Compute external state variables
varc_curr = phys_pb.getExternalStateVariables(timeToTest)

# Compute nodal forces
forc_noda = disc_comp.getMechanicalNodalForces(sief_elga, disp, varc_curr=varc_curr)


# Compute nodal forces with given map
behaviourMap = resultTrac.getField("COMPORTEMENT", para="INST", value=timeToTest)
forc_noda2 = disc_comp.getMechanicalNodalForces(
    sief_elga, disp, varc_curr=varc_curr, behaviourMap=behaviourMap
)


TEST_RESU(
    CHAM_NO=(
        _F(
            GROUP_NO="NO8",
            CHAM_GD=forc_noda,
            NOM_CMP="DX",
            VALE_CALC=valeCalc,
            VALE_REFE=valeRefe,
            REFERENCE="ANALYTIQUE",
            PRECISION=1.0e-3,
        ),
    )
)

TEST_RESU(
    CHAM_NO=(
        _F(
            GROUP_NO="NO8",
            CHAM_GD=forc_noda2,
            NOM_CMP="DX",
            VALE_CALC=valeCalc,
            VALE_REFE=valeRefe,
            REFERENCE="ANALYTIQUE",
            PRECISION=1.0e-3,
        ),
    )
)

# ==================================================================================================
#
# Compute FORC_EXTE
#
# ==================================================================================================

# Add loads
phys_pb.addLoad(CHR_LIAI)
phys_pb.addLoadFromDict({"CHARGE": CHR_FORC, "FONC_MULT": F_CHAR, "TYPE_CHARGE": "FIXE_CSTE"})

# Compute new numbering
phys_pb.computeDOFNumbering()

# Compute forces
forc_exte = disc_comp.getMechanicalForces(time_curr=timeToTest, varc_curr=varc_curr)

TEST_RESU(
    CHAM_NO=(
        _F(
            GROUP_NO="NO8",
            CHAM_GD=forc_exte,
            NOM_CMP="DX",
            VALE_CALC=valeCalc,
            VALE_REFE=valeRefe,
            REFERENCE="ANALYTIQUE",
            PRECISION=1.0e-3,
        ),
    )
)

# ==================================================================================================
#
# Compute REAC_NODA
#
# ==================================================================================================

# Compute forces
reac_noda = disc_comp.getMechanicalReactionForces(
    sief_elga, disp, time_curr=timeToTest, varc_curr=varc_curr
)

# # ORDRE_GRANDEUR n'existe pas pour TEST_RESU/CHAM_NO.
# # Bien que le résultat soit ZERO, je mets -1.1920928955078125e-06 :(
# TEST_RESU(
#     CHAM_NO=(
#         _F(
#             GROUP_NO="NO8",
#             CHAM_GD=reac_noda,
#             NOM_CMP="DX",
#             VALE_CALC=3e-7,
#             PRECISION=1.0e-3,
#             VALE_REFE=0.0,
#             CRITERE="ABSOLU",
#             REFERENCE="ANALYTIQUE",
#         ),
#     )
# )


########## tests interpolation
test = CA.TestCase()

t0 = 0.0
t1 = 1.0
v0 = 0.0
v1 = 2.0

depl_t0 = CA.FieldOnNodesReal(model)
depl_t0.setValues(v0)
depl_t1 = CA.FieldOnNodesReal(model)
depl_t1.setValues(v1)

sief_t0 = CA.FieldOnCellsReal(model, "ELGA", "SIEF_R", phys_pb.getBehaviourProperty())
sief_t0.setValues(v0)
sief_t1 = CA.FieldOnCellsReal(model, "ELGA", "SIEF_R", phys_pb.getBehaviourProperty())
sief_t1.setValues(v1)

resu = CREA_RESU(
    OPERATION="AFFE",
    TYPE_RESU="EVOL_ELAS",
    AFFE=(
        _F(NOM_CHAM="DEPL", MODELE=model, CHAM_GD=depl_t0, INST=t0),
        _F(NOM_CHAM="DEPL", MODELE=model, CHAM_GD=depl_t1, INST=t1),
        _F(NOM_CHAM="SIEF_ELGA", MODELE=model, CHAM_GD=sief_t0, INST=t0),
        _F(NOM_CHAM="SIEF_ELGA", MODELE=model, CHAM_GD=sief_t1, INST=t1),
    ),
)
t_interp = 0.5
v_interp_resu = 1.0
depl_interp = resu.interpolateField("DEPL", t_interp)
sief_interp = resu.interpolateField("SIEF_ELGA", t_interp)

test.assertTrue(all(abs(i - v_interp_resu) < 1.0e-12 for i in depl_interp.getValues()))
test.assertTrue(all(abs(i - v_interp_resu) < 1.0e-12 for i in sief_interp.getValues()))
with test.assertRaises(ValueError):
    err = resu.interpolateField("DEEEEEPL", t_interp)
with test.assertRaises(ValueError):
    err = resu.interpolateField("DEPL", 1.5)
with test.assertRaises(ValueError):
    err = resu.interpolateField("DEPL", -0.5)
with test.assertRaises(ValueError):
    err = resu.interpolateField("DEPL", t_interp, prec=-1.0)
with test.assertRaises(ValueError):
    err = resu.interpolateField("SIEF_ELGA", t_interp, crit="ABSOREL")

test.printSummary()
FIN()
#
