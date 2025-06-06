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

# Mesh and model
mesh = LIRE_MAILLAGE(FORMAT="MED", UNITE=20)
model = AFFE_MODELE(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))

# Material
acier = DEFI_MATERIAU(ELAS=_F(E=200000.0, NU=0.3), ECRO_LINE=_F(D_SIGM_EPSI=2000.0, SY=200000.0))
mater = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=acier))

# The physical problem
phys_pb = CA.PhysicalProblem(model, mater)

# Numbering of equations
phys_pb.computeDOFNumbering()

# Load
value = 0.0
zero = CA.FieldOnNodesReal(phys_pb.getDOFNumbering())
kinematic = AFFE_CHAR_CINE(
    MODELE=model, MECA_IMPO=(_F(GROUP_MA="BAS", DX=0, DY=0.0, DZ=0.0), _F(GROUP_MA="HAUT", DZ=1.0))
)
kinematicField = CALC_CHAR_CINE(NUME_DDL=phys_pb.getDOFNumbering(), CHAR_CINE=kinematic)

# Create behaviour
phys_pb.computeBehaviourProperty(COMPORTEMENT=(_F(RELATION="VMIS_ISOT_LINE", TOUT="OUI"),))

# Create displacement field
value = 0.0
disp = CA.FieldOnNodesReal(phys_pb.getDOFNumbering())
disp.setValues(value)
disp_incr = CA.FieldOnNodesReal(phys_pb.getDOFNumbering())
disp_incr.setValues(value)

# Create stress field
value = 0.0
stress = CA.FieldOnCellsReal(phys_pb.getModel(), "ELGA", "SIEF_R")
stress.setValues(value)

# Create internal state variable field
value = 0.0
internVar = CA.FieldOnCellsReal(
    phys_pb.getModel(), "ELGA", "VARI_R", phys_pb.getBehaviourProperty()
)
internVar.setValues(value)

# Object for discrete computation
disc_comp = CA.DiscreteComputation(phys_pb)

# Time
timeBeginStep = 0.0
timeEndStep = 1.0

# Compute tangent prediction matrix
res = disc_comp.getTangentStiffnessMatrix(
    disp, disp_incr, stress, internVar, internVar, timeBeginStep, timeEndStep
)
matrElem = res[2]

# Assemblying
matrAsseRef = ASSE_MATRICE(
    MATR_ELEM=matrElem, NUME_DDL=phys_pb.getDOFNumbering(), CHAR_CINE=kinematic
)

# Solve
matrAsseRef = FACTORISER(reuse=matrAsseRef, METHODE="MULT_FRONT", MATR_ASSE=matrAsseRef)
fieldResuRefe = RESOUDRE(MATR=matrAsseRef, CHAM_NO=zero, CHAM_CINE=kinematicField)

# Compute tangent prediction matrices
res1 = disc_comp.getTangentStiffnessMatrix(
    disp,
    disp_incr,
    stress,
    internVar,
    internVar,
    timeBeginStep,
    timeEndStep,
    groupOfCells=["CoucheHaut"],
)
matrElem1 = res1[2]
res2 = disc_comp.getTangentStiffnessMatrix(
    disp,
    disp_incr,
    stress,
    internVar,
    internVar,
    timeBeginStep,
    timeEndStep,
    groupOfCells=["CoucheBas"],
)
matrElem2 = res2[2]

# Assemblying
matrAsseComb = ASSE_MATRICE(
    MATR_ELEM=(matrElem1, matrElem2), NUME_DDL=phys_pb.getDOFNumbering(), CHAR_CINE=kinematic
)

# Solve
matrAsseComb = FACTORISER(reuse=matrAsseComb, METHODE="MULT_FRONT", MATR_ASSE=matrAsseComb)
fieldResuComb = RESOUDRE(MATR=matrAsseComb, CHAM_NO=zero, CHAM_CINE=kinematicField)

TEST_RESU(
    CHAM_NO=_F(
        REFERENCE="ANALYTIQUE",
        NOM_CMP="DZ",
        TYPE_TEST="MAX",
        CHAM_GD=fieldResuRefe,
        VALE_CALC=1.0,
        VALE_REFE=1.0,
    )
)

TEST_RESU(
    CHAM_NO=_F(
        REFERENCE="ANALYTIQUE",
        NOM_CMP="DZ",
        TYPE_TEST="MAX",
        CHAM_GD=fieldResuComb,
        VALE_CALC=1.0,
        VALE_REFE=1.0,
    )
)

FIN()
