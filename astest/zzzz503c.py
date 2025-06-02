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


from code_aster.CA import MPI
from code_aster.Commands import *
from code_aster import CA
from code_aster.Utilities import SharedTmpdir

from code_aster.MedUtils import readMedFileToResults

CA.init("--test")

test = CA.TestCase()

rank = MPI.ASTER_COMM_WORLD.Get_rank()
pMesh = CA.ParallelMesh()
pMesh.readMedFile("mesh004a/%d.med" % rank, partitioned=True)
DEFI_GROUP(
    reuse=pMesh,
    MAILLAGE=pMesh,
    CREA_GROUP_MA=(_F(NOM="BLABLA", OPTION="SPHERE", POINT=(0.2, 0.2, 0.2), RAYON=0.2),),
)

list_cells = pMesh.getCells("BLABLA")
nb_cells = [72, 4, 0, 4]
exi_grp = [True, True, False, True]
test.assertEqual(len(list_cells), nb_cells[rank])
test.assertEqual(pMesh.hasGroupOfCells("BLABLA", True), exi_grp[rank])
test.assertEqual(pMesh.hasGroupOfCells("BLABLA", False), True)

monModel = CA.Model(pMesh)
monModel.addModelingOnMesh(CA.Physics.Mechanics, CA.Modelings.Tridimensional)
monModel.build()

testMesh = monModel.getMesh()
test.assertEqual(testMesh.getType(), "MAILLAGE_P")

acier = DEFI_MATERIAU(ELAS=_F(E=2.0e11, NU=0.3))

affectMat = CA.MaterialField(pMesh)
affectMat.addMaterialOnMesh(acier)
affectMat.build()

testMesh2 = affectMat.getMesh()
test.assertEqual(testMesh2.getType(), "MAILLAGE_P")

charCine = CA.MechanicalDirichletBC(monModel)
charCine.addBCOnNodes(CA.PhysicalQuantityComponent.Dx, 0.0, "COTE_B")
charCine.addBCOnNodes(CA.PhysicalQuantityComponent.Dy, 0.0, "COTE_B")
charCine.addBCOnNodes(CA.PhysicalQuantityComponent.Dz, 0.0, "COTE_B")
charCine.build()

charCine2 = CA.MechanicalDirichletBC(monModel)
charCine2.addBCOnNodes(CA.PhysicalQuantityComponent.Dz, 1.0, "COTE_H")
charCine2.build()

resu = MECA_STATIQUE(
    MODELE=monModel,
    CHAM_MATER=affectMat,
    EXCIT=(_F(CHARGE=charCine), _F(CHARGE=charCine2)),
    SOLVEUR=_F(METHODE="PETSC", RENUM="SANS", PRE_COND="SOR"),
)


test.assertFalse(resu.hasElementaryCharacteristics())
test.assertFalse(resu.hasElementaryCharacteristics(1))

resu = CALC_CHAMP(RESULTAT=resu, reuse=resu, CONTRAINTE=("SIEF_ELGA"))

DEPL = resu.getField("DEPL", 1)
sfon = DEPL.toSimpleFieldOnNodes()
sfon.build()

SIEF = resu.getField("SIEF_ELGA", 1)


val = [0.134228076192, 0.134176297047, 0.154099687654, 0.154189676715]
test.assertAlmostEqual(sfon[4, 1], val[rank])


with SharedTmpdir("zzzz503c_") as tmpdir:
    medfile = osp.join(tmpdir.path, "resu.zzzz503c.no.med")
    DEFI_FICHIER(UNITE=81, FICHIER=medfile, TYPE="BINARY")

    IMPR_RESU(
        FICHIER_UNIQUE="OUI",
        FORMAT="MED",
        RESU=_F(RESULTAT=resu, GROUP_NO="COTE_H"),
        VERSION_MED="4.0.0",
        UNITE=81,
    )

    DEFI_FICHIER(ACTION="LIBERER", UNITE=81)

with SharedTmpdir("zzzz503c_") as tmpdir:
    medfile = osp.join(tmpdir.path, "resu" + str(rank) + ".zzzz503c.no.med")
    DEFI_FICHIER(UNITE=81, FICHIER=medfile, TYPE="BINARY")

    IMPR_RESU(FORMAT="MED", RESU=_F(RESULTAT=resu), VERSION_MED="4.0.0", UNITE=81)

    DEFI_FICHIER(ACTION="LIBERER", UNITE=81)

    myResu = readMedFileToResults(
        medfile,
        {"resu____DEPL": "DEPL", "resu____SIEF_ELGA": "SIEF_ELGA"},
        CA.ElasticResult,
        monModel,
    )
    depl1 = myResu.getField("DEPL", 1).getValues()
    depl2 = resu.getField("DEPL", 1).getValues()
    test.assertEqual(depl1, depl2)

with SharedTmpdir("zzzz503c_") as tmpdir:
    medfile = osp.join(tmpdir.path, "resu.zzzz503c.ma.med")
    DEFI_FICHIER(UNITE=82, FICHIER=medfile, TYPE="BINARY")

    IMPR_RESU(
        FICHIER_UNIQUE="OUI",
        FORMAT="MED",
        RESU=_F(RESULTAT=resu, GROUP_MA="BLABLA"),
        VERSION_MED="4.0.0",
        UNITE=82,
    )

    DEFI_FICHIER(ACTION="LIBERER", UNITE=82)

# load result in sequential
with SharedTmpdir("zzzz503c_") as tmpdir:
    medfile = osp.join(tmpdir.path, "resu.zzzz503c.med")
    DEFI_FICHIER(UNITE=80, FICHIER=medfile, TYPE="BINARY")

    IMPR_RESU(
        FICHIER_UNIQUE="OUI", FORMAT="MED", UNITE=80, RESU=_F(RESULTAT=resu), VERSION_MED="4.0.0"
    )

    mesh_std = LIRE_MAILLAGE(UNITE=80, FORMAT="MED")

    model_std = AFFE_MODELE(
        AFFE=_F(MODELISATION=("3D",), PHENOMENE="MECANIQUE", TOUT="OUI"), MAILLAGE=mesh_std
    )

    affectMat_std = CA.MaterialField(mesh_std)
    affectMat_std.addMaterialOnMesh(acier)
    affectMat_std.build()

    charCine_std = CA.MechanicalDirichletBC(model_std)
    charCine_std.addBCOnNodes(CA.PhysicalQuantityComponent.Dx, 0.0, "COTE_B")
    charCine_std.addBCOnNodes(CA.PhysicalQuantityComponent.Dy, 0.0, "COTE_B")
    charCine_std.addBCOnNodes(CA.PhysicalQuantityComponent.Dz, 0.0, "COTE_B")
    charCine_std.build()

    charCine2_std = CA.MechanicalDirichletBC(model_std)
    charCine2_std.addBCOnNodes(CA.PhysicalQuantityComponent.Dz, 1.0, "COTE_H")
    charCine2_std.build()

    resu_std = LIRE_RESU(
        MODELE=model_std,
        FORMAT="MED",
        UNITE=80,
        TYPE_RESU="EVOL_ELAS",
        CHAM_MATER=affectMat_std,
        EXCIT=(_F(CHARGE=charCine_std), _F(CHARGE=charCine2_std)),
        FORMAT_MED=(
            _F(NOM_RESU="resu", NOM_CHAM="DEPL"),
            _F(NOM_RESU="resu", NOM_CHAM="SIEF_ELGA"),
        ),
        TOUT_ORDRE="OUI",
    )

    DEFI_FICHIER(ACTION="LIBERER", UNITE=80)


SIEF_std = resu_std.getField("SIEF_ELGA", 1)
DEPL_std = resu_std.getField("DEPL", 1)

rela = abs(DEPL.norm("NORM_2") - DEPL_std.norm("NORM_2")) / DEPL_std.norm("NORM_2")
test.assertAlmostEqual(rela, 0.0, delta=1e-12)
rela = abs(SIEF.norm("NORM_2") - SIEF_std.norm("NORM_2")) / SIEF_std.norm("NORM_2")
test.assertAlmostEqual(rela, 0.0, delta=1e-12)


test.printSummary()

FIN()
