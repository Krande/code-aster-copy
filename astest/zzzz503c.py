# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

import code_aster
from code_aster import MPI
from code_aster.Commands import *
from code_aster.Utilities import shared_tmpdir

code_aster.init("--test")

test = code_aster.TestCase()

rank = MPI.COMM_WORLD.Get_rank()
pMesh = code_aster.ParallelMesh()
pMesh.readMedFile("mesh004a/%d.med"%rank, True)
DEFI_GROUP(reuse=pMesh,
           MAILLAGE=pMesh,
           CREA_GROUP_MA=(_F(NOM='BLABLA',
                             OPTION='SPHERE',
                             POINT=(0.2, 0.2, 0.2),
                             RAYON = 0.2),),)

pMesh.debugPrint(30+rank)

list_cells = pMesh.getCells( 'BLABLA' )
nb_cells = [ 72 , 4 , 0, 4]
exi_grp = [True, True, False, True]
test.assertEqual(len(list_cells), nb_cells[rank])
test.assertEqual(pMesh.hasGroupOfCells('BLABLA', True), exi_grp[rank])
test.assertEqual(pMesh.hasGroupOfCells('BLABLA', False), True)

monModel = code_aster.Model(pMesh)
monModel.addModelingOnMesh(code_aster.Physics.Mechanics,
                              code_aster.Modelings.Tridimensional)
monModel.build()

testMesh = monModel.getMesh()
test.assertEqual(testMesh.getType(), "MAILLAGE_P")

acier = DEFI_MATERIAU(ELAS = _F(E = 2.e11,
                                NU = 0.3,),)

affectMat = code_aster.MaterialField(pMesh)
affectMat.addMaterialsOnMesh( acier )
affectMat.buildWithoutExternalStateVariables()

testMesh2 = affectMat.getMesh()
test.assertEqual(testMesh2.getType(), "MAILLAGE_P")

charCine = code_aster.MechanicalDirichletBC(monModel)
charCine.addBCOnNodes(code_aster.PhysicalQuantityComponent.Dx, 0., "COTE_B")
charCine.addBCOnNodes(code_aster.PhysicalQuantityComponent.Dy, 0., "COTE_B")
charCine.addBCOnNodes(code_aster.PhysicalQuantityComponent.Dz, 0., "COTE_B")
charCine.build()

charCine2 = code_aster.MechanicalDirichletBC(monModel)
charCine2.addBCOnNodes(code_aster.PhysicalQuantityComponent.Dz, 1., "COTE_H")
charCine2.build()

monSolver = code_aster.PetscSolver(code_aster.Renumbering.Sans)
monSolver.setPreconditioning(code_aster.Preconditioning.Sor)

mecaStatique = code_aster.LinearStaticAnalysis(monModel, affectMat)
mecaStatique.addDirichletBC(charCine)
mecaStatique.addDirichletBC(charCine2)
mecaStatique.setLinearSolver(monSolver)

resu = mecaStatique.execute()

test.assertFalse(resu.hasElementaryCharacteristics())
test.assertFalse(resu.hasElementaryCharacteristics(1))

resu=CALC_CHAMP(RESULTAT=resu, reuse=resu, CONTRAINTE=('SIEF_ELGA'))

DEPL = resu.getFieldOnNodesReal("DEPL", 0)
sfon = DEPL.exportToSimpleFieldOnNodes()
sfon.updateValuePointers()

SIEF = resu.getFieldOnCellsReal("SIEF_ELGA", 0)


val = [0.134228076192 , 0.134176297047, 0.154099687654, 0.154189676715]
test.assertAlmostEqual(sfon.getValue(4, 1), val[rank])


with shared_tmpdir("zzzz503c_") as tmpdir:
    medfile = osp.join(tmpdir, "resu.zzzz503c.no.med")
    DEFI_FICHIER( UNITE=81, FICHIER=medfile, TYPE='BINARY')

    IMPR_RESU(FICHIER_UNIQUE='OUI',
            FORMAT='MED',
            RESU=_F(RESULTAT=resu, GROUP_NO="COTE_H"),
            VERSION_MED='4.0.0',
            UNITE=81)

    DEFI_FICHIER(ACTION='LIBERER',UNITE=81)

with shared_tmpdir("zzzz503c_") as tmpdir:
    medfile = osp.join(tmpdir, "resu.zzzz503c.ma.med")
    DEFI_FICHIER( UNITE=82, FICHIER=medfile, TYPE='BINARY')

    IMPR_RESU(FICHIER_UNIQUE='OUI',
            FORMAT='MED',
            RESU=_F(RESULTAT=resu, GROUP_MA="BLABLA"),
            VERSION_MED='4.0.0',
            UNITE=82)

    DEFI_FICHIER(ACTION='LIBERER',UNITE=82)

# load result in sequential
with shared_tmpdir("zzzz503c_") as tmpdir:
    medfile = osp.join(tmpdir, "resu.zzzz503c.med")
    DEFI_FICHIER( UNITE=80, FICHIER=medfile, TYPE='BINARY')

    IMPR_RESU(FICHIER_UNIQUE='OUI',
            FORMAT='MED',
            UNITE=80,
            RESU=_F(RESULTAT=resu,),
            VERSION_MED='4.0.0')

    mesh_std = LIRE_MAILLAGE(UNITE=80, FORMAT='MED',)

    model_std = AFFE_MODELE(
        AFFE=_F(MODELISATION=('3D', ),
                PHENOMENE='MECANIQUE', TOUT='OUI'),
        MAILLAGE=mesh_std
    )

    affectMat_std = code_aster.MaterialField(mesh_std)
    affectMat_std.addMaterialsOnMesh( acier )
    affectMat_std.buildWithoutExternalStateVariables()

    charCine_std = code_aster.MechanicalDirichletBC(model_std)
    charCine_std.addBCOnNodes(code_aster.PhysicalQuantityComponent.Dx, 0., "COTE_B")
    charCine_std.addBCOnNodes(code_aster.PhysicalQuantityComponent.Dy, 0., "COTE_B")
    charCine_std.addBCOnNodes(code_aster.PhysicalQuantityComponent.Dz, 0., "COTE_B")
    charCine_std.build()

    charCine2_std = code_aster.MechanicalDirichletBC(model_std)
    charCine2_std.addBCOnNodes(code_aster.PhysicalQuantityComponent.Dz, 1., "COTE_H")
    charCine2_std.build()

    resu_std = LIRE_RESU(MODELE=model_std,
                    FORMAT='MED',
                    UNITE=80,
                    TYPE_RESU='EVOL_ELAS',
                    CHAM_MATER=affectMat_std,
                    EXCIT=(
                            _F(CHARGE=charCine_std, ),
                            _F(CHARGE=charCine2_std, ),
                        ),
                    FORMAT_MED=(_F(NOM_RESU="resu", NOM_CHAM='DEPL'),
                                _F(NOM_RESU="resu", NOM_CHAM='SIEF_ELGA'),),
                    TOUT_ORDRE="OUI")

    DEFI_FICHIER(ACTION='LIBERER',UNITE=80)


SIEF_std = resu_std.getFieldOnCellsReal("SIEF_ELGA", 0)
DEPL_std = resu_std.getFieldOnNodesReal("DEPL", 0)

rela = abs(DEPL.norm("NORM_2")- DEPL_std.norm("NORM_2"))/DEPL_std.norm("NORM_2")
test.assertAlmostEqual(rela, 0.0, delta=1e-12)
rela = abs(SIEF.norm("NORM_2")- SIEF_std.norm("NORM_2"))/SIEF_std.norm("NORM_2")
test.assertAlmostEqual(rela, 0.0, delta=1e-12)


test.printSummary()

FIN()
