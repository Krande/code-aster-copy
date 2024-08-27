# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
from code_aster.CA import MPI
from code_aster.Utilities import petscInitialize, PETSc
import operator
from code_aster.Utilities.ExecutionParameter import ExecutionParameter, Options

test = CA.TestCase()

CA.init("--test")

petscInitialize("-on_error_abort ")

rank = MPI.ASTER_COMM_WORLD.Get_rank()

pMesh2 = CA.ParallelMesh()
pMesh2.readMedFile("zzzz504b/%d.med" % rank, partitioned=True)

model = AFFE_MODELE(
    MAILLAGE=pMesh2, AFFE=_F(MODELISATION="D_PLAN_INCO_UPG", PHENOMENE="MECANIQUE", TOUT="OUI")
)

char_cin = AFFE_CHAR_CINE(
    MODELE=model,
    MECA_IMPO=(
        _F(GROUP_NO="N2", DX=0.0, DY=0.0, DZ=0.0),
        _F(GROUP_NO="N4", DX=0.0, DY=0.0, DZ=0.0),
    ),
)

char_meca = AFFE_CHAR_MECA(
    MODELE=model,
    LIAISON_DDL=_F(GROUP_NO=("N1", "N3"), DDL=("PRES", "PRES"), COEF_MULT=(1.0, 1.0), COEF_IMPO=0),
    DDL_IMPO=(_F(GROUP_NO="N1", PRES=200000.0),),
)
char_meca.debugPrint(10 + rank)

MATER1 = DEFI_MATERIAU(ELAS=_F(E=200000.0, NU=0.499999, RHO=7.8))

AFFMAT = AFFE_MATERIAU(MAILLAGE=pMesh2, AFFE=_F(TOUT="OUI", MATER=MATER1))

LINSTC = DEFI_LIST_REEL(VALE=(1.0, 2.0))


resu1 = MECA_STATIQUE(
    CHAM_MATER=AFFMAT,
    MODELE=model,
    LIST_INST=LINSTC,
    INST_FIN=1.0,
    EXCIT=(_F(CHARGE=char_cin), _F(CHARGE=char_meca)),
)

ranks = resu1.getIndexes()
test.assertEqual(len(ranks), 1)
test.assertAlmostEqual(ranks[0], 1.0)
resu1.debugPrint(10 + rank)

test.assertEqual("resu1", resu1.userName)
test.assertFalse(resu1.hasElementaryCharacteristics())
test.assertFalse(resu1.hasElementaryCharacteristics(1))

resu1.printMedFile("test" + str(rank) + ".med")
# from shutil import copyfile
# copyfile("test"+str(rank)+".med", "/home/siavelis/test"+str(rank)+".med")

field1 = resu1.getField("DEPL", 1)
sfon = field1.toSimpleFieldOnNodes()
sfon.debugPrint(10 + rank)
sfon.build()


# DX displacement on nodes "N1" and "N3", comparison with sequential results
if rank == 0:
    test.assertAlmostEqual(sfon[1, 0], 1.14977255749554, 6)
elif rank == 1:
    test.assertAlmostEqual(sfon[1, 0], 1.14977255749554, 6)

myOptions = "-pc_type lu -pc_factor_mat_solver_type mumps -ksp_type fgmres -snes_linesearch_type basic  -snes_max_it 10"
resu2 = MECA_NON_LINE(
    MODELE=model,
    CHAM_MATER=AFFMAT,
    EXCIT=(_F(CHARGE=char_cin), _F(CHARGE=char_meca)),
    COMPORTEMENT=_F(RELATION="ELAS"),
    NEWTON=_F(REAC_INCR=1, PREDICTION="ELASTIQUE", MATRICE="TANGENTE", REAC_ITER=1),
    METHODE="SNES",
    CONVERGENCE=_F(RESI_GLOB_RELA=1e-8),
    INCREMENT=_F(LIST_INST=LINSTC),
    SOLVEUR=_F(METHODE="PETSC", OPTION_PETSC=myOptions),
    INFO=1,
)

field2 = resu2.getField("DEPL", 1)

test.assertAlmostEqual(field1.norm("NORM_2"), field2.norm("NORM_2"), 6)


# -------------------------------------------------------------------------

# Export of the local portions of a vector and a matrix

ExecutionParameter().disable(Options.UseLegacyMode)

model = AFFE_MODELE(
    MAILLAGE=pMesh2, AFFE=_F(MODELISATION="D_PLAN", PHENOMENE="MECANIQUE", TOUT="OUI")
)

char_cin = AFFE_CHAR_CINE(
    MODELE=model,
    MECA_IMPO=(
        _F(GROUP_NO="N2", DX=0.0, DY=0.0, DZ=0.0),
        _F(GROUP_NO="N4", DX=0.0, DY=0.0, DZ=0.0),
    ),
)

acier = DEFI_MATERIAU(ELAS=_F(E=200000.0, NU=0.3, RHO=38.0))

chmat = AFFE_MATERIAU(MAILLAGE=pMesh2, AFFE=_F(TOUT="OUI", MATER=acier))

pesa = AFFE_CHAR_MECA(MODELE=model, PESANTEUR=_F(GRAVITE=9.81, DIRECTION=(0.0, -1.0, 0.0)))

asse = ASSEMBLAGE(
    MODELE=model,
    CHAM_MATER=chmat,
    CHARGE=pesa,
    CHAR_CINE=char_cin,
    NUME_DDL=CO("numeddl"),
    MATR_ASSE=_F(MATRICE=CO("K"), OPTION="RIGI_MECA"),
    VECT_ASSE=(_F(VECTEUR=CO("F"), OPTION="CHAR_MECA")),
)

K, F = asse.K, asse.F


# print helper
def Print(*args):
    print(*args, flush=True)


# Export the parallel matrix
pK = K.toPetsc()
Print(f"pK.size={pK.size}")
test.assertEqual(pK.size, (36, 36))


# Export the local matrix
pK_loc = K.toPetsc(local=True)
Print(f"pK_loc.size={pK_loc.size}")
test.assertEqual(pK_loc.size, (26, 26))

# Change values of F and export to PETSc
F.setValues(1)
pF = F.toPetsc()
Print(f"pF.size={pF.size}")
test.assertEqual(pF.size, 36)

pF_loc = F.toPetsc(local=True)
Print(f"pF_loc.size={pF_loc.size}")
Print(f"pF_loc.comm={pF_loc.comm}")
Print(f"MPI.ASTER_COMM_SELF={MPI.ASTER_COMM_SELF}")
test.assertEqual(pF_loc.size, 26)
if rank == 0:
    pF *= 2
    pF_loc *= 2
pF_loc.view()
Print(f"F_init={F.getValues()}")
# Import fro PETSc - use local since pF_loc is a seq vec
F.fromPetsc(pF_loc, local=True)
Print(f"F_post={F.getValues()}")

# Assemble the local in a parallel one and check it is the same
pF2 = PETSc.Vec().create(comm=MPI.ASTER_COMM_WORLD)
pF2.setType("mpi")
globNume = F.getDescription()
ownedRows = globNume.getNoGhostDOFs()
neql = len(ownedRows)
neqg = globNume.getNumberOfDOFs(local=False)
pF2.setSizes((neql, neqg))
val = pF_loc.getArray()
l2g = globNume.getLocalToGlobalMapping()
extract = operator.itemgetter(*ownedRows)
ll2g = extract(l2g)
lval = extract(val)
pF2.setValues(ll2g, lval, PETSc.InsertMode.INSERT_VALUES)
pF2.assemble()

test.assertAlmostEqual((pF - pF2).norm(), 0.0)

# We only loop on the local rows (no ghost row!)
for i in ownedRows:
    for j in range(neql):
        test.assertAlmostEqual(pK_loc[i, j], pK[l2g[i], l2g[j]])


FIN()
