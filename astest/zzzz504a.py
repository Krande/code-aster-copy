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

test = CA.TestCase()

CA.init("--test")

nProc = MPI.ASTER_COMM_WORLD.Get_size()
parallel = nProc > 1


rank = MPI.ASTER_COMM_WORLD.Get_rank()
comm = MPI.ASTER_COMM_WORLD

if parallel:
    rank = MPI.ASTER_COMM_WORLD.Get_rank()
    pMesh2 = CA.ParallelMesh()
    pMesh2.readMedFile(f"mesh004c/{rank}.med", partitioned=True)
    # os.system('echo "-mat_view :/tmp/par.txt:ascii_matlab " > ~/.petscrc')
    # os.system('echo "-ksp_view_rhs ascii:/tmp/rhs_par.txt " >> ~/.petscrc')
    # os.system('echo "-ksp_view_solution ascii:/tmp/sol_par.txt  " >> ~/.petscrc')
else:
    pMesh2 = CA.Mesh()
    pMesh2.readMedFile("zzzz504a.med")
    # os.system('echo "-mat_view :/tmp/seq.txt:ascii_matlab " > ~/.petscrc')
    # os.system('echo "-ksp_view_rhs ascii:/tmp/rhs_seq.txt  " >> ~/.petscrc')
    # os.system('echo "-ksp_view_solution ascii:/tmp/sol_seq.txt  " >> ~/.petscrc')

# print a unique file
with SharedTmpdir("zzzz504a_") as tmpdir:
    medfile = osp.join(tmpdir.path, "zzzz504a.med")
    DEFI_FICHIER(UNITE=80, FICHIER=medfile, TYPE="BINARY")

    IMPR_RESU(
        FICHIER_UNIQUE="OUI", FORMAT="MED", UNITE=80, RESU=_F(MAILLAGE=pMesh2), VERSION_MED="4.0.0"
    )

    DEFI_FICHIER(ACTION="LIBERER", UNITE=80)

with SharedTmpdir("zzzz504a_") as tmpdir:
    medfile = osp.join(tmpdir.path, "zzzz504a_new_1.med")
    pMesh2.printMedFile(medfile, False)

# print separated file
with SharedTmpdir("zzzz504a_") as tmpdir:
    medfile = osp.join(tmpdir.path, f"zzzz504a_new_{rank}.med")
    pMesh2.printMedFile(medfile, True)

with SharedTmpdir("zzzz504a_") as tmpdir:
    medfile = osp.join(tmpdir.path, f"zzzz504a_{rank}_0.med")
    pMesh2.printMedFile(medfile)

model = AFFE_MODELE(
    MAILLAGE=pMesh2,
    AFFE=_F(MODELISATION="D_PLAN", PHENOMENE="MECANIQUE", TOUT="OUI"),
    DISTRIBUTION=_F(METHODE="CENTRALISE"),
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
    LIAISON_DDL=_F(GROUP_NO=("N1", "N3"), DDL=("DX", "DX"), COEF_MULT=(1.0, -1.0), COEF_IMPO=0),
    DDL_IMPO=_F(GROUP_NO="N1", DX=1.0),
)

MATER1 = DEFI_MATERIAU(ELAS=_F(E=200000.0, NU=0.3))

AFFMAT = AFFE_MATERIAU(MAILLAGE=pMesh2, AFFE=_F(TOUT="OUI", MATER=MATER1))

LI = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=1.0, NOMBRE=1))

resu = STAT_NON_LINE(
    CHAM_MATER=AFFMAT,
    METHODE="NEWTON",
    COMPORTEMENT=_F(RELATION="ELAS"),
    CONVERGENCE=_F(RESI_GLOB_RELA=1.0e-1),
    EXCIT=(_F(CHARGE=char_cin), _F(CHARGE=char_meca)),
    INCREMENT=_F(LIST_INST=LI),
    MODELE=model,
    NEWTON=_F(MATRICE="TANGENTE", REAC_ITER=1),
    SOLVEUR=_F(METHODE="PETSC", RESI_RELA=1.0e-10, PRE_COND="LDLT_SP"),
)

# print single file
with SharedTmpdir("zzzz504a_") as tmpdir:
    medfile = osp.join(tmpdir.path, "resu_new.med")
    resu.printMedFile(medfile, local=False)

# print multiple files
with SharedTmpdir("zzzz504a_") as tmpdir:
    medfile = osp.join(tmpdir.path, f"resu_new_{rank}.med")
    resu.printMedFile(medfile, local=True)

# print single files
depl = resu.getField("DEPL", 1.0, "INST")
with SharedTmpdir("zzzz504a_") as tmpdir:
    medfile = osp.join(tmpdir.path, "depl.med")
    depl.printMedFile(medfile, local=False)

# if (parallel):
# rank = MPI.ASTER_COMM_WORLD.Get_rank()
# myFile='par.txt'
# os.system("sed 's/Mat_.*\=/par\ \=/g' /tmp/par.txt > /tmp/par_clean.txt && mv /tmp/par_clean.txt /tmp/par.txt")
# if (rank==0): os.system( """grep -v %% /tmp/%s | grep -v zzz | grep -v \] | grep -v Mat | awk '{print $3}' | LANG=en_US.UTF-8  sort -g > /tmp/%s_sorted"""%(myFile,myFile) )
# if (rank==0): os.system( """grep -v Object /tmp/rhs_par.txt | grep -v type | grep -v Process | LANG=en_US.UTF-8  sort -g > /tmp/rhs_par_sorted.txt  """)
# if (rank==0): os.system( """grep -v Object /tmp/sol_par.txt | grep -v type | grep -v Process | LANG=en_US.UTF-8  sort -g > /tmp/sol_par_sorted.txt  """)
# else:
# myFile='seq.txt'
# os.system("sed 's/Mat_.*\=/seq\ \=/g' /tmp/seq.txt > /tmp/seq_clean.txt && mv /tmp/seq_clean.txt /tmp/seq.txt")
# os.system("grep -v %% /tmp/%s | grep -v zzz | grep -v \] | grep -v Mat | awk '{print $3}' | LANG=en_US.UTF-8  sort -g > /tmp/%s_sorted"%(myFile,myFile))
# os.system( """grep -v Object /tmp/rhs_seq.txt | grep -v type | grep -v Process | LANG=en_US.UTF-8  sort -g > /tmp/rhs_seq_sorted.txt  """)
# os.system( """grep -v Object /tmp/sol_seq.txt | grep -v type | grep -v Process | LANG=en_US.UTF-8  sort -g > /tmp/sol_seq_sorted.txt  """)

# if (parallel):
# rank = MPI.ASTER_COMM_WORLD.Get_rank()
# if (rank==0):
# os.system( """cp fort.11 /tmp/ddl0.txt """ )
# os.system( """cp fort.30 /tmp/sol_petsc_par_0.txt """ )
# if (rank==1):
# os.system( """cp fort.12 /tmp/ddl1.txt """ )
# os.system( """cp fort.31 /tmp/sol_petsc_par_1.txt """ )
# else:
# os.system( """cp fort.11 /tmp/ddl_seq.txt """ )
# os.system( """cp fort.19 /tmp/sol_petsc_seq.txt """ )

# if parallel:
# rank = MPI.ASTER_COMM_WORLD.Get_rank()
# resu.printMedFile('/tmp/par_%d.resu.med'%rank)
# else:
# resu.printMedFile('/tmp/seq.resu.med')

MyFieldOnNodes = resu.getField("DEPL", 1)
sfon = MyFieldOnNodes.toSimpleFieldOnNodes()
if parallel:
    test.assertAlmostEqual(sfon[2, 0], 0.5175556367605225)
else:
    test.assertAlmostEqual(sfon[6, 0], 0.0)

pMesh2 = CA.Mesh()
pMesh2.readMedFile("zzzz504a.med")

# Test if med file print by #0 is not overwrite by others procs
with SharedTmpdir("zzzz504a_") as tmpdir:
    medfile = osp.join(tmpdir.path, f"resu_new{rank}.med")
    resu.printMedFile(medfile, local=True)

model = AFFE_MODELE(
    MAILLAGE=pMesh2,
    AFFE=_F(MODELISATION="D_PLAN", PHENOMENE="MECANIQUE", TOUT="OUI"),
    DISTRIBUTION=_F(METHODE="CENTRALISE"),
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
    LIAISON_DDL=_F(GROUP_NO=("N1", "N3"), DDL=("DX", "DX"), COEF_MULT=(1.0, -1.0), COEF_IMPO=0),
    DDL_IMPO=_F(GROUP_NO="N1", DX=1.0),
)

MATER1 = DEFI_MATERIAU(ELAS=_F(E=200000.0, NU=0.3))

AFFMAT = AFFE_MATERIAU(MAILLAGE=pMesh2, AFFE=_F(TOUT="OUI", MATER=MATER1))

LI = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=_F(JUSQU_A=1.0, NOMBRE=1))

resu = STAT_NON_LINE(
    CHAM_MATER=AFFMAT,
    METHODE="NEWTON",
    COMPORTEMENT=_F(RELATION="ELAS"),
    CONVERGENCE=_F(RESI_GLOB_RELA=1.0e-1),
    EXCIT=(_F(CHARGE=char_cin), _F(CHARGE=char_meca)),
    INCREMENT=_F(LIST_INST=LI),
    MODELE=model,
    NEWTON=_F(MATRICE="TANGENTE", REAC_ITER=1),
    SOLVEUR=_F(METHODE="PETSC", RESI_RELA=1.0e-10, PRE_COND="LDLT_SP"),
)

with SharedTmpdir("zzzz504a_") as tmpdir:
    medfile = osp.join(tmpdir.path, f"resu_new{rank}.med")
    resu.printMedFile(medfile, local=True)

FIN()
