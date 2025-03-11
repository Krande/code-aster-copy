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
from code_aster.CA import MPI

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()


def printRank(mesh, filename):
    fon = CREA_CHAMP(
        OPERATION="AFFE",
        TYPE_CHAM="NOEU_TEMP_R",
        MAILLAGE=mesh,
        AFFE=_F(TOUT="OUI", NOM_CMP="TEMP", VALE=-1.0),
    )
    nodes_owner = mesh.getNodesOwner()
    assert len(nodes_owner) == fon.size()
    fon.updateValuePointers()
    for i in range(fon.size()):
        fon[i] = nodes_owner[i]
    fon.printMedFile(filename)


def printNumGlob(mesh, filename):
    fon = CREA_CHAMP(
        OPERATION="AFFE",
        TYPE_CHAM="NOEU_TEMP_R",
        MAILLAGE=mesh,
        AFFE=_F(TOUT="OUI", NOM_CMP="TEMP", VALE=-1.0),
    )
    global_num = mesh.getNodes(False)
    assert len(global_num) == fon.size()
    fon.updateValuePointers()
    for i in range(fon.size()):
        fon[i] = global_num[i]
    fon.printMedFile(filename)


def transfo(mesh):
    mesh_line = mesh.convertToLinear()

    mesh_raf = mesh_line.refine(2)

    mesh_quad = mesh_raf.convertToQuadratic()

    return mesh_quad


def computation(mesh, solv):
    mesh = MODI_MAILLAGE(reuse=mesh, MAILLAGE=mesh, ORIE_PEAU=_F(GROUP_MA_PEAU="L2"))

    MODE = AFFE_MODELE(
        MAILLAGE=mesh, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="C_PLAN")
    )

    ACIER = DEFI_MATERIAU(
        ELAS=_F(E=2.000000000e12, NU=0.3e00, RHO=0.000000000e03, ALPHA=0.000000000e00)
    )

    MATE = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=ACIER))

    R = FORMULE(VALE="(X-1.5)*(X-1.5) + (Y-0.5)*(Y-0.5)", NOM_PARA=["X", "Y"])

    DIRI = AFFE_CHAR_CINE_F(MODELE=MODE, MECA_IMPO=(_F(GROUP_MA="L1", DX=R, DY=R),))

    CHAR = AFFE_CHAR_MECA(MODELE=MODE, PRES_REP=_F(GROUP_MA="L2", PRES=-100.0))

    CHXN = CREA_CHAMP(
        OPERATION="EXTR", TYPE_CHAM="NOEU_GEOM_R", NOM_CHAM="GEOMETRIE", MAILLAGE=mesh, INFO=1
    )

    TEMP1 = CREA_CHAMP(
        OPERATION="AFFE",
        TYPE_CHAM="NOEU_NEUT_F",
        MAILLAGE=mesh,
        AFFE=(_F(TOUT="OUI", NOM_CMP="X1", VALE_F=R),),
    )
    TEMP2 = CREA_CHAMP(OPERATION="EVAL", TYPE_CHAM="NOEU_NEUT_R", CHAM_F=TEMP1, CHAM_PARA=CHXN)

    ASSEMBLAGE(
        MODELE=MODE,
        CHARGE=(CHAR,),
        CHAR_CINE=DIRI,
        CHAM_MATER=MATE,
        NUME_DDL=CO("NDDL"),
        MATR_ASSE=(_F(MATRICE=CO("MATK"), OPTION="RIGI_MECA"),),
    )

    Uana = CREA_CHAMP(
        OPERATION="ASSE",
        TYPE_CHAM="NOEU_DEPL_R",
        NUME_DDL=NDDL,
        MAILLAGE=mesh,
        PROL_ZERO="OUI",
        ASSE=(
            _F(TOUT="OUI", CHAM_GD=TEMP2, CUMUL="OUI", NOM_CMP="X1", NOM_CMP_RESU="DX"),
            _F(TOUT="OUI", CHAM_GD=TEMP2, CUMUL="OUI", NOM_CMP="X1", NOM_CMP_RESU="DY"),
        ),
    )

    MDEP = MATK * Uana

    CHAR = AFFE_CHAR_MECA(MODELE=MODE, VECT_ASSE=MDEP)

    RESU = MECA_STATIQUE(
        MODELE=MODE,
        CHAM_MATER=MATE,
        OPTION="SANS",
        EXCIT=(_F(CHARGE=CHAR), _F(CHARGE=DIRI)),
        **solv
    )

    return RESU


# test case


rank = MPI.ASTER_COMM_WORLD.Get_rank()
nbproc = MPI.ASTER_COMM_WORLD.Get_size()

mesh_p = LIRE_MAILLAGE(PARTITIONNEUR="PTSCOTCH", UNITE=20)
mesh = LIRE_MAILLAGE(UNITE=20)

mesh_p2 = transfo(mesh_p)
mesh2 = transfo(mesh)

solvers = [
    {"SOLVEUR": {"METHODE": "PETSC", "RESI_RELA": 1.0e-14, "PRE_COND": "ML"}},
    {"SOLVEUR": {"METHODE": "MUMPS", "TYPE_RESOL": "NONSYM"}},
    {"SOLVEUR": {"METHODE": "MUMPS"}},
]

for solver in solvers:
    resu_p = computation(mesh_p2, solver)
    resu = computation(mesh2, solver)

    depl = resu.getField("DEPL", 1)
    depl_p = resu_p.getField("DEPL", 1)

    test.assertAlmostEqual(depl.norm("NORM_2"), depl_p.norm("NORM_2"))


test.printSummary()

FIN()
