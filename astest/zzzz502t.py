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


from code_aster.CA import MPI
from code_aster.Commands import *
from code_aster import CA
from code_aster.Utilities import ExecutionParameter, Options
from code_aster.Utilities import petscInitialize

test = CA.TestCase()

CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

petscInitialize()
rank = MPI.ASTER_COMM_WORLD.Get_rank()

# modeling = "D_PLAN_INCO_UPG"
modeling = "D_PLAN_INCO_UP"
part = "PTSCOTCH"  #'SANS'

# -----------------------------------------------------------------------------
# ----------------------------------- Modèle  ---------------------------------
# -----------------------------------------------------------------------------

Mesh = LIRE_MAILLAGE(UNITE=20, FORMAT="MED", PARTITIONNEUR=part)

Model = AFFE_MODELE(
    MAILLAGE=Mesh,
    AFFE=_F(MODELISATION=modeling, PHENOMENE="MECANIQUE", TOUT="OUI"),
    DISTRIBUTION=_F(METHODE="CENTRALISE"),
)

CharCin = AFFE_CHAR_CINE(
    MODELE=Model,
    MECA_IMPO=(
        _F(GROUP_NO="N2", DX=0.0, DY=0.0, DZ=0.0),
        _F(GROUP_NO="N4", DX=0.0, DY=0.0, DZ=0.0),
    ),
)

CharMeca = AFFE_CHAR_MECA(
    MODELE=Model,
    LIAISON_DDL=_F(GROUP_NO=("N1", "N3"), DDL=("PRES", "PRES"), COEF_MULT=(1.0, 1.0), COEF_IMPO=0),
    DDL_IMPO=(_F(GROUP_NO="N1", PRES=200000.0),),
    INFO=2,
)

MATER1 = DEFI_MATERIAU(ELAS=_F(E=200000.0, NU=0.499999))

AffMat = AFFE_MATERIAU(MAILLAGE=Mesh, AFFE=_F(TOUT="OUI", MATER=MATER1))

# -----------------------------------------------------------------------------
# --------------------------------- Assemblage --------------------------------
# -----------------------------------------------------------------------------

ExecutionParameter().disable(Options.UseLegacyMode)

AssemblyObj = ASSEMBLAGE(
    MODELE=Model,
    CHAM_MATER=AffMat,
    CHARGE=CharMeca,
    CHAR_CINE=CharCin,
    NUME_DDL=CO("asterNume"),
    MATR_ASSE=(_F(MATRICE=CO("asterRigi"), OPTION="RIGI_MECA"),),
)

nume_ddl = AssemblyObj.asterNume
lagr_components = [comp for comp in nume_ddl.getComponents() if comp.startswith("LAGR")]

local_lagr_rows = set()
for comp in lagr_components:
    local_lagr_rows.update(nume_ddl.getDOFsAssociatedToComponent(comp))

global_lagr_rows = set()
for comp in lagr_components:
    global_lagr_rows.update(nume_ddl.getDOFsAssociatedToComponent(comp, local=False))

test.assertListEqual(sorted(local_lagr_rows), sorted(nume_ddl.getLagrangeDOFs(local=True)))
test.assertListEqual(sorted(global_lagr_rows), sorted(nume_ddl.getLagrangeDOFs(local=False)))
dicLag = nume_ddl.getDictOfLagrangeDOFs(local=True)
ref1 = {0: [3], 1: [5, 6]}
ref2 = {0: [7], 1: [10, 11]}
test.assertListEqual(dicLag[1], ref1[rank])
test.assertListEqual(dicLag[2], ref2[rank])
dicLag = nume_ddl.getDictOfLagrangeDOFs(local=False)
ref1 = {0: [22], 1: [46, 22]}
ref2 = {0: [23], 1: [47, 23]}
test.assertListEqual(dicLag[1], ref1[rank])
test.assertListEqual(dicLag[2], ref2[rank])

nume_eq = nume_ddl.getEquationNumbering()
localDX = nume_eq.getDOFsWithDescription(["DX"], local=True, same_rank=True)
ref1 = {0: [0, 1, 2, 3, 4, 5, 6, 7, 8], 1: [0, 1, 2, 3, 4, 5, 6, 7, 8]}
ref2 = {0: [0, 4, 8, 11, 14, 16, 18, 20, 22], 1: [0, 2, 7, 12, 15, 18, 20, 22, 24]}
test.assertListEqual(localDX[0][0], ref1[rank])
test.assertListEqual(localDX[-1], ref2[rank])
globalDX = nume_eq.getDOFsWithDescription(["DX"], local=False, same_rank=True)
ref1 = {0: [0, 1, 2, 3, 4, 5, 6, 7, 8], 1: [9, 10, 11, 12, 13, 14, 15, 16, 17]}
ref2 = {0: [0, 3, 6, 9, 12, 14, 16, 18, 20], 1: [24, 26, 29, 32, 35, 38, 40, 42, 44]}
test.assertListEqual(globalDX[0][0], ref1[rank])
test.assertListEqual(globalDX[-1], ref2[rank])


petscMat = AssemblyObj.asterRigi.toPetsc()
print("Norm: ", petscMat.getSizes())
ref = 1823496.3881588143
test.assertAlmostEqual(petscMat.norm(), ref, delta=ref * 1.0e-6)
test.assertSequenceEqual(petscMat.getSizes(), ((24, 48), (24, 48)))

# petscMat.view()

# -----------------------------------------------------------------------------
# --------------------------------- Fin Aster ---------------------------------
# -----------------------------------------------------------------------------

FIN()
