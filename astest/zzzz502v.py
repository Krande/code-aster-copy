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

import numpy as np
from code_aster.CA import MPI
from code_aster.Utilities import PETSc


CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))

test = CA.TestCase()

if MPI.ASTER_COMM_WORLD.Get_size() > 1:
    is_parallel = True
else:
    is_parallel = False

if is_parallel:
    MA = CA.ParallelMesh()
else:
    MA = CA.Mesh()

MA.readMedFile("zzzz502s.mmed")

# MA=MA.refine(1)

MO = AFFE_MODELE(
    MAILLAGE=MA,
    AFFE=(
        _F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"),
        _F(GROUP_MA=("COTE_H"), PHENOMENE="MECANIQUE", MODELISATION="3D_ABSO"),
    ),
)

MAT = DEFI_MATERIAU(ELAS=_F(E=1000.0, NU=0.3, RHO=1000, AMOR_ALPHA=0.1, AMOR_BETA=0.1))

CHMAT = AFFE_MATERIAU(MAILLAGE=MA, AFFE=_F(TOUT="OUI", MATER=MAT))

sinus = FORMULE(VALE="sin(4*INST*pi*2.)", NOM_PARA="INST")

x0 = 0.0
y0 = 0.0
z0 = 0.0

x1 = 0.0
y1 = 0.0
z1 = 0.5

# TODO When issue33111 is fixed, set the plane wave loading
ONDE = AFFE_CHAR_MECA_F(
    MODELE=MO,
    ONDE_PLANE=_F(
        DIRECTION=(0.0, 0.0, 1.0),
        TYPE_ONDE="P",
        FONC_SIGNAL=sinus,
        COOR_SOURCE=(x0, y0, z0),
        COOR_REFLECHI=(x1, y1, z1),
        GROUP_MA=("COTE_H",),
    ),
)
BLOQ = AFFE_CHAR_CINE(MODELE=MO, MECA_IMPO=_F(GROUP_MA=("COTE_B",), DX=0, DY=0, DZ=0))
# ONDE = AFFE_CHAR_MECA_F(MODELE=MO, PRES_REP=_F(GROUP_MA=("COTE_H",), PRES=sinus), VERI_NORM="NON")

KEL = CALC_MATR_ELEM(OPTION="RIGI_MECA", MODELE=MO, CHAM_MATER=CHMAT)
MEL = CALC_MATR_ELEM(OPTION="MASS_MECA", MODELE=MO, CHAM_MATER=CHMAT)
CEL = CALC_MATR_ELEM(OPTION="AMOR_MECA", MODELE=MO, CHAM_MATER=CHMAT, RIGI_MECA=KEL, MASS_MECA=MEL)
FELEM = CALC_VECT_ELEM(OPTION="CHAR_MECA", CHARGE=ONDE, CHAM_MATER=CHMAT, INST=0.12)

NUMEDDL = NUME_DDL(MATR_RIGI=KEL)

STIFFNESS = ASSE_MATRICE(MATR_ELEM=KEL, NUME_DDL=NUMEDDL, CHAR_CINE=BLOQ)

DAMPING = ASSE_MATRICE(MATR_ELEM=CEL, NUME_DDL=NUMEDDL, CHAR_CINE=BLOQ)

MASS = ASSE_MATRICE(MATR_ELEM=MEL, NUME_DDL=NUMEDDL, CHAR_CINE=BLOQ)
study = CA.PhysicalProblem(MO, CHMAT)
study.addDirichletBC(BLOQ)
study.setDOFNumbering(NUMEDDL)

FASS = ASSE_VECTEUR(VECT_ELEM=FELEM, NUME_DDL=NUMEDDL)

study.zeroDirichletBCDOFs(FASS)

freq = 1.0
REFE = DYNA_VIBRA(
    TYPE_CALCUL="HARM",
    BASE_CALCUL="PHYS",
    MODELE=MO,
    MATR_MASS=MASS,
    MATR_RIGI=STIFFNESS,
    MATR_AMOR=DAMPING,
    FREQ=freq,
    EXCIT=(_F(VECT_ASSE=FASS, COEF_MULT=1.0), _F(CHARGE=BLOQ, COEF_MULT=1.0)),
    INFO=1,
)
ref = REFE.getField("DEPL", 1)
# sequential comparaison
seq_value = 27.702679151680293
mpi_value = ref.norm("NORM_1")
test.assertAlmostEqual(mpi_value, seq_value)


# ==========================================================================
# here an equivalent real system is set-up and solved

nref = ref.getValues()
ref_r = CREA_CHAMP(OPERATION="C2R", TYPE_CHAM="NOEU_DEPL_R", CHAM_GD=ref, PARTIE="REEL", INFO=1)
ref_i = CREA_CHAMP(OPERATION="C2R", TYPE_CHAM="NOEU_DEPL_R", CHAM_GD=ref, PARTIE="IMAG", INFO=1)
pRef_r = ref_r.toPetsc()
pRef_i = ref_i.toPetsc()

ref_rr = ref.getRealPart()
norm_r = (ref_rr - ref_r).norm()
test.assertAlmostEqual(norm_r, 0)

ref_ii = ref.getImaginaryPart()
norm_i = (ref_ii - ref_i).norm()
test.assertAlmostEqual(norm_i, 0)

ref_rr = ref_r.getRealPart()
norm_r = (ref_rr - ref_r).norm()
test.assertAlmostEqual(norm_r, 0)

ref_ii = ref_r.getImaginaryPart()
norm_i = ref_ii.norm()
test.assertAlmostEqual(norm_i, 0)
# ----------------------------------------

om = 2 * np.pi * freq

mat_r = (STIFFNESS - om * om * MASS).toPetsc()
mat_r.setBlockSize(3)
mat_i = (om * DAMPING).toPetsc()
mat_i.setBlockSize(3)

rhs = FASS.toPetsc()

comm = mat_r.getComm()

# Decomposition in 2x2 blocks
A = PETSc.Mat().create(comm=comm)
A.createNest([[mat_i, -mat_r], [mat_r, mat_i]])
A.assemble()

# Decomposition in 2x2 blocks
y = PETSc.Vec().create(comm=comm)
y.createNest([0 * rhs, rhs])

yy = A.getVecRight()
yy[...] = y.getArray()[...]

optDB = PETSc.Options()
optDB.setValue("-ksp_monitor_true_residual", "")
# optDB.setValue("-ksp_view", "")

nested_IS = A.getNestISs()
rIS = nested_IS[0][0]
iIS = nested_IS[0][1]


# ------------------------------------------------------------------------------------
# Using MatNest
# "Solving Complex-Valued Linear Systems via Equivalent Real Formulations"
# David Day and Michael A. Heroux, SIAM Journal on Scientific ComputingVol. 23 (2) (2001)
# -----
ksp = PETSc.KSP().create(comm=comm)

ksp.setOperators(A, A)
ksp.setType("fgmres")

ksp.getPC().setType("fieldsplit")
ksp.getPC().setFieldSplitType(PETSc.PC.CompositeType.MULTIPLICATIVE)

ksp.getPC().setFieldSplitIS(("r", rIS), ("i", iIS))

ksp_r, ksp_i = ksp.getPC().getFieldSplitSubKSP()
ksp_r.setType("preonly")
ksp_r.getPC().setType("hypre")
ksp_i.setType("preonly")
ksp_i.getPC().setType("hypre")
ksp.setTolerances(rtol=1.0e-6)
ksp.setFromOptions()

x = y.duplicate()
ksp.solve(y, x)

xr = x.getSubVector(rIS)
xi = -x.getSubVector(iIS)
print(f"Error = {(xr-pRef_r).norm() / pRef_r.norm()*100:.2E} %")
print(f"Error = {(xi-pRef_i).norm() / pRef_i.norm()*100:.2E} %")
test.assertAlmostEqual(xr.norm(), pRef_r.norm(), places=6)
test.assertAlmostEqual(xi.norm(), pRef_i.norm(), places=6)


# ------------------------------------------------------------------------------------
# Using renumbering 2x2 blocks
# "Solving Complex-Valued Linear Systems via Equivalent Real Formulations"
# David Day and Michael A. Heroux, SIAM Journal on Scientific ComputingVol. 23 (2) (2001)
# -----
r = rIS.getIndices()
i = iIS.getIndices()

z = np.ravel(np.column_stack((r, i)))  # interlace r and i indices

z = PETSc.IS().createGeneral(z.tolist(), comm=comm)

Aij = A.copy()
Aij.convert("aij")
Aij = Aij.permute(z, z)
Aij.setBlockSize(2 * 3)

ksp = PETSc.KSP().create(comm=comm)

ksp.setOperators(Aij, Aij)
ksp.setType("fgmres")
ksp.getPC().setType("bjacobi")
optDB.setValue("-sub_pc_type", "lu")
ksp.setFromOptions()

yy.permute(z)
x = yy.duplicate()
ksp.solve(yy, x)

zz = z.invertPermutation()
x.permute(zz)

xr = x.getSubVector(rIS)
xi = -x.getSubVector(iIS)
print(f"Error = {(xr-pRef_r).norm() / pRef_r.norm()*100:.2E} %")
print(f"Error = {(xi-pRef_i).norm() / pRef_i.norm()*100:.2E} %")
test.assertAlmostEqual(xr.norm(), pRef_r.norm(), places=5)
test.assertAlmostEqual(xi.norm(), pRef_i.norm(), places=5)


# ------------------------------------------------------------------------------------
# Using C-to-R
# "A comparison of iterative methods to solve complex valued linear algebraic systems"
# Owe Axelsson · Maya Neytcheva · Bashir Ahmad, Numer Algor (2014)
# -----
class C2R(object):
    def create(self, pc, mat_r=mat_r, K=None):
        self.mat_r = mat_r
        self.ksp = PETSc.KSP().create(comm=mat_r.getComm())
        K_, _ = pc.getOperators()
        A = K_.getNestSubMatrix(0, 0)
        B = K_.getNestSubMatrix(1, 0)
        H = A + B
        self.ksp.setOperators(H, H)
        self.ksp.setType("fgmres")
        self.ksp.getPC().setType("jacobi")
        self.ksp.setTolerances(max_it=5, rtol=1.0e-3)
        self.ksp.setTabLevel(1)
        self.ksp.setFromOptions()
        l = K_.getNestISs()
        self.i1, self.i2 = l[0][0], l[0][1]
        self.A = A

    def apply(self, pc, x, y):
        f1 = x.getSubVector(self.i1)
        f2 = x.getSubVector(self.i2)
        r1 = f1 + f2
        g = r1.duplicate()
        self.ksp.solve(r1, g)
        r2 = f1 - self.A(g)
        h = r2.duplicate()
        self.ksp.solve(r2, h)
        w = PETSc.Vec().create(comm=comm)
        w.createNest([g + h, -h])
        y[...] = w


ksp = PETSc.KSP().create(comm=comm)

ksp.setOperators(A, A)
ksp.setType("fgmres")
ksp.setFromOptions()

ksp.getPC().setType(PETSc.PC.Type.PYTHON)
ksp.getPC().setPythonContext(C2R())

x = y.duplicate()
ksp.solve(y, x)

xr = x.getSubVector(rIS)
xi = -x.getSubVector(iIS)
print(f"Error = {(xr-pRef_r).norm() / pRef_r.norm()*100:.2E} %")
print(f"Error = {(xi-pRef_i).norm() / pRef_i.norm()*100:.2E} %")
test.assertAlmostEqual(xr.norm(), pRef_r.norm(), places=6)
test.assertAlmostEqual(xi.norm(), pRef_i.norm(), places=6)


FIN()
