# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

import code_aster

# Pour des raisons de performances (x5), on force SDVERI='NON'.
DEBUT(CODE=_F(NIV_PUB_WEB="INTERNET"), ERREUR=_F(ALARME="EXCEPTION"), DEBUG=_F(SDVERI="NON"))

import numpy as NP

# Read mesh
MAIL = LIRE_MAILLAGE(FORMAT="MED")

# Functions


CONDUC = DEFI_FONCTION(
    NOM_PARA="TEMP",
    NOM_RESU="LAMBDA",
    VALE=(0.0, 2.0e2, 500.0, 7.0e2),
    PROL_DROITE="LINEAIRE",
    PROL_GAUCHE="LINEAIRE",
)

ENTHAL = DEFI_FONCTION(
    NOM_PARA="TEMP",
    NOM_RESU="CP",
    VALE=(0.0, 0.0, 500.0, 4000.0e6),
    PROL_DROITE="LINEAIRE",
    PROL_GAUCHE="LINEAIRE",
)


# Material properties
MATER = DEFI_MATERIAU(THER_NL=_F(LAMBDA=CONDUC, BETA=ENTHAL))

CHMAT = AFFE_MATERIAU(MAILLAGE=MAIL, AFFE=_F(TOUT="OUI", MATER=MATER))

# Finite elements
MOTH = AFFE_MODELE(MAILLAGE=MAIL, AFFE=_F(TOUT="OUI", MODELISATION="PLAN", PHENOMENE="THERMIQUE"))

# Loads
TGAUCHE = DEFI_FONCTION(
    NOM_RESU="TEMP",
    NOM_PARA="INST",
    VALE=(0.0e0, 200.0e0, 10.0e0, 200.0e0, 10.001e0, 100.0e0, 100.0e0, 100.0e0),
)

TDROITE = DEFI_FONCTION(NOM_RESU="TEMP", NOM_PARA="INST", VALE=(0.0e0, 100.0e0, 100.0e0, 100.0e0))
CHTH = AFFE_CHAR_THER_F(
    MODELE=MOTH,
    TEMP_IMPO=(_F(GROUP_NO="NOE_GAU", TEMP=TGAUCHE), _F(GROUP_NO="NOE_DRO", TEMP=TDROITE)),
)

#
# DEFINITION DE LA STATEGIE DE CALCUL -----------------------------
#

LIST = DEFI_LIST_REEL(
    DEBUT=0.0,
    INTERVALLE=(
        _F(JUSQU_A=1.0e-3, NOMBRE=10),
        _F(JUSQU_A=1.0e-2, NOMBRE=9),
        _F(JUSQU_A=1.0e-1, NOMBRE=9),
        _F(JUSQU_A=1.0e0, NOMBRE=9),
        _F(JUSQU_A=10.0e0, NOMBRE=9),
        _F(JUSQU_A=13.0e0, NOMBRE=3),
    ),
)

temp_init = 100.0
theta = 0.57
NbIterNewtonMax = 20
ResiGlobRela = 1.0e-6
inst_fin = 13.0

l_inst = LIST.getValues()

phys_pb = code_aster.PhysicalProblem(MOTH, CHMAT)
phys_pb.addLoad(CHTH)

phys_pb.computeDOFNumbering()
phys_pb.computeBehaviourProperty({"RELATION": "THER_NL", "TOUT": "OUI"})

disc_comp = code_aster.DiscreteComputation(phys_pb)

# Create numbering
NU = phys_pb.getDOFNumbering()

# Fields
T = code_aster.FieldOnNodesReal(NU)
T.setValues({"TEMP": temp_init})

DELTA_T = T.copy() * 0.0

# initial residuals
_, _, RESI_THER_PREV = disc_comp.getInternalThermalForces(T, DELTA_T)
RESI_MASS_PREV = disc_comp.getNonLinearCapacityForces(T, DELTA_T)

linear_solver = code_aster.MultFrontSolver()

for i_inst in range(1, len(l_inst)):
    inst_prev = l_inst[i_inst - 1]
    inst = l_inst[i_inst]
    if inst >= inst_fin:
        break
    d_inst = inst - inst_prev

    print(" ")
    print("        ####################################")
    print("             instant de calcul ", inst)
    print("        ####################################")
    print(" ")

    print("IterNewton | Resi_Glob_rela   | Resi_Glob_Maxi  | Convergence")

    ####################################################################################################
    # Routine NXPRED
    ####################################################################################################

    # Compute non linear quantities
    mass_ther = disc_comp.getTangentCapacityMatrix(T, DELTA_T)

    rigi_ther = disc_comp.getTangentConductivityMatrix(T, DELTA_T, with_dual=False)
    rigi_ther_dual = disc_comp.getDualLinearConductivityMatrix()

    MT_AS = code_aster.AssemblyMatrixTemperatureReal()
    MT_AS.setDOFNumbering(NU)
    MT_AS.addElementaryMatrix(rigi_ther, theta)
    MT_AS.addElementaryMatrix(rigi_ther_dual)
    MT_AS.addElementaryMatrix(mass_ther, 1.0 / d_inst)
    MT_AS.assemble()

    # CHAR_THER_EVNL
    RESI_PREV = RESI_MASS_PREV / d_inst - (1.0 - theta) * RESI_THER_PREV

    # Linear loads - B * u
    CHAR_AS = disc_comp.getImposedDualBC(inst)

    # Factor matrix
    scaling = MT_AS.getLagrangeScaling()
    isFacto = linear_solver.factorize(MT_AS)
    assert isFacto

    # =========================================================
    #               BOUCLE DE NEWTON
    # =========================================================

    Conv = False

    for IterNewton in range(NbIterNewtonMax + 1):

        ####################################################################################################
        # Routine NXNEWT
        ####################################################################################################

        # Compute non linear quantities

        # Residual
        _, _, RESI_THER = disc_comp.getInternalThermalForces(T, DELTA_T)
        RESI_MASS = disc_comp.getNonLinearCapacityForces(T, DELTA_T)

        RESI_AS = RESI_MASS / d_inst + theta * RESI_THER

        # Reaction
        BTLA_AS = disc_comp.getDualForces(T + DELTA_T)

        DUAL_BC = disc_comp.getDualPrimal(T + DELTA_T, scaling)

        DUAL_AS = BTLA_AS + DUAL_BC - CHAR_AS

        ####################################################################################################
        # Routine NXRESI
        ####################################################################################################
        # Evaluate equilibrium
        CN2MBR = RESI_PREV - RESI_AS - DUAL_AS

        print("EVNL_AS: ", RESI_PREV.norm("NORM_2"), flush=True)
        print("RESI_AS: ", RESI_AS.norm("NORM_2"), flush=True)
        print("DUAL_AS: ", BTLA_AS.norm("NORM_2"), flush=True)

        # Evaluate residual
        resi_maxi = CN2MBR.norm("NORM_INFINITY")
        resi_rela = CN2MBR.norm("NORM_2")
        vnorm = (RESI_PREV - DUAL_AS).norm("NORM_2")

        Residu = resi_rela / vnorm
        ResiduX = resi_maxi

        # Estimation de la convergence
        if IterNewton > 0:
            if Residu <= ResiGlobRela:
                Conv = True

            print("     %d     |   %e   |   %e  |    %d " % (IterNewton, Residu, ResiduX, Conv))

        if Conv:
            RESI_THER_PREV = RESI_THER
            RESI_MASS_PREV = RESI_MASS
            break

        # Solve
        DDT = linear_solver.solve(CN2MBR)

        # Update TEMP increment
        DELTA_T += DDT

    if not Conv:
        raise code_aster.ConvergenceError("echec de la convergence des iterations de Newton")

    # Somme en python T+ = DELTA_t + T-
    T += DELTA_T
    DELTA_T.setValues(0)


T.printMedFile("/home/C00976/tmp/pynl03a.med")

FIN()
#
