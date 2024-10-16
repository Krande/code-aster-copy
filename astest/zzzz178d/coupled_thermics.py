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

import os.path as osp

from code_aster.Commands import *
from code_aster.Coupling import ExternalCoupling, MEDProj
from code_aster import CA
import medcoupling as MEDC


################################################################################
# definition of the coupled study
################################################################################
BASE = osp.dirname(__file__)
ther_med = osp.join(BASE, "mesh", "zzzz178b_ther.mmed")
interf_med = "fort.81"

# prepare MEDCoupling meshes
mfm = MEDC.MEDFileUMesh(ther_med)
mc_solid = mfm.getMeshAtLevel(0)
interf_ids = mfm.getGroupArr(0, "Volume")
interf_ids.setName("interf")
mc_interf = mc_solid[interf_ids]
mc_interf.setName("interface")

MEDC.WriteUMesh(interf_med, mc_interf, True)

cpl = ExternalCoupling("thermics", starter=True, debug=True)
cpl.setup(
    "mechanics",
    mc_interf,
    input_fields=[("DEPL", ["DX", "DY", "DZ"], "NODES")],
    output_fields=[("TEMP", ["TEMP"], "NODES")],
)


################################################################################
# setup the simulation
################################################################################
# send signal 6 (abort) to produce a traceback
test = CA.TestCase()

CA.init("--test", comm=cpl.comm, debug=False, ERREUR=_F(ERREUR_F="ABORT"))


ML = CA.Mesh()
ML.readMedFile(ther_med)

interf = CA.Mesh()
interf.readMedFile(interf_med)

# Assign thermal model on linear model
MODE_TL = AFFE_MODELE(MAILLAGE=ML, AFFE=_F(TOUT="OUI", MODELISATION="3D", PHENOMENE="THERMIQUE"))

model_meca = AFFE_MODELE(AFFE=_F(MODELISATION="3D", PHENOMENE="MECANIQUE", TOUT="OUI"), MAILLAGE=ML)

modinterf = AFFE_MODELE(
    AFFE=_F(MODELISATION="3D", PHENOMENE="MECANIQUE", TOUT="OUI"), MAILLAGE=interf
)

medp = MEDProj(mc_interf, interf_ids, 0, modinterf, model_meca)

# Assign thermal loading
CLIMT = AFFE_CHAR_THER(MODELE=MODE_TL, FLUX_REP=(_F(GROUP_MA="L1", FLUN=-400),))

# Assign thermal Dirichlet BC
BLOQT = AFFE_CHAR_CINE(MODELE=MODE_TL, THER_IMPO=_F(GROUP_MA="L2", TEMP=10))

# Define the thermal material
ACIER_T = DEFI_MATERIAU(THER=_F(LAMBDA=33.5, RHO_CP=526.0e4))

# Assign the thermal material
MATE_T = AFFE_MATERIAU(MAILLAGE=ML, AFFE=_F(TOUT="OUI", MATER=ACIER_T))

L_INST = DEFI_LIST_REEL(VALE=(-1.0, 0.0))


################################################################################
# define one iteration
################################################################################


cpl.ctxt["timedone"] = []


def exec_iteration(i_iter, current_time, delta_t, data, ctxt):
    """Execute one iteration.

    Arguments:
        i_iter (int): Iteration number.
        current_time (float): Current time.
        delta_t (float): Time step.
        data (list[*MEDCouplingField*]): List of input fields, on cells.

    Returns:
        list[*MEDCouplingField*]: Output fields, on nodes.
    """

    assert len(data) == 1, "expecting one field"

    # no coupling with displacement

    # Solve the thermal problem
    ctxt["result"] = THER_LINEAIRE(
        MODELE=MODE_TL,
        CHAM_MATER=MATE_T,
        EXCIT=(_F(CHARGE=CLIMT), _F(CHARGE=BLOQT)),
        # SOLVEUR=_F(METHODE="PETSC", PRE_COND="HPDDM", RESI_RELA=1.0e-9),
        SOLVEUR=_F(METHODE="PETSC", PRE_COND="LDLT_DP", RESI_RELA=1.0e-9),
        TYPE_CALCUL="STAT",
        INCREMENT=_F(LIST_INST=DEFI_LIST_REEL(VALE=current_time)),
        INFO=1,
    )

    temp = ctxt["result"].getField("TEMP", ctxt["result"].getLastIndex())
    mc_temp = medp.exportMEDCTemperature(temp, "TEMP")
    print("[Convert] Temperature field info:")
    print(mc_temp.simpleRepr(), flush=True)

    return True, {"TEMP": mc_temp}


################################################################################
# loop on time steps
################################################################################


cpl.run(exec_iteration, time_list=L_INST.getValues())

# analyse result

# Create the quadratic mesh
MQ = CREA_MAILLAGE(MAILLAGE=ML, LINE_QUAD=_F(TOUT="OUI"), INFO=1)

# Assign thermal model on quadratic model (to be used in the projection)
MODE_TQ = AFFE_MODELE(MAILLAGE=MQ, AFFE=_F(TOUT="OUI", MODELISATION="3D", PHENOMENE="THERMIQUE"))


# Project the field from the linear to the quadratic thermal model
TEMP2 = PROJ_CHAMP(
    METHODE="COLLOCATION", RESULTAT=self.ctxt["result"], MODELE_1=MODE_TL, MODELE_2=MODE_TQ
)


# Validate the result
temp2 = TEMP2.getField("TEMP", 1)
test.assertAlmostEqual(temp2.norm("NORM_INFINITY"), temp.norm("NORM_INFINITY"))

test.assertAlmostEqual(temp2.norm("NORM_2"), 412.46150094775635)

################################################################################
# Finalize the coupled study
################################################################################
cpl.finalize()

test.printSummary()

CA.close()
