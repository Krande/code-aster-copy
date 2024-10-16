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
# definition of the coupled cpl
################################################################################
BASE = osp.dirname(__file__)
meca_med = osp.join(BASE, "mesh", "zzzz178b_meca.mmed")
interf_med = "fort.88"

# prepare MEDCoupling MQes
mfm = MEDC.MEDFileUMesh(meca_med)
mc_solid = mfm.getMeshAtLevel(0)
interf_ids = mfm.getGroupArr(0, "Volume")
interf_ids.setName("interf")
mc_interf = mc_solid[interf_ids]
mc_interf.setName("interface")

MEDC.WriteUMesh(interf_med, mc_interf, True)


cpl = ExternalCoupling("mechanics", starter=False, debug=True)
cpl.setup(
    "thermics",
    mc_interf,
    input_fields=[("TEMP", ["TEMP"], "NODES")],
    output_fields=[("DEPL", ["DX", "DY", "DZ"], "NODES")],
)

################################################################################
# setup the simulation
################################################################################
# send signal 6 (abort) to produce a traceback
test = CA.TestCase()

CA.init("--test", comm=cpl.comm, debug=False, ERREUR=_F(ERREUR_F="ABORT"))


interf = CA.Mesh()
interf.readMedFile(interf_med)

MQ = CA.Mesh()
MQ.readMedFile(meca_med)

# Check the orientation of the boundary
MQ = MODI_MAILLAGE(reuse=MQ, MAILLAGE=MQ, ORIE_PEAU=_F(GROUP_MA_PEAU="L2"))

# Assign Mechanical model
MODE_MQ = AFFE_MODELE(MAILLAGE=MQ, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))


modinterf = AFFE_MODELE(
    AFFE=_F(MODELISATION="3D", PHENOMENE="THERMIQUE", TOUT="OUI"), MAILLAGE=interf
)

medp = MEDProj(mc_interf, interf_ids, 0, modinterf, MODE_MQ)

# Define the mechanical material
ACIER_M = DEFI_MATERIAU(ELAS=_F(E=2.0e12, NU=0.3e00, RHO=1.0e03, ALPHA=1.0e-4))


# Assign mechanical Dirichlet BC
DIRI = AFFE_CHAR_CINE(MODELE=MODE_MQ, MECA_IMPO=(_F(GROUP_MA="L1", DX=0, DY=0.0, DZ=0.0),))

# Assign mechanical loading
CHAR = AFFE_CHAR_MECA(MODELE=MODE_MQ, PRES_REP=_F(GROUP_MA="L2", PRES=-100.0))


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
        data (dict[*MEDCouplingField*]): dict of input fields, on cells.

    Returns:
        dict[*MEDCouplingField*]: Output fields, on nodes.
    """

    assert len(data) == 1, "expecting one field"
    mc_ther = data["TEMP"]

    # MEDC field => .med => code_aster field
    TEMPE = medp.importMEDCTemperature(mc_ther)

    ctxt["evol_ther"] = CREA_RESU(
        TYPE_RESU="EVOL_THER",
        OPERATION="AFFE",
        AFFE=_F(NOM_CHAM="TEMP", CHAM_GD=TEMPE, INST=current_time),
    )

    # Assign the mechanical material with the thermal field
    MATE_M = AFFE_MATERIAU(
        MAILLAGE=MQ,
        AFFE=_F(TOUT="OUI", MATER=ACIER_M),
        AFFE_VARC=_F(
            TOUT="OUI", EVOL=ctxt["evol_ther"], VALE_REF=0.0, NOM_VARC="TEMP", NOM_CHAM="TEMP"
        ),
    )

    # Solve the mechanical problem
    ctxt["result"] = MECA_STATIQUE(
        MODELE=MODE_MQ,
        CHAM_MATER=MATE_M,
        EXCIT=(_F(CHARGE=CHAR), _F(CHARGE=DIRI)),
        SOLVEUR=_F(METHODE="PETSC", PRE_COND="HPDDM"),
    )

    displ = ctxt["result"].getField("DEPL", ctxt["result"].getLastIndex())
    mc_displ = medp.exportMEDCDisplacement(displ, "Displ")
    print("[Convert] Displacement field info:")
    print(mc_displ.simpleRepr(), flush=True)

    return True, {"DEPL": mc_displ}


################################################################################
# loop on time steps
################################################################################

cpl.run(exec_iteration)


# Extract the field from the result
displ = cpl.ctxt["result"].getField("DEPL", 1)

# Validate the result againt sequential run
# norm = displ.norm("NORM_2")
norm = 0.10070829146943601
test.assertAlmostEqual(displ.norm("NORM_2"), norm)

################################################################################
# Finalize the coupled study
################################################################################
cpl.finalize()

test.printSummary()

CA.close()
