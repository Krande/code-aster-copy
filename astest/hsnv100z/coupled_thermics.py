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
solid_med = osp.join(BASE, "mesh", "hsnv100a.mmed")
interf_med = "fort.81"

# prepare MEDCoupling meshes
mfm = MEDC.MEDFileUMesh(solid_med)
mc_solid = mfm.getMeshAtLevel(0)
interf_ids = mfm.getGroupArr(0, "M1")
interf_ids.setName("interf")
mc_interf = mc_solid[interf_ids]
mc_interf.setName("interface")

MEDC.WriteUMesh(interf_med, mc_interf, True)

cpl = ExternalCoupling("thermics", starter=True, debug=True)
cpl.setup(
    "mechanics",
    mc_interf,
    input_fields=[("DEPL", ["DX", "DY"], "NODES")],
    output_fields=[("TEMP", ["TEMP"], "NODES")],
)


################################################################################
# setup the simulation
################################################################################
# send signal 6 (abort) to produce a traceback
CA.init(comm=cpl.comm, debug=False, ERREUR=_F(ERREUR_F="ABORT"))


MAIL = CA.Mesh()
MAIL.readMedFile(solid_med)

interf = CA.Mesh()
interf.readMedFile(interf_med)

model = AFFE_MODELE(AFFE=_F(MODELISATION="AXIS", PHENOMENE="THERMIQUE", TOUT="OUI"), MAILLAGE=MAIL)

model_meca = AFFE_MODELE(
    AFFE=_F(MODELISATION="AXIS", PHENOMENE="MECANIQUE", TOUT="OUI"), MAILLAGE=MAIL
)
modinterf = AFFE_MODELE(
    AFFE=_F(MODELISATION="AXIS", PHENOMENE="MECANIQUE", TOUT="OUI"), MAILLAGE=interf
)

medp = MEDProj(mc_interf, interf_ids, 0, modinterf, model_meca)

MAT = DEFI_MATERIAU(THER=_F(LAMBDA=1.0e-3, RHO_CP=0.0e-3))

CM = AFFE_MATERIAU(MAILLAGE=MAIL, AFFE=_F(TOUT="OUI", MATER=MAT))

TIMPVAR = DEFI_FONCTION(NOM_PARA="INST", NOM_RESU="TEMP", VALE=(0.0e0, 0.0e0, 100.0e0, 100.0e0))


CHTHER = AFFE_CHAR_THER_F(
    MODELE=model,
    TEMP_IMPO=(
        _F(GROUP_NO="GRNO1", TEMP=TIMPVAR),
        _F(GROUP_NO="GRNO2", TEMP=TIMPVAR),
        _F(GROUP_NO="GRNO3", TEMP=TIMPVAR),
        _F(GROUP_NO="GRNO4", TEMP=TIMPVAR),
    ),
)

t1 = 66.666

t3 = 80.0

t5 = 90.0

t4 = 85.0

t2 = t1 + ((t3 - t1) / 2.0)


L_INST = DEFI_LIST_REEL(
    DEBUT=0.0,
    INTERVALLE=(_F(JUSQU_A=t1, NOMBRE=1), _F(JUSQU_A=t3, NOMBRE=2), _F(JUSQU_A=t5, NOMBRE=2)),
)


################################################################################
# define one iteration
################################################################################


class Context:
    """Context to be saved between iterations"""


ctxt = Context()
ctxt.timedone = [cpl.params.init_time]


def exec_iteration(i_iter, current_time, delta_t, data):
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
    mc_displ = data["DEPL"]

    # MEDC field => .med => code_aster field
    if i_iter > 1:
        warp = medp.importMEDCDisplacement(mc_displ)
        # distorsion of the mesh (very low impact on the thermal calculation!)
        MAIL = MODI_MAILLAGE(reuse=MAIL, MAILLAGE=MAIL, DEFORME=_F(OPTION="TRAN", DEPL=warp))

    ctxt.timedone.append(current_time)
    listr = DEFI_LIST_REEL(VALE=ctxt.timedone)

    opts = {}
    if i_iter > 1:
        opts["reuse"] = ctxt.result
        opts["RESULTAT"] = ctxt.result
        opts["ETAT_INIT"] = _F(EVOL_THER=ctxt.result)
    else:
        opts["ETAT_INIT"] = _F(CHAM_NO=T0)

    ctxt.result = THER_LINEAIRE(
        MODELE=model, CHAM_MATER=CM, EXCIT=_F(CHARGE=CHTHER), INCREMENT=_F(LIST_INST=listr), **opts
    )

    if i_iter > 1:
        # distorsion of the mesh (very low impact on the thermal calculation!)
        MAIL = MODI_MAILLAGE(reuse=MAIL, MAILLAGE=MAIL, DEFORME=_F(OPTION="TRAN", DEPL=-warp))

    temp = ctxt.result.getField("TEMP", ctxt.result.getLastIndex())
    mc_temp = medp.exportMEDCTemperature(temp, "TEMP")
    print("[Convert] Temperature field info:")
    print(mc_temp.simpleRepr(), flush=True)

    return True, {"TEMP": mc_temp}


################################################################################
# Initialization
################################################################################

T0 = CA.FieldOnNodesReal(model)
T0.setValues(0.0)

mc_temp = medp.exportMEDCTemperature(T0, "TEMP")
cpl.send_output_data({"TEMP": mc_temp})

################################################################################
# loop on time steps
################################################################################


cpl.run(exec_iteration, time_list=L_INST.getValues())

TEST_RESU(
    RESU=(
        _F(
            INST=90.0,
            RESULTAT=ctxt.result,
            NOM_CHAM="TEMP",
            GROUP_NO="N1",
            NOM_CMP="TEMP",
            VALE_CALC=90.0,
        ),
    )
)


################################################################################
# Finalize the coupled study
################################################################################
cpl.finalize()
