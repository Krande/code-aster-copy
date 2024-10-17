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
import medcoupling as MEDC

from code_aster.Commands import *
from code_aster.Coupling import ExternalCoupling, MEDProj
from code_aster import CA


def coupled_mechanics(cpl):
    """Run mechanical coupling.

    Arguments:
        cpl (ExternalCoupling): Mechanical coupler
    """

    ################################################################################
    # setup the simulation
    ################################################################################
    # send signal 6 (abort) to produce a traceback
    CA.init("--test", comm=cpl.comm, debug=False, ERREUR=_F(ERREUR_F="ABORT"))

    # Study directory (containing this file and the datafiles)
    STUDY = osp.dirname(__file__)

    solid_med = "hsnv100z.mmed"
    interf_med = "fort.88"

    # prepare MEDCoupling meshes
    mfm = MEDC.MEDFileUMesh(osp.join(STUDY, "hsnv100z.mmed"))
    mc_solid = mfm.getMeshAtLevel(0)
    interf_ids = mfm.getGroupArr(0, "M1")
    interf_ids.setName("interf")
    mc_interf = mc_solid[interf_ids]
    mc_interf.setName("interface")

    MEDC.WriteUMesh(interf_med, mc_interf, True)

    MAIL = CA.Mesh()
    MAIL.readMedFile(solid_med)

    interf = CA.Mesh()
    interf.readMedFile(interf_med)

    model = AFFE_MODELE(
        AFFE=_F(MODELISATION="AXIS", PHENOMENE="MECANIQUE", TOUT="OUI"), MAILLAGE=MAIL
    )

    modinterf = AFFE_MODELE(
        AFFE=_F(MODELISATION="AXIS", PHENOMENE="THERMIQUE", TOUT="OUI"), MAILLAGE=interf
    )

    cpl.setup(
        interface=(MAIL, ["M1"]),
        input_fields=[("TEMP", ["TEMP"], "NODES")],
        output_fields=[("DEPL", ["DX", "DY"], "NODES")],
    )

    medp = MEDProj(mc_interf, interf_ids, 0, modinterf, model)

    # DONNEES DE MODELISATION

    FCT1 = DEFI_FONCTION(
        NOM_PARA="EPSI",
        VALE=(0.200e-2, 400.0, 0.400e-2, 500.0),
        PROL_DROITE="LINEAIRE",
        PROL_GAUCHE="LINEAIRE",
    )

    #

    FCT2 = DEFI_FONCTION(
        NOM_PARA="EPSI",
        VALE=(0.100e-2, 200.0, 0.300e-2, 300.0),
        PROL_DROITE="LINEAIRE",
        PROL_GAUCHE="LINEAIRE",
    )

    #

    CTRACB = DEFI_NAPPE(
        NOM_PARA="TEMP",
        PARA=(0.0, 50.0),
        FONCTION=(FCT1, FCT2),
        PROL_DROITE="LINEAIRE",
        PROL_GAUCHE="LINEAIRE",
    )

    MAT = DEFI_MATERIAU(ELAS=_F(E=200.0e3, NU=0.3, ALPHA=10.0e-6), TRACTION=_F(SIGM=CTRACB))

    CHMECA = AFFE_CHAR_CINE(
        MODELE=model, MECA_IMPO=(_F(GROUP_NO="GRNO1", DY=0.0), _F(GROUP_NO="GRNO3", DY=0.0))
    )

    ################################################################################
    # define one iteration
    ################################################################################

    cpl.ctxt["timedone"] = [0.0]

    def exec_iteration(i_iter, current_time, delta_t, data, ctxt):
        """Execute one iteration.

        Arguments:
            i_iter (int): Iteration number.
            current_time (float): Current time.
            delta_t (float): Time step.
            data (dict[*MEDCouplingField*]): dict of input fields, on cells.
            ctxt (object): context of the computation

        Returns:
            dict[*MEDCouplingField*]: Output fields, on nodes.
        """

        assert len(data) == 1, "expecting one field"
        mc_ther = data["TEMP"]

        # MEDC field => .med => code_aster field
        TEMPE = medp.importMEDCTemperature(mc_ther)

        ctxt["evol_ther"] = CREA_RESU(
            reuse=ctxt["evol_ther"],
            RESULTAT=ctxt["evol_ther"],
            TYPE_RESU="EVOL_THER",
            OPERATION="AFFE",
            AFFE=_F(NOM_CHAM="TEMP", CHAM_GD=TEMPE, INST=current_time),
        )

        CTM = AFFE_MATERIAU(
            MAILLAGE=MAIL,
            AFFE=_F(TOUT="OUI", MATER=MAT),
            AFFE_VARC=_F(TOUT="OUI", NOM_VARC="TEMP", EVOL=ctxt["evol_ther"], VALE_REF=0.0),
        )

        ctxt["timedone"].append(current_time)
        listr = DEFI_LIST_REEL(VALE=ctxt["timedone"])

        opts = {}
        if i_iter > 1:
            opts["reuse"] = ctxt["result"]
            opts["RESULTAT"] = ctxt["result"]
            opts["ETAT_INIT"] = _F(EVOL_NOLI=ctxt["result"])

        ctxt["result"] = STAT_NON_LINE(
            MODELE=model,
            CHAM_MATER=CTM,
            COMPORTEMENT=_F(RELATION="VMIS_ISOT_TRAC"),
            EXCIT=_F(CHARGE=CHMECA),
            INCREMENT=_F(LIST_INST=listr),
            **opts,
        )

        displ = ctxt["result"].getField("DEPL", ctxt["result"].getLastIndex())
        mc_displ = medp.exportMEDCDisplacement(displ, "Displ")
        print("[Convert] Displacement field info:")
        print(mc_displ.simpleRepr(), flush=True)

        return True, {"DEPL": mc_displ}

    ################################################################################
    # Initialization
    ################################################################################

    input_data = cpl.recv_input_data()
    TEMPE = medp.importMEDCTemperature(input_data["TEMP"])

    cpl.ctxt["evol_ther"] = CREA_RESU(
        TYPE_RESU="EVOL_THER", OPERATION="AFFE", AFFE=_F(NOM_CHAM="TEMP", CHAM_GD=TEMPE, INST=0.0)
    )

    ################################################################################
    # loop on time steps
    ################################################################################

    cpl.run(exec_iteration)

    cpl.ctxt["result"] = CALC_CHAMP(
        reuse=cpl.ctxt["result"],
        RESULTAT=cpl.ctxt["result"],
        CONTRAINTE=("SIGM_ELNO",),
        DEFORMATION=("EPSI_ELNO",),
        VARI_INTERNE=("VARI_ELNO",),
    )

    TEST_RESU(
        RESU=(
            _F(
                INST=66.665999999999997,
                RESULTAT=cpl.ctxt["result"],
                NOM_CHAM="EPSI_ELNO",
                GROUP_NO="N1",
                NOM_CMP="EPXX",
                VALE_CALC=8.66658e-4,
                GROUP_MA="M1",
            ),
            _F(
                INST=66.665999999999997,
                RESULTAT=cpl.ctxt["result"],
                NOM_CHAM="SIGM_ELNO",
                GROUP_NO="N1",
                NOM_CMP="SIYY",
                VALE_CALC=-133.332,
                GROUP_MA="M1",
            ),
            _F(
                INST=80.0,
                RESULTAT=cpl.ctxt["result"],
                NOM_CHAM="EPSI_ELNO",
                GROUP_NO="N2",
                NOM_CMP="EPZZ",
                VALE_CALC=1.1000000000000001e-3,
                GROUP_MA="M1",
            ),
            _F(
                INST=80.0,
                RESULTAT=cpl.ctxt["result"],
                NOM_CHAM="VARI_ELNO",
                GROUP_NO="N2",
                NOM_CMP="V1",
                VALE_CALC=2.9999999999999997e-4,
                GROUP_MA="M1",
            ),
            _F(
                INST=80.0,
                RESULTAT=cpl.ctxt["result"],
                NOM_CHAM="SIGM_ELNO",
                GROUP_NO="N2",
                NOM_CMP="SIYY",
                VALE_CALC=-100.0,
                GROUP_MA="M1",
            ),
            _F(
                INST=90.0,
                RESULTAT=cpl.ctxt["result"],
                NOM_CHAM="EPSI_ELNO",
                GROUP_NO="N3",
                NOM_CMP="EPZZ",
                VALE_CALC=1.2750000000000001e-3,
                GROUP_MA="M1",
            ),
            _F(
                INST=90.0,
                RESULTAT=cpl.ctxt["result"],
                NOM_CHAM="VARI_ELNO",
                GROUP_NO="N3",
                NOM_CMP="V1",
                VALE_CALC=5.2499999999999997e-4,
                GROUP_MA="M1",
            ),
            _F(
                INST=90.0,
                RESULTAT=cpl.ctxt["result"],
                NOM_CHAM="SIGM_ELNO",
                GROUP_NO="N3",
                NOM_CMP="SIYY",
                VALE_CALC=-75.0,
                GROUP_MA="M1",
            ),
        )
    )

    ################################################################################
    # Finalize the coupled study
    ################################################################################
    cpl.finalize()

    FIN()
