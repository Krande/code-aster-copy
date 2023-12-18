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

# person_in_charge: francesco.bettonte at edf.fr

import os.path as osp

from ...Cata.Syntax import _F
from ...Objects import PhysicalProblem, FieldOnCellsReal, FieldOnNodesReal
from ...Commands import CREA_RESU
from ...Messages import UTMESS
from ...Utilities import ExecutionParameter
from .mac3coeur_ac_permute import MACRO_AC_PERMUTE
from .mac3coeur_coeur import CoeurFactory


def perm_mac3coeur_ops(self, **args):
    """Corps principal de la macro pour la permutation des AC dans MAC3COEUR"""

    _type_coeur_N = args.get("TYPE_COEUR_N")
    _type_coeur_P = args.get("TYPE_COEUR_NP1")
    _long_ligne_N = None
    _long_ligne_P = None
    if _type_coeur_N[:5] == "LIGNE":
        _long_ligne_N = args.get("NB_ASSEMBLAGE_N")
    if _type_coeur_P[:5] == "LIGNE":
        _long_ligne_P = args.get("NB_ASSEMBLAGE_NP1")
    _TAB_N = args.get("TABLE_N")

    l_RESUI = args.get("RESU_N")
    assert len(_TAB_N) == len(l_RESUI)
    l_last_i = []
    _l_MA_N = []
    for RESUI in l_RESUI:
        l_last_i.append(RESUI.LIST_PARA()["INST"][-1])
        # on recupere le concept maillage
        _MA_N = RESUI.getModel().getMesh()
        _l_MA_N.append(_MA_N)

    _l_coeur = []
    for _tabn1 in _TAB_N:
        _coeur = CoeurFactory.build(_type_coeur_N, _tabn1, _long_ligne_N)
        _l_coeur.append(_coeur)

    _TAB_NP1 = args.get("TABLE_NP1")
    _coeurp1 = CoeurFactory.build(_type_coeur_P, _TAB_NP1, _long_ligne_P)

    _MA1 = args.get("MAILLAGE_NP1")
    _MA_NP1 = _coeurp1.affectation_maillage(_MA1)
    _MO_NP1 = _coeurp1.affectation_modele(_MA_NP1)
    _coeurp1.init_from_mesh(_MA_NP1)
    _GFF_NP1 = _coeurp1.definition_geom_fibre()
    _CARANP1 = _coeurp1.definition_cara_coeur(_MO_NP1, _GFF_NP1)

    _fluence = 0.0
    _timep1 = _coeurp1.definition_time(_fluence, 1.0)
    _FLU_NP1 = _coeurp1.definition_fluence(_fluence, _MA_NP1, 0.0)
    _CHTHNP1 = _coeurp1.definition_champ_temperature(_MA_NP1)
    _AFSCNP1 = _coeurp1.definition_materiau(_MA_NP1, _GFF_NP1, _FLU_NP1, _CHTHNP1, CONTACT="NON")

    compor = (
        _F(
            RELATION="MULTIFIBRE",
            GROUP_MA=("CRAYON", "T_GUIDE"),
            PARM_THETA=0.5,
            DEFORMATION="PETIT",
        ),
        _F(RELATION="DIS_GRICRA", GROUP_MA="ELA"),
        _F(RELATION="DIS_CHOC", GROUP_MA=("CREIC", "RES_TOT")),
        _F(RELATION="ELAS", GROUP_MA=("CREI", "EBOINF", "EBOSUP", "RIG", "DIL")),
        _F(RELATION="VMIS_ISOT_TRAC", GROUP_MA="MAINTIEN", DEFORMATION="PETIT"),
    )

    phys_pb = PhysicalProblem(_MO_NP1, _AFSCNP1, _CARANP1)
    phys_pb.computeBehaviourProperty(compor)
    bhv_prop = phys_pb.getBehaviourProperty()

    __ASSDEP = FieldOnNodesReal(_MO_NP1)
    __ASSDEP.setValues(0.0)

    __ASSSIE = FieldOnCellsReal(_MO_NP1, "ELGA", "SIEF_R", bhv_prop, _CARANP1)
    __ASSSIE.setValues(0.0)

    __ASSVAR = FieldOnCellsReal(_MO_NP1, "ELGA", "VARI_R", bhv_prop, _CARANP1)
    __ASSVAR.setValues(0.0)

    BIDON = CREA_RESU(
        OPERATION="AFFE",
        TYPE_RESU="EVOL_NOLI",
        COMPORTEMENT=compor,
        NOM_CHAM="DEPL",
        AFFE=_F(
            CHAM_GD=__ASSDEP, INST=0.0, CHAM_MATER=_AFSCNP1, CARA_ELEM=_CARANP1, MODELE=_MO_NP1
        ),
    )

    CREA_RESU(
        reuse=BIDON,
        RESULTAT=BIDON,
        COMPORTEMENT=compor,
        OPERATION="AFFE",
        TYPE_RESU="EVOL_NOLI",
        NOM_CHAM="SIEF_ELGA",
        AFFE=_F(
            CHAM_GD=__ASSSIE, INST=0.0, CHAM_MATER=_AFSCNP1, CARA_ELEM=_CARANP1, MODELE=_MO_NP1
        ),
    )

    CREA_RESU(
        reuse=BIDON,
        RESULTAT=BIDON,
        COMPORTEMENT=compor,
        OPERATION="AFFE",
        TYPE_RESU="EVOL_NOLI",
        NOM_CHAM="VARI_ELGA",
        AFFE=_F(
            CHAM_GD=__ASSVAR, CHAM_MATER=_AFSCNP1, CARA_ELEM=_CARANP1, INST=0.0, MODELE=_MO_NP1
        ),
    )

    nbresu = len(l_RESUI)
    assert len(_l_coeur) == nbresu
    assert len(l_last_i) == nbresu
    assert len(_l_MA_N) == nbresu

    tran_x = 0.0
    indice = 1

    for nom in list(_coeurp1.nameAC.keys()):
        for i, _coeur in enumerate(_l_coeur):

            last_i = l_last_i[i]
            RESUI = l_RESUI[i]
            _MA_N = _l_MA_N[i]
            if nom in _coeur.nameAC:
                idx_z_p1, idx_z_p = _coeurp1.get_index(_coeurp1.nameAC[nom][0]), _coeurp1.get_index(
                    _coeur.nameAC[nom][0]
                )
                idx_y_p1, idx_y_p = _coeurp1.get_index(_coeurp1.nameAC[nom][2]), _coeurp1.get_index(
                    _coeur.nameAC[nom][2]
                )
                # print('index z : ', idx_z_p1, idx_z_p)
                # print('index y : ', idx_y_p1, idx_y_p)
                tran_z = _coeurp1.pas_assemblage * (idx_z_p1 - idx_z_p)
                tran_y = _coeurp1.pas_assemblage * (idx_y_p1 - idx_y_p)
                # print('tran_z, tran_y, tran_x = ', tran_z, tran_y, tran_x)
                # print('AC init, AC_fin = ', _coeur.nameAC[nom], _coeurp1.nameAC[nom])

                BIDON = MACRO_AC_PERMUTE(
                    POS_INIT=_coeur.nameAC[nom],
                    POS_FIN=_coeurp1.nameAC[nom],
                    RESU_INI=RESUI,
                    RESU_FIN=BIDON,
                    MAILLAGE_INIT=_MA_N,
                    INSTANT=last_i,
                    MAILLAGE_FINAL=_MA_NP1,
                    MODELE_FINAL=_MO_NP1,
                    TRAN=(tran_x, tran_y, tran_z),
                )
                valk = (
                    _coeur.position_todamac(_coeur.nameAC[nom]),
                    _coeurp1.position_todamac(_coeurp1.nameAC[nom]),
                )
                UTMESS("I", "COEUR0_3", valk=valk)
                indice += 1
                break

    UTMESS("I", "COEUR0_2", vali=(indice))

    return BIDON
