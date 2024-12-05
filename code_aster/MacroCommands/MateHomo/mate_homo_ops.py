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

"""
Calcul de proprietés homo
"""
from ...Cata.Syntax import _F
from ...CodeCommands import IMPR_RESU
from ...Messages import ASSERT, UTMESS

from .mate_homo_utilities import setup_calcul
from .mate_homo_mesh import prepare_mesh_syme
from .mate_homo_massif import calc_tabpara_massif, calc_corr_massif_syme
from .mate_homo_plaque import calc_tabpara_plaque, calc_corr_plaque_syme


def mate_homo_ops(self, **kwargs):
    meshin = kwargs.get("MAILLAGE")
    ls_affe = kwargs.get("AFFE")
    ls_varc = kwargs.get("VARC")
    unit = kwargs.get("UNITE")
    verbose = kwargs.get("INFO") == 2
    type_homo = kwargs.get("TYPE_HOMO")

    dir_plaque = kwargs.get("VECT_NORM")

    affe_all = any("TOUT" in i for i in ls_affe)
    affe_groups = list(
        set([j for i in [k for k in ls_affe if "GROUP_MA" in k] for j in i["GROUP_MA"]])
    )

    varc_name = ls_varc["NOM_VARC"]
    varc_values = sorted(ls_varc["VALE"])

    mesh, group_tout, volume_ver, dirthick = prepare_mesh_syme(meshin, affe_groups, affe_all)

    DEPLMATE, MODME, CHMATME, MODTH, CHMATTH, L_INST, alpha_calc = setup_calcul(
        type_homo, mesh, (group_tout,), ls_affe, varc_name, varc_values
    )

    if type_homo in ("MASSIF",):
        elas_fields, ther_fields = calc_corr_massif_syme(
            MODME, CHMATME, MODTH, CHMATTH, L_INST, alpha_calc, (group_tout,)
        )

        A_hom, K_hom, tabpara = calc_tabpara_massif(
            DEPLMATE,
            volume_ver,
            (group_tout,),
            varc_name,
            varc_values,
            **elas_fields,
            **ther_fields
        )

    elif type_homo in ("PLAQUE",):
        elas_fields, ther_fields = calc_corr_plaque_syme(
            MODME, CHMATME, MODTH, CHMATTH, L_INST, (group_tout,), dir_plaque
        )

        C_hom, D_hom, G_hom, tabpara = calc_tabpara_plaque(
            DEPLMATE,
            volume_ver,
            (group_tout,),
            varc_name,
            varc_values,
            dir_plaque,
            dirthick,
            **elas_fields,
            **ther_fields
        )
    else:
        ASSERT(False)

    # Return correctors
    if kwargs.get("CORR_MECA") is not None:
        self.register_result(elas_fields, kwargs.get("CORR_MECA"))

    if kwargs.get("CORR_THER") is not None:
        self.register_result(ther_fields, kwargs.get("CORR_THER"))

    # Save MED
    med_name = (
        lambda s: s.replace("MECA", "ME")
        .replace("THER", "TH")
        .replace("DILA", "DIL")
        .replace("CORR_", "")
        .replace("MEMB", "ME")
        .replace("FLEX", "FL")
        .replace("PINT", "PINT")
    )

    save_fields = dict(**elas_fields, **ther_fields)
    if unit is not None:
        IMPR_RESU(
            FORMAT="MED",
            RESU=[
                _F(RESULTAT=fld, NOM_RESU_MED=med_name(fldname))
                for fldname, fld in save_fields.items()
            ],
            UNITE=unit,
        )

    return tabpara
