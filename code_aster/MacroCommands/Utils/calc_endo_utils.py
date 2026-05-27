# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
Utilitaires pour CALC_ENDO
"""

from ...Messages import UTMESS
from ...CodeCommands import CALC_TABLE, FORMULE
from ...Behaviours import regu_visc_elas, endo_loca_tc, endo_fiss_tc


def get_obs_values(_ctrl_resu, _name_obs):
    """Get values of a given OBSERVATION

    Args:
        _ctrl_resu (*Table*): observation table
        _name_obs (str): name of the quantity to extract from the table

    Returns:
        _ctrl_obs_values (list): list of obervation table values
    """

    _ctrl_obs = CALC_TABLE(
        TABLE=_ctrl_resu,
        ACTION=(
            _F(OPERATION="FILTRE", NOM_PARA="NOM_OBSERVATION", VALE_K=_name_obs),
            _F(OPERATION="EXTR", NOM_PARA=("INST", "VALE")),
        ),
    )

    _ctrl_obs_values = _ctrl_obs.EXTR_TABLE().values()["VALE"]

    return _ctrl_obs_values


def concatenate_table(_table_out, _obs_table, _t_comp):
    """Create output table for CALC_ENDO

    Args:
        _table_out (*Table*): output table of CALC_ENDO
        _obs_table (*Table*): observation table
        _t_comp (float): physical inst for _obs_table

    Returns:
        _obs_table (*Table*): output table of CALC_ENDO
    """

    _obs_table = CALC_TABLE(
        TABLE=_obs_table,
        ACTION=(
            _F(OPERATION="RENOMME", NOM_PARA=("INST", "INST_VISC")),
            _F(OPERATION="SUPPRIME", NOM_PARA=("NUME_ORDRE")),
        ),
    )

    _obs_table = CALC_TABLE(
        reuse=_obs_table,
        TABLE=_obs_table,
        ACTION=_F(OPERATION="AJOUT_COLONNE", NOM_PARA="INST", VALE=_t_comp),
    )

    if _table_out:
        _table_out = CALC_TABLE(
            reuse=_table_out,
            TABLE=_table_out,
            ACTION=_F(OPERATION="COMB", TABLE=_obs_table, RESTREINT="NON", NOM_PARA="INST"),
        )
        return _table_out
    else:
        return _obs_table


def Get_Vari_By_Name(bhv, vi_name):
    """
    Return the Vxx name of an internal variable vi_name with respect to the constituive relation bhv

    Args:
        bhv (*Behaviour*): behaviour law
        vi_name (str): name of the internal variable

    Returns:
        num_name (str): number of internal variable

    """
    pos = bhv.loi.get_nom_vari().index(vi_name)
    num = pos + 1
    num_name = "V" + repr(num)
    return num_name


def set_default_observation(kwds):
    """Create default OBSERVATION syntax
    SIGMVISC for ENDO_LOCA_TC and ENDO_FISS_TC
    Viscous stress for REGU_VISC

    Args:
        kwds (*dict*): arguments for CALC_ENDO

    Returns:
        obs_stab_visc (tuple): OBSERVATION syntax for STAT_NON_LINE. Used for evaluation of stabilisation criteria
        crit_stab_visc (list): values of stabilisation criteria
        other_obs (tuple): OBSERVATION syntax for STAT_NON_LINE. NOT used for evaluation of stabilisation criteria

    """

    obs_stab_visc = []
    crit_stab_visc = []
    other_obs = []

    cham_mater = kwds["CHAM_MATER"]
    _mat_by_mesh = cham_mater.getMaterialsOnMeshEntities()
    _names_mesh_ent = [x[1].getNames() for x in _mat_by_mesh]

    for _compor in kwds["COMPORTEMENT"]:
        if (
            _compor["REGU_VISC"] == "OUI"
            or _compor["RELATION"] == "ENDO_LOCA_TC"
            or _compor["RELATION"] == "ENDO_FISS_TC"
        ):

            _ldc_name = _compor["RELATION"]
            if _ldc_name == "ENDO_LOCA_TC":
                _ldc = endo_loca_tc
            elif _ldc_name == "ENDO_FISS_TC":
                _ldc = endo_fiss_tc
            else:
                assert False

            assert len(_compor["GROUP_MA"]) == 1
            _grp_ma = _compor["GROUP_MA"][0]
            _index_mesh_ent = [_names_mesh_ent.index(x) for x in _names_mesh_ent if _grp_ma in x]
            assert len(_index_mesh_ent) == 1
            _mat_on_grp_ma = _mat_by_mesh[_index_mesh_ent[0]][0][0]

            if "RESI_REFE_RELA" not in kwds["CONVERGENCE"]:
                UTMESS("F", "CALCENDO_11")
            _resi_refe_rela = kwds["CONVERGENCE"]["RESI_REFE_RELA"]
            _ft = _mat_on_grp_ma.getValueReal(_ldc_name, "FT")

            if _compor["REGU_VISC"] == "OUI":

                _kvisc_elas = _mat_on_grp_ma.getValueReal("VISC_ELAS", "K")
                _nbvi_ldc = _ldc.loi.get_nb_vari()
                ##ATTENTION : il faudrait s'assurer qu'il n'y a jamais d'autres variables internes insérées entre celles de la LDC et REGU_VISC
                _nom_vari_ener_elas = "V%d" % (
                    _nbvi_ldc + int(Get_Vari_By_Name(regu_visc_elas, "VISCELAS")[1:])
                )
                _f_sigm_visc_elas = FORMULE(
                    NOM_PARA=_nom_vari_ener_elas,
                    VALE="(2*k*%s)**0.5" % _nom_vari_ener_elas,
                    k=_kvisc_elas,
                )

                _obs_visc = _F(
                    TITRE=_grp_ma + "_VISCELAS",
                    PAS_OBSE=1,
                    NOM_CHAM="VARI_ELGA",
                    NOM_CMP=_nom_vari_ener_elas,
                    GROUP_MA=_grp_ma,
                    EVAL_CMP="FORMULE",
                    FORMULE=_f_sigm_visc_elas,
                    EVAL_ELGA="MAX",
                    EVAL_CHAM="MAX",
                )
                obs_stab_visc.append(_obs_visc)
                crit_stab_visc.append(_resi_refe_rela * _ft)

            if _compor["RELATION"] == "ENDO_LOCA_TC" or _compor["RELATION"] == "ENDO_FISS_TC":
                _vari_endotot = Get_Vari_By_Name(_ldc, "ENDOTOT")
                _obs_endomoy = _F(
                    TITRE=_grp_ma + "_ENDOTOT",
                    PAS_OBSE=1,
                    NOM_CHAM="VARI_ELGA",
                    NOM_CMP=_vari_endotot,
                    GROUP_MA=_grp_ma,
                    EVAL_ELGA="MAX",
                    EVAL_CHAM="MOY",
                )
                other_obs.append(_obs_endomoy)

                _vari_visc = Get_Vari_By_Name(_ldc, "SIGMVISC")
                _obs_visc = _F(
                    TITRE=_grp_ma + "_VISCENDO",
                    PAS_OBSE=1,
                    NOM_CHAM="VARI_ELGA",
                    NOM_CMP=_vari_visc,
                    GROUP_MA=_grp_ma,
                    EVAL_ELGA="MAX",
                    EVAL_CHAM="MAX",
                )
                obs_stab_visc.append(_obs_visc)
                crit_stab_visc.append(_resi_refe_rela * _ft)

    return tuple(obs_stab_visc), crit_stab_visc, tuple(other_obs)
