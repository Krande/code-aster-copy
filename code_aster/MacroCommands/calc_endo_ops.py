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

from ..Utilities import logger, no_new_attributes
from ..Messages import UTMESS, MasquerAlarme, RetablirAlarme
from ..Objects import EntityType, ExternalVariableTraits, ListOfFloats, TimesList
from libaster import MechanicalLoadFunction, Function

from ..CodeCommands import (
    DEFI_FONCTION,
    STAT_NON_LINE,
    RECU_TABLE,
    CREA_CHAMP,
    CREA_RESU,
    AFFE_MATERIAU,
)
from .Utils.calc_endo_utils import (
    get_obs_values,
    set_default_observation,
    concatenate_table,
    set_default_listinst,
)

from numpy import isclose, where, array


class CalcEndo:

    _kwds = None
    _tau = _visc_list_inst = _user_list_inst = None
    _t_init_ramp = _dt_stab = _t_fin = None
    _fixed_loads = _visc_loads = None
    _user_time_loads = _user_fixed_didi_loads = _user_time_didi_loads = None
    _fixed_varc = _user_time_varc = None
    _crit_stab_visc = _stab = _arret = None
    _obs_stab_visc = _other_obs = None
    _arch = _arch_visc = None
    _resu = _resu_visc = _tab_out = None
    _info_prefix = None

    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, args):
        """Initialization and definition of main arguments.

        Args:
            result (dict): User's keywords.
        """

        self._kwds = args.copy()
        self._visc_list_inst = self._kwds["ENDO_VISC"]["LIST_INST_VISC"]
        self._crit_stab_visc = self._kwds["ENDO_VISC"]["CRIT_STAB_VISC"]
        self._arret = self._kwds["ENDO_VISC"]["ARRET"]

        ##Initialisation of OBSERVATION for stabilisation criteria
        if self._crit_stab_visc:
            if len(self._kwds["ENDO_VISC"]["CRIT_STAB_VISC"]) != len(
                self._kwds["ENDO_VISC"]["OBSERVATION_VISC"]
            ):
                UTMESS("F", "CALCENDO_10")

            self._obs_stab_visc = ()
            self._other_obs = ()
            _titre_obs_visc = self._kwds["ENDO_VISC"]["OBSERVATION_VISC"]
            for _titre in _titre_obs_visc:
                _l_obs = [x for x in self._kwds["OBSERVATION"] if _titre == x["TITRE"]]
                if len(_l_obs) != 1:
                    UTMESS("F", "CALCENDO_3", valk=(_titre))
                self._obs_stab_visc = self._obs_stab_visc + tuple(_l_obs)
            self._kwds["OBSERVATION"] = tuple(
                [x for x in args["OBSERVATION"] if x["TITRE"] not in _titre_obs_visc]
            )
        else:
            self._obs_stab_visc, self._crit_stab_visc, self._other_obs = set_default_observation(
                self._kwds
            )

        self._resu_visc = []
        if "ARCHIVAGE_VISC" in self._kwds["ENDO_VISC"]:
            self._arch_visc = self._kwds["ENDO_VISC"]["ARCHIVAGE_VISC"]

        self._info_prefix = "Info CALC_ENDO : "
        if "OBSERVATION" not in args:
            self._kwds["OBSERVATION"] = ()

    def check_consistency(self):
        """Check consistency of viscosity material parameters."""

        _cham_mater = self._kwds["CHAM_MATER"]

        _tau_name = {
            "ENDO_LOCA_TC": "TAU_REGU_VISC",
            "ENDO_FISS_TC": "TAU_REGU_VISC",
            "VISC_ELAS": "TAU",
        }
        self._tau = 0.0

        for _mater in _cham_mater.getVectorOfMaterial():
            for _ldc in _tau_name:
                if _ldc in _mater.getMaterialNames():
                    _tau = _mater.getValueReal(_ldc, _tau_name[_ldc])
                    if _tau > self._tau:
                        if self._tau > 0.0:
                            UTMESS("A", "CALCENDO_1")
                        self._tau = _tau

        logger.info(
            self._info_prefix
            + "On a retenu TAU = "
            + str(self._tau)
            + " pour la création de la liste d'instants."
        )

        if not self._tau > 0.0:
            ##TODO : supprimer la valeur par défaut et repasser en message d'erreur ?
            self._tau = 1.0
            UTMESS("I", "CALCENDO_2")

    def create_user_list_inst(self):
        """Create list of physical load sequences and list of physical timestep for archivage"""

        full_user_list = self._kwds["INCREMENT"]["LIST_INST"].getValues()
        self._user_list_inst = full_user_list

        if "NUME_INST_FIN" in self._kwds["INCREMENT"] or "INST_FIN" in self._kwds["INCREMENT"]:
            if "NUME_INST_FIN" in self._kwds["INCREMENT"]:
                _nume_inst_fin = self._kwds["INCREMENT"]["NUME_INST_FIN"]
            else:
                _tol = self._kwds["INCREMENT"]["PRECISION"]
                _inst_fin = self._kwds["INCREMENT"]["INST_FIN"]
                _nume_inst_fin = (
                    where(isclose(array(full_user_list), _inst_fin, rtol=_tol))[0][0] + 1
                )
            self._user_list_inst = self._user_list_inst[:_nume_inst_fin]

        if "NUME_INST_INIT" in self._kwds["INCREMENT"] or "INST_INIT" in self._kwds["INCREMENT"]:
            if "NUME_INST_INIT" in self._kwds["INCREMENT"]:
                _nume_inst_init = self._kwds["INCREMENT"]["NUME_INST_INIT"]
            else:
                _tol = self._kwds["INCREMENT"]["PRECISION"]
                _inst_init = self._kwds["INCREMENT"]["INST_INIT"]
                _nume_inst_init = where(isclose(array(full_user_list), _inst_init, rtol=_tol))[0][0]
            self._user_list_inst = self._user_list_inst[_nume_inst_init:]

        logger.info(
            self._info_prefix
            + "Une séquence de chargement fictive (rampe et stabilisation) sera calculée pour chacun des instants physiques suivants : "
            + str(self._user_list_inst)
        )

        if "PAS_ARCH" in self._kwds["ARCHIVAGE"]:
            self._arch = self._user_list_inst[:: self._kwds["ARCHIVAGE"]["PAS_ARCH"]]
        if "INST" in self._kwds["ARCHIVAGE"]:
            self._arch = self._kwds["ARCHIVAGE"]["INST"]
        if "LIST_INST" in self._kwds["ARCHIVAGE"]:
            self._arch = self._kwds["ARCHIVAGE"]["LIST_INST"].getValues()

    def create_endo_list_inst(self):
        """Create list of timestep for a load sequence (ramp and stabilisation)"""

        _values_visc = self._visc_list_inst.getValues()
        if len(_values_visc) != 3:
            UTMESS("F", "CALCENDO_9")

        self._t_init_ramp = -self._tau * _values_visc[0]

        if abs(self._t_init_ramp) > self._tau:
            _rampe = [self._t_init_ramp, self._t_init_ramp + self._tau, 0.0]
        else:
            _rampe = [self._t_init_ramp, 0.0]

        self._dt_stab = self._tau * _values_visc[1]
        _nb_stab_max = int(_values_visc[2] / _values_visc[1])

        _l_fict_endo = _rampe + [
            _val
            for _i in range(_nb_stab_max)
            for _val in [_i * self._dt_stab + self._tau, _i * self._dt_stab + self._dt_stab]
        ]
        self._t_fin = _l_fict_endo[-1]

        logger.info(
            self._info_prefix
            + "Une séquence de chargement fictive (rampe et stabilisation) sera discrétisée par la liste d'instants suivante : "
            + str(_l_fict_endo)
        )

        if isinstance(self._kwds["ENDO_VISC"]["LIST_INST_VISC"], ListOfFloats):
            self._visc_list_inst = set_default_listinst(self._kwds, self._tau, _l_fict_endo)
        elif isinstance(self._kwds["ENDO_VISC"]["LIST_INST_VISC"], TimesList):
            ##TODO : copier le concept. Pour le moment, on modifie le DEFI_LIST_INST donné en entrée
            self._visc_list_inst = self._kwds["ENDO_VISC"]["LIST_INST_VISC"]  # .copy()
            self._visc_list_inst.setValues(_l_fict_endo)  # bug : issue35770

    def sort_loads(self):
        """Sort loads in two lists : time dependant loads and other loads"""

        self._fixed_loads = []
        self._user_time_loads = []
        self._user_time_didi_loads = []
        self._user_fixed_didi_loads = []

        for _load in self._kwds["EXCIT"]:
            ##Dépendance au temps via FONC_MULT
            if "FONC_MULT" in _load:
                if _load["TYPE_CHARGE"] == "DIDI":
                    self._user_time_didi_loads.append(_load)
                else:
                    self._user_time_loads.append(_load)
            ##Dépendance au temps via AFFE_CHAR_MECA_F ?
            elif isinstance(_load["CHARGE"], MechanicalLoadFunction):
                _find_function = False
                for _dependency in _load["CHARGE"].getDependencies():
                    if isinstance(_dependency, Function):
                        _param_function = _dependency.Parametres()["NOM_PARA"]
                        if _param_function == "INST":
                            UTMESS("F", "CALCENDO_4")
                        else:
                            self._fixed_loads.append(_load)
                            _find_function = True
                assert _find_function
            else:
                if _load["TYPE_CHARGE"] == "DIDI":
                    self._user_fixed_didi_loads.append(_load)
                else:
                    self._fixed_loads.append(_load)

        logger.info("")
        logger.info(
            self._info_prefix
            + "On a trouvé "
            + str(len(self._user_time_loads))
            + " chargement(s) fonction du temps, et "
            + str(len(self._fixed_loads))
            + " chargement(s) indépendant(s) du temps"
        )
        logger.info("")

    def sort_varc(self):
        """Sort varc in two lists : time dependant varc and other varc"""

        self._user_time_varc = []
        self._fixed_varc = []
        _user_mat = self._kwds["CHAM_MATER"]

        if _user_mat.hasExternalStateVariable():
            for _varc in _user_mat.getExtStateVariablesOnMeshEntities():
                _name = ExternalVariableTraits.getExternVarTypeStr(_varc[0].getType())
                logger.info(self._info_prefix + "Variable de commande : " + str(_name))
                _field = _varc[0].getField()
                _transient = _varc[0].getTransientResult()
                assert _field or _transient

                if _field:
                    logger.info(self._info_prefix + _name + " indépendant du temps.")
                    self._fixed_varc.append(_varc)

                if _transient:
                    logger.info(self._info_prefix + _name + " fonction du temps.")
                    self._user_time_varc.append(_varc)

            logger.info("")

        else:
            logger.info(self._info_prefix + "Aucune variable de commande détectée.")
            logger.info("")

        if len(self._user_time_loads) == 0 and len(self._user_time_varc) == 0:
            UTMESS("A", "CALCENDO_5")

    def init_sequence(self, _nume_ordre, _t_init, _t_comp):
        """Initialisation of a new load sequence

        Args:
            _nume_ordre (int): nume_ordre of last computed load sequence
            _t_init (float): initial time of last computed load sequence
            _t_comp (float): end time of new load sequence

        Returns:
            _nume_ordre (int): nume_ordre of new load sequence
            _t_init (float): initial time of new load sequence

        """

        if _nume_ordre is None:
            _nume_ordre = 0
            _t_init = self._user_list_inst[_nume_ordre]
        else:
            _nume_ordre += 1
            if _t_init is None:
                _t_init = self._user_list_inst[_nume_ordre]
        self._stab = False

        logger.info(self._info_prefix + "Séquence de chargement : " + str(_nume_ordre))
        logger.info(
            self._info_prefix
            + "Instant physique initial de cette séquence de chargement : "
            + str(_t_init)
        )
        logger.info(
            self._info_prefix
            + "Instant physique final de cette séquence de chargement : "
            + str(_t_comp)
        )

        return _nume_ordre, _t_init

    def set_init_state(self):
        """Prepare initial state of non linear computation

        Returns:
            _nume_ordre (int): nume_ordre of new load sequence
            _t_init (float): initial time of new load sequence
            _depl_init (*FieldOnNodes*): depl field initial state
            _sief_init (*FieldOnCells*): stress field initial state
            _vari_init (*FieldOnCells*): internal variable field initial state
            _strx_init (*FieldOnCells*): initial state for a few structural elements

        """

        _nume_ordre = _t_init = _depl_init = _sief_init = _vari_init = _strx_init = None

        if "ETAT_INIT" in self._kwds:
            if "EVOL_NOLI" in self._kwds["ETAT_INIT"]:
                self._resu = self._kwds["ETAT_INIT"]["EVOL_NOLI"]

                if (
                    "NUME_ORDRE" in self._kwds["ETAT_INIT"]
                    or "INST" in self._kwds["ETAT_INIT"]
                    or "NUME_DIDI" in self._kwds["ETAT_INIT"]
                ):
                    UTMESS("F", "CALCENDO_8")

                if (
                    "NUME_ORDRE" not in self._kwds["ETAT_INIT"]
                    and "INST" not in self._kwds["ETAT_INIT"]
                ):
                    _idx = self._resu.getAccessParameters()["NUME_ORDRE"][-1]
                    _t_init = self._resu.getAccessParameters()["INST"][-1]
                    _depl_init = self._resu.getField("DEPL", _idx)
                    _sief_init = self._resu.getField("SIEF_ELGA", _idx)
                    _vari_init = self._resu.getField("VARI_ELGA", _idx)
                    if "STRX_ELGA" in self._resu.getFieldsNames():
                        _strx_init = self._resu.getField("STRX_ELGA", _idx)
                    else:
                        _strx_init = None
                    _index_init = where(isclose(array(self._user_list_inst), _t_init))[0][0]
                    self._user_list_inst = self._user_list_inst[_index_init:]
            else:
                if "DEPL" in self._kwds["ETAT_INIT"]:
                    _depl_init = self._kwds["ETAT_INIT"]["DEPL"]
                if "SIGMA" in self._kwds["ETAT_INIT"]:
                    _sief_init = self._kwds["ETAT_INIT"]["SIGMA"]
                if "VARI" in self._kwds["ETAT_INIT"]:
                    _vari_init = self._kwds["ETAT_INIT"]["VARI"]
                if "STRX" in self._kwds["ETAT_INIT"]:
                    _strx_init = self._kwds["ETAT_INIT"]["STRX"]

        return _nume_ordre, _t_init, _depl_init, _sief_init, _vari_init, _strx_init

    def eval_loads(self, _stab_seq, _t_init, _t_comp, _nume_ordre):
        """Evaluate time dependant loads for a given load sequence

        Args:
            _stab_seq (bool): False if ramp, True if stabilisation
            _t_init (float): initial time of current load sequence
            _t_comp (float): end time of current load sequence
            _nume_ordre (int): local value of NUME_ORDRE

        """

        self._visc_loads = []
        ##TYPE_CHARGE = "FIXE_CSTE"or "SUIV" + FONC_MULT
        for _load in self._user_time_loads:
            _ramp_load = DEFI_FONCTION(
                NOM_PARA="INST",
                ABSCISSE=(self._t_init_ramp, 0),
                ORDONNEE=(_load["FONC_MULT"](_t_init), _load["FONC_MULT"](_t_comp)),
                PROL_DROITE="CONSTANT",
            )

            _endo_load = {key: item for key, item in _load.items() if key != "FONC_MULT"}
            _endo_load["FONC_MULT"] = _ramp_load

            self._visc_loads.append(_endo_load)

        ##TYPE_CHARGE = "DIDI" (no FONC_MULT)
        for _load in self._user_fixed_didi_loads:
            _endo_load = {key: item for key, item in _load.items()}
            if _nume_ordre > 0 or _stab_seq:
                ##Stabilisation : load set to zero
                _ramp_load = DEFI_FONCTION(
                    NOM_PARA="INST",
                    ABSCISSE=(self._t_init_ramp, 0),
                    ORDONNEE=(0, 0),
                    PROL_DROITE="CONSTANT",
                )
                _endo_load["FONC_MULT"] = _ramp_load
            self._visc_loads.append(_endo_load)

        ##TYPE_CHARGE = "DIDI" + FONC_MULT
        for _load in self._user_time_didi_loads:
            _endo_load = {key: item for key, item in _load.items() if key != "FONC_MULT"}
            if _stab_seq:
                ##Stabilisation : load set to zero
                _ramp_load = DEFI_FONCTION(
                    NOM_PARA="INST",
                    ABSCISSE=(self._t_init_ramp, 0),
                    ORDONNEE=(0, 0),
                    PROL_DROITE="CONSTANT",
                )
            else:
                _ramp_load = DEFI_FONCTION(
                    NOM_PARA="INST",
                    ABSCISSE=(self._t_init_ramp, 0),
                    ORDONNEE=(0, _load["FONC_MULT"](_t_comp) - _load["FONC_MULT"](_t_init)),
                    PROL_DROITE="CONSTANT",
                )
            _endo_load["FONC_MULT"] = _ramp_load
            self._visc_loads.append(_endo_load)

    def eval_varc(self, _t_init, _t_comp):
        """Evaluate external state variables for a given load sequence

        Args:
            _t_init (float): initial time of current load sequence
            _t_comp (float): end time of current load sequence

        Returns:
            _visc_mat_field (*MaterialField*): material field with time dependant varc evaluated at _t_init and _t_comp

        """

        _l_affe_varc = []

        for _varc in self._fixed_varc + self._user_time_varc:
            _affe_varc = self.get_affe_varc_syntax(_varc, _t_init, _t_comp)
            _l_affe_varc.append(_affe_varc)
        MasquerAlarme("MATERIAL2_61")
        _visc_mat_field = AFFE_MATERIAU(
            MODELE=self._kwds["MODELE"], CHAM_MATER=self._kwds["CHAM_MATER"], AFFE_VARC=_l_affe_varc
        )
        RetablirAlarme("MATERIAL2_61")

        return _visc_mat_field

    def get_affe_varc_syntax(self, _varc, _t_init, _t_comp):
        """Get syntax for an external state variable

        Args:
            _varc (list): external state variable given by getExtStateVariablesOnMeshEntities
            _t_init (float): initial time of current load sequence
            _t_comp (float): end time of current load sequence

        Returns:
            _affe_varc : input for AFFE_VARC in AFFE_MATERIAU

        """
        _dict_varc = {}
        _extevarc = _varc[0]
        _meshvarc = _varc[1]

        _name = ExternalVariableTraits.getExternVarTypeStr(_extevarc.getType())
        _dict_varc["NOM_VARC"] = _name

        if _meshvarc.getType() == EntityType.GroupOfCellsType:
            _dict_varc["GROUP_MA"] = _meshvarc.getNames()
        else:
            _dict_varc["TOUT"] = "OUI"

        if _extevarc.isSetRefe():
            _dict_varc["VALE_REF"] = _extevarc.getReferenceValue()

        _field = _extevarc.getField()
        _evol = _extevarc.getEvolutionParameter()
        if _field:
            _dict_varc["CHAM_GD"] = _field
        if _evol:
            _transient = _evol.getTransientResult()

            ##FONC_INST interdit
            if _evol.getTimeFormula() or _evol.getTimeFunction():
                UTMESS("F", "CALCENDO_7")

            _field_init = _transient.interpolateField(
                _name, _t_init, left=_evol.getLeftExtension(), right=_evol.getRightExtension()
            )
            _field_comp = _transient.interpolateField(
                _name, _t_comp, left=_evol.getLeftExtension(), right=_evol.getRightExtension()
            )

            _ramp_varc = CREA_RESU(
                OPERATION="AFFE",
                TYPE_RESU=_transient.getType(),
                AFFE=(
                    _F(
                        NOM_CHAM=_evol.getFieldName(),
                        CHAM_GD=_field_init,
                        INST=self._t_init_ramp,
                        MODELE=self._kwds["MODELE"],
                    ),
                    _F(
                        NOM_CHAM=_evol.getFieldName(),
                        CHAM_GD=_field_comp,
                        INST=0,
                        MODELE=self._kwds["MODELE"],
                    ),
                    _F(
                        NOM_CHAM=_evol.getFieldName(),
                        CHAM_GD=_field_comp,
                        INST=self._t_fin,
                        MODELE=self._kwds["MODELE"],
                    ),
                ),
            )
            _dict_varc["EVOL"] = _ramp_varc

        assert _field or _evol

        _affe_varc = _F(_dict_varc)
        return _affe_varc

    def eval_stab_crit(self, _evol_endo):
        """Evaluate user-defined stabilisation criteria at the last timestep of a given result

        Args:
            _evol_endo (*evol_noli*): result of the current load sequence

        """

        _current_inst = _evol_endo.getAccessParameters()["INST"][-1]
        logger.info("")
        logger.info(
            self._info_prefix
            + "Instant courant de la séquence de stabilisation : "
            + str(_current_inst)
        )
        _ctrl_resu = RECU_TABLE(CO=_evol_endo, NOM_TABLE="OBSERVATION")

        for _obs in self._kwds["OBSERVATION"]:
            if ("EVAL_CHAM" in _obs) and (_obs["EVAL_CHAM"] != "VALE"):
                _name_obs = _obs["TITRE"]
                _v_obs = get_obs_values(_ctrl_resu, _name_obs)

                logger.info(self._info_prefix + _name_obs + " courant : " + str(_v_obs[-1]))
                logger.info(
                    self._info_prefix
                    + "Variation de "
                    + _name_obs
                    + " courante : "
                    + str(_v_obs[-1] - _v_obs[0])
                )

        _crit = False
        for _num_obs, _obs_visc in enumerate(self._obs_stab_visc):
            _name_obs_visc = _obs_visc["TITRE"]
            _v_obs_visc = get_obs_values(_ctrl_resu, _name_obs_visc)
            _crit += _v_obs_visc[-1] < self._crit_stab_visc[_num_obs]

            logger.info(
                self._info_prefix + "Evaluation du critère de stabilisation pour " + _name_obs_visc
            )
            logger.info(self._info_prefix + _name_obs_visc + " courant : " + str(_v_obs_visc[-1]))
            logger.info(
                self._info_prefix
                + "Ratio "
                + _name_obs_visc
                + " sur seuil de stabilité : "
                + str(_v_obs_visc[-1] / self._crit_stab_visc[_num_obs])
            )

        self._stab = _crit == len(self._obs_stab_visc)
        logger.info("")
        if self._stab:
            logger.info(self._info_prefix + "Le critère de stabilisation global est atteint.")
        else:
            logger.info(self._info_prefix + "Le critère de stabilisation global n'est pas atteint.")

        if self._arret == "NON" and _current_inst == self._t_fin:
            self._stab = True
            UTMESS("A", "CALCENDO_6")

    def compute_ramp(self, _depl_init, _sief_init, _vari_init, _strx_init, _visc_mat_field):
        """Non linear computation during the ramp part of the load sequence

        Args:
            _depl_init (*FieldOnNodes*): depl field initial state
            _sief_init (*FieldOnCells*): stress field initial state
            _vari_init (*FieldOnCells*): internal variable field initial state
            _strx_init (*FieldOnCells*): initial state for a few structural elements
            _visc_mat_field (*MaterialField*) : material field with time dependant varc evaluated for load sequence

        Returns:
            _evol_endo (*evol_noli*): result of the current load sequence (ramp only)

        """

        _params_snl = {
            key: item
            for key, item in self._kwds.items()
            if key
            not in [
                "ENDO_VISC",
                "EXCIT",
                "INCREMENT",
                "ARCHIVAGE",
                "ETAT_INIT",
                "OBSERVATION",
                "CHAM_MATER",
                "reuse",
            ]
        }

        _params_snl["CHAM_MATER"] = _visc_mat_field

        if ("ETAT_INIT" in self._kwds and "EVOL_NOLI" in self._kwds["ETAT_INIT"]) or (
            _depl_init or _sief_init or _vari_init or _strx_init
        ):
            l_ETAT_INIT = _F(DEPL=_depl_init, SIGM=_sief_init, VARI=_vari_init, STRX=_strx_init)
        else:
            l_ETAT_INIT = []

        _params_snl["ETAT_INIT"] = l_ETAT_INIT
        _params_snl["INCREMENT"] = _F(
            LIST_INST=self._visc_list_inst, INST_FIN=0.0, NUME_INST_INIT=0
        )
        _params_snl["OBSERVATION"] = (
            self._kwds["OBSERVATION"] + self._obs_stab_visc + self._other_obs
        )
        _params_snl["EXCIT"] = self._fixed_loads + self._visc_loads
        _params_snl["ARCHIVAGE"] = self._arch_visc

        _evol_endo = STAT_NON_LINE(**_params_snl)
        self.eval_stab_crit(_evol_endo)

        return _evol_endo

    def compute_stab(self, _evol_endo, _visc_mat_field):
        """Non linear computation of one stabilisation sequence

        Args:
            _evol_endo (*evol_noli*): result of the current load sequence (ramp only)
            _visc_mat_field (*MaterialField*) : material field with time dependant varc evaluated for load sequence

        Returns:
            _evol_endo (*evol_noli*): result of the current load sequence (ramp and stab)

        """
        _current_inst = _evol_endo.getAccessParameters()["INST"][-1]
        _next_inst = _current_inst + self._dt_stab

        _params_snl = {
            key: item
            for key, item in self._kwds.items()
            if key
            not in [
                "ENDO_VISC",
                "EXCIT",
                "INCREMENT",
                "ARCHIVAGE",
                "ETAT_INIT",
                "OBSERVATION",
                "CHAM_MATER",
                "reuse",
            ]
        }

        _params_snl["CHAM_MATER"] = _visc_mat_field
        _params_snl["ETAT_INIT"] = _F(EVOL_NOLI=_evol_endo)
        _params_snl["INCREMENT"] = _F(LIST_INST=self._visc_list_inst, INST_FIN=_next_inst)
        _params_snl["OBSERVATION"] = (
            self._kwds["OBSERVATION"] + self._obs_stab_visc + self._other_obs
        )
        _params_snl["EXCIT"] = self._fixed_loads + self._visc_loads
        _params_snl["ARCHIVAGE"] = self._arch_visc

        _evol_endo = STAT_NON_LINE(reuse=_evol_endo, **_params_snl)

        self.eval_stab_crit(_evol_endo)

        return _evol_endo

    def arch_resu(self, _evol_endo, _t_comp):
        """Save results for output

        Args:
            _evol_endo (*evol_noli*): result of the current load sequence
            _t_comp (float): physical time at the end of current load sequence

        Returns:
            _t_init (float): reset _t_init to None
            _depl_arch (*FieldOnNodes*): depl field at the end of load sequence
            _sief_arch (*FieldOnCells*): stress field at the end of load sequence
            _vari_arch (*FieldOnCells*): internal variable field at the end of load sequence
            _strx_arch (*FieldOnCells*): strx field at the end of load sequence

        """

        _nume = _evol_endo.getAccessParameters()["NUME_ORDRE"][-1]
        _depl_arch = CREA_CHAMP(
            TYPE_CHAM="NOEU_DEPL_R",
            OPERATION="EXTR",
            NOM_CHAM="DEPL",
            RESULTAT=_evol_endo,
            NUME_ORDRE=_nume,
        )
        _sief_arch = CREA_CHAMP(
            TYPE_CHAM="ELGA_SIEF_R",
            OPERATION="EXTR",
            NOM_CHAM="SIEF_ELGA",
            RESULTAT=_evol_endo,
            NUME_ORDRE=_nume,
        )
        _vari_arch = CREA_CHAMP(
            TYPE_CHAM="ELGA_VARI_R",
            OPERATION="EXTR",
            NOM_CHAM="VARI_ELGA",
            RESULTAT=_evol_endo,
            NUME_ORDRE=_nume,
        )

        if "STRX_ELGA" in _evol_endo.getFieldsNames():
            _strx_arch = CREA_CHAMP(
                TYPE_CHAM="ELGA_STRX_R",
                OPERATION="EXTR",
                NOM_CHAM="STRX_ELGA",
                RESULTAT=_evol_endo,
                NUME_ORDRE=_nume,
            )
        else:
            _strx_arch = None

        _ctrl_obs = RECU_TABLE(CO=_evol_endo, NOM_TABLE="OBSERVATION")
        self._tab_out = concatenate_table(self._tab_out, _ctrl_obs, _t_comp)

        _arch_t_comp = (_t_comp == self._user_list_inst[-1]) or (_t_comp in self._arch)
        if _arch_t_comp:
            _d_affe = {"MODELE": self._kwds["MODELE"], "CHAM_MATER": self._kwds["CHAM_MATER"]}
            if "CARA_ELEM" in self._kwds:
                _d_affe["CARA_ELEM"] = self._kwds["CARA_ELEM"]

            l_affe = (
                _F(NOM_CHAM="DEPL", CHAM_GD=_depl_arch, INST=_t_comp, **_d_affe),
                _F(NOM_CHAM="SIEF_ELGA", CHAM_GD=_sief_arch, INST=_t_comp, **_d_affe),
                _F(NOM_CHAM="VARI_ELGA", CHAM_GD=_vari_arch, INST=_t_comp, **_d_affe),
            )

            if _strx_arch:
                l_affe += _F(_F(NOM_CHAM="STRX_ELGA", CHAM_GD=_strx_arch, INST=_t_comp, **_d_affe))

            self._resu = CREA_RESU(
                reuse=self._resu,
                OPERATION="AFFE",
                TYPE_RESU="EVOL_NOLI",
                AFFE=l_affe,
                COMPORTEMENT=self._kwds["COMPORTEMENT"],
            )

        if self._arch_visc:
            self._resu_visc.append(_evol_endo)

        return None, _depl_arch, _sief_arch, _vari_arch, _strx_arch


def calc_endo_ops(self, **args):
    """Execute the command.

    Arguments:
        **args (dict): User's keywords.

    Returns:
        _calc_endo._resu (*evol_noli*): Non linear result (physical time)
        _calc_endo._tab_out (*Table*): Observation table
        _calc_endo._resu_visc (list): List of evol_noli for viscous time discretisation. Optionnal.

    """

    _calc_endo = CalcEndo(args)
    _calc_endo.check_consistency()
    _calc_endo.create_user_list_inst()
    _calc_endo.create_endo_list_inst()
    _calc_endo.sort_loads()
    _calc_endo.sort_varc()

    _nume_ordre, _t_init, _depl_init, _sief_init, _vari_init, _strx_init = (
        _calc_endo.set_init_state()
    )

    for _t_comp in _calc_endo._user_list_inst[1:]:
        _nume_ordre, _t_init = _calc_endo.init_sequence(_nume_ordre, _t_init, _t_comp)
        _calc_endo.eval_loads(False, _t_init, _t_comp, _nume_ordre)
        _visc_mat_field = _calc_endo.eval_varc(_t_init, _t_comp)

        _evol_endo = _calc_endo.compute_ramp(
            _depl_init, _sief_init, _vari_init, _strx_init, _visc_mat_field
        )

        if not _calc_endo._stab:
            _calc_endo.eval_loads(True, _t_init, _t_comp, _nume_ordre)

        while not _calc_endo._stab:
            _evol_endo = _calc_endo.compute_stab(_evol_endo, _visc_mat_field)

        _t_init, _depl_init, _sief_init, _vari_init, _strx_init = _calc_endo.arch_resu(
            _evol_endo, _t_comp
        )

    if _calc_endo._arch_visc:
        return _calc_endo._resu, _calc_endo._tab_out, _calc_endo._resu_visc
    else:
        return _calc_endo._resu, _calc_endo._tab_out  # , None
