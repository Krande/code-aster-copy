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

from ...Cata.Syntax import _F
from ...CodeCommands import EXTR_MODE, IMPR_RESU, NORM_MODE


def calc_modes_post(self, modes, lmatphys, norme_mode, filtre_mode, impression):
    """
    Macro-command CALC_MODES, post-treatment
    """

    # import the definitions of the commands to use in the macro-command
    # The name of the variable has to be the name of the command

    motscles = {}

    # for all the types
    if norme_mode is not None:
        modes = NORM_MODE(
            reuse=modes, MODE=modes, NORME=norme_mode["NORME"], INFO=norme_mode["INFO"]
        )

    # ONLY for (TYPE_RESU == 'DYNAMIQUE' + the matrices are physical)
    if lmatphys:
        # copy the modes concept in a temporary concept, in order to free its name
        __modes_temp = EXTR_MODE(FILTRE_MODE=_F(MODE=modes, TOUT_ORDRE="OUI"))

        impr_tout = False
        if filtre_mode is not None:
            motscles["FILTRE_MODE"] = _F(
                MODE=__modes_temp, CRIT_EXTR=filtre_mode["CRIT_EXTR"], SEUIL=filtre_mode["SEUIL"]
            )
        else:
            motscles["FILTRE_MODE"] = _F(MODE=__modes_temp, TOUT_ORDRE="OUI")
        if impression is not None:
            motscles["IMPRESSION"] = _F(
                CUMUL=impression["CUMUL"], CRIT_EXTR=impression["CRIT_EXTR"]
            )

        modes = EXTR_MODE(**motscles)

        if impression is not None:
            if impression["TOUT_PARA"] == "OUI":
                IMPR_RESU(
                    FORMAT="RESULTAT",
                    RESU=_F(RESULTAT=modes, TOUT_ORDRE="OUI", TOUT_CHAM="NON", TOUT_PARA="OUI"),
                )

    return modes
