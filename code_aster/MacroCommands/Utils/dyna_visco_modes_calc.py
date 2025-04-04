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

import aster
import numpy as NP

from ...Cata.Syntax import _F
from ...CodeCommands import CALC_MODES, COMB_MATR_ASSE, CREA_CHAMP, CREA_RESU


def dyna_visco_modes_calc(
    self,
    TYPE_MODE,
    freq1,
    nmode,
    RESI_RELA,
    i,
    j,
    MATER_ELAS_FO,
    e0,
    eta0,
    __asseMg,
    __asseKgr,
    __asseKg,
    __listKv,
    trKg,
    ltrv,
    TYPE_RESU,
    reuse="non",
    **args
):
    """
    Macro-command DYNA_VISCO,
    function to compute with iterations one eigenmode,
    and store it
    """

    dfreq = freq1

    while abs(dfreq) >= RESI_RELA * freq1:
        if i > 10:
            nmode = nmode + 5
            i = 0

        if TYPE_MODE == "REEL":
            __asseKw = __asseKgr
        elif TYPE_MODE == "BETA_REEL":
            __asseKw = __asseKg
            betab = NP.real(trKg)
            betah = NP.imag(trKg)
        elif TYPE_MODE == "COMPLEXE":
            __asseKw = __asseKg
        else:
            assert False

        ny = 0
        for y in MATER_ELAS_FO:
            e = float(y["E"](freq1))
            eta = float(y["AMOR_HYST"](freq1))
            if TYPE_MODE == "REEL":
                __asseKw = COMB_MATR_ASSE(
                    COMB_R=(
                        _F(MATR_ASSE=__asseKw, COEF_R=1.0),
                        _F(MATR_ASSE=__listKv[ny], COEF_R=e / e0[ny] - 1.0),
                    )
                )

            if TYPE_MODE in ["BETA_REEL", "COMPLEXE"]:
                __asseKw = COMB_MATR_ASSE(
                    COMB_C=(
                        _F(MATR_ASSE=__asseKw, COEF_R=1.0),
                        _F(
                            MATR_ASSE=__listKv[ny],
                            COEF_C=(complex(e / e0[ny] - 1.0, eta * e / e0[ny] - eta0[ny])),
                        ),
                    )
                )
                if TYPE_MODE == "BETA_REEL":
                    betab = betab + (e / e0[ny] - 1.0) * ltrv[ny]
                    betah = betah + (eta * e / e0[ny] - eta0[ny]) * ltrv[ny]

            ny = ny + 1

        if TYPE_MODE == "BETA_REEL":
            __asseKw = COMB_MATR_ASSE(
                COMB_R=(
                    _F(MATR_ASSE=__asseKw, PARTIE="REEL", COEF_R=1.0),
                    _F(MATR_ASSE=__asseKw, PARTIE="IMAG", COEF_R=betah / betab),
                )
            )

        # IMPR_CO(CONCEPT=_F(NOM=__asseKw))

        __modtmp = CALC_MODES(
            MATR_RIGI=__asseKw,
            MATR_MASS=__asseMg,
            OPTION="CENTRE",
            CALC_FREQ=_F(FREQ=freq1, NMAX_FREQ=nmode),
            VERI_MODE=_F(STOP_ERREUR="OUI", SEUIL=1.0e-3, STURM="NON"),
        )

        freq2 = aster.GetResu(__modtmp.getName(), "VARI_ACCES")["FREQ"]
        dfreq = abs(freq1 - freq2[0])
        __numod = 0

        for ii in range(1, nmode):
            __new_dfreq = abs(freq1 - freq2[ii])
            if __new_dfreq < dfreq:
                dfreq = __new_dfreq
                __numod = ii

        freq1 = freq2[__numod]
        if TYPE_MODE == "COMPLEXE":
            amor_red1 = aster.GetResu(__modtmp.getName(), "PARAMETRES")["AMOR_REDUIT"][__numod]

        if __numod + 1 == nmode:
            nmode = nmode + 5
            dfreq = freq1

        i = i + 1

    if TYPE_MODE in ["REEL", "BETA_REEL"]:
        type_cham = "NOEU_DEPL_R"
    elif TYPE_MODE == "COMPLEXE":
        type_cham = "NOEU_DEPL_C"
    else:
        assert False

    # extract the modal shape
    __unmod = CREA_CHAMP(
        OPERATION="EXTR",
        NOM_CHAM="DEPL",
        TYPE_CHAM=type_cham,
        RESULTAT=__modtmp,
        NUME_ORDRE=__numod + 1,
    )

    motcles = {}

    if TYPE_MODE in ["REEL", "BETA_REEL"]:
        type_resu = "MODE_MECA"
        motcles["AFFE"] = _F(NOM_CHAM="DEPL", CHAM_GD=__unmod, NUME_MODE=j + 1, FREQ=freq1)
    elif TYPE_MODE == "COMPLEXE":
        type_resu = "MODE_MECA_C"
        motcles["AFFE"] = _F(
            NOM_CHAM="DEPL", CHAM_GD=__unmod, NUME_MODE=j + 1, FREQ=freq1, AMOR_REDUIT=amor_red1
        )
    else:
        assert False

    if reuse == "oui":
        motcles["reuse"] = args["co_reuse"]
        motcles["RESULTAT"] = args["co_reuse"]

    # fill the concept containing the eigenmodes
    _modes = CREA_RESU(OPERATION="AFFE", TYPE_RESU=type_resu, MATR_MASS=__asseMg, **motcles)

    freq1 = freq2[__numod + 1]
    return _modes, freq1, nmode
