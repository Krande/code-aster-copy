# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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

import math

import numpy as NP

from ..Objects.function_py import t_fonction, t_fonction_c, t_nappe
from ..Cata.Syntax import _F
from ..CodeCommands import CALC_FONCTION, CREA_TABLE, IMPR_TABLE
from ..Objects.table_py import Table
from ..Messages import UTMESS


def info_fonction_ops(self, INFO, **args):
    """
    Ecriture de la macro INFO_FONCTION
    """

    # On importe les definitions des commandes a utiliser dans la macro

    # type de traitement

    #
    if "MAX" in args:
        MAX = args["MAX"]
        # liste des t_fonction
        l_cofonc = MAX["FONCTION"]
        if type(l_cofonc) not in (list, tuple):
            l_cofonc = [l_cofonc]
        l_fonc = [concept.convert() for concept in l_cofonc]

        # intervalles
        mc_interv = MAX["INTERVALLE"]
        with_intervalle = mc_interv is not None
        interv = []
        if with_intervalle:
            nbv = len(mc_interv)
            if nbv % 2 != 0:
                UTMESS("F", "FONCT0_55")
            tint = NP.array(mc_interv)
            tint.shape = (nbv // 2, 2)
            dx = tint[:, 1] - tint[:, 0]
            if min(dx) < 0.0:
                UTMESS("F", "FONCT0_56")
            interv = tint.tolist()

        # vérifications de cohérence
        typobj = set()
        npara = set()
        nparf = set()
        nresu = set()
        l_nom = []
        for tf in l_fonc:
            typobj.add(tf.__class__)
            npara.add(tf.para["NOM_PARA"])
            nparf.add(tf.para.get("NOM_PARA_FONC"))
            nresu.add(tf.para["NOM_RESU"])
            l_nom.append(tf.nom)
        if len(typobj) > 1:
            # types (fonction, fonction_c, nappe) non homogènes
            UTMESS("F", "FONCT0_37")
        is_nappe = typobj.pop() is t_nappe
        if len(npara) > 1:
            # NOM_PARA non homogènes
            UTMESS("F", "FONCT0_38", valk=" ".join(npara))
        if len(nparf) > 1:
            # NOM_PARA_FONC non homogènes
            UTMESS("F", "FONCT0_38", valk=" ".join(nparf))
        if len(nresu) > 1:
            # NOM_RESU non homogènes
            UTMESS("F", "FONCT0_39", valk=" ".join(nresu))

        # nom des paramètres et leurs types
        k_para = npara.pop()
        k_parf = nparf.pop()
        k_ordo = nresu.pop()
        k_min = k_para + "_MIN"
        k_max = k_para + "_MAX"
        ordered_params = ["FONCTION", "TYPE"]
        ordered_type = ["K8", "K8"]
        if with_intervalle:
            ordered_params.extend(["INTERVALLE", k_min, k_max])
            ordered_type.extend(["I", "R", "R"])
        ordered_params.append(k_para)
        ordered_type.append("R")
        if is_nappe:
            ordered_params.append(k_parf)
            ordered_type.append("R")
        ordered_params.append(k_ordo)
        ordered_type.append("R")

        # boucle sur les fonctions, intervalles, min/max, extrema
        _type = {"min": "MINI", "max": "MAXI"}
        _PREC = 1.0e-6
        tab = Table(para=ordered_params, typ=ordered_type)
        for tf in l_fonc:
            if not with_intervalle:
                if not is_nappe:
                    interv = [[float(min(tf.vale_x)), float(max(tf.vale_x))]]
                else:
                    interv = [[-1.0e-300, 1.0e300]]
            for num_int, bornes in enumerate(interv):
                x1, x2 = bornes
                if not is_nappe:
                    stf = tf.cut(x1, x2, _PREC, nom=tf.nom)
                else:
                    stf = tf
                extrema = stf.extreme()
                for key in ("min", "max"):
                    nb = len(extrema[key])
                    for i in range(nb):
                        line = {
                            "FONCTION": tf.nom.strip(),
                            "TYPE": _type[key],
                            k_para: extrema[key][i][0],
                        }
                        if is_nappe:
                            line.update({k_parf: extrema[key][i][1], k_ordo: extrema[key][i][2]})
                        else:
                            line.update({k_ordo: extrema[key][i][1]})
                        if with_intervalle:
                            line.update({"INTERVALLE": num_int + 1, k_min: x1, k_max: x2})
                        tab.append(line)
        tab.titr = "Extrema de " + ", ".join(l_nom)
        # table résultat
        dprod = tab.dict_CREA_TABLE()
        C_out = CREA_TABLE(**dprod)

    #
    if "ECART_TYPE" in args:
        ECART_TYPE = args["ECART_TYPE"]
        __ff = ECART_TYPE["FONCTION"].convert()
        if ECART_TYPE["INST_INIT"] is not None:
            tini = ECART_TYPE["INST_INIT"]
        else:
            tini = __ff.vale_x[0]
        if ECART_TYPE["INST_FIN"] is not None:
            tfin = ECART_TYPE["INST_FIN"]
        else:
            tfin = __ff.vale_x[-1]
        __ff = __ff.cut(tini, __ff.vale_x[-1], ECART_TYPE["PRECISION"], ECART_TYPE["CRITERE"])
        __ff = __ff.cut(__ff.vale_x[0], tfin, ECART_TYPE["PRECISION"], ECART_TYPE["CRITERE"])
        if ECART_TYPE["METHODE"] == "SIMPSON":
            __ex = __ff.simpson(0.0)
        if ECART_TYPE["METHODE"] == "TRAPEZE":
            __ex = __ff.trapeze(0.0)
        fmoy = __ex.vale_y[-1] / (__ff.vale_x[-1] - __ff.vale_x[0])
        __ff = __ff + (-1 * fmoy)
        __ff = __ff * __ff
        if ECART_TYPE["METHODE"] == "SIMPSON":
            __ez = __ff.simpson(0.0)
        if ECART_TYPE["METHODE"] == "TRAPEZE":
            __ez = __ff.trapeze(0.0)
        sigma = math.sqrt(__ez.vale_y[-1] / (__ff.vale_x[-1] - __ff.vale_x[0]))
        C_out = CREA_TABLE(
            LISTE=(
                _F(LISTE_K=ECART_TYPE["FONCTION"].getName(), PARA="FONCTION"),
                _F(LISTE_K=ECART_TYPE["METHODE"], PARA="METHODE"),
                _F(LISTE_R=fmoy, PARA="MOYENNE"),
                _F(LISTE_R=sigma, PARA="ECART_TYPE"),
                _F(LISTE_R=tini, PARA="INST_INIT"),
                _F(LISTE_R=tfin, PARA="INST_FIN"),
            )
        )

    #
    if "RMS" in args:
        RMS = list(args["RMS"])
        sigm = []
        tmpi = []
        tmpf = []
        nomf = []
        meth = []
        for i_rms in RMS:
            __ff = i_rms["FONCTION"].convert()
            if i_rms["INST_INIT"] is not None:
                tini = i_rms["INST_INIT"]
            else:
                tini = __ff.vale_x[0]
            if i_rms["INST_FIN"] is not None:
                tfin = i_rms["INST_FIN"]
            else:
                tfin = __ff.vale_x[-1]
            __ff = __ff.cut(tini, __ff.vale_x[-1], i_rms["PRECISION"], i_rms["CRITERE"])
            __ff = __ff.cut(__ff.vale_x[0], tfin, i_rms["PRECISION"], i_rms["CRITERE"])
            __ff = __ff * __ff
            if i_rms["METHODE"] == "SIMPSON":
                __ez = __ff.simpson(0.0)
            if i_rms["METHODE"] == "TRAPEZE":
                __ez = __ff.trapeze(0.0)
            sigm.append(math.sqrt(__ez.vale_y[-1] / (__ff.vale_x[-1] - __ff.vale_x[0])))
            tmpi.append(tini)
            tmpf.append(tfin)
            nomf.append(i_rms["FONCTION"].getName())
            meth.append(i_rms["METHODE"])
        C_out = CREA_TABLE(
            LISTE=(
                _F(LISTE_K=nomf, PARA="FONCTION"),
                _F(LISTE_K=meth, PARA="METHODE"),
                _F(LISTE_R=tmpi, PARA="INST_INIT"),
                _F(LISTE_R=tmpf, PARA="INST_FIN"),
                _F(LISTE_R=sigm, PARA="RMS"),
            )
        )

    #
    if "NORME" in args:
        NORME = args["NORME"]
        __ff = NORME["FONCTION"].convert()
        norme = []
        for __fi in __ff.l_fonc:
            norme.append(__fi.normel2())
        nom = [NORME["FONCTION"].getName()] * len(norme)
        C_out = CREA_TABLE(
            LISTE=(_F(LISTE_R=norme, PARA="NORME"), _F(LISTE_K=nom, PARA="FONCTION"))
        )

    #
    if "NOCI_SEISME" in args:
        NOCI_SEISME = args["NOCI_SEISME"]
        l_table = []
        if NOCI_SEISME["SPEC_OSCI"] is not None:
            # cas intensité spectrale d'une nappe de SRO
            # la seule option licite est INTE_SPEC
            # intensite spectrale, il est prudent de verifier la norme de la nappe sur laquelle \
            # porte le calcul, ceci peut etre une source d erreurs.''')
            UTMESS("I", "FONCT0_40")
            amor = NOCI_SEISME["AMOR_REDUIT"]
            fini = NOCI_SEISME["FREQ_INIT"]
            ffin = NOCI_SEISME["FREQ_FIN"]
            __sp = NOCI_SEISME["SPEC_OSCI"].convert()
            vale_x = __sp.l_fonc[0].vale_x
            vale_y = [__sp(amor, f) for f in vale_x]
            para = __sp.l_fonc[0].para
            __srov = t_fonction(vale_x, vale_y, para)
            if NOCI_SEISME["NATURE"] == "DEPL":
                __srov.vale_y = (__srov.vale_y / __srov.vale_x) * 2.0 * math.pi
            elif NOCI_SEISME["NATURE"] == "VITE":
                __srov.vale_y = __srov.vale_y / __srov.vale_x / __srov.vale_x
            elif NOCI_SEISME["NATURE"] == "ACCE":
                __srov.vale_y = __srov.vale_y / __srov.vale_x / __srov.vale_x
                __srov.vale_y = __srov.vale_y / __srov.vale_x / 2.0 / math.pi
            __srov = __srov.cut(fini, ffin, NOCI_SEISME["PRECISION"], NOCI_SEISME["CRITERE"])
            insp = __srov.trapeze(0.0).vale_y[-1]
            l_table.append(_F(LISTE_R=fini, PARA="FREQ_INIT"))
            l_table.append(_F(LISTE_R=ffin, PARA="FREQ_FIN"))
            l_table.append(_F(LISTE_R=amor, PARA="AMOR_REDUIT"))
            l_table.append(_F(LISTE_R=insp, PARA="INTE_SPECT"))
        if NOCI_SEISME["FONCTION"] is not None:
            # cas fonction
            l_table.append(_F(LISTE_K=NOCI_SEISME["FONCTION"].getName(), PARA="FONCTION"))
            __ac = NOCI_SEISME["FONCTION"].convert()
            option = NOCI_SEISME["OPTION"]
            if NOCI_SEISME["INST_INIT"] is not None:
                tdeb = NOCI_SEISME["INST_INIT"]
            else:
                tdeb = __ac.vale_x[0]
            if NOCI_SEISME["INST_FIN"] is not None:
                tfin = NOCI_SEISME["INST_FIN"]
            else:
                tfin = __ac.vale_x[-1]
            # calcul de la vitesse :
            __vi = __ac.trapeze(NOCI_SEISME["COEF"])
            # calcul du déplacement :
            __de = __vi.trapeze(NOCI_SEISME["COEF"])
            # calcul de |acceleration| :
            __aa = __ac.abs()
            # calcul de integrale(|acceleration|) :
            # on "coupe" la fonction entre tdeb et tfin
            __ac = __ac.cut(tdeb, tfin, NOCI_SEISME["PRECISION"], NOCI_SEISME["CRITERE"])
            __vi = __vi.cut(tdeb, tfin, NOCI_SEISME["PRECISION"], NOCI_SEISME["CRITERE"])
            __de = __de.cut(tdeb, tfin, NOCI_SEISME["PRECISION"], NOCI_SEISME["CRITERE"])
            __aa = __aa.cut(tdeb, tfin, NOCI_SEISME["PRECISION"], NOCI_SEISME["CRITERE"])
            if NOCI_SEISME["FREQ"] is not None:
                l_freq = NOCI_SEISME["FREQ"]
            elif NOCI_SEISME["LIST_FREQ"] is not None:
                l_freq = NOCI_SEISME["LIST_FREQ"].getValues()
            else:
                # fréquences par défaut
                l_freq = []
                for i in range(56):
                    l_freq.append(0.2 + 0.050 * i)
                for i in range(8):
                    l_freq.append(3.0 + 0.075 * i)
                for i in range(14):
                    l_freq.append(3.6 + 0.100 * i)
                for i in range(24):
                    l_freq.append(5.0 + 0.125 * i)
                for i in range(28):
                    l_freq.append(8.0 + 0.250 * i)
                for i in range(6):
                    l_freq.append(15.0 + 0.500 * i)
                for i in range(4):
                    l_freq.append(18.0 + 1.000 * i)
                for i in range(10):
                    l_freq.append(22.0 + 1.500 * i)
            if option in ("TOUT", "MAXI", "ACCE_SUR_VITE"):
                #   calcul du max des valeurs absolues
                maxa_ac = __ac.abs().extreme()["max"][0][1]
                maxa_vi = __vi.abs().extreme()["max"][0][1]
                maxa_de = __de.abs().extreme()["max"][0][1]
                l_table.append(_F(LISTE_R=maxa_ac, PARA="ACCE_MAX"))
                l_table.append(_F(LISTE_R=maxa_vi, PARA="VITE_MAX"))
                l_table.append(_F(LISTE_R=maxa_de, PARA="DEPL_MAX"))
                l_table.append(_F(LISTE_R=maxa_ac / maxa_vi, PARA="ACCE_SUR_VITE"))
            if option in ("TOUT", "INTE_ARIAS"):
                __a2 = __ac * __ac
                inte_arias = __a2.trapeze(0.0).vale_y[-1]
                inte_arias = inte_arias * math.pi / NOCI_SEISME["PESANTEUR"] / 2.0
                l_table.append(_F(LISTE_R=inte_arias, PARA="INTE_ARIAS"))
            if option in ("TOUT", "VITE_ABSO_CUMU"):
                __vc = __aa.trapeze(0.0)
                vite_abso = __vc.vale_y[-1]
                l_table.append(_F(LISTE_R=vite_abso, PARA="VITE_ABSO_CUMU"))
            if option in ("TOUT", "ASA"):
                amor = NOCI_SEISME["AMOR_REDUIT"]
                freq_osci = NOCI_SEISME["FREQ_FOND"]
                ratio = NOCI_SEISME["RATIO"]
                freq_pas = NOCI_SEISME["FREQ_PAS"]
                f_ini = (1 - ratio) * freq_osci
                liste_freq = NP.arange(f_ini, freq_osci + freq_pas, freq_pas)
                __so = CALC_FONCTION(
                    SPEC_OSCI=_F(
                        NATURE="ACCE",
                        NATURE_FONC="ACCE",
                        FONCTION=NOCI_SEISME["FONCTION"],
                        METHODE="NIGAM",
                        NORME=NOCI_SEISME["NORME"],
                        FREQ=liste_freq,
                        AMOR_REDUIT=(amor,),
                    )
                )
                __srov = __so.convert().l_fonc[0]
                ASA_R = 1.0 / (ratio * freq_osci) * NP.trapz(__srov.vale_y, __srov.vale_x)
                l_table.append(_F(LISTE_R=ASA_R, PARA="ASA"))
                l_table.append(_F(LISTE_R=ratio, PARA="RATIO"))
                if option == "ASA":
                    l_table.append(_F(LISTE_R=amor, PARA="AMOR_REDUIT"))
            if option in ("TOUT", "INTE_SPEC"):
                amor = NOCI_SEISME["AMOR_REDUIT"]
                fini = NOCI_SEISME["FREQ_INIT"]
                ffin = NOCI_SEISME["FREQ_FIN"]
                __so = CALC_FONCTION(
                    SPEC_OSCI=_F(
                        NATURE="VITE",
                        NATURE_FONC="ACCE",
                        FONCTION=NOCI_SEISME["FONCTION"],
                        METHODE="NIGAM",
                        NORME=NOCI_SEISME["NORME"],
                        FREQ=l_freq,
                        AMOR_REDUIT=(amor,),
                    )
                )
                __srov = __so.convert().l_fonc[0]
                __srov = __srov.cut(fini, ffin, NOCI_SEISME["PRECISION"], NOCI_SEISME["CRITERE"])
                __srov.vale_y = __srov.vale_y / __srov.vale_x / __srov.vale_x
                insp = __srov.trapeze(0.0).vale_y[-1]
                l_table.append(_F(LISTE_R=fini, PARA="FREQ_INIT"))
                l_table.append(_F(LISTE_R=ffin, PARA="FREQ_FIN"))
                l_table.append(_F(LISTE_R=amor, PARA="AMOR_REDUIT"))
                l_table.append(_F(LISTE_R=insp, PARA="INTE_SPECT"))
            if option in ("TOUT", "DUREE_PHAS_FORT"):
                __a2 = __ac * __ac
                __i2 = __a2.trapeze(0.0)
                arias = __i2.vale_y[-1] * math.pi / NOCI_SEISME["PESANTEUR"] / 2.0
                valinf = arias * NOCI_SEISME["BORNE_INF"]
                valsup = arias * NOCI_SEISME["BORNE_SUP"]
                for i in range(len(__i2.vale_x)):
                    ariask = __i2.vale_y[i] * math.pi / NOCI_SEISME["PESANTEUR"] / 2.0
                    if ariask >= valinf:
                        break
                for j in range(len(__i2.vale_x) - 1, -1, -1):
                    ariask = __i2.vale_y[j] * math.pi / NOCI_SEISME["PESANTEUR"] / 2.0
                    if ariask <= valsup:
                        break
                dphfor = __i2.vale_x[j] - __i2.vale_x[i]
                l_table.append(_F(LISTE_R=__i2.vale_x[i], PARA="DEBUT_PHAS_FORT"))
                l_table.append(_F(LISTE_R=dphfor, PARA="DUREE_PHAS_FORT"))
        C_out = CREA_TABLE(LISTE=l_table)

    if INFO > 1:
        IMPR_TABLE(UNITE=6, TABLE=C_out)
    return C_out
