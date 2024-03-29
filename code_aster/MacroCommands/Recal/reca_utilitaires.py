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

import copy
import glob
import os
import sys

import numpy as NP

from ...Messages import UTMESS

from ...Cata.Syntax import _F
from ...CodeCommands import DEFI_FICHIER, IMPR_FONCTION, INFO_EXEC_ASTER

try:
    import Gnuplot
except ImportError:
    pass


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------

# _____________________________________________
#
# DIVERS UTILITAIRES POUR LA MACRO
# _____________________________________________


def transforme_list_Num(parametres, res_exp):
    """
    Transforme les données entrées par l'utilisateur en tableau numpy
    """

    dim_para = len(parametres)  # donne le nb de parametres
    val_para = NP.zeros(dim_para)
    borne_inf = NP.zeros(dim_para)
    borne_sup = NP.zeros(dim_para)
    para = []
    for i in range(dim_para):
        para.append(parametres[i][0])
        val_para[i] = parametres[i][1]
        borne_inf[i] = parametres[i][2]
        borne_sup[i] = parametres[i][3]
    return para, val_para, borne_inf, borne_sup


# ------------------------------------------------------------------------
def Random_Tmp_Name(prefix=None):
    crit = False
    while crit == False:
        nombre = int(random.random() * 10000000)
        if prefix:
            fic = prefix + str(nombre)
        else:
            if "TEMP" in os.environ:
                fic = os.path.join(os.environ["TEMP"], "file%s" % str(nombre))
            else:
                fic = "/tmp/file" + str(nombre)
        if not os.path.isfile(fic):
            crit = True
    return fic


# _____________________________________________
#
# CALCUL DU TEMPS CPU RESTANT
# _____________________________________________
def temps_CPU(restant_old, temps_iter_old):
    """
    Fonction controlant le temps restant
    """
    __cpu = INFO_EXEC_ASTER(LISTE_INFO=("TEMPS_RESTANT",))
    TEMPS = __cpu["TEMPS_RESTANT", 1]
    err = 0
    # Indique une execution interactive
    if TEMPS > 1.0e9:
        return 0.0, 0.0, 0
    # Indique une execution en batch
    else:
        restant = TEMPS
        # Initialisation
        if restant_old == 0.0:
            temps_iter = -1.0
        else:
            # Première mesure
            if temps_iter_old == -1.0:
                temps_iter = restant_old - restant
            # Mesure courante
            else:
                temps_iter = (temps_iter_old + (restant_old - restant)) / 2.0
            if (temps_iter > 0.96 * restant) or (restant < 0.0):
                err = 1
                UTMESS("F", "RECAL0_53")

    return restant, temps_iter, err


# _____________________________________________
#
# IMPRESSIONS GRAPHIQUES
# _____________________________________________
def graphique(FORMAT, L_F, res_exp, reponses, iter, UL_out, pilote, fichier=None, INFO=0):
    if iter:
        txt_iter = "Iteration : " + str(iter)
    else:
        txt_iter = ""

    # Le try/except est la pour eviter de planter betement dans un trace de
    # courbes (DISPLAY non defini, etc...)
    try:
        if FORMAT == "XMGRACE":
            for i in range(len(L_F)):
                _tmp = []
                courbe1 = res_exp[i]
                _tmp.append(
                    {
                        "ABSCISSE": courbe1[:, 0].tolist(),
                        "ORDONNEE": courbe1[:, 1].tolist(),
                        "COULEUR": 1,
                        "LEGENDE": "Expérience",
                    }
                )
                courbe2 = L_F[i]
                _tmp.append(
                    {
                        "ABSCISSE": courbe2[:, 0].tolist(),
                        "ORDONNEE": courbe2[:, 1].tolist(),
                        "COULEUR": 2,
                        "LEGENDE": "Calcul",
                    }
                )

                motscle2 = {"COURBE": _tmp}
                motscle2["PILOTE"] = pilote

                IMPR_FONCTION(
                    FORMAT="XMGRACE",
                    UNITE=int(UL_out),
                    TITRE="Courbe : " + reponses[i][0],
                    SOUS_TITRE=txt_iter,
                    LEGENDE_X=reponses[i][1],
                    LEGENDE_Y=reponses[i][2],
                    **motscle2
                )
                dic = {
                    "": "",
                    "POSTSCRIPT": ".ps",
                    "EPS": ".eps",
                    "MIF": ".mif",
                    "SVG": ".svg",
                    "PNM": ".pnm",
                    "PNG": ".png",
                    "JPEG": ".jpg",
                    "PDF": ".pdf",
                    "INTERACTIF": ".agr",
                }
                ext = dic[pilote]
                if ext != "":
                    os.system(
                        "mv ./fort.%s ./REPE_OUT/courbes_%s_iter_%s%s"
                        % (str(UL_out), reponses[i][0], str(iter), ext)
                    )

        elif FORMAT == "GNUPLOT":
            if fichier:
                if INFO >= 2:
                    UTMESS("I", "RECAL0_41", valk=fichier)
                # On efface les anciens graphes
                liste = glob.glob(fichier + "*.ps")
                for fic in liste:
                    try:
                        os.remove(fic)
                    except:
                        pass

            graphe = []
            impr = Gnuplot.Gnuplot()
            Gnuplot.GnuplotOpts.prefer_inline_data = 1
            impr("set data style linespoints")
            impr("set grid")
            impr("set pointsize 1.")
            impr("set terminal postscript color")
            impr('set output "fort.' + str(UL_out) + '"')

            for i in range(len(L_F)):
                graphe.append(Gnuplot.Gnuplot(persist=0))
                graphe[i]("set data style linespoints")
                graphe[i]("set grid")
                graphe[i]("set pointsize 1.")
                graphe[i].xlabel(reponses[i][1])
                graphe[i].ylabel(reponses[i][2])
                graphe[i].title(reponses[i][0] + "  " + txt_iter)
                graphe[i].plot(
                    Gnuplot.Data(L_F[i], title="Calcul"),
                    Gnuplot.Data(res_exp[i], title="Experimental"),
                )
                if pilote == "INTERACTIF":
                    graphe[i]("pause 5")
                else:
                    if fichier:
                        if INFO >= 2:
                            UTMESS("I", "RECAL0_41", valk=fichier + "_" + str(i) + ".ps")
                        graphe[i].hardcopy(fichier + "_" + str(i) + ".ps", enhanced=1, color=1)

                impr.xlabel(reponses[i][1])
                impr.ylabel(reponses[i][2])
                impr.title(reponses[i][0] + "  Iteration " + str(iter))
                impr.plot(
                    Gnuplot.Data(L_F[i], title="Calcul"),
                    Gnuplot.Data(res_exp[i], title="Experimental"),
                )

    except Exception as err:
        UTMESS("A", "RECAL0_42", valk=str(err))
