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

# person_in_charge: jerome.laverne at edf.fr
# ----------------------------------------------------------------------
#  POST_CZM_FISS :
#  ---------------
#  OPTION = 'LONGUEUR'
#    - CALCUL DE LA LONGUEUR DES FISSURES COHESIVES 2D
#    - PRODUIT UNE TABLE
#  OPTION = 'TRIAXIALITE'
#    - CALCUL DE LA TRIAXIALITE DANS LES ELEMENTS MASSIFS CONNECTES A
#      L'INTERFACE COHESIVE
#    - PRODUIT UNE CARTE
# ----------------------------------------------------------------------

from math import acos, pi

from numpy import *

import aster
from ..Messages import UTMESS, MasquerAlarme

from ..Cata.Syntax import _F
from ..CodeCommands import CALC_CHAM_ELEM, CREA_CHAMP, CREA_TABLE
from .Fracture.post_voisin_czm import POST_VOISIN_CZM


def distance(x, y, xref, yref, xdir, ydir):
    xvect = x - xref
    yvect = y - yref
    si = sign(xdir * xvect + ydir * yvect)
    di = (((xvect) ** 2 + (yvect) ** 2) ** 0.5) * si
    return di


def post_czm_fiss_ops(self, OPTION, RESULTAT, **args):
    """Corps de POST_CZM_FISS"""

    #
    # calcul de la longueur d'une fissure cohesive 2D
    #
    if OPTION == "LONGUEUR":
        # Mots cles specifiques au bloc "LONGUEUR"
        GROUP_MA = args["GROUP_MA"]
        POINT_ORIG = args["POINT_ORIG"]
        VECT_TANG = args["VECT_TANG"]

        # On importe les definitions des commandes a utiliser dans la macro

        # Recuperation du nom du modele
        __MODEL = RESULTAT.getModel()

        if __MODEL is None:
            UTMESS("F", "RUPTURE0_18")

        # Calcul des coordonnees des points de Gauss
        __CHAMEL = CALC_CHAM_ELEM(MODELE=__MODEL, GROUP_MA=GROUP_MA, OPTION="COOR_ELGA")
        __CORX = __CHAMEL.getValuesWithDescription("X", list(GROUP_MA))
        __CORY = __CHAMEL.getValuesWithDescription("Y", list(GROUP_MA))

        xg = array(__CORX[0])
        yg = array(__CORY[0])
        nbpg = len(xg)

        xmin = min(xg)
        ymin = min(yg)
        xmax = max(xg)
        ymax = max(yg)

        # A = coef dir et B=ordo orig de la droite des pg
        A = (ymax - ymin) / (xmax - xmin)
        B = ymax - A * xmax

        # Vecteur des points de gauss (qui va du point min au point max) et sa
        # norme
        vx = xmax - xmin
        vy = ymax - ymin
        nv = (vx**2 + vy**2) ** 0.5

        # vecteur directeur dans la direction ou l'on calcule la longueur
        xdir = VECT_TANG[0]
        ydir = VECT_TANG[1]
        ndir = (xdir**2 + ydir**2) ** 0.5
        # point de reference a partir duquel on calcule la longueur
        xref = POINT_ORIG[0]
        yref = POINT_ORIG[1]

        # angle entre le vecteur des points de gauss et le vecteur donne par
        # l'utilisateur
        alpha = acos((xdir * vx + ydir * vy) / (ndir * nv))

        # petit parametre de tolerence
        eps = 0.0001
        # cas ou le point de reference n'est pas aligne avec les points de
        # Gauss

        if abs(yref - A * xref - B) >= eps:
            UTMESS("F", "POST0_45", valk=list(GROUP_MA), valr=(xmin, xmax, ymin, ymax))

        # cas ou le vecteur n'est pas colineaire a la droite des points de
        # Gauss
        if (abs(alpha) >= eps) and (abs(alpha - pi) >= eps):
            UTMESS("F", "POST0_46", valk=list(GROUP_MA), valr=(xmin, xmax, ymin, ymax))

        # Calcul de la distance signee des points de Gauss au point de
        # reference
        disg = distance(xg, yg, xref, yref, xdir, ydir)
        ming = min(disg)
        maxg = max(disg)

        __INST = aster.GetResu(RESULTAT.getName(), "VARI_ACCES")["INST"]
        nbinst = len(__INST)

        Lfis = [0] * (nbinst)
        Ltot = [0] * (nbinst)
        Lcoh = [0] * (nbinst)
        __VI = [0] * (nbinst)

        for j in range(0, nbinst):
            __VI[j] = CREA_CHAMP(
                TYPE_CHAM="ELGA_VARI_R",
                OPERATION="EXTR",
                RESULTAT=RESULTAT,
                NOM_CHAM="VARI_ELGA",
                NUME_ORDRE=j,
            )

            __VI3 = __VI[j].getValuesWithDescription("V3", list(GROUP_MA))

            mat_v3 = __VI3[0]
            nbpg = len(mat_v3)

            # Evaluation du nombre de points de gauss dans chaque etat
            cpt0 = 0
            cpt1 = 0
            cpt2 = 0

            max0 = ming
            min0 = maxg
            max1 = ming
            min1 = maxg
            max2 = ming
            min2 = maxg

            for i in range(0, nbpg):
                if disg[i] >= 0.0:
                    if mat_v3[i] == 1.0:  # si c'est un pdg en zone cohesive
                        cpt1 = cpt1 + 1
                        max1 = max(max1, disg[i])
                        min1 = min(min1, disg[i])
                    else:
                        if mat_v3[i] == 2.0:  # si c'est un pdg en fissure
                            cpt2 = cpt2 + 1
                            max2 = max(max2, disg[i])
                            min2 = min(min2, disg[i])
                        else:
                            if mat_v3[i] == 0.0:  # si c'est un pdg sain
                                cpt0 = cpt0 + 1
                                max0 = max(max0, disg[i])
                                min0 = min(min0, disg[i])

            # verification qu'entre min1 et max1 on a que des mailles 1
            for i in range(0, nbpg):
                if cpt1 != 0:
                    if (disg[i] >= min1) and (disg[i] <= max1):
                        if mat_v3[i] != 1.0:
                            UTMESS("A", "POST0_48")

            # Verification qu'il y a bien des points de Gauss sur la
            # demi-droite definie par l'utilisateur
            if cpt0 + cpt1 + cpt2 == 0:
                UTMESS("F", "POST0_47", valk=list(GROUP_MA), valr=(xmin, xmax, ymin, ymax))

            # Verification de la taille de la zone cohesive
            if (cpt2 != 0) and (cpt1 <= 3):
                UTMESS("A", "POST0_49")

            # Evaluation des longueurs
            if (cpt1 == 0) and (cpt2 == 0):
                Ltot[j] = min0
                Lfis[j] = min0
            else:
                if (cpt1 != 0) and (cpt2 == 0):
                    Ltot[j] = (min0 + max1) / 2.0
                    Lfis[j] = min1
                else:
                    if (cpt1 == 0) and (cpt2 != 0):
                        Ltot[j] = (min0 + max2) / 2.0
                        Lfis[j] = (min0 + max2) / 2.0
                    else:
                        Ltot[j] = (min0 + max1) / 2.0
                        Lfis[j] = (min1 + max2) / 2.0

            Lcoh[j] = Ltot[j] - Lfis[j]

        TABLE_OUT = CREA_TABLE(
            LISTE=(
                _F(LISTE_R=__INST, PARA="INST"),
                _F(LISTE_R=Lfis, PARA="LONG_FIS"),
                _F(LISTE_R=Ltot, PARA="LONG_TOT"),
                _F(LISTE_R=Lcoh, PARA="LONG_COH"),
            )
        )

        return TABLE_OUT

    #
    # calcul de la triaxialite dans les elements massifs voisins de l'interface cohesive
    #
    elif OPTION == "TRIAXIALITE":
        CARTE_OUT = POST_VOISIN_CZM(RESULTAT=RESULTAT)

        return CARTE_OUT
