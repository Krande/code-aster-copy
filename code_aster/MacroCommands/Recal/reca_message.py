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

# person_in_charge: mathieu.courtois@edf.fr

import os

import numpy as NP

from ...Messages import UTMESS, MessageLog

from .recal import Affiche_Param

# =========================================================================


# AFFICHAGE DES MESSAGES


class Message:
    """
    classe gérant l'affichage des messages concernant le déroulement de l'optmisation
    """

    # ------------------------------------------------------------------------
    def __init__(self, para, val_init, resu_exp, ul_out):
        self.nom_para = para
        self.resu_exp = resu_exp
        self.val_init = val_init
        self.resu_exp = resu_exp
        self.ul_out = ul_out

    # ------------------------------------------------------------------------
    def get_filename(self):
        return os.getcwd() + "/fort." + str(self.ul_out)

    # ------------------------------------------------------------------------
    def initialise(self):
        """Initialisation du fichier"""
        UTMESS("I", "RECAL0_1", files=self.get_filename())

    # ------------------------------------------------------------------------
    def affiche_valeurs(self, val):
        """Affichage de la valeur des parametres"""
        txt = Affiche_Param(self.nom_para, val)
        UTMESS("I", "RECAL0_32", valk=txt, files=self.get_filename())

    # ------------------------------------------------------------------------
    def affiche_fonctionnelle(self, J):
        """Affichage de la fonctionnelle"""
        UTMESS("I", "RECAL0_33", valr=J, files=self.get_filename())

    # ------------------------------------------------------------------------
    def affiche_result_iter(self, iter, J, val, residu, Act=[], ecart_para=None, ecart_fonc=None):
        """Affichage du message recapitulatif de l'iteration"""
        UTMESS("I", "RECAL0_30")
        UTMESS("I", "RECAL0_79", files=self.get_filename())
        UTMESS("I", "RECAL0_31", vali=iter, files=self.get_filename())
        self.affiche_fonctionnelle(J)
        UTMESS("I", "RECAL0_34", valr=residu, files=self.get_filename())
        if ecart_para:
            UTMESS("I", "RECAL0_37", valr=ecart_para, files=self.get_filename())
        if ecart_fonc:
            UTMESS("I", "RECAL0_38", valr=ecart_fonc, files=self.get_filename())

        # Affichage des parametres
        self.affiche_valeurs(val)

        # Si les parametres sont en butee
        if len(Act) != 0:
            lpara = " ".join([self.nom_para[i] for i in Act])
            if len(Act) == 1:
                UTMESS("I", "RECAL0_46", valk=lpara, files=self.get_filename())
            else:
                UTMESS("I", "RECAL0_47", valk=lpara, files=self.get_filename())

        UTMESS("I", "RECAL0_80", files=self.get_filename())

    # ------------------------------------------------------------------------
    def affiche_etat_final_convergence(
        self, iter, max_iter, iter_fonc, max_iter_fonc, prec, residu, Act=[]
    ):
        """Affichage du message recapitulatif a la fin du processus d'optimisation"""
        if (iter < max_iter) and (residu <= prec) and (iter_fonc < max_iter_fonc):
            UTMESS("I", "RECAL0_56", files=self.get_filename())
            if len(Act) != 0:
                UTMESS("I", "RECAL0_58", files=self.get_filename())
        else:
            UTMESS("I", "RECAL0_57", files=self.get_filename())
            if iter >= max_iter:
                UTMESS("I", "RECAL0_55", files=self.get_filename())
            if iter_fonc >= max_iter_fonc:
                UTMESS("I", "RECAL0_54", files=self.get_filename())

        UTMESS("I", "RECAL0_80", files=self.get_filename())

    # ------------------------------------------------------------------------
    def affiche_calcul_etat_final(
        self, para, Hessien, valeurs_propres, vecteurs_propres, sensible, insensible
    ):
        """Affichage des informations de l'optimisation (valeurs propres, vecteurs propres, etc.)"""
        UTMESS("I", "RECAL0_60", valk=str(valeurs_propres), files=self.get_filename())
        UTMESS("I", "RECAL0_61", valk=str(vecteurs_propres), files=self.get_filename())
        UTMESS("I", "RECAL0_62", files=self.get_filename())

        if len(sensible) != 0 or len(insensible) != 0:
            UTMESS("I", "RECAL0_63", files=self.get_filename())

        # Parametres sensibles
        if len(sensible) != 0:
            UTMESS("I", "RECAL0_64", files=self.get_filename())
            k = 0
            for i in sensible:
                k = k + 1
                colonne = vecteurs_propres[:, i]
                numero = NP.nonzero(NP.greater(abs(colonne / max(abs(colonne))), 1.0e-1))[0]
                txt = "\n   " + str(k) + ") "
                for j in numero:
                    txt += "%+3.1E " % colonne[j] + "* " + para[j] + " "
                dict_args = dict(valk=(txt, str(valeurs_propres[i])), files=self.get_filename())
                UTMESS("I", "RECAL0_65", **dict_args)

        # Parametres insensibles
        if len(insensible) != 0:
            UTMESS("I", "RECAL0_66", files=self.get_filename())
            k = 0
            for i in insensible:
                k = k + 1
                colonne = vecteurs_propres[:, i]
                numero = NP.nonzero(NP.greater(abs(colonne / max(abs(colonne))), 1.0e-1))[0]
                txt = "\n   " + str(k) + ") "
                for j in numero:
                    txt += "%+3.1E " % colonne[j] + "* " + para[j] + " "
                dict_args = dict(valk=(txt, str(valeurs_propres[i])), files=self.get_filename())
                UTMESS("I", "RECAL0_65", **dict_args)

        if len(sensible) != 0 or len(insensible) != 0:
            UTMESS("I", "RECAL0_62", files=self.get_filename())
