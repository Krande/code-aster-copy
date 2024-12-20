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

import pickle

import numpy as NP
from numpy import linalg as LA

from ...Cata.Commons import *
from ...Cata.DataStructure import *
from ...Cata.Syntax import *
from ...Cata.Syntax import _F
from ...CodeCommands import CREA_TABLE
from ...Supervis.ExecuteCommand import UserMacro

from ...Objects import NonLinearResult


class TANGENT:

    """Vérification sur les matrices tangentes

    Attributs publics :
        mat       : matrice tangente
        ddl       : nom des degres de liberte
        nddl      : nombre de ddl
        norme     : norme de la matrice tangente
        prec_zero : en-dessous de prec_zero, on ne compare pas les matrices
    """

    def __init__(self, ddl="", prec_zero=1.0e-12):
        """
        ddl       : chaine de caracteres designant les ddl (ex: 'UUP')
        prec_zero : en-dessous de prec_zero, on ne compare pas les matrices
        """
        self.ddl = ddl
        self.prec_zero = prec_zero

    def Load(self, nom_fichier):
        """lit la matrice depuis un fichier"""
        with open(nom_fichier, "rb") as pick:
            self.__dict__ = pickle.load(pick)

    def Save(self, nom_fichier):
        """sauvegarde la matrice dans un fichier"""
        with open(nom_fichier, "wb") as pick:
            pickle.dump(self.__dict__, pick)

    def Aster(self, suffix):
        """lit la matrice depuis l'espace Aster.
        nom : suffixe de l'objet jeveux
        """
        values = NonLinearResult.getTangentMatrix(suffix)
        if not values:
            raise RuntimeError("TANGENT : OBJET JEVEUX DE SUFFIXE " + suffix + " INEXISTANT")
        self.Matrice(values)

    def Eigen(self):
        """Retourne les valeurs propres de la matrice"""
        self.vp = NP.sort(LA.eigvals(self.mat))

    def Matrice(self, matrice):
        """range la matrice.
        matrice   : la matrice tangente (rangement carre)
        """
        if type(matrice) == type((1,)):
            matrice = NP.array(list(matrice))
        elif type(matrice) == type([]):
            matrice = NP.array(matrice)
        matrice = matrice.astype(float)
        nddl = int(len(matrice) ** 0.5 + 0.5)
        matrice.shape = (nddl, nddl)
        self.mat = matrice
        self.nddl = nddl
        if not self.ddl:
            self.ddl = "D" * nddl
        elif len(self.ddl) != nddl:
            raise RuntimeError("Nommage des DDL incoherents avec la taille de la matrice")
        self.norme = NP.trace(NP.dot(NP.transpose(self.mat), self.mat))

    def Difference(self, matp, affi_ok=0, prec_diff=1.0e-4):
        """Comparaison relative de la matrice tangente avec une autre matrice
        matp      : matrice avec laquelle self.mat est comparee
        affi_ok   : si oui, on affiche egalement les valeurs qui collent bien
        prec_diff : ecart au-dessus duquel on considere que ce n'est pas OK
        """
        if type(matp) is tuple:
            matp = NP.array(list(matp))
        elif type(matp) is list:
            matp = NP.array(matp)
        elif type(matp) == type(self):
            matp = matp.mat
        elif type(matp) is NP.ndarray:
            pass
        else:
            raise RuntimeError(
                "1er argument doit etre une matrice (tuple,liste,TANGENT " "ou tableau numpy)"
            )
        matp = NP.ravel(matp)
        matp = matp.astype(float)
        if len(matp) != self.nddl * self.nddl:
            raise RuntimeError("Matrices de tailles differentes")
        matp.shape = (self.nddl, self.nddl)
        refe = NP.abs(self.mat) + NP.abs(matp)
        diff = NP.where(refe > self.prec_zero, NP.abs(self.mat - matp) / (refe + self.prec_zero), 0)
        nook = (diff.ravel() > prec_diff).nonzero()[0]
        ok = (diff.ravel() <= prec_diff).nonzero()[0]
        nook = nook.astype(NP.int32)
        ok = ok.astype(NP.int32)
        if affi_ok:
            affi = [ok, nook]
        else:
            affi = [nook]
        liste_i = []
        liste_j = []
        liste_matt = []
        liste_matp = []
        liste_diff = []
        for ind in affi:
            for pos in ind:
                i = int(pos // self.nddl)
                j = int(pos % self.nddl)
                liste_i.append(i + 1)
                liste_j.append(j + 1)
                liste_matt.append(self.mat[i, j])
                liste_matp.append(matp[i, j])
                if NP.abs(self.mat[i, j]) + NP.abs(matp[i, j]) >= self.prec_zero:
                    liste_diff.append(
                        NP.abs(self.mat[i, j] - matp[i, j])
                        / (NP.abs(self.mat[i, j]) + NP.abs(matp[i, j]) + self.prec_zero)
                    )
                else:
                    liste_diff.append(0.0)
        if self.norme > self.prec_zero:
            ecart = (self.mat - matp) / 2.0
            nor_ecart = NP.trace(NP.dot(NP.transpose(ecart), ecart))
            nor_diff = nor_ecart / self.norme
        else:
            nor_diff = 0.0
        max_diff = 0.0
        if len(liste_diff) > 0:
            max_diff = NP.max(liste_diff)
        return liste_i, liste_j, liste_matt, liste_matp, liste_diff, nor_diff, max_diff

    def Symetrie(self, prec_diff=1.0e-4):
        """Vérification que la matrice tangente est symétrique
        On retourne la norme relative de l'ecart a la symetrie : || (A-At)/2|| / ||A||
        On affiche les termes qui s'ecartent de la symetrie
        prec_diff : ecart au-dessus duquel on considere que ce n'est pas OK
        """
        tran = NP.transpose(self.mat)
        return self.Difference(tran, affi_ok=0, prec_diff=prec_diff)


def veri_matr_tang_ops(self, **args):
    """
    Ecriture de la macro verif_matrice_tangente_ops
    """

    # On importe les definitions des commandes a utiliser dans la macro

    # Le concept sortant (de type fonction) est nomme ROTGD dans
    # le contexte de la macro

    prec_zero = args["PREC_ZERO"]
    tgt = TANGENT(prec_zero=prec_zero)
    tgt.Aster("MATA")
    matp = TANGENT(prec_zero=prec_zero)
    matp.Aster("MATC")
    prec_diff = args["PRECISION"]
    if args["SYMETRIE"] == "OUI":
        symetgt = tgt.Symetrie(prec_diff)[-2]
        symeper = matp.Symetrie(prec_diff)[-2]
        print("Symetrie de la matrice tangente", symetgt)
        print("Symetrie de la matrice pr pertubation", symeper)
    if args["DIFFERENCE"] == "OUI":
        liste_i, liste_j, liste_matt, liste_matp, liste_diff, nor_diff, max_diff = tgt.Difference(
            matp, prec_diff
        )
        print(
            "difference entre matrice tangente et matrice par pertubation : norme=",
            nor_diff,
            " max=",
            max_diff,
        )
        TAB_MAT = CREA_TABLE(
            LISTE=(
                _F(PARA="I", LISTE_I=liste_i),
                _F(PARA="J", LISTE_I=liste_j),
                _F(PARA="MAT_TGTE", LISTE_R=liste_matt),
                _F(PARA="MAT_PERT", LISTE_R=liste_matp),
                _F(PARA="MAT_DIFF", LISTE_R=liste_diff),
            )
        )
    return TAB_MAT


VERI_MATR_TANG_cata = MACRO(
    nom="VERI_MATR_TANG",
    op=None,
    sd_prod=table_sdaster,
    docu="",
    reentrant="n",
    fr="verification de la matrice tangente : symetrie et difference "
    "par rapport a la matrice calculee par perturbation",
    regles=(AU_MOINS_UN("SYMETRIE", "DIFFERENCE")),
    SYMETRIE=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
    DIFFERENCE=SIMP(statut="f", typ="TXM", defaut="OUI", into=("OUI", "NON")),
    PRECISION=SIMP(statut="f", typ="R", defaut=1.0e-4),
    PREC_ZERO=SIMP(statut="f", typ="R", defaut=1.0e-12),
)


VERI_MATR_TANG = UserMacro("VERI_MATR_TANG", VERI_MATR_TANG_cata, veri_matr_tang_ops)
