# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

from ..Utilities import _

cata_msg = {
    1: {
        "message": _(
            """
 L'option %(k1)s est déjà calculée pour le numéro d'ordre %(k2)s.
 On la recalcule car les données peuvent être différentes.
"""
        ),
        "flags": "DECORATED",
    },
    2: _(
        """
Champ inexistant SIEF_ELGA ou SIEF_ELGA numéro d'ordre %(k1)s pour le calcul de l'option %(k2)s
"""
    ),
    3: _(
        """
Champ inexistant DEPL numéro d'ordre %(k1)s pour le calcul de l'option %(k2)s
"""
    ),
    8: _(
        """

 la taille mémoire   nécessaire au vecteur de travail dans   lequel nous stockons les composantes   u et v du vecteur TAU est trop importante   par rapport a la place disponible.
 taille disponible :  %(i1)d
 taille nécessaire :  %(i2)d
"""
    ),
    10: _(
        """
 le noeud traité n'est associé à aucune maille volumique.
 numéro du noeud =  %(i1)d
 nombre de mailles attachées au noeud =  %(i2)d
"""
    ),
    16: _(
        """
 appel erroné numéro d'ordre %(i1)d code retour  : %(i2)d
 Problème CHAM_NO %(k1)s
Séchage moins %(r2)f  Séchage plus %(r3)f
"""
    ),
    19: _(
        """
 nombre de noeud(s) éliminé(s) du maillage  %(i1)d
"""
    ),
    20: _(
        """
 nombre de maille(s) éliminée(s) du maillage  %(i1)d
"""
    ),
    41: _(
        """
  le paramètre  %(k1)s n'existe pas %(k2)s
"""
    ),
    45: _(
        """
 noeud inconnu dans le fichier  IDEAS  noeud numéro :  %(i1)d
"""
    ),
    46: _(
        """
 élément inconnu dans le fichier IDEAS élément numéro :  %(i1)d
"""
    ),
    48: _(
        """
Erreur utilisateur dans la commande EXTR_RESU / RESTREINT :
  La structure de donnée contient des champs par éléments.
  Le mot clé MODELE est obligatoire.
"""
    ),
    49: _(
        """
Erreur utilisateur dans la commande EXTR_RESU / RESTREINT :
  Le maillage associé à certains champs de %(k1)s est %(k2)s.
  Il devrait être %(k3)s.
"""
    ),
    50: _(
        """
 Attention : la valeur d'amortissement associée au mode propre %(i1)d est négative ou nulle : %(r1)f.
"""
    ),
    51: _(
        """
 La valeur d'amortissement associée au mode propre %(i1)d est négative ou nulle : %(r1)f.
 Vous avez demandé qu'elle soit corrigée. Cet amortissement est mis à %(r2)f.
"""
    ),
    52: _(
        """
 La valeur d'amortissement associée au mode propre %(i1)d est négative ou nulle : %(r1)f.
"""
    ),
    57: _(
        """
 problème dans  le traitement de l'instant  %(r1)f
  récupération de  %(k1)s
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    58: _(
        """
 problème dans  le traitement de l'instant  %(r1)f
  récupération  pour  %(k1)s
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    59: _(
        """
 problème dans  le traitement de l'instant  %(r1)f
  récupération pour le secteur  %(i1)d
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    64: _(
        """
 la valeur d'amortissement réduit est trop grande
 la valeur d'amortissement :  %(r1)f
  du mode propre  %(i1)d
  est tronquée au seuil :  %(r2)f
"""
    ),
    67: _(
        """

 la taille mémoire   nécessaire au vecteur de travail   est trop importante   par rapport a la place disponible.
 taille disponible :  %(i1)d
 taille nécessaire :  %(i2)d
"""
    ),
    68: _(
        """

 la taille du vecteur  contenant les caractéristiques des   paquets de mailles est trop petite.
 nombre de paquets max :  %(i1)d
 nombre de paquets réels:  %(i2)d
"""
    ),
    70: _(
        """

 la taille du vecteur  contenant les caractéristiques des   paquets de noeuds est trop petite.
 nombre de paquets max :  %(i1)d
 nombre de paquets réels:  %(i2)d
"""
    ),
    73: _(
        """
 appel erroné  résultat :  %(k1)s   archivage numéro :  %(i1)d
   code retour :  %(i2)d
   problème champ :  %(k2)s
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    74: _(
        """
 on ne trouve pas l'instant  %(r1)f  dans la table  %(k1)s
"""
    ),
    75: _(
        """
 on trouve   plusieurs instants  %(r1)f  dans la table  %(k1)s
"""
    ),
    76: _(
        """
 noeud non contenu dans une  maille sachant calculer l" option
 noeud numéro :  %(i1)d
"""
    ),
    77: _(
        """
banque de données pour le type de géométrie  %(k1)s
  le couple de matériaux  %(k2)s
  ne se trouve pas dans la banque. %(k3)s
"""
    ),
    78: _(
        """
 le calcul du rayon n'est pas assez précis.
    %(r1)f
    %(r2)f
    %(r3)f
    %(r4)f
    %(r5)f
    %(r6)f
    %(i1)d
    %(r7)f
    %(r8)f
    %(r9)f
    %(r10)f
    %(r11)f
    %(r12)f
    %(r13)f
    %(r14)f
    (%(r15)f
    %(r16)f

Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
}
