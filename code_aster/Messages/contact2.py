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

from ..Utilities import _

cata_msg = {
    3: _(
        """
  Toutes vos zones de contact sont en mode RESOLUTION='NON'.
  Le mode REAC_GEOM = 'SANS' est forcé.
"""
    ),
    4: _(
        """
  Toutes vos zones de contact sont en mode RESOLUTION='NON'.
  Le mode ALGO_RESO_GEOM = 'POINT_FIXE' est forcé.
"""
    ),
    10: _(
        """
Au moins une des mailles de contact que vous avez définies est de type < %(k1)s > , or les types attendus sont POI1,
SEG2, SEG3, TRIA3, TRIA6, TRIA7, QUAD4, QUAD8, QUAD9
Conseil :
Vérifiez votre AFFE_MODELE et le type de vos mailles dans la définition des surfaces de contact.
"""
    ),
    11: _(
        """
Au moins une des mailles de contact que vous avez définies est de dimension %(i1)i, or la dimension de votre problème est : %(i2)i.
Cette maille n'est donc pas une maille de bord. Il doit y avoir une erreur dans votre mise en données.

Conseil :
Vérifiez votre AFFE_MODELE et le type de vos mailles dans la définition des surfaces de contact.
"""
    ),
    12: _(
        """
Contact avec formulation continue.
Votre modèle contient des surfaces de contact qui s'appuient sur un mélange d'éléments axisymétriques et non axisymétriques.
Cela n'a pas de sens. Toute la modélisation doit être axisymétrique.

Conseil :
Vérifiez votre AFFE_MODELE et le type de vos mailles dans la définition des surfaces de contact.
"""
    ),
    13: _(
        """
Contact méthodes maillées.
La zone de contact numéro %(i1)i contient %(i2)i noeuds communs aux surfaces maîtres et esclaves.
Vérifiez la définition de vos surfaces de contact ou bien renseignez un des mots-clés SANS_NOEUD/SANS_GROUP_NO/SANS_GROUP_MA.
"""
    ),
    14: _(
        """
Contact méthode continue.
  -> Une zone de contact est définie sur une modélisation axisymétrique. Le Jacobien
     est nul car un noeud de la surface de contact esclave appartient à l'axe.
     La pression de contact (degré de liberté LAGS_C) risque d'être erronée.
  -> Conseil :
     Il faut changer de schéma d'intégration et utiliser 'GAUSS'.
"""
    ),
    15: _(
        """
Contact formulation discrète.
Les zones de contact numéro %(i1)i et numéro %(i2)i ont %(i3)i noeuds communs à leurs surfaces esclaves. Cela peut parfois conduire à une matrice de contact singulière.

Si le calcul venait à échouer, vérifiez la définition de vos surfaces de contact ou bien renseignez un des mots-clés SANS_NOEUD/SANS_GROUP_NO/SANS_GROUP_MA.
"""
    ),
    16: _(
        """
Les zones de contact numéro %(i1)i et numéro %(i2)i ont %(i3)i noeuds communs à leurs surfaces esclaves : c'est interdit.
Conseil :
 - changez vos surfaces de contact,
 - définissez alternativement les zones possédant un noeud commun comme maître et esclave,
 - pour la méthode LAC, il faut désactiver le lissage.
"""
    ),
    17: _(
        """
L'option DECOUPE_LAC de CREA_MAILLAGE n'a pas traité le même nombre de zones que celles définies dans DEFI_CONTACT
Conseil :
 - Assurez vous d'avoir renseigné le bon nombre de zones esclaves dans CREA_MAILLAGE/DECOUPE_LAC et DEFI_CONTACT
"""
    ),
    18: _(
        """
L'option DECOUPE_LAC de CREA_MAILLAGE n'a pas traité la zone %(k1)s définie dans DEFI_CONTACT.
Conseil :
 - Assurez vous d'avoir renseigné cette zone dans dans CREA_MAILLAGE/DECOUPE_LAC.
"""
    ),
    19: _(
        """
La méthode de gestion du contact LAC n'est pas compatible avec l'option MATR_DISTRIBUEE='OUI'
Conseil :
 - Assurez vous de ne pas utiliser MATR_DISTRIBUEE.
"""
    ),
    20: _(
        """
Au moins une maille de peau sans volume a été détecté. L'orientation ne peut pas être vérifiée.
Si l'orientation est incorrecte vous pouvez avoir des problèmes de détection de contact avec la méthode de gestion du contact LAC.

Pour désactiver cette alarme, utilisez VERI_NORM='NON' dans DEFI_CONTACT.
"""
    ),
    21: _(
        """La méthode de gestion du contact LAC n'est pas compatible avec le critère de convergence RESI_REFE_RELA."""
    ),
}
