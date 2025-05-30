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

from ..Utilities import _

cata_msg = {
    1: _(
        """   Liste des comportements
"""
    ),
    2: _(
        """
Comportement - Occurrence  %(i1)d"""
    ),
    4: _("""   Affecté sur %(i1)d éléments"""),
    5: _("""     Relation                             : %(k1)s"""),
    6: _("""     Déformation                          : %(k1)s"""),
    7: _("""     Régularisation visqueuse             : %(k1)s"""),
    8: _("""     Algorithme contraintes planes (ou 1D): Deborst"""),
    9: _("""     Nombre total de variables internes   : %(i1)d"""),
    10: _("""     Relation                             : prototype (MFRONT ou UMAT)"""),
    11: _("""     Relation                             : défini par DEFI_COMPOR"""),
    12: _(
        """Une relation de type MULTIFIBRE (définie par DEFI_COMPOR) a été affectée sur %(i1)d éléments."""
    ),
    13: _(
        """
  Il y a deux types de modélisations différents pour une même affectation du comportement MFRONT.
  Chaque occurrence de COMPORTEMENT avec relation MFRONT doit contenir un seul type de modélisation.
  Par exemple, une occurrence pour les groupes de mailles contenant des éléments 3D (et comportement MFRONT),
  puis une occurrence pour les groupes de mailles contenant des éléments discrets (et comportement ELAS).

  Les types de modélisation rencontrés sont : %(i1)d et %(i2)d.
  Voir la signification de l'index dans '%(k1)s'.
"""
    ),
    14: _(
        """
  La modélisation %(i1)d dans le modèle n'est pas utilisable avec MFront.
  S'il s'agit d'éléments de type TUYAU ou PMF et que vous êtes en mode prototypage, vous devez activer ALGO_CPLAN='DEBORST' pour pouvoir l'utiliser.

  Voir la signification de l'index dans '%(k1)s'.
"""
    ),
    15: _(
        """
  Pour les comportements multiples (définis par DEFI_COMPOR), on n'a pas d'informations sur le nom des variables internes.
"""
    ),
    16: _(
        """
  Pour les comportements externes prototypes (MFRONT ou UMAT), on n'a pas d'informations sur le nom des variables internes.
"""
    ),
    18: _("""     Pas de régularisation visqueuse"""),
    19: _("""     Pas de post-traitement à la fin du pas de temps"""),
    20: _("""            V%(i1)d : %(k1)s"""),
    21: _(
        """
Le comportement %(k1)s contient %(i1)d variables d'états externes (variables de commande):
"""
    ),
    22: _("""            Variable externe %(i1)d : %(k1)s"""),
    23: _(
        """Le comportement utilise la variable d'état externe %(k1)s (variables de commande). Or elle n'a pas été définie dans AFFE_MATERIAU.
"""
    ),
    24: _(
        """
Le comportement MFront utilise les phases du Zircaloy.
Ce n'est actuellement pas possible.
"""
    ),
    25: _(
        """Le comportement utilise la variable d'état externe %(k1)s. Ce n'est actuellement pas possible avec un solveur externe (MFRONT ou UMAT)."""
    ),
    26: _(
        """Le comportement utilise la variable d'état externe %(k1)s, or la température n'est pas définie dans AFFE_VARC."""
    ),
    27: _("""     Post-traitement fin du pas de temps  : %(k1)s"""),
    32: _(
        """Le comportement %(k1)s n'est pas possible pour les éléments de type barre, grille ou autres éléments unidimensionnels."""
    ),
    35: _("""Le modèle de déformation utilisé dans MFront est incompatible avec %(k1)s."""),
    53: _(
        """
Comportement POLYCRISTAL
      Nombre de grains  %(i1)d : localisation %(k1)s
      Nombre d'occurrences de MONOCRISTAL différentes : %(i2)d - nombre de variables internes : %(i3)d
      Noms des variables internes: """
    ),
    54: _(""" A partir de la variable interne %(i1)d : pour chaque grain : """),
    55: _(""" Dernière variable interne V%(i1)d : %(k1)s"""),
    56: _(""" ... jusqu'à V%(i1)d """),
    62: _(
        """
  -> Le critère de convergence pour intégrer le comportement 'RESI_INTE'
     est lâche (très supérieur à la valeur par défaut).
  -> Risque & Conseil :
     Cela peut nuire à la qualité de la solution et à la convergence.
"""
    ),
    63: _(
        """
La définition explicite du comportement est obligatoire.
"""
    ),
    64: _(
        """
Comme vous n'avez pas défini explicitement le comportement, tout le modèle est supposé élastique en petites perturbations.
"""
    ),
    65: _(
        """
Il y a trop d'occurrences du mot-clef facteur COMPORTEMENT. On n'affichera aucune information sur les comportements."""
    ),
    71: _(
        """
La valeur propre numéro %(i1)d du module tangent local est négative et vaut %(r1)f.
L'énergie libre n'est donc pas convexe ce qui peut mener à des problèmes de convergence.
"""
    ),
    72: _(
        """
L'occurrence %(i1)d du mot-clef COMPORTEMENT n'affecte aucune maille du modèle.
Par défaut, on affecte le comportement élastique en petites déformations sur les mailles du modèle non affectées par l'utilisateur
Conseils: vérifier que ce comportement est voulu (pas d'oubli dans AFFE_MODELE).
"""
    ),
    73: _(
        """
Utilisation de la loi CZM_TURON : il est recommandé de commencer le calcul par un pas de chargement purement élastique.
"""
    ),
    74: _(
        """
RUPT_TURON : Au moins l'un des paramètres matériau (%(k1)s, %(k2)s, %(k3)s, %(k4)s, %(k5)s) vaut 0.
Ils doivent être non-nuls.
"""
    ),
    75: _(
        """
La loi de comportement MFront que vous utilisez nécessite la variable externe %(k2)s qui correspond à la variable de commande %(k1)s.
Veuillez renseigner les valeurs de cette variable de commande dans le mot-clé AFFE_VARC de AFFE_MATERIAU.
"""
    ),
    76: _(
        """
La loi de comportement réclame une température mais elle n'est pas définie par AFFE_VARC de AFFE_MATERIAU.
"""
    ),
    77: _(
        """
La loi de comportement demande à calculer une déformation non-mécanique externe mais le
nombre de composantes fournies par AFFE_VARC n'est pas cohérent avec le modèle.
"""
    ),
    79: _(
        """
La loi de comportement demande à calculer une déformation non-mécanique externe sous forme volumique
mais le problème n'est pas isotrope.
Cette erreur survient quant on tente d'utiliser de la thermique non-isotrope avec MFront.
"""
    ),
}
