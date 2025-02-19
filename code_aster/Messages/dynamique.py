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
    1: _(
        """
 schéma inconnu
"""
    ),
    2: _(
        """
 la liste d'instants fournie ne respecte pas la condition de stabilité.
"""
    ),
    3: _(
        """
 la condition de stabilité n'a pas pu être calculée pour tous les éléments. elle peut être trop grande.
"""
    ),
    4: _(
        """
  -> La condition de stabilité n'a pu être calculée pour aucun élément.
  -> Risque & Conseil :
     Vous prenez le risque de sortir du cadre de la stabilité conditionnelle du schéma de temps explicite. Vérifiez bien
     que vos éléments finis ont une taille et un matériau (module de Young) compatibles avec le respect de la condition
     de Courant vis-à-vis du pas de temps que vous avez imposé (temps de propagation des ondes dans la maille, voir
     documentation). Si c'est le cas, lever l'arrêt fatal en utilisant l'option "STOP_CFL", à vos risques et périls
     (risques de résultats faux).
"""
    ),
    5: _(
        """
 Pas de temps maximal (condition CFL) pour le schéma des différences centrées : %(r1)g s, sur la maille : %(k1)s
"""
    ),
    6: _(
        """
  Pas de temps maximal (condition CFL) pour le schéma de Tchamwa-Wilgosz : %(r1)g s, sur la maille : %(k1)s
"""
    ),
    7: _(
        """
 Pas de temps maximal (condition CFL) pour le schéma des différences centrées : %(r1)g s
"""
    ),
    8: _(
        """
  Pas de temps maximal (condition CFL) pour le schéma de Tchamwa-Wilgosz : %(r1)g s
"""
    ),
    9: _(
        """
  On ne peut pas avoir plus d'une charge de type FORCE_SOL.
"""
    ),
    10: _(
        """
   Arrêt par manque de temps CPU au groupe de pas de temps : %(i1)d
                                 au "petit" pas de temps   : %(i2)d
      - Temps moyen par "petit" pas : %(r1)f
      - Temps restant               : %(r2)f

   La base globale est sauvegardée. Elle contient les pas archivés avant l'arrêt.
"""
    ),
    11: _(
        """
   Arrêt par manque de temps CPU après le calcul de %(i1)d pas.
      - Dernier instant archivé : %(r1)f
      - Numéro d'ordre correspondant : %(i2)d
      - Temps moyen pour les %(i3)d pas de temps : %(r2)f
      - Temps restant                : %(r3)f

   La base globale est sauvegardée. Elle contient les pas archivés avant l'arrêt.
"""
    ),
    12: _(
        """
 Dans l'intervalle : %(i2)d
 Le pas de temps est trop grand : %(r1)f
 le pas de temps maximal est    : %(r2)f

 Avec le pas de temps maximal, le nombre de pas de calcul est %(i1)d
"""
    ),
    13: _(
        """
   Arrêt par manque de temps CPU à la fréquence : %(i1)d
      - Temps moyen par pas fréquence : %(r1)f
      - Temps restant                 : %(r2)f

   La base globale est sauvegardée. Elle contient les pas archivés avant l'arrêt.
"""
    ),
    14: _(
        """
   La matrice est presque singulière à la fréquence : %(r1)f
   Cette fréquence est probablement une fréquence propre du système.
"""
    ),
    15: _(
        """
 Pas de temps maximal (mot-clé PAS_MAXI) demandé : %(r1)f plus petit que
 le pas de temps initial demandé par l'utilisateur (mot-clé PAS) : %(r2)f
 Il faut s'assurer que PAS est bien inférieur ou égal à PAS_MAXI
"""
    ),
    16: _(
        """
 Pas de temps maximal calculé pour le schéma ADAPT : %(r1)f

 Risque & Conseil : la méthode de calcul automatique de ce pas maximal semble être prise en défaut.
 On recommande donc de définir explicitement cette valeur avec le mot-clé PAS_MAXI (sous INCREMENT).
"""
    ),
    17: _(
        """
 Pas de temps maximal (mot-clé PAS_MAXI) demandé trop grand :   %(r1)f
 Pas de temps nécessaire pour le calcul: %(r2)f
 Risques de problèmes de précision
"""
    ),
    18: _(
        """
 Le nombre maximal de sous division du pas : %(i1)d est atteint à l'instant : %(r1)f
 Le pas de temps vaut alors : %(r2)f
 On continue cependant la résolution en passant au pas suivant.

 Risque & Conseil : la solution calculée risque d'être imprécise.
 Il faudrait relancer la calcul en autorisant le schéma ADAPT à utiliser un pas de temps plus petit.
 Pour cela on peut jouer sur au moins un des trois paramètres suivants :
 - diminuer le pas de temps initial (mot-clé PAS),
 - augmenter le nombre maximal de sous découpages du pas (mot-clé NMAX_ITER_PAS),
 - augmenter le facteur de division du pas (mot-clé COEF_DIVI_PAS)
"""
    ),
    19: _(
        """
 Le chargement contient plus d'une charge répartie.
 Le calcul n'est pas possible pour les modèles de poutre.
"""
    ),
    20: _(
        """
 La fréquence d'actualisation de FORCE_SOL est prise dans le fichier des raideurs.
"""
    ),
    21: _(
        """
 La fréquence d'actualisation de FORCE_SOL est prise dans le fichier des masses.
"""
    ),
    22: _(
        """
 La fréquence d'actualisation de FORCE_SOL est prise dans le fichier des amortissements.
"""
    ),
    23: _(
        """
    Nombre de fréquences: %(i1)d
    Intervalle des fréquences: %(r1)f
"""
    ),
    25: _(
        """
 La fréquence d'actualisation de FORCE_SOL n'est pas cohérente avec la fréquence d'archivage des résultats dans
 DYNA_NON_LINE.
"""
    ),
    26: _(
        """
 Deux des fréquences %(r1)f Hz et %(r2)f HZ de la liste LIST_RAFFINE sont proches.
Les intervalles de raffinement entourant ces deux valeurs se chevauchent.
Si une valeur du premier intervalle est trop proche d'une valeur du deuxième
intervalle (écart inférieur à PAS_MINI), l'une des deux sera supprimée de la liste.
"""
    ),
    28: _(
        """
 Le nombre d'obstacles de choc est limité à %(i1)d en cas de traitement implicite des non-linéarités dans
 DYNA_VIBRA.
 Conseil : si la modélisation ne permet pas de réduire le nombre de lieux de choc, il faudrait
 repasser le calcul mais avec un traitement explicite des chocs.
"""
    ),
    29: _(
        """
 La matrice d'amortissement n'est pas diagonale. Or, la méthode ITMI ne permet pas de modéliser
 de couplage dynamique par l'amortissement. Les termes diagonaux seront alors extraits pour la suite
 du calcul.
"""
    ),
    30: _(
        """
 La fréquence d'actualisation de FORCE_SOL dans le fichier des masses est incohérente avec
celle choisie précédemment.
"""
    ),
    31: _(
        """
 La fréquence d'actualisation de FORCE_SOL dans le fichier des amortissements est incohérente avec
celle choisie précédemment.
"""
    ),
    32: _(
        """
La condition de stabilité n'a pas pu être calculée car il s'agit d'élasticité non-isotrope.
"""
    ),
    33: _(
        """
Il y a une incohérence dans les type de résultats, le résultat selon X n'est pas le même que celui selon Y.
"""
    ),
    34: _(
        """
Il y a une incohérence dans les type de résultats, le résultat selon X et Y n'est pas le même que celui selon Z.
"""
    ),
    35: _(
        """
Il semble que vos calculs dynamiques ont été réalisés sur des listes de fréquences ou d'instants différentes.
L'utilisation de la macro-commande nécessite d'avoir réalisé les calculs dynamiques sur une unique liste de fréquences ou d'instants.
"""
    ),
    36: _(
        """
Les signaux d'entraînements ont une fréquence finale inférieure à celle du calcul dynamique.
"""
    ),
    37: _(
        """
Les signaux d'entraînements ont un instant final inférieur à celui du calcul dynamique.
"""
    ),
    38: _(
        """
Les signaux servant de supports à la détermination des signaux d'entraînement ont une fréquence finale inférieure à celle du calcul dynamique.
"""
    ),
    39: _(
        """
Les signaux servant de supports à la détermination des signaux d'entraînement ont un instant final inférieur à  celui du calcul dynamique.
"""
    ),
    40: _(
        """
Les signaux d'entraînements ne sont pas discrétisés de la même manière. Vérifier le pas de chaque signaux ainsi que leur longueur.
"""
    ),
    41: _(
        """
Les signaux servant de supports à la détermination des signaux d'entraînement ne sont pas discrétisés de la même manière. Vérifier le pas de chaque signaux ainsi que leur longueur.
"""
    ),
    42: _(
        """
 Le pas de temps est négatif ou nul : %(r1)f

 Vérifiez votre mise en données
"""
    ),
    43: _(
        """
 L'instant de fin de calcul (t=%(r1)f) est inférieur à l'instant initial (t=%(r2)f)

 Vérifiez votre mise en données
"""
    ),
    50: _(
        """
Schéma multi-pas
On n'a pas trouvé l'instant précédent dans la structure de données résultat du mot-clef ETAT_INIT.
C'est probablement parce qu'il n'y a pas assez d'instants archivés.
On ignore donc le calcul du second membre pour cet instant.
"""
    ),
    51: _(
        """
Schéma multi-pas
On n'a pas trouvé l'instant précédent dans la structure de données résultat du mot-clef ETAT_INIT.
C'est probablement parce que la structure de données vient d'un calcul statique (STAT_NON_LINE) ou d'une lecture directe (LIRE_RESU).
On ignore donc le calcul du second membre pour cet instant.
"""
    ),
    52: _(
        """
Schéma multi-pas
L'instant précédent et l'instant initial sont presque confondus.
On ignore donc le calcul du second membre pour cet instant.
"""
    ),
    53: _(
        """
Schéma multi-pas
On n'a pas de structure de données résultat dans le mot-clef ETAT_INIT parce que l'état initial est entré champ par champ.
On ignore donc le calcul du second membre pour cet instant.
"""
    ),
    54: _(
        """
  -> La numérotation des modes pour l'amortissement modal est incohérente avec la numérotation utilisée pour la dynamique transitoire
  -> Risque & Conseil :
     Soit le modèle n'est pas le même, soit, plus probablement, vous n'utilisez pas les mêmes conditions limites de Dirichlet dualisées (AFFE_CHAR_MECA) entre le calcul modal et la dynamique transitoire.
"""
    ),
    55: _(
        """
--------------------------------------------------------------------------------------

%(k1)s
=====================================================================================
                     Calcul %(k2)s sur base %(k3)s
=====================================================================================
"""
    ),
    56: _(
        """
Superposition modale classique
--------------------------------------------------------------------------------------
    >> base modale de projection : %(k1)s
    >> nombre de DDL avant projection (physiques) : %(i1)d"""
    ),
    57: _(
        """
Modèle de sous-structuration dynamique
--------------------------------------------------------------------------------------
    >> modèle généralisé : %(k1)s
    >> numérotation généralisée : %(k2)s"""
    ),
    58: _(
        """
Modèle sous interaction fluide-structure
--------------------------------------------------------------------------------------
    >> base de couplage fluide-élastique : %(k1)s
    >> vitesse du fluide  :%(r1)12.5e"""
    ),
    59: _(
        """    >> nombre de modes    : %(i1)d
    >> fréquence minimale :%(r1)12.5e
    >> fréquence maximale :%(r2)12.5e
"""
    ),
    60: _(
        """
Matrices dynamiques pour la résolution
--------------------------------------------------------------------------------------"""
    ),
    61: _(
        """    >> matrice de masse        : %(k1)s
    >> matrice de rigidité     : %(k2)s"""
    ),
    62: _(
        """    >> matrice d'amortissement : %(k1)s
"""
    ),
    63: _(
        """    >> amortissement modal diagonal
"""
    ),
    64: _(
        """    >> système conservatif, sans amortissement.
"""
    ),
    65: _(
        """    >> masse diagonale extraite de la base de couplage fluide-élastique
    >> rigidité diagonale extraite de la base de couplage fluide-élastique
    >> amortissement modal diagonal, extrait de la base de couplage fluide-élastique
"""
    ),
    66: _(
        """
Schéma d'intégration %(k1)s à pas adaptatif
--------------------------------------------------------------------------------------
    >> type de schéma                 : explicite
    >> pas d'intégration initial      :%(r1)12.5e
    >> pas d'intégration minimal      :%(r2)12.5e  (arrêt du calcul si inférieur)
    >> pas d'intégration maximal      :%(r3)12.5e  (plafond maximal du pas d'intégration)"""
    ),
    67: _("""    >> tolérance                      :%(r1)12.5e"""),
    68: _("""    >> coefficient de division du pas :%(r1)12.5e"""),
    69: _(
        """    >> nombre minimum de pas calculés : %(i1)d
    >> nombre maximum de pas calculés : %(i2)d"""
    ),
    70: _(
        """
Schéma d'intégration %(k1)s à pas constant
--------------------------------------------------------------------------------------
    >> type de schéma         : %(k2)s
    >> pas d'intégration      :%(r1)12.5e
    >> nombre de pas calculés : %(i1)d"""
    ),
    71: _(
        """
Non-linéarités localisées
--------------------------------------------------------------------------------------"""
    ),
    72: _(
        """    >> nombre de lieux de choc                  : %(i1)d
    >> méthode de traitement de chocs           : %(k1)s"""
    ),
    73: _("""    >> nombre de dispositifs anti-sismique      : %(i1)d"""),
    74: _("""    >> nombre de lieux de choc avec flambement  : %(i1)d"""),
    75: _("""    >> nombre de relations effort-déplacement   : %(i1)d"""),
    76: _("""    >> nombre de relations effort-vitesse       : %(i1)d"""),
    77: _("""    >> nombre de couplages de type %(k1)s     : %(i1)d"""),
    78: _(
        """
État initial
--------------------------------------------------------------------------------------
    >> extrait à partir du résultat : %(k1)s
"""
    ),
    79: _(
        """
État initial
--------------------------------------------------------------------------------------
    >> déplacement  : %(k1)s
    >> vitesse      : %(k2)s"""
    ),
    80: _(
        """
Durée de la simulation
--------------------------------------------------------------------------------------
    >> instant initial :%(r1)12.5e
    >> instant final   :%(r2)12.5e"""
    ),
    81: _(
        """
Archivage
--------------------------------------------------------------------------------------
    >> fréquence d'archivage         : tous les %(i1)d pas calculés
    >> instants d'archivage forcés   : instants initial et final, et %(i2)d autres instants
    >> champs archivés (généralisés) : DEPL, VITE, ACCE

======================================================================================

Avancement du calcul
--------------------------------------------------------------------------------------
"""
    ),
    82: _(
        """
Modèle physique
--------------------------------------------------------------------------------------
    >> modèle mécanique : %(k1)s
    >> nombre de DDL physiques : %(i1)d
"""
    ),
    83: _(
        """
La méthode intégrale (ITMI) est uniquement disponible si les matrices dynamiques sont
diagonales. Vérifiez que le stockage diagonal a été choisi lors de la numérotation des
DDL généralisés.
"""
    ),
    84: _("""    >> accélération : %(k1)s"""),
    85: _(
        """
Archivage
--------------------------------------------------------------------------------------
    >> fréquence d'archivage : tous les %(i1)d pas calculés"""
    ),
    86: _(
        """
Archivage
--------------------------------------------------------------------------------------
    >> nombre d'instants d'archivage : %(i1)d instants"""
    ),
    87: _("""    >> matrice d'impédance     : %(k1)s"""),
    88: _(
        """
Résolution
--------------------------------------------------------------------------------------
    >> fréquence minimale :%(r1)12.5e
    >> fréquence maximale :%(r2)12.5e
    >> nombre de fréquences calculées : %(i1)d

Archivage
--------------------------------------------------------------------------------------
    >> fréquences archivées : toutes (%(i1)d fréquences)
    >> champs archivés      : %(k1)s

======================================================================================

Avancement du calcul
--------------------------------------------------------------------------------------
"""
    ),
    90: _(
        """
Entrée/Changement d'état de choc détecté à l'instant %(r1)12.5e
---------------------------------------------------------------------------------------------
    >> Descriptif de l'état :"""
    ),
    91: _(
        """           %(k1)s
          | Choc numéro | %(k2)s"""
    ),
    92: _(
        """          | État        %(k1)s
           %(k2)s"""
    ),
    93: _(
        """    >> Repassage en état de vol libre détecté à l'instant %(r1)12.5e

"""
    ),
    94: _(
        """    >> Premier passage dans cet état, mise à jour des matrices dynamiques et calcul d'une
       nouvelle base de modes propres :
       """
    ),
    96: _(
        """    >> champs archivés       : %(k1)s

======================================================================================

Avancement du calcul
--------------------------------------------------------------------------------------
"""
    ),
    97: _(
        """
======================================================================================
                                                                        Fin du calcul

"""
    ),
    98: _("""    >> nombre de modélisations de rotor fissuré : %(i1)d"""),
    99: _("""    >> nombre de modélisations de palier : %(i1)d"""),
}
