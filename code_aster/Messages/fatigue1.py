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
    1: _(
        """
 Aucun élément du modèle ne sait calculer l'option
 de fatigue multiaxiale,
 Il se peut que la modélisation affectée au groupe de mailles
 sur lequel vous faites un calcul de fatigue ne soit pas "3D".

 Le critère de fatigue que vous utilisez n'est utilisable qu'en 3D.

"""
    ),
    2: _(
        """
 La modélisation affectée au groupe de mailles sur lequel vous
 faites un calcul de fatigue n'est probablement pas "3D".
 La composante %(i1)d du tenseur des contraintes n'existe pas.

 Le critère de fatigue que vous utilisez n'est utilisable qu'en 3D.

"""
    ),
    3: _(
        """
 La modélisation affectée au groupe de mailles sur lequel vous
 faites un calcul de fatigue n'est probablement pas "3D".
 La composante %(i1)d du tenseur des déformations n'existe pas.

 Le critère de fatigue que vous utilisez n'est utilisable qu'en 3D.

"""
    ),
    4: _(
        """
 le coefficient de GOODMAN n'est pas calculable
"""
    ),
    5: _(
        """
 le coefficient de Gerber n'est pas calculable
"""
    ),
    6: _(
        """
 pour calculer le dommage de Lemaitre-Sermage,
 il faut définir le comportement DOMMA_LEMAITRE dans DEFI_MATERIAU
"""
    ),
    7: _(
        """
 pour calculer le dommage de Lemaitre-Sermage,
 il faut définir le comportement ELAS_FO dans DEFI_MATERIAU
"""
    ),
    8: _(
        """
 le matériau est obligatoire pour le calcul du dommage par TAHERI_MANSON
"""
    ),
    9: _(
        """
 une fonction doit être introduite sous le mot clé TAHERI_FONC
"""
    ),
    10: _(
        """
 une nappe doit être introduite sous le mot clé TAHERI_NAPPE
"""
    ),
    11: _(
        """
 la courbe de MANSON_COFFIN est nécessaire pour le calcul du dommage TAHERI_MANSON_COFFIN
"""
    ),
    12: _(
        """
 le matériau est obligatoire pour le calcul du dommage par TAHERI_MIXTE
"""
    ),
    13: _(
        """
 la courbe de MANSON_COFFIN est nécessaire pour le calcul du dommage TAHERI_MIXTE
"""
    ),
    14: _(
        """
 la courbe de WOHLER est nécessaire pour le calcul du dommage TAHERI_MIXTE
"""
    ),
    15: _(
        """
 méthode de comptage inconnue
"""
    ),
    16: _(
        """
 nombre de cycles nul
"""
    ),
    17: _(
        """
 l'utilisation de MANSON_COFFIN est réservé à des histoires de chargements en déformations
"""
    ),
    18: _(
        """
 la courbe de MANSON_COFFIN doit être donnée dans DEFI_MATERIAU
"""
    ),
    19: _(
        """
 les lois de TAHERI sont réservées pour des chargements en déformations
"""
    ),
    20: _(
        """
 loi de dommage non compatible
"""
    ),
    21: _(
        """
 l'histoire de chargement doit avoir même discrétisation pour toutes les composantes
"""
    ),
    22: _(
        """
 l'histoire de la déformation plastique cumulée doit avoir même discrétisation que l'histoire des contraintes
"""
    ),
    23: _(
        """
 l'histoire de la température doit avoir même discrétisation que l'histoire des contraintes
"""
    ),
    24: _(
        """
 pour calculer le dommage, il faut définir le comportement "FATIGUE" dans DEFI_MATERIAU
"""
    ),
    25: _(
        """
 la méthode 'TAHERI_MANSON' ne peut pas être utilisée avec l'option %(k1)s
"""
    ),
    26: _(
        """
 le nom de la fonction nappe doit être présent sous le mot clé 'TAHERI_NAPPE'
"""
    ),
    27: _(
        """
 le nom de la fonction doit être présent sous le mot clé 'TAHERI_FONC'
"""
    ),
    28: _(
        """
 la méthode 'TAHERI_MIXTE' ne peut pas être utilisée avec l'option %(k1)s
"""
    ),
    29: _(
        """
 la méthode 'WOHLER' ne peut pas être utilisée avec l'option %(k1)s
"""
    ),
    30: _(
        """
 une courbe de WOHLER doit être définie dans DEFI_MATERIAU
"""
    ),
    31: _(
        """
 la méthode 'MANSON_COFFIN' ne peut pas être utilisée avec l'option %(k1)s
"""
    ),
    32: _(
        """
 une courbe de MANSON_COFFIN doit être définie dans DEFI_MATERIAU
"""
    ),
    33: _(
        """
 les mailles attachées au noeud traité ne sont pas affectées du même matériau.
"""
    ),
    34: _(
        """
 la donnée d'une courbe de WOHLER est obligatoire
"""
    ),
    35: _(
        """
 la donnée du moment spectral d'ordre 4 est obligatoire pour le comptage des pics de contraintes
"""
    ),
    36: _(
        """
 la valeur du moment spectral d'ordre 0 (lambda_0) est certainement nulle
"""
    ),
    37: _(
        """
 la valeur du moment spectral d'ordre 2 (lambda_2) est nulle
"""
    ),
    38: _(
        """
 la valeur du moment spectral d'ordre 4 (lambda_4) est nulle
"""
    ),
    39: _(
        """
Le chargement à compter est un chargement constant. On considère tous les chargements comme un cycle
avec valeur_max = valeur_min = valeur du chargement, i.e., amplitude = 0.
"""
    ),
    63: _(
        """
 pour calculer le dommage max,
 il faut renseigner CISA_PLAN_CRIT dans la commande DEFI_MATERIAU
"""
    ),
    64: _(
        """
 nous ne pouvons pas récupérer la valeur du paramètre A du critère de MATAKE
 commande: DEFI_MATERIAU
 opérande: CISA_PLAN_CRIT
"""
    ),
    65: _(
        """
 nous ne pouvons pas récupérer la valeur du paramètre B du critère de MATAKE
 commande: DEFI_MATERIAU
 opérande: CISA_PLAN_CRIT
"""
    ),
    66: _(
        """
 nous ne pouvons pas récupérer la valeur du coefficient de passage flexion-torsion
 commande: DEFI_MATERIAU
 opérande: CISA_PLAN_CRIT
"""
    ),
    67: _(
        """
 nous ne pouvons pas récupérer la valeur du paramètre A du critère de DANG_VAN_MODI_AC
 commande: DEFI_MATERIAU
 opérande: CISA_PLAN_CRIT
"""
    ),
    68: _(
        """
 nous ne pouvons  pas récupérer la valeur du paramètre B du critère de DANG_VAN_MODI_AC
 commande: DEFI_MATERIAU
 opérande: CISA_PLAN_CRIT
"""
    ),
    69: _(
        """
 nous ne pouvons  pas récupérer la valeur du coefficient de passage cisaillement-traction
 commande: DEFI_MATERIAU
 opérande: CISA_PLAN_CRIT
"""
    ),
    70: _(
        """
 nous ne pouvons pas récupérer la valeur du paramètre A de MATAKE
 commande: DEFI_MATERIAU
 opérande: CISA_PLAN_CRIT
"""
    ),
    71: _(
        """
 nous ne pouvons pas récupérer la valeur du paramètre B de MATAKE
 commande: DEFI_MATERIAU
 opérande: CISA_PLAN_CRIT
"""
    ),
    72: _(
        """
 nous ne pouvons pas récupérer la valeur du coefficient de passage cisaillement-traction
 commande: DEFI_MATERIAU
 opérande: CISA_PLAN_CRIT
"""
    ),
    73: _(
        """
 nous ne pouvons pas récupérer la valeur du paramètre A du critère DANG_VAN_MODI_AV
 commande: DEFI_MATERIAU
 opérande: CISA_PLAN_CRIT
"""
    ),
    74: _(
        """
 nous ne pouvons pas récupérer la valeur du paramètre B du critère DANG_VAN_MODI_AV
 commande: DEFI_MATERIAU
 opérande: CISA_PLAN_CRIT
"""
    ),
    75: _(
        """
 nous ne pouvons pas récupérer la valeur du paramètre A du critère FATESOCI
 commande: DEFI_MATERIAU
 opérande: CISA_PLAN_CRIT
"""
    ),
    76: _(
        """
 Le nombre d'instants calculés est égal à %(i1)d.

 Il faut que l'histoire du chargement comporte au moins 2 instants
 pour calculer un dommage.

"""
    ),
    77: _(
        """
 Le nombre de cycles extraits est égal à %(i1)d.

 Votre chargement est constant.
 On ne peut donc pas extraire de cycles.

"""
    ),
    78: _(
        """
 Le nombre de points à traiter n'est pas correcte.

 Soit les mailles comportent des sous-points, or ce cas n'est pas prévu.

 Soit le nombre total de composantes de l'option : %(k1)s a changé et n'est plus égal à %(i1)d.

"""
    ),
    79: _(
        """
   *** Point  %(i1)d
   Contrainte statique        =  %(r1)f
   Contrainte dynamique       =  %(r2)f
   Amplitude maximale admissible en ce point  =  %(r3)f
"""
    ),
    80: _(
        """
   Attention, la contrainte statique en ce point est supérieure à la contrainte
   à la rupture du matériau.
"""
    ),
    81: _(
        """
 Calcul du dommage en %(k1)s (composante grandeur équivalente %(k3)s)
 Points de calcul du dommage : %(k2)s
 Nombre de points de calcul : %(i1)d
 Nombre de modes considérés : %(i2)d
"""
    ),
    82: _(
        """

 Amplitude de vibration maximale admissible par la structure :  %(r1)f

"""
    ),
    83: _(
        """
   Attention, la contrainte statique en un ou plusieurs points est supérieure à
   la contrainte à la rupture du matériau.
   Contrainte statique maximale = %(r1)f
   Contrainte à la rupture      = %(r2)f
"""
    ),
    84: _(
        """
   Le résultat correspondant à la contrainte statique (mot clé
   RESULTAT_STATIQUE) comporte %(i1)d instants.
   Le calcul en fatigue vibratoire n'est possible que si le résultat statique
   comporte un et un seul instant. Vérifiez les données.
"""
    ),
    85: _(
        """
   Le nombre de points de calcul est différent entre la contrainte statique
   ( %(i1)d points) et la contrainte modale (%(i2)d points).
   Vérifiez la cohérence des données.
"""
    ),
    86: _(
        """
   La longueur de la liste des coefficients modaux COEF_MODE est différente
   du nombre de modes retenus pour le calcul NUME_MODE.
   Vérifiez la cohérence des données.
"""
    ),
    87: _(
        """
 Contrainte à la rupture : %(r1)f
 Limite d'endurance : %(r2)f
"""
    ),
    88: _(
        """
 Le paramètre %(k1)s est absent dans la définition du matériau. Le calcul est impossible.
 Risques et conseils : pour l'option FATIGUE_VIBR de CALC_FATIGUE, il est obligatoire de
 définir les propriétés matériaux suivantes (dans DEFI_MATERIAU) :
 - la contrainte à la rupture (opérande SU, mot clé facteur RCCM) ;
 - la limite d'endurance (première abscisse de la courbe de fatigue définie dans l'opérande
 WOHLER du mot clé facteur FATIGUE).

"""
    ),
    89: _(
        """
 Dans la commande DEFI_MATERIAU, l'opérande WOHLER du mot clé facteur
 FATIGUE est incompatible avec le critère %(k1)s requis dans le mot
 clé facteur CISA_PLAN_CRIT.

"""
    ),
    90: _(
        """
 Dans la commande DEFI_MATERIAU, l'opérande MANSON_COFFIN du mot clé
 facteur FATIGUE est incompatible avec le critère %(k1)s requis dans
 le mot clé facteur CISA_PLAN_CRIT.

"""
    ),
    91: _(
        """
 Dans CALC_FATIGUE/POST_FATIGUE  pour le critère d'amorçage fournis par la formule, le calcul
 de la grandeur %(k1)s n'est pas disponible. Merci de fournir les noms des grandeurs
 disponibles ou contacter le développeur.
"""
    ),
    92: _(
        """
 Dans CALC_FATIGUE/POST_FATIGUE, pour le critère d'amorçage fournis par la formule, pour déterminer
 le plan de dommage maximal, il n'est pas possible de projeter simultanément la contrainte
 et la déformation. Les grandeurs sont incompatibles avec le critère requis.
"""
    ),
    93: _(
        """
 Dans CALC_FATIGUE/POST_FATIGUE, pour le critère d'amorçage fournis par la formule, le mot-clé
 FORMULE_VIE est fournis par une formule le seul paramètre accepté est N_f,
 c'est-à-dire, N_f, car le critère formule est pour GRDEQ = f(N_f).
 Changez le nom et vérifiez bien que la fonction est de type: GRDEQ = f(N_f).
"""
    ),
    94: _(
        """
 Dans CALC_FATIGUE/POST_FATIGUE, pour le critère d'amorçage fournis par la formule et le mot-clé
 FORMULE_VIE est fournis par la formule, la grandeur équivalente pour l'instant est
 plus grande que f(N_f =1 ). Vérifiez la formule de la courbe FORMULE_VIE.
"""
    ),
    95: _(
        """
 Pour le critère utilisé, au moins d'une histoire de la tenseur de déformation est demandée.
"""
    ),
    96: _(
        """
 On note que DEFORMATION ELASTIQUE  = DEFORMATION TOTALE - DEFORMATION PLASTIQUE. Si les déformations
 totale ou plastique ne sont pas fournies, on prendre la valeur zéro.
"""
    ),
    97: _(
        """
 Pour le critère utilisé, l'histoire de la tenseur de contrainte(SIGM_XX...) est demandée.
"""
    ),
    98: _(
        """
 Pour le critère utilisé, l'histoire de la tenseur de déformation totale (EPS_XX...)est demandée.
"""
    ),
    99: _(
        """
 Pour le critère utilisée, l'histoire de la tenseur de déformation plastique (EPSP_XX...) est demandée.
"""
    ),
    53: _(
        """
 Dans CALC_FATIGUE/POST_FATIGUE  pour le critère d'amorçage du type plan critique, aucune grandeur
 de la FORMULE_CRITIQUE dépend de l'orientation du plan. Veuillez vérifier et/ou enlever cet mot-clé.
"""
    ),
}
