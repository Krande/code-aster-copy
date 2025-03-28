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
 Erreur d'utilisation de POST_CHAMP :
   Dans la structure de données %(k2)s,
   vous avez demandé l'extraction du champ %(k1)s pour le numéro d'ordre %(i1)d.
   Mais ce champ n'existe pas.
"""
    ),
    4: _(
        """
Problème création CHAM_ELEM nul
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    6: _(
        """
Erreur utilisateur :
 Vous utilisez le mot clé NOM_CMP, mais l'une (au moins) des composantes indiquées
 n'appartient pas à la grandeur : %(k1)s
"""
    ),
    9: _(
        """
Erreur utilisateur :
  Vous ne pouvez pas utiliser la méthode ECLA_PG avec le mot-clé RESULTAT.
Conseil :
   Extrayez le champ aux ELGA que contient votre résultat puis utilisez la méthode ECLA_PG avec le mot-clé CHAM_GD.
"""
    ),
    20: _(
        """
 le GROUP_NO  %(k1)s  contient  %(k2)s  noeuds
"""
    ),
    21: _(
        """
 le GROUP_MA  %(k1)s  contient  %(k2)s  mailles
"""
    ),
    28: _(
        """
PROJ_CHAMP :
  La méthode SOUS_POINT accepte uniquement des évolutions de type TEMP HYDR NEUT
  pour définir des évolutions de variables de commandes EVOL_VARC"""
    ),
    29: _(
        """
PROJ_CHAMP :
  La méthode SOUS_POINT accepte uniquement des champs de type TEMP HYDR NEUT SIEF
  Un champ de type %(k1)s n'est actuellement pas géré."""
    ),
    30: _(
        """
PROJ_CHAMP :
  La méthode SOUS_POINT accepte uniquement les résultats de type
  EVOL_THER."""
    ),
    31: _(
        """
PROJ_CHAMP :
  Le mot-clé %(k1)s est interdit avec la méthode SOUS_POINT."""
    ),
    32: _(
        """
PROJ_CHAMP (ou LIAISON_MAILLE) :
  La méthode %(k1)s est incompatible avec les champs aux noeuds."""
    ),
    33: _(
        """
PROJ_CHAMP (ou LIAISON_MAILLE) :
  La méthode %(k1)s est incompatible avec les champs par élément de type %(k2)s."""
    ),
    34: _(
        """
   Maillage quadratique obligatoire avec terme source non nul."""
    ),
    35: _(
        """
PROJ_CHAMP (ou LIAISON_MAILLE) :
  Vous cherchez à projeter un champ par élément (ELGA).
  Pour cela, il vous faut renseigner le mot-clé MODELE_1."""
    ),
    36: _(
        """
PROJ_CHAMP (ou LIAISON_MAILLE) :
  Le mot-clé TYPE_CHAM est incompatible avec le mot-clé CHAM_GD.
  Il n'est utilisable qu'avec le mot-clé RESULTAT."""
    ),
    37: _(
        """
PROJ_CHAMP (ou LIAISON_MAILLE) :
  Vous cherchez à projeter un champ par élément (ELNO, ELEM ou ELGA).
  Pour cela, il vous faut renseigner le mot-clé MODELE_2."""
    ),
    38: _(
        """
  il faut définir un champ de vitesse
"""
    ),
    39: _(
        """
 la grandeur pour la variable:  %(k1)s  doit être:  %(k2)s  mais elle est:  %(k3)s
"""
    ),
    40: _(
        """
PROJ_CHAMP  :
  Vous utilisez la méthode SOUS_POINT.
  Pour cela, il vous faut renseigner le mot-clé  %(k1)s."""
    ),
    41: _(
        """
PROJ_CHAMP  :
  Vous utilisez la méthode SOUS_POINT. Dans ce cas :
  PROJECTION = "OUI" ET le mot clef MATR_PROJECTION ne doit pas être renseigné."""
    ),
    43: _(
        """
PROJ_CHAMP (ou LIAISON_MAILLE) :
  Le noeud %(k1)s de coordonnées (%(r1)e,%(r2)e,%(r3)e) est projeté à la distance %(r4)e"""
    ),
    44: _(
        """
Le champ doit être un CHAM_ELEM
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    45: _(
        """
Longueurs des modes locaux incompatibles entre eux
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    46: _(
        """
Terme normalisation global nul
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    47: _(
        """
 PROJ_CHAMP, occurrence %(i1)d de VIS_A_VIS :

 Le noeud %(k1)s est jugé distant de la maille sur laquelle il devrait être
 projeté. Comme la projection de ce noeud a déjà été faite lors d'une occurrence
 précédente, on ne tient pas compte de cette projection distante.

 Conseil : il est cependant recommandé de vérifier l'affectation des groupes dans VIS_A_VIS.
"""
    ),
    48: _(
        """
 Vous utilisez la commande PROJ_CHAMP ou un mot clé nécessitant de "projeter"
 des noeuds sur des mailles (par exemple LIAISON_MAIL).
 Il y a %(i1)d noeuds qui ont été projetés sur des mailles jugées distantes.

 Les noeuds sont jugés distants si :
  * la distance à la maille la plus proche est supérieure à DISTANCE_ALARME
    (si ce mot clé est utilisé).
  * ou si la distance à la maille la plus proche est supérieure à 1/10ème
    de la taille de cette maille (si le mot clé DISTANCE_ALARME n'est pas utilisé).

 Les %(i2)d noeuds les plus éloignés ont été imprimés ci-dessus.

Risques et conseils :
  * Un maillage constitué des noeuds distants a été imprimé au format MED.
    La visualisation de ce maillage pourra vous rassurer (ou non).
    Le nom du fichier MED est : %(k1)s
  * Le mot clé DISTANCE_MAX permet d'éviter que les noeuds trop distants ne soient
    projetés (ou "liés" quand on utilise le mot clé LIAISON_MAIL).
  * Le mot clé DISTANCE_ALARME permet d'éviter cette alarme.

"""
    ),
    52: _(
        """
 Le calcul du champ SIGM_ELNO n'a pas été fait sur la maille volumique %(k1)s qui borde
 la maille surfacique %(k2)s.

 Conseils :
  Il faut faire le calcul du champ SIGM_ELNO sur les éléments volumiques de l'autre "coté"
  de la face en choisissant le bon groupe de mailles soit en faisant le calcul sur tout
  le volume.
  Il est aussi possible de supprimer le calcul préalable de SIGM_ELNO, le calcul sera fait
  automatiquement sur les bonnes mailles volumiques.
"""
    ),
    53: _(
        """
 La SUPER_MAILLE %(k1)s n'existe pas dans le maillage %(k2)s.
"""
    ),
    54: _(
        """
 Aucune maille de peau n'a été fournie.

 Vous devez renseigner le mot-clé MAILLE/GROUP_MA en donnant une liste de mailles ou
 un groupe de maille contenant des mailles de peau.
 Si vous avez renseigné le mot-clé TOUT='OUI', cela signifie qu'il n'y a pas de mailles
 de peau dans votre modèle ; il faut revoir le maillage.
"""
    ),
    55: _(
        """
Alarme utilisateur :
  Vous avez utilisé le mot clé LIAISON_SOLIDE pour solidifier un ensemble
  de noeuds.
  Le nuage formé par ces noeuds est volumique mais il est très aplati.
  Le rapport entre les dimensions 3 et 1 est faible : %(r1)f
  Les relations cinématiques engendrées peuvent être proches de la
  redondance et provoquer des problèmes de type "pivot nul".

Risques et Conseils :
  En utilisant le mot clé DIST_MIN, vous pouvez faire en sorte que le
  programme considère le nuage de points comme surfacique.
  Pour cela, vous devez choisir un DIST_MIN > %(r2)f
"""
    ),
    56: _(
        """
Alarme utilisateur :
  Vous avez utilisé le mot clé LIAISON_SOLIDE pour solidifier un ensemble
  de noeuds.
  Le nuage formé par ces noeuds est surfacique mais il est très allongé.
  Le rapport entre les dimensions 2 et 1 est faible : %(r1)f
  Les relations cinématiques engendrées peuvent être proches de la
  redondance et provoquer des problèmes de type "pivot nul".

Risques et Conseils :
  En utilisant le mot clé DIST_MIN, vous pouvez faire en sorte que le
  programme considère le nuage de points comme linéique.
  Pour cela, vous devez choisir un DIST_MIN > %(r2)f
"""
    ),
    57: _(
        """
Erreur utilisateur dans la commande PROJ_CHAMP :
  La structure de données résultat à projeter ne contient que des champs 'ELGA'.
  La méthode de projection adaptée à ces champs est la méthode 'ECLA_PG' mais
  elle ne fonctionne qu'avec un champ isolé (mot clé CHAM_GD).
"""
    ),
    58: _(
        """
Erreur utilisateur dans la commande PROJ_CHAMP :
  On cherche à projeter un champ par éléments 'ELGA' isolé (mot clé CHAM_GD).
  La méthode de projection doit être 'ECLA_PG' et non pas %(k1)s
"""
    ),
    65: _(
        """
 composante non définie dans  la grandeur composante:  %(k1)s
"""
    ),
    66: _(
        """

 le nombre de composantes affectées n'est pas égal  au nombre de composantes à affecter
 occurrence de AFFE numéro %(i1)d
 nombre de composante affectées :  %(i2)d
 nombre de composante a affecter :  %(i3)d
"""
    ),
    67: _(
        """
 erreurs données le GROUP_MA  %(k1)s
  n'a pas le même nombre de mailles  que le GROUP_MA  %(k2)s
"""
    ),
    68: _(
        """
 erreurs données le GROUP_MA  %(k1)s
  n'a pas les mêmes types de maille  que le GROUP_MA  %(k2)s
"""
    ),
    69: _(
        """
 Problème lors de la vérification de correspondance entre le
 GROUP_MA_INIT '%(k1)s' et le GROUP_MA_FINAL '%(k2)s' :

 - Aucune maille du GROUP_MA_FINAL '%(k2)s' ne semble être en correspondance
 avec la maille '%(k3)s' du GROUP_MA_INIT '%(k1)s'.

 - Vérifiez que le groupe '%(k2)s' (GROUP_MA_FINAL) du maillage %(k5)s
 (MAILLAGE_FINAL) est bien la translation du groupe '%(k1)s' (GROUP_MA_INIT)
 du maillage  %(k4)s (MAILLAGE_INIT), suivant le vecteur de translation
   ( %(r1)f , %(r2)f , %(r3)f ), avec la précision donnée.

 - Une erreur fréquente lors de l'utilisation de la commande PERM_MAC3COEUR
 est que les 2 assemblages à permuter ne sont pas du même type de conception.
"""
    ),
    70: _(
        """
 l'instant  de calcul  %(r1)f  n'existe pas dans  %(k1)s
"""
    ),
    71: _(
        """
 plusieurs numéros d'ordre trouves pour l'instant  %(r1)f
"""
    ),
    72: _(
        """
 Cette commande est réentrante :   structure de données résultat en sortie     %(k1)s
 Structure de données RESU_FINAL  %(k2)s
"""
    ),
    73: _(
        """
La structure de données résultat en sortie  %(k1)s doit contenir qu'un seul NUME_ORDRE %(k2)s
"""
    ),
    76: _(
        """
 Il n'est pas encore possible de découper le type_élément :  %(k1)s  en sous-éléments
    élément fini :  %(k2)s
    famille      :  %(k3)s.
 Faites une demande d'évolution.
"""
    ),
    78: _(
        """
 Il n'est pas encore possible de découper le type_élément :  %(k1)s  en sous-éléments
    élément fini :  %(k2)s.
 Faites une demande d'évolution.
"""
    ),
    79: _(
        """
Vous utilisez la commande PROJ_CHAMP (ou une fonctionnalité utilisant la
projection d'un maillage sur un autre).
Pour des raisons de performance, les mailles du maillage à projeter sont
placées dans une grille cartésienne régulière.

Le maillage est tel que certaines cellules de la grille contiennent beaucoup
plus de mailles que d'autres.
Dans ces conditions, les performances CPU de la projection sont dégradées.

Nombre moyen de mailles dans une cellule   : %(i1)d
Nombre maximal de mailles dans une cellule : %(i2)d

Conseil : dans ce type de cas de figure, l'utilisation de VIS_A_VIS peut
améliorer très significativement les performances de PROJ_CHAMP.

"""
    ),
    85: _(
        """
 Problème liste de mailles carte : %(k1)s  numéro entité : %(i1)d
  position dans liste : %(i2)d
  numéro de maille  : %(i3)d
"""
    ),
    86: _(
        """
 Le champ %(k1)s  n'appartient pas au résultat qu'on essaie d'imposer comme
  condition limite. Rajouter ce champ au préalable à travers un CALC_CHAMP.
"""
    ),
}
