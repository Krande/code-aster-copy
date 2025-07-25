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
    2: _(
        """
 le CHAMP_S:  %(k1)s  est a la fois CHAM_ELEM_S et CHAM_NO_S.
"""
    ),
    3: _(
        """
 le CHAMP_S:  %(k1)s n'existe pas.
"""
    ),
    4: _(
        """
Erreur de programmation :
 On essaye de calculer l'intégrale d'un CHAM_ELEM / ELGA.
 Malheureusement, la famille de points de Gauss : %(k1)s n'est pas autorisée dans la programmation.

Conseil :
  Si nécessaire, il faut demander une évolution du code.
"""
    ),
    5: _(
        """
Erreur dans CREA_RESU :
  Quand on utilise la commande CREA_RESU avec le mot clé AFFE / CHAM_GD, les
  composantes du champ de fonctions %(k2)s (de la géométrie et/ou du temps),
  doivent être au même rang que celles du champ de réels %(k1)s.
Incohérence dans les grandeurs %(k1)s et %(k2)s :
  Le rang de la composante %(k3)s de %(k1)s correspond
  au rang de la composante %(k4)s de %(k2)s.

Conseil :
  Si nécessaire, il faut demander une évolution du code.
"""
    ),
    6: _(
        """
Erreur utilisateur dans PROJ_CHAMP :
  Le champ utilisé dans le mot clé CHAM_NO_REFE (%(k1)s) est associé au maillage %(k2)s
  Il doit être associé au maillage %(k3)s
"""
    ),
    7: _(
        """
 trop d'antécédents
 vérifiez si le maillage de l'interface ne contient pas de noeuds coïncidents ou diminuez DIST_REFE.
"""
    ),
    8: _(
        """
  %(k1)s  valeurs de CHAMNO de déplacement n'ont pas été recopiées sur  %(k2)s  noeuds
  a affecter  ce qui peut entraîner des erreurs de calcul sur la masse ajoutée des sous structures
  déduites par rotation et translation définies dans le modèle  généralisé. augmentez DIST_REFE
  ou assurez vous de l' invariance du maillage de structure par la translation et la rotation
  définies dans le modèle généralisé.
"""
    ),
    9: _(
        """
  -> plus de 50 %% des valeurs de CHAM_NO de déplacement n'ont pas été recopiées
     ce qui peut entraîner des erreurs graves de calcul sur la masse ajoutée des
     sous structures déduites par rotation et translation définies dans le modèle généralisé
  -> Risque & Conseil :
     augmentez DIST_REFE
"""
    ),
    10: _(
        """
 trop de noeuds affectés
"""
    ),
    11: _(
        """
 Erreur d'utilisation :
   Le maillage associé au modèle : %(k1)s
   n'est pas le même que celui du champ de matériaux : %(k2)s
"""
    ),
    12: _(
        """
Erreur lors de la transformation du CHAM_NO_S (%(k1)s) en CHAM_NO (%(k2)s):
Le CHAM_NO_S est vide (i.e. il n'a aucune valeur).
"""
    ),
    13: _(
        """
Erreur lors d'une transformation d'un CHAM_NO_S en CHAM_NO :
  Il manque la composante: %(k1)s  sur le noeud: %(k2)s pour le CHAM_NO: %(k3)s

Risques & conseils :
  Si cette erreur se produit dans la commande CREA_CHAMP, il est possible de
  mettre à zéro les composantes manquantes en utilisant le mot-clé PROL_ZERO.
"""
    ),
    14: _(
        """
Erreur utilisateur dans la commande POST_CHAMP :
 On demande l'extraction des champs sur une couche de numéro supérieur au nombre de couches.
"""
    ),
    15: _(
        """
Erreur utilisateur dans la commande POST_CHAMP :
 On demande l'extraction pour des champs n'ayant pas de "sous-points".
"""
    ),
    16: _(
        """
Erreur utilisateur dans la commande POST_CHAMP :
 On demande l'extraction des champs sur une fibre de numéro supérieur au nombre de fibres.
"""
    ),
    17: _(
        """
Erreur utilisateur dans la commande CREA_CHAMP :
 Incohérence entre le champ %(k1)s associé au maillage %(k2)s
 et le maillage %(k3)s
"""
    ),
    18: _(
        """
Erreur utilisateur dans la commande POST_RCCM / MOMENT_EQUIVALENT :
Le champ calculé est vide pour le numéro d'ordre %(i1)d .
"""
    ),
    19: _(
        """
Erreur utilisateur dans la commande POST_CHAMP / COQUE_EXCENT :
 Pour l'occurrence %(i1)d du mot clé COQUE_EXCENT,
 et pour le numéro d'ordre %(i2)d le champ calculé est vide.
"""
    ),
    20: _(
        """
Erreur utilisateur dans la commande POST_CHAMP / COQUE_EXCENT :
 La structure de donnée produite est vide.
"""
    ),
    21: _(
        """
 grandeur :  %(k1)s  inexistante au catalogue
"""
    ),
    22: _(
        """
 composante :  %(k1)s  inexistante au catalogue pour la grandeur : %(k2)s
"""
    ),
    23: _(
        """
 la grandeur : %(k1)s  n est pas de type réel.
"""
    ),
    24: _(
        """
Erreur utilisateur dans la commande MACR_ECLA_PG
  On ne sait pas où sont situés les points de Gauss du CHAM_ELEM (ELGA) %(k1)s.
  Cela arrive par exemple pour les champs correspondants à
  NOM_CHAM= ('VARI_ELGA', 'UT01_ELGA',...)
  car il n'existe pas d'option de calcul pour ces NOM_CHAM.

Conseil :
  Pour pouvoir post-traiter ce CHAM_ELEM, il faut lui associer
  un champ qui ne pose pas problème et qui partage la même localisation de ses
  points de Gauss.
  Par exemple, pour post-traiter VARI_ELGA, on fera :
  MACR_ECLA_PG(... NOM_CHAM=('SIEF_ELGA', 'VARI_ELGA'),
"""
    ),
    25: _(
        """
  Attention : %(i1)d mailles de type %(k1)s n'ont pas été projetées car la famille
  de points de Gauss sur le champ en question a une maille support 1D.

  Il s'agit certainement d'éléments de joint 2D.
"""
    ),
    26: _(
        """
On ne trouve pas le champ de nom %(k1)s pour le numéro d'ordre %(i1)d dans
le concept résultat fourni.
"""
    ),
    27: _(
        """
On ne peut pas transférer le champ sur le nouveau support car on n'a pas le modèle.
"""
    ),
    31: _(
        """
Erreur utilisateur dans la commande AFFE_CARA_ELEM :
  On a affecté un excentrement non nul (mot clé COQUE / EXCENTREMENT)
  sur un élément qui ne sait pas traiter l'excentrement (maille %(k1)s).
"""
    ),
    32: _(
        """
Erreur utilisateur :
  On cherche à déterminer de quel coté se situe une maille
  par rapport à une maille de "peau".
  Mais cette maille semble dégénérée (sans épaisseur)
  ce qui empêche de répondre à la question.
  La maille "coupable" est : %(k1)s

Risques et conseils :
  Si le problème concerne le mots clé ORIE_PEAU de la commande MODI_MAILLAGE,
  vous pouvez utiliser le mot clé GROUP_MA_INTERNE pour éviter de traiter
  certaines mailles.
"""
    ),
    33: _(
        """
Erreur Utilisateur :
 Pour le modèle  %(k1)s  on ne peut pas visualiser ensemble plusieurs champs ELGA (%(k2)s,  ...)
 car les familles de points de Gauss sont différentes
"""
    ),
    35: _(
        """
Erreur Utilisateur :
 Aucun élément du modèle n'est visualisable avec ECLA_PG
"""
    ),
    36: _(
        """
 On ne trouve aucun point de Gauss
"""
    ),
    39: _(
        """
 les seuls champs autorisés pour ECLA_PG sont les champs réels.
"""
    ),
    40: _(
        """
Erreur :
 Après avoir retiré tous les éléments à sous-points du champ %(k1)s (grandeur: %(k2)s), celui-ci est vide.
"""
    ),
    41: _(
        """
 les seuls champs autorises sont ELGA.
"""
    ),
    42: _(
        """
%(k1)s  n'a pas le nombre de points de Gauss déclaré dans la routine
Ce message est un message d'erreur développeur.
Contactez le support technique.
"""
    ),
    43: _(
        """
 nombre de noeuds > 27
"""
    ),
    44: _(
        """
   Le modèle n'a pas été trouvé. Le calcul n'est pas possible.
"""
    ),
    45: _(
        """
   Le CARA_ELEM n'a pas été trouvé. Le calcul n'est pas possible.
"""
    ),
    46: _(
        """
  mode ligne  %(k1)s  /= mode colonne  %(k2)s
"""
    ),
    47: _(
        """
  le mode  %(k1)s  de code  %(k2)s  référence le mode  %(k3)s  dont le code :  %(k4)s  > 3
"""
    ),
    48: _(
        """
  pour le mode  %(k1)s  nombre de points  %(k2)s  < argument k :  %(k3)s
"""
    ),
    49: _(
        """
 carte inexistante.
"""
    ),
    50: _(
        """
Erreur utilisateur :
  Le champ %(k1)s n'est pas associé au maillage %(k2)s.
"""
    ),
    51: _(
        """
  La méthode ECLA_PG de PROJ_CHAMP ne peut pas traiter les champs de plus de %(i1)d composantes.
  Le champ à traiter en comporte %(i2)d.
"""
    ),
    52: _(
        """
Pour le champ %(k1)s
Vous utilisez la commande CREA_TABLE sur un champ incomplet issu d'éléments à sous-points :
    poutres multifibres, plaques multicouches, tuyaux , ...
Ce champ a peut-être été crée par POST_CHAMP/EXTR_**

Le champ extrait  a %(i1)4d sous-point(s).
L'élément support a %(i2)4d sous-point(s).

Conseil : Il ne faut pas modifier le nombre de composantes d'un champ pour ensuite
          le réutiliser dans une autre commande.
"""
    ),
    57: _(
        """
  Erreur utilisateur dans la commande CREA_CHAMP :
    On demande l'affectation du champ sur les sous-points mais le champ a déjà des sous-points.

  Conseil :
    Il ne faut pas utiliser le mot clé AFFE_SP dans cette commande.
"""
    ),
    58: _(
        """
Erreur Utilisateur :
  Vous utilisez le mot clé PILOTAGE et votre calcul utilise des variables de
  commande qui dépendent du temps (mot clé AFFE_VARC / EVOL).
  C'est interdit.
  Champ de matériau : %(k1)s
Conseil :
  Lorsque l'on utilise le mot clé PILOTAGE, les variables de commande ne doivent
  pas dépendre du temps.
  Il faut utiliser le mot clé AFFE_VARC / CHAM_GD
"""
    ),
    62: _(
        """
 Erreur utilisateur POST_CHAMP /MIN_MAX_SP :
  Il n'y a rien à calculer car le champ  %(k1)s n'existe pas pour les numéros d'ordre indiqués.
"""
    ),
    65: _(
        """
Erreur d'utilisation :
  -> Le modèle %(k1)s n'a pas d'éléments sachant calculer la rigidité.

  -> Risque & Conseil :
     Ce modèle ne peut pas être utilisé pour faire des calculs.
     Vérifier la définition du modèle (AFFE_MODELE) et assurez-vous que les
     types de mailles du maillage (SEG2, TRIA3, QUAD4, ...) sont compatibles avec votre
     modélisation.
     Exemples d'erreur :
       * affecter une modélisation "3D" sur un maillage formé de facettes.
       * affecter une modélisation qui ne sait pas traiter tous les types de mailles du maillage
         (par exemple 'PLAN_DIAG' en thermique, 'AXIS_SI' en mécanique)
"""
    ),
    66: _(
        """Erreur d'utilisation :
 On ne peut pas utiliser plus de 50 paramètres pour évaluer une fonction.
 Ici, les différents champs du mot-clé CHAM_PARA possèdent au total plus de 50 composantes.
"""
    ),
    67: _(
        """Erreur d'utilisation :
 On ne peut pas filtrer les mailles de type %(k1)s car ce n'est pas un type de maille connu.
"""
    ),
    78: _(
        """
Erreur utilisateur :
  Pour évaluer une fonction, Il ne faut pas fournir plusieurs fois le même paramètre.
  Ici, le paramètre %(k1)s a été fourni plus d'une fois.

Conseil :
  Si les champs utilisés avec le mot clé CHAM_PARA ont été obtenus avec la commande CREA_CHAMP,
  on peut voir leur contenu avec le mot clé INFO=2.
"""
    ),
    79: _(
        """
Erreur dans CREA_RESU :
  Quand on utilise la commande CREA_RESU / EVOL_VARC avec le mot clé AFFE / CHAM_GD et
  que le champ de fonctions est %(k2)s_F (de la géométrie et/ou du temps),
  NOM_CHAM doit être %(k2)s et pas %(k1)s.
"""
    ),
    80: _(
        """
Erreur dans CREA_RESU :
  Quand on utilise la commande CREA_RESU avec le mot clé AFFE / CHAM_GD :
  le champ de fonctions %(k2)s doit avoir le même nombre d'entier codés que le
  champ de réels %(k1)s.
    %(k2)s : %(i2)d entier codé
    %(k1)s : %(i1)d entier codé

Conseil :
  Si nécessaire, il faut demander une évolution du code.
"""
    ),
    81: _(
        """
Erreur utilisateur :
  Calcul de la déformation thermique d'un élément de grille.
  On ne trouve pas de température sur le maille %(k1)s.
"""
    ),
    82: _(
        """
 il faut un MODELE
"""
    ),
    86: _(
        """
 La carte de COMPORTEMENT est absente.
 Votre résultat a peut-être été produit par LIRE_RESU ou CREA_RESU.
 Si votre résultat a été produit par LIRE_RESU, il faut renseigner le mot-clé COMPORTEMENT.
"""
    ),
    88: _(
        """
 L'option %(k1)s  n'est disponible pour aucun des éléments de votre modèle.
 Le calcul d'indicateur d'erreur est donc impossible.
"""
    ),
    89: _(
        """
 Alarme utilisateur :
   Le champ  %(k1)s  n'a pas pu être calculé.

 Risques & conseils :
   * Si le champ est un champ par éléments, c'est que le calcul élémentaire n'est pas disponible
     pour les éléments finis utilisés. Cela peut se produire soit parce que ce
     calcul n'a pas été encore programmé, soit parce que ce calcul n'a pas de sens.
     Par exemple, le champ EFGE_ELNO n'a pas de sens pour les éléments de la modélisation '3D'.
   * Si le champ est un champ aux noeuds (XXXX_NOEU), cela veut dire que le champ XXXX_ELNO
     n'existe pas sur les éléments spécifiés.
     Par exemple, le calcul de SIGM_NOEU sur les éléments de bord est impossible.

"""
    ),
    90: _(
        """
Erreur dans CREA_RESU :
  Quand on utilise la commande CREA_RESU avec le mot clé AFFE / CHAM_GD et que
  le champ de fonctions est %(k2)s_F (de la géométrie et/ou du temps),
  le champ de réels doit être %(k1)s_R.

Conseil :
  Vérifier la construction du champ de fonctions est %(k2)s_F
"""
    ),
    92: _(
        """
 votre chargement contient plus d'une charge répartie
 le calcul n'est pas possible pour les modèles de poutre.
"""
    ),
}
