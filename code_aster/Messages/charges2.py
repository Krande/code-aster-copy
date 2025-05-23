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
        """
Le modèle sur lequel est défini le chargement n'est pas cohérent avec le modèle existant.
"""
    ),
    3: _(
        """
 La fonction %(k1)s fournie à l'opérande FLUN de FLUX_NL doit avoir le paramètre TEMP.
"""
    ),
    4: _(
        """
  -> Le modèle contient un mélange d'éléments finis 2D (plan Oxy) et 3D

  -> Risque & Conseil :
     Sur ce genre de modèle, on ne sait pas déterminer s'il est 2D ou 3D.
     Certains chargements ne seront pas possibles.
"""
    ),
    5: _(
        """
 Le chargement de type %(k1)s est interdit en 2D.
"""
    ),
    6: _(
        """
Erreur d'utilisation :
 Le modèle contient un mélange d'éléments 2D (vivant dans le plan Oxy) et 3D.
 Le code n'a pas prévu ce cas de figure dans l'application du chargement demandé.

Risques et conseils :
 Il faut peut être émettre une demande d'évolution pour pouvoir traiter ce problème.
"""
    ),
    7: _(
        """
Le modèle est de dimension %(i1)d . ARETE_IMPO s'applique sur des arêtes d'éléments 3D,
donc un modèle de dimension 3. Pour les arêtes d'éléments 2D utiliser FACE_IMPO.
"""
    ),
    8: _(
        """
Pour un chargement de type LIAISON_COQUE, il faut avoir autant de noeuds dans les deux listes.
"""
    ),
    9: _(
        """
Échec de l'appariement des deux listes de noeuds pour le chargement de type LIAISON_COQUE.

 Conseils :
   - Si la distance entre les deux surfaces à apparier est grande devant leurs dimensions, précisez l'isométrie qui permet de les superposer par l'intermédiaire des mots-clés CENTRE, ANGL_NAUT et TRAN.
"""
    ),
    10: _(
        """
Pour le chargement de type LIAISON_CHAMNO, on doit utiliser le mot clé CHAM_NO pour donner le CHAM_NO dont les composantes seront les coefficients de la relation linéaire.
"""
    ),
    11: _(
        """
Pour le chargement de type LIAISON_CHAMNO, il faut que le CHAM_NO dont les termes servent de coefficients à la relation linéaire à écrire ait été défini.
"""
    ),
    12: _(
        """
Pour le chargement de type LIAISON_CHAMNO, tous les coefficients donnés par le mot-clef CHAM_NO sont nuls.
"""
    ),
    17: _(
        """
 Pour le chargement courant, la liste des noeuds donnée est réduite à un seul terme et l'on ne fait aucun traitement.
"""
    ),
    23: _(
        """
 Il est impossible de calculer la tangente de la maille %(k1)s, des noeuds doivent être confondus.
"""
    ),
    24: _(
        """
 Il est impossible de calculer la normale de la maille %(k1)s, des noeuds doivent être confondus.
"""
    ),
    25: _(
        """
 Il n'est pas possible de calculer la normale d'un segment en 3d.
 Il ne doit pas y avoir de segment dans le groupe sur lequel vous appliquez la condition limite.
"""
    ),
    26: _(
        """
 Il est impossible de calculer la normale de la maille %(k1)s .
 Il y a un problème dans la définition de vos mailles: des arêtes doivent être confondues.
"""
    ),
    29: _(
        """
 L'angle formé par le vecteur normal courant à une face et le vecteur normal moyen, au noeud %(k1)s, est supérieur a 10 degrés et vaut %(k2)s degrés.
"""
    ),
    30: _(
        """
Erreur d'utilisation :
 La norme du vecteur normal (moyenne des normales des éléments concourants) est presque nulle.
 Les facettes concourantes au noeud  %(k1)s ne définissent pas une normale fiable.
 Il y a un problème dans la définition des mailles de bord .

Suggestion :
 Pensez à réorienter les mailles de bord avec l'opérateur MODI_MAILLAGE.
"""
    ),
    31: _(
        """
 Erreur utilisateur:
    On cherche à imposer une condition aux limites sur le ddl %(k1)s du noeud %(k2)s.
    Mais ce noeud ne porte pas ce ddl.

    Conseils :
     - vérifiez le modèle et les conditions aux limites :
        - le noeud incriminé fait-il partie du modèle ?
        - le noeud porte-t-il le ddl que l'on cherche à contraindre ?
"""
    ),
    32: _(
        """
 Il y a un problème sur une relation linéaire car les coefficients sont trop petits.
 Ce problème survient dans le cas de liaisons linéaires automatiques. Vérifiez les éventuelles
 alarmes émises précédemment dans des mots clefs comme LIAISON_MAILLE par exemple.
"""
    ),
    33: _(
        """
 Le noeud <%(k1)s> ne fait pas partie du modèle.
"""
    ),
    34: _(
        """
 Les relations suivantes sont redondantes et donc supprimées en appliquant le principe de
surcharge.
"""
    ),
    35: _(
        """
 Liste des noeuds en vis-à-vis de l'occurrence %(i1)d de LIAISON_GROUP :"""
    ),
    36: _("""    <%(k1)s> en  vis-à-vis de <%(k2)s>"""),
    37: _(
        """
Le noeud de nom <%(k1)s> est connecté à plus d'une maille.
Vous n'avez pas défini de repère local.
Il est donc impossible de définir le repère local pour appliquer DDL_POUTRE.
Revoyez la définition de votre repère local.
"""
    ),
    38: _(
        """
Le noeud de nom <%(k1)s> est connecté à plus d'une maille.
Il y a plusieurs repères locaux possibles, pour la définition de l'axe X.
Il faut définir la maille support donnant l'axe X.
Il est donc impossible de définir le repère local pour appliquer DDL_POUTRE.
Revoyez la définition de votre repère local.
"""
    ),
    39: _(
        """
Le repère local que vous avez défini ne contient pas de maille attachée au noeud <%(k1)s>.
Il est donc impossible de définir le repère local pour appliquer DDL_POUTRE.
Revoyez la définition de votre repère local.
"""
    ),
    40: _(
        """
La maille <%(k1)s>, attaché au noeud <%(k2)s> n'est pas de type "SEG".
Il est donc impossible de définir le repère local pour appliquer DDL_POUTRE.
"""
    ),
    41: _(
        """
La maille <%(k1)s>, attaché au noeud <%(k2)s> est de longueur nulle.
Il est donc impossible de définir le repère local pour appliquer DDL_POUTRE.
"""
    ),
    42: _(
        """
 Pour définir le mot-clef %(k1)s, on doit utiliser %(i1)d valeurs car nous sommes dans le
cas d'un modèle à %(i2)d dimensions.
"""
    ),
    45: _(
        """
Erreur utilisateur :
  Vous voulez contraindre le ddl %(k1)s sur un ensemble de noeuds,
  Mais ce ddl n'existe sur aucun de ces noeuds.
"""
    ),
    46: _(
        """
Erreur utilisateur :
  Vous voulez imposer une force selon la composante %(k1)s sur un ensemble de noeuds,
  mais cette composante de force n'existe sur aucun de ces noeuds.
"""
    ),
    48: _(
        """
 Le concept CABLE_BP de nom %(k1)s ne contient pas de relations linéaires. L'option RELA_CINE est donc inutilisable.
"""
    ),
    50: _(
        """
 On ne trouve pas de noeud assez près du noeud %(k1)s .
"""
    ),
    51: _(
        """
 Il y a un conflit dans les vis-à-vis des noeuds.
 Le noeud  %(k1)s est à la fois le vis-à-vis du noeud %(k2)s et du noeud %(k3)s.
"""
    ),
    52: _(
        """
 Il y a un conflit dans les vis-à-vis des noeuds.
 Le noeud  %(k1)s n'est l'image d'aucun noeud par la correspondance inverse.
"""
    ),
    53: _(
        """
 Vous avez donné une direction nulle pour l'axe de rotation.
"""
    ),
    61: _(
        """
 Le type d'onde S est interdit en 3D pour le chargement ONDE_PLANE, précisez SV ou SH.
"""
    ),
    62: _(
        """
 Les types d'onde SV et SH sont interdits en 2D pour le chargement ONDE_PLANE, on a pris le type S.
"""
    ),
    63: _(
        """
 Vous ne pouvez pas bloquer le déplacement tangent sur des faces d'éléments 3D.
 Utiliser DDL_IMPO ou LIAISON_DDL.
"""
    ),
    64: _(
        """
  Vous définissez une charge avec %(k1)s sur un modèle de type %(k2)s, ce n'est pas possible.
"""
    ),
    65: _(
        """
Le vecteur définissant l'axe de rotation a une composante non nulle suivant Ox ou Oz,
ce qui induit un chargement non axisymétrique. Avec une modélisation AXIS ou AXIS_FOURIER,
l'axe de rotation doit être dirigé suivant Oy.
"""
    ),
    66: _(
        """
Les coordonnées du centre de rotation ont au moins une composante non nulle, ce qui induit
un chargement non axisymétrique. Avec une modélisation AXIS ou AXIS_FOURIER,
le centre de rotation doit être confondu avec l'origine.
"""
    ),
    67: _(
        """
Le vecteur définissant l'axe de rotation a une composante non nulle suivant Ox ou Oy,
ce qui induit des forces centrifuges hors plan. Avec une modélisation C_PLAN ou D_PLAN,
l'axe de rotation doit être dirigé suivant Oz.
"""
    ),
    68: _(
        """
Vous ne pouvez bloquer DRNOR qu'en dimension 3.
"""
    ),
    69: _(
        """
L'utilisation de DRNOR n'est pas compatible avec la présence d'une modélisation X-FEM.
"""
    ),
    82: _(
        """
 Il faut au moins deux noeuds pour LIAISON_UNIF.
"""
    ),
    85: _(
        """
 La maille de nom %(k1)s n'est pas de bon ordre d'interpolation.
 Elle ne sera pas affectée par %(k2)s .
"""
    ),
    86: _(
        """
 La maille de nom %(k1)s n'est pas de type linéique (segments).
 Elle ne sera pas affectée par %(k2)s .
"""
    ),
    87: _(
        """
 La maille de nom %(k1)s n'est pas de type surfacique (triangles ou quadrangles).
 Elle ne sera pas affectée par %(k2)s .
"""
    ),
    88: _(
        """
 La maille de nom %(k1)s n'est pas de type volumique.
 Elle ne sera pas affectée par %(k2)s .
"""
    ),
    89: _(
        """
Erreur dans les mailles du mot-clé facteur %(k1)s.
Aucune maille n'est du bon type. Elles sont toutes ignorées.
"""
    ),
    92: _(
        """
Un noeud dans DDL_IMPO ne porte aucun DDL de type %(k1)s.
"""
    ),
    97: _(
        """
 Tous les coefficients de la relation linéaire sont strictement nuls.
 Cette erreur peut survenir si le maillage est incorrect (par exemple des noeuds confondus) ou si
vous affectez des coefficients nuls.
"""
    ),
    99: _(
        """
Problème :
  Une relation linéaire entre ddls a un second membre de type "fonction".
  On ne peut pas la normaliser (afin que son plus grand coefficient soit 1.) car on ne
  sait pas "diviser" une fonction par un réel.

  Le plus grand coefficient de la relation est très différent de 1.  (<1.d-3 ou > 1.d3).
  Cette équation (non normalisée) peut conduire à des difficultés numériques lors de
  la résolution des systèmes linéaires.

Conseil :
  Utilisez le solveur MUMPS afin de contrôler la qualité de la résolutions des systèmes linéaires.
"""
    ),
}
