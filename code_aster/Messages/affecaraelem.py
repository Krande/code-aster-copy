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

# person_in_charge: jean-luc.flejou at edf.fr

from ..Utilities import _

cata_msg = {
    # Messages dans OP0019
    2: _(
        """AFFE_CARA_ELEM
Au moins un des mot-clés facteur d'AFFE_CARA_ELEM n'a aucune affection sur des GROUP_MA,
des MAILLES ou sur TOUT='OUI'.

Vérifiez vos données.
"""
    ),
    # Messages dans ace_mass_rep
    10: _(
        """AFFE_CARA_ELEM / <%(k1)s> / occurrence %(i1)d
Une maille dans le groupe <%(k2)s> n'est pas de la bonne topologique.
Vous êtes en 2D, la topologie des mailles doit être 1.
Les mailles de <%(k2)s> doivent être des segments.

Pour information :
 - nom de la maille incriminée                       : %(k3)s
 - dimension topologique de la maille incriminée     : %(i3)d
 - type de la maille incriminée                      : %(k4)s

Vérifiez vos données.
"""
    ),
    11: _(
        """AFFE_CARA_ELEM / <%(k1)s> / occurrence %(i1)d
Une maille dans le groupe  <%(k2)s> n'a pas la bonne dimension topologique.

Les mailles de <%(k2)s> doivent être des segments ou des éléments de surfaces.
Pour information :
 - maille concernée         : %(k3)s
 - sa dimension topologique : %(i2)d
 - son type                 : %(k4)s

Vérifiez vos données.
"""
    ),
    12: _(
        """AFFE_CARA_ELEM / <%(k1)s> / occurrence %(i1)d
La maille %(k3)s est présente au moins 2 fois dans la définition de la surface.

Les mailles ne doivent être présente qu'une seule fois, en cas de doublon :
- La surface totale de répartition ne peut pas être correctement déterminée.
- La contribution de cette maille n'est pas correcte.

Vérifiez vos données.
"""
    ),
    13: _(
        """AFFE_CARA_ELEM / <%(k1)s> / occurrence %(i1)d
La maille %(k3)s présente dans <%(k2)s> n'a pas le bon nombre de noeuds.
Les mailles gérées doivent avoir %(k5)s noeuds

Pour information :
 - maille concernée              : %(k3)s
 - son type                      : %(k4)s
 - nombre de noeuds de la maille : %(i2)d

Vérifiez vos données.
"""
    ),
    14: _(
        """AFFE_CARA_ELEM / <%(k1)s> / occurrence %(i1)d
Au moins un élément appartenant au GROUP_MA_POI1 <%(k2)s> n'est pas connecté aux mailles
surfaciques du GROUP_MA.

Pour pouvoir calculer la contribution de la maille, il faut que tous les noeuds de <%(k2)s> soient
en relation avec un seul noeud d'une maille surfacique du GROUP_MA.

Pour information :
 - maille détectée : %(k3)s

Vérifiez vos données.
"""
    ),
    15: _(
        """AFFE_CARA_ELEM / <%(k1)s> / occurrence %(i1)d
Au moins un noeud appartenant à un élément du GROUP_MA surfacique n'est pas connecté à une maille du
GROUP_MA_POI1.

Pour pouvoir répartir la contribution des mailles surfaciques, il faut que tous les noeuds
soient en relation avec une et une seule maille de GROUP_MA_POI1.

Pour information :
 - noeud concerné : %(k2)s

Vérifiez vos données.
"""
    ),
    16: _(
        """AFFE_CARA_ELEM / <%(k1)s> / occurrence %(i1)d
Des éléments appartenant au GROUP_MA_POI1 <%(k2)s> n'ont pas le bon nombre de noeuds.
Les mailles doivent être de type POI1 et posséder 1 seul noeud.

Pour information :
 - maille détectée : %(k3)s

Vérifiez vos données.
"""
    ),
    17: _(
        """AFFE_CARA_ELEM / <%(k1)s> / occurrence %(i1)d
Les mailles appartenant au GROUP_MA sont surfaciques. Vous ne devez pas utiliser
TYPE=LINEIQUE mais TYPE=SURFACIQUE ou TOTALE.

Vérifiez vos données.
"""
    ),
    18: _(
        """AFFE_CARA_ELEM / <%(k1)s> / occurrence %(i1)d
Les mailles appartenant au GROUP_MA sont linéiques. Vous ne devez pas utiliser
TYPE=SURFACIQUE mais TYPE=LINEIQUE ou TOTALE.

Vérifiez vos données.
"""
    ),
    19: _(
        """AFFE_CARA_ELEM / <%(k1)s> / occurrence %(i1)d
Le noeud <%(k5)s> de la maille <%(k3)s> appartenant au GROUP_MA_POI1 <%(k2)s> est connecté
avec un noeud d'une maille surfacique du GROUP_MA qui est déjà connecté à la maille
<%(k4)s> du GROUP_MA_POI1.

Les mailles de <%(k2)s> ne doivent être connecté qu'a un seul noeud des mailles de surface.
Les noeuds des mailles de surface ne doivent être connecté qu'a un seul noeud des mailles POI1.
La relation doit être bijective.

Vérifiez vos données.
"""
    ),
    20: _(
        """
AFFE_CARA_ELEM. Il y a %(i1)d occurrences du mot clef facteur <%(k1)s>.
Entre ces différentes occurrences les GROUP_MA_POI1 ont %(i2)d mailles en communs.
La règle de surcharge est donc appliquée %(i2)d fois.
La liste des mailles surchargées est affichée avec INFO=2.
"""
    ),
    21: _(
        """
AFFE_CARA_ELEM / MASS_REP. Occurrence %(i1)d.

Le nombre des mailles POI1 affectées lors des différentes occurrences de MASS_REP dépasse
le nombre de DIS_T/DIS_TR présent dans le modèle : %(i3)d

Pour information :
- le nombre de maille de type POI1 précédemment affecté par MASS_REP est de %(i4)d
- le nombre de maille de type POI1 dans %(k1)s est %(i2)d

Conseils : Vérifiez que lors de votre AFFE_MODELE, vous n'avez pas oublié d'affecter des
           DIS_T ou DIS_TR.
"""
    ),
    22: _(
        """AFFE_CARA_ELEM / <%(k1)s> / occurrence %(i1)d
Des éléments appartenant au GROUP_MA_POI1 <%(k2)s> ne sont pas dans le modèle.

Pour information :
 - maille détectée : %(k3)s

Vérifiez vos données.
"""
    ),
    23: _(
        """AFFE_CARA_ELEM / <%(k1)s> / occurrence %(i1)d
Une maille dans le groupe <%(k2)s> n'est pas de la bonne topologique.
Vous êtes en 3D, la topologie des mailles dans le groupe doit être identique.

Les mailles de <%(k2)s> doivent être des segments ou des éléments de surfaces.

Pour information :
 - nom de la maille incriminée                       : %(k3)s
 - dimension topologique de la maille incriminée     : %(i3)d
 - type de la maille incriminée                      : %(k4)s

La topologie du groupe est déterminée par la 1ère maille du groupe.
 - nom de la 1ère maille du groupe                   : %(k5)s
 - dimension topologique du groupe                   : %(i2)d
 - type de la 1ère maille du groupe                  : %(k6)s

Vérifiez vos données.
"""
    ),
    98: _(
        """AFFE_CARA_ELEM / RIGI_PARASOL
Il y a %(i1)d occurrences de RIGI_PARASOL.
Pour chaque occurrence, les discrets sont définis par GROUP_MA_[SEG2|POI1] et la
surface est définie par GROUP_MA.

Lorsqu'il y a plusieurs occurrences, si un même discret est présent dans plusieurs occurrences
de GROUP_MA_[SEG2|POI1] les raideurs qui lui sont associées seront CUMULÉES.

Pour éviter cela et ne pas cumuler les raideurs, le discret ne doit être présent que dans une
seule occurrence, dans un seul GROUP_MA_[SEG2|POI1]
"""
    ),
    99: _(
        """AFFE_CARA_ELEM / RIGI_PARASOL
Pour l'occurrence %(i1)d de RIGI_PARASOL vous donnez %(i2)d groupes sous GROUP_MA.
Les discrets sont définis par GROUP_MA_[SEG2|POI1], la surface est définie par GROUP_MA.

Si des mailles sont présentes dans plusieurs groupes donnés sous GROUP_MA, elles seront
décomptées plusieurs fois. Par conséquent :
- la surface totale sera la somme de la surface de chaque groupe de mailles. La surface
  associées aux mailles présentes plusieurs fois sera donc prise en compte plusieurs fois.
- la raideur associée aux discrets connectés à ces mailles, sera pondérée par la surface
  de ces mailles, donc prise en compte plusieurs fois.

Pour éviter cela et ne pas pondérer les discrets associés par une surface comptée plusieurs
fois, il faut que l’intersection des groupes de mailles soit vide.
"""
    ),
}
