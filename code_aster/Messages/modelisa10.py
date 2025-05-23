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
AFFE_CARA_ELEM :
Vous avez activé le mot-clé MODI_METRIQUE='OUI' sur des éléments pour lesquels cette fonctionnalité n'est
pas disponible.

Cette fonctionnalité est disponible avec les modélisations suivantes uniquement :
COQUE_AXIS.
"""
    ),
    2: _(
        """
AFFE_CARA_ELEM :
Vous souhaitez affecter une orientation sur des éléments massifs 3D à l'aide de ORIG_AXE du
mot-clé facteur MASSIF. En 3D, il faut également fournir le mot-clé ANGL_AXE pour définir le repère cylindrique.
"""
    ),
    3: _(
        """
AFFE_CARA_ELEM :
Vous souhaitez affecter une orientation de type repère cylindrique sur des éléments massifs 2D.
En 2D, ANGL_AXE ne doit pas être fourni, seul ORIG_AXE est demandé.
"""
    ),
    4: _(
        """
Les mailles affectées à la modélisation TUYAU ne semblent pas former des lignes continues.
Il y a probablement un problème dans le maillage (superposition d'éléments par exemple).
Pour obtenir le détail des mailles affectées, utilisez INFO=2.
"""
    ),
    5: _(
        """
Le quadrangle de nom %(k1)s est dégénéré : les cotés 1-2 et 1-3 sont colinéaires.
Reprenez votre maillage.
"""
    ),
    7: _(
        """
Il faut au moins un noeud esclave.
"""
    ),
    8: _(
        """
Le groupe d'esclaves %(k1)s est vide.
"""
    ),
    9: _(
        """
Le groupe du noeud maître %(k1)s contient %(i1)d noeuds alors qu'il en faut un seul.
"""
    ),
    10: _(
        """
Arguments incompatibles : il y a %(i1)d degrés de liberté esclaves mais %(i2)d noeuds esclaves.
"""
    ),
    11: _(
        """
Le degré de liberté  %(k1)s est invalide.
"""
    ),
    12: _(
        """
Arguments incompatibles : il y a %(i1)d degrés de liberté esclaves mais %(i2)d coefficients esclaves.
"""
    ),
    13: _(
        """
Arguments incompatibles : il y a %(i1)d coefficients mais %(i2)d noeuds esclaves.
"""
    ),
    14: _(
        """
Pour un spectre de type SPEC_CORR_CONV_3, il faut donner le nom du
MODELE_INTERFACE dans PROJ_SPEC_BASE
"""
    ),
    15: _(
        """
La géométrie de la section utilisée n'est pas prévue par l'opérande SECTION = 'RECTANGLE' de AFFE_CARA_ELEM.
L'un des bords est trop fin.
Utilisez l'opérande SECTION = 'GENERALE'.
"""
    ),
    17: _(
        """
La valeur du mot clé DEFORMATION='%(k1)s' et incompatible avec la modélisation.
"""
    ),
    18: _(
        """
Certains éléments à interpolation linéaires du modèle ne sont pas compatibles avec le modèle de déformation DEFORMATION='%(k1)s'.
"""
    ),
    19: _(
        """
AFFE_CARA_ELEM Pour l'occurrence n° %(i1)d de BARRE le nombre de caractéristiques et de valeurs doivent être identiques.
"""
    ),
}
