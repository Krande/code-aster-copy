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
        """Pour affecter une liste de modélisations, il faut qu'elles soient de même dimension topologique."""
    ),
    2: _("""La maille %(k1)s de type %(k2)s n'a pas pu être affectée."""),
    3: _(
        """La modélisation choisie est axisymétrique ou plane. Il faut que le maillage soit dans le plan z=0 et ce n'est pas le cas."""
    ),
    4: _(
        """Sur les %(i1)d mailles du maillage %(k1)s, on a demandé l'affectation de %(i2)d, on a pu en affecter %(i3)d."""
    ),
    5: _(
        """Votre modèle doit contenir au moins un élément fini car il n'est pas possible de n'avoir que des macro-éléments si le maillage est lu par LIRE_MAILLAGE. Si vous voulez définir entièrement un modèle avec des macro-éléments, il faut définir ces derniers avec DEFI_MAILLAGE."""
    ),
    6: _("""Aucune maille du maillage %(k1)s n'a été affectée par des éléments finis."""),
    7: _(
        """Attention l'élément HEXA8 en 3D_SI ne fonctionne correctement que sur les parallélépipèdes. Sur les éléments quelconques on peut obtenir des résultats faux."""
    ),
    8: _("""Liste des noeuds affectés pour la modélisation:"""),
    9: _("""Liste des mailles affectées pour la modélisation:"""),
    10: _("""Le modèle contient un mélange d'éléments HHO et non HHO. Ce n'est pas possible."""),
    11: _("""On ne peut utiliser qu'un seul phénomène."""),
    12: _("""La formulation %(k1)s n'est pas disponible pour les modélisations choisies."""),
    13: _(
        """Les méthodes de partition MAIL_CONTIGU et MAIL_DISPERSE sont interdites avec des
modélisations nécessitant des informations sur les mailles voisines (3D_JOINT, 3D_INTERFACE, ...)"""
    ),
    14: _(
        """Le modèle contient un mélange d'éléments finis 2D (plan Oxy) et 3D. Sur ce genre de modèle, on ne sait pas déterminer s'il est 2D ou 3D.
Parfois, cela empêche de faire le "bon choix". Cette alarme risque de se transformer en erreur fatale dans la suite du calcul.
"""
    ),
    15: {
        "message": _(
            """Attention l'élément PYRAM13 donne des résultats imprécis à cause de l'intégration."""
        ),
        "flags": "DECORATED",
    },
    20: _(""" Modélisation     Formulation      Type maille  Élément fini     Nombre"""),
    21: _(""" %(k1)-16s %(k2)-16s %(k3)-12s %(k4)-16s %(i1)d"""),
    38: _("""%(k1)-8s %(k2)-8s %(k3)-8s %(k4)-8s %(k5)-8s %(k6)-8s %(k7)-8s %(k8)-8s"""),
    53: _(
        """Le maillage est 3D (tous les noeuds ne sont pas dans le même plan Z = constante), mais les éléments du modèle sont de dimension 2.
Si les facettes supportant les éléments ne sont pas dans un plan Z = constante, les résultats seront faux.
Assurez-vous de la cohérence entre les mailles à affecter et la modélisation souhaitée dans la commande AFFE_MODELE.
"""
    ),
    54: _("""Il est interdit de mélanger des éléments discrets 2D et 3D dans le même modèle."""),
    58: _(
        """Les éléments du modèle sont de dimension 2. Mais les noeuds du maillage sont un même plan Z qui n'est pas à l'origine.
Il est d'usage d'utiliser un maillage dans le plan Z passant par l'origine pour les modélisations planes ou axisymétriques
"""
    ),
    63: _(
        """La maille %(k1)s porte un élément fini de bord, mais elle ne borde aucun élément ayant une "rigidité".

  -> Risque & Conseil :
     Cela peut entraîner des problèmes de "pivot nul" lors de la résolution.
     Si la résolution des systèmes linéaires ne pose pas de problèmes, vous
     pouvez ignorer ce message.
     Sinon, vérifier la définition du modèle (AFFE_MODELE) en évitant l'utilisation
     de l'opérande TOUT='OUI'.
"""
    ),
    64: _(
        """Le modèle n'a pas d'éléments sachant calculer la rigidité.

  -> Risque & Conseil :
     Ce modèle ne pourra donc pas (en général) être utilisé pour faire des calculs.
     Vérifier la définition du modèle (AFFE_MODELE) et assurez-vous que les
     types de mailles du maillage (SEG2, TRIA3, QUAD4, ...) sont compatibles avec votre
     modélisation.
     Exemples d'erreur :
       * affecter une modélisation "3D" sur un maillage formé de facettes.
       * affecter une modélisation qui ne sait pas traiter tous les types de mailles du maillage
         (par exemple 'PLAN_DIAG' en thermique, 'AXIS_SI' en mécanique)
"""
    ),
    70: _(
        """
 Possible erreur utilisateur dans la commande AFFE_MODELE :
   Un problème a été détecté lors de l'affectation des éléments finis.
   Pour l'occurrence AFFE de numéro %(i1)d, certaines mailles de même dimension topologique
   que la modélisation (ici dimension = %(i3)d) n'ont pas pu être affectées et
   ne l'étaient pas déjà.

   Cela veut dire que la modélisation que l'on cherche à affecter
   ne supporte pas tous les types des mailles présentes dans le maillage et
   non encore affectées d'une modélisation.

   Le nombre de mailles non encore affectées que l'on n'a pas pu affecter
   (lors de cette occurrence de AFFE) est :  %(i2)d

 Risques & conseils :
   * Comme certaines mailles n'ont peut-être pas été affectées, il y a un risque
     de résultats faux (présence de "trous" dans la modélisation).
     Pour connaître les mailles non affectées (à la fin de l'opérateur), on peut utiliser INFO=2.
   * Ce problème est fréquent quand on souhaite une modélisation "sous intégrée"
     (par exemple AXIS_SI). Pour l'éviter, il faut donner, dans une occurrence précédent l'occurrence
     en question, une modélisation de "substitution" pour les mailles qui n'existent pas dans
     la modélisation désirée (ici 'AXIS_SI').
     On fera par exemple :
        MO=AFFE_MODELE( MAILLAGE=MA,  INFO=2, AFFE=(
                        _F(TOUT='OUI', PHENOMENE='MECANIQUE', MODELISATION=('AXIS')),
                        _F(TOUT='OUI', PHENOMENE='MECANIQUE', MODELISATION=('AXIS_SI')),
                        ))
"""
    ),
    80: _(
        """On a détecté %(i1)d mailles qui ont une modélisation de type interaction fluide-structure dont on va vérifier l'orientation des normales."""
    ),
    81: _("""Les mailles concernées sont de dimension topologique 1 (segments)."""),
    82: _("""Les mailles concernées sont de dimension topologique 2 (triangles ou quadrangles)."""),
    99: _(
        """Le modèle est basé sur le maillage parallèle %(k2)s. La méthode de distribution des éléments doit être CENTRALISE (et non pas %(k1)s)"""
    ),
}
