# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

    1 : _("""La relation linéaire destinée à éliminer un des noeuds esclaves est une tautologie car la maille maître en vis à vis de ce noeud possède ce même noeud dans sa connectivité. On ne l'écrit donc pas."""),

    2 : _("""La composante normale (DNOR) doit être la seule des composantes de la liste."""),

    3 : _("""On ne trouve pas de noeud assez près du noeud %(k1)s."""),

    4 : _("""Un des éléments esclave n'est pas du bon type.
 Pour le calcul de la normale, il faut que les éléments soient de la bonne dimension: des segments en 2D ou des faces en 3D."""),

    6 : _(""" Le modèle contient un mélange d'éléments 2D (vivant dans le plan Oxy) et 3D.
 Il n'est pas possible de réaliser une liaison dans cette configuration."""),

    7 : _("""Les erreurs d'appariement précédentes sont fatales."""),

    9 : _("""Il est interdit d'avoir deux mailles de type POI1 simultanément sur les deux surfaces en vis-à-vis."""),

   48 : _("""Il n'y a aucun noeud esclave à lier pour l'occurrence %(i1)d du mot clé LIAISON_MAIL.
   Peut-être que tous les noeuds esclaves ont déjà été éliminés dans des occurrences précédentes."""),

   49 : _(""" Pour le calcul de la normale sur le côté esclave, il faut donner des éléments de facette."""),

   77 : _("""Il y a un conflit dans les vis-à-vis des noeuds. Un noeud est apparié deux fois.
 Conseils :
   - Si la distance entre les deux surfaces à apparier est grande devant leurs dimensions, précisez l'isométrie qui permet de les superposer par l'intermédiaire des mots-clés CENTRE, ANGL_NAUT et TRAN.
   - Si les maillages sont incompatibles, utilisez plutôt le chargement LIAISON_MAIL.
"""),

   88 : _("""Il y a un conflit dans les vis-à-vis des noeuds. Certains noeuds ne sont pas appariés."""),

}
