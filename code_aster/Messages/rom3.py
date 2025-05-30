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

# person_in_charge: mickael.abbas at edf.fr

from ..Utilities import _

cata_msg = {
    29: _("""Le système global résultant utilise %(i1)d équations."""),
    30: _("""Le système est une combinaison de %(i1)d matrices."""),
    31: _("""Matrice réelle %(i1)d de nom %(k1)s."""),
    32: _("""Matrice complexe %(i1)d de nom %(k1)s."""),
    33: _("""Le vecteur second membre de nom %(k1)s est de type réel."""),
    34: _("""Le vecteur second membre de nom %(k1)s est de type complexe ."""),
    35: _("""Le système global résultant est de type réel."""),
    36: _("""Le système global résultant est de type complexe."""),
    37: _("""Paramètres du système."""),
    39: _("""Paramètres du problème multi-paramétrique."""),
    41: _("""Coefficient scalaire complexe: (%(r1)19.12e,%(r2)19.12e)."""),
    42: _("""Coefficient scalaire réel: %(r1)19.12e."""),
    43: _("""Coefficient fonction complexe: %(k1)s."""),
    44: _("""Coefficient fonction réel: %(k1)s."""),
    45: _("""Paramètres pour la variation des coefficients."""),
    46: _("""Type de définition pour la variation des coefficients: %(k1)s."""),
    47: _("""Nombre de variations des coefficients: %(i1)d par mode calculé."""),
    48: _("""Nombre de paramètres définis pour la variation des coefficients:  %(i1)d."""),
    50: _("""Nom du paramètre: %(k1)s."""),
    51: _("""Nombre total de valeurs pour le paramètre: %(i1)d."""),
    52: _("""Valeur initiale du paramètre : %(r1)19.12e."""),
    60: _("""Le second membre est une combinaison de %(i1)d vecteurs."""),
}
