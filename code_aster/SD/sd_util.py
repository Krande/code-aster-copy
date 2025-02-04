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

# person_in_charge: mathieu.courtois at edf.fr

"""
   Utilitaires pour la vérification des SD
"""

import copy

import numpy

from ..Objects import PhysicalQuantityManager

#  1) Utilitaires pour vérifier certaines propriétés.
#     Ces utilitaires ne provoquent pas d'arret mais écrivent des messages dans un "checker"
#  -----------------------------------------------------------------------

#   1.1 Utilitaires pour des scalaires :
#   ------------------------------------
def sdu_assert(ojb, checker, bool, comment=""):
    """Vérifie que le booléen (bool) est vrai"""
    if not bool:
        checker.err(ojb, "condition non respectée :  (%s)" % (comment,))


def sdu_compare(ojb, checker, val1, comp, val2, comment=""):
    """Vérifie que la relation de comparaison entre val1 et val2 est respectée :
    comp= '==' / '!=' / '>=' / '>' / '<=' / '<'"""
    comp = comp.strip()
    ok = 0
    if comp == "==":
        ok = val1 == val2
        comp1 = "n'est pas égale au"
    elif comp == "!=":
        ok = val1 != val2
        comp1 = "est égale au"
    elif comp == ">=":
        ok = val1 >= val2
        comp1 = "est inférieure strictement au"
    elif comp == "<=":
        ok = val1 <= val2
        comp1 = "est supérieure strictement au"
    elif comp == ">":
        ok = val1 > val2
        comp1 = "est inférieure ou égale au"
    elif comp == "<":
        ok = val1 < val2
        comp1 = "est supérieure ou égale au"
    else:
        sdu_assert(ojb, checker, 0, "sdu_compare: opérateur de comparaison interdit: " + comp)

    if not ok:
        checker.err(
            ojb,
            "condition non respectée pour le test suivant : longueur séquence (%s) %s nombre d'éléments différents dans la séquence (%s) (%s)"
            % (val1, comp1, val2, comment),
        )


#   1.2 Utilitaires pour des séquences :
#   ------------------------------------
def sdu_tous_differents(ojb, checker, sequence=None, comment=""):
    """Vérifie que les éléments de la séquence sont tous différents.
    Si l'argument sequence est None, on prend l'ensemble de l'ojb."""
    if sequence:
        seq = sequence
    else:
        seq = ojb.get()
    sdu_compare(
        ojb,
        checker,
        len(seq),
        "==",
        len(set(seq)),
        comment="Tous les éléments de la séquence devraient être différents, mais ils ne le sont pas"
        + comment,
    )


def sdu_tous_non_blancs(ojb, checker, sequence=None, comment=""):
    """Vérifie que les éléments (chaines) de la séquence sont tous "non blancs".
    Si l'argument sequence est None, on prend l'ensemble de l'ojb."""
    if sequence:
        seq = sequence
    else:
        seq = ojb.get()
    for elem in seq:
        assert len(elem.strip()) > 0, (seq, self, 'tous "non blancs" ' + comment)


def sdu_tous_compris(ojb, checker, sequence=None, vmin=None, vmax=None, comment=""):
    """Vérifie que toutes les valeurs de la sequence sont comprises entre vmin et vmax
    Les bornes vmin et vmax sont autorisées
    Si l'argument sequence est None, on prend l'ensemble de l'ojb."""
    assert (not vmin is None) or (
        not vmax is None
    ), "Il faut fournir au moins une des valeurs vmin ou vmax"
    if sequence:
        seq = sequence
    else:
        seq = ojb.get()
    ier = 0
    for v in seq:
        if vmin and v < vmin:
            ier = 1
        if vmax and v > vmax:
            ier = 1
    if ier == 1:
        checker.err(
            ojb, "L'objet doit contenir des valeurs dans l'intervalle : [%s, %s] " % (vmin, vmax)
        )


def sdu_monotone(seqini):
    """vérifie qu'une séquence est triée par ordre croissant (ou décroissant)
    retourne :
       3 : ni croissant ni décroissant  (désordre)
       1 : croissant
      -1 : décroissant
       0 : constant"""
    if len(seqini) < 2:
        return 0
    tv = numpy.array(seqini)
    diff = tv[1:] - tv[:-1]
    croiss = min(diff) >= 0
    decroiss = max(diff) <= 0
    if croiss and decroiss:
        return 0
    elif croiss and not decroiss:
        return 1
    elif not croiss and decroiss:
        return -1
    else:
        return 3


#  2) Utilitaires de questionnement :
#  -----------------------------------------------------------------------


def sdu_verif_nom_gd(nomgd):
    """vérifie que nomgd est bien un nom de grandeur"""
    if not PhysicalQuantityManager.hasQuantityOfName(nomgd.strip()):
        checker.err(ojb, "condition non respectée : " + nomgd + " n'est pas un nom de grandeur.")


def sdu_nom_gd(numgd):
    """retourne le nom de la grandeur de numéro (numgd)"""
    assert numgd > 0 and numgd < 1000, numgd
    return PhysicalQuantityManager.getPhysicalQuantityName(numgd).strip()


def sdu_licmp_gd(numgd):
    """retourne la liste des cmps de la grandeur de numéro (numgd)"""
    assert numgd > 0 and numgd < 1000, numgd
    return PhysicalQuantityManager.getComponentNames(numgd)


def sdu_nb_ec(numgd):
    """retourne le nombre d'entiers codés pour décrire les composantes de la grandeur (numgd)"""
    assert numgd > 0 and numgd < 1000, numgd
    return PhysicalQuantityManager.getNumberOfEncodedInteger(numgd)


#  3) Utilitaires pour la vérification de l'existence des objets :
#  -----------------------------------------------------------------------


def sdu_ensemble(lojb):
    """vérifie que les objets JEVEUX de lojb existent simultanément :"""
    assert len(lojb) > 1, lojb
    lexi = []
    for obj1 in lojb:
        lexi.append(obj1.exists)
    for x in lexi[1:]:
        assert x == lexi[0], (lojb, lexi)
