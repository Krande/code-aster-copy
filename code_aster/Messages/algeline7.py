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
        """                 Participation du mode : %(i1)2d --> %(r1)12.5E
"""
    ),
    2: _(
        """ =======================================================================
                Calcul modal par %(k1)s

"""
    ),
    3: _(
        """
numéro    itération      erreur              valeur propre
"""
    ),
    4: _(
        """%(i1)4d        %(i2)4d       %(r1)10.3E      ( %(r2)9.2E, %(r3)9.2E )
"""
    ),
    5: _(
        """ Normalisation des modes : %(k1)s

"""
    ),
    6: _(
        """ La bande de fréquence est vide.
"""
    ),
    7: _(
        """
 NUME_ORDRE                    norme
"""
    ),
    8: _(
        """ %(i1)5d        %(k1)24s
"""
    ),
    9: _(
        """
 NUME_ORDRE             ancienne norme              nouvelle norme
"""
    ),
    10: _(
        """ %(i1)5d        %(k1)24s    %(k2)24s
"""
    ),
    11: _(
        """ On saute la valeur propre numéro %(i1)d .
"""
    ),
    12: _(
        """ ALPHA = %(r1)12.5E
 BETA = %(r2)12.5E
"""
    ),
    13: _(
        """ Erreur directe LAPACK %(r1)12.5E
"""
    ),
    14: _(
        """ Elle correspond soit à un Lagrange soit à un DDL physique bloqué.
"""
    ),
    15: _(
        """ Fréquence     = %(r1)12.5E
 Amortissement = %(r2)12.5E
"""
    ),
    16: _(
        """ LAMBDA = %(r1)12.5E
"""
    ),
    17: _(
        """  Le nombre total de DDL est       : %(i1)10d
  Le nombre de DDL de Lagrange est : %(i2)10d
  Le nombre de DDL actifs est      : %(i3)10d
"""
    ),
    18: _(
        """  Le nombre total de DDL est               : %(i1)10d
  Le nombre de DDL de Lagrange est         : %(i2)10d
  Le nombre de DDL bloqués cinématiquement : %(i3)10d
  Le nombre de DDL actifs est              : %(i4)10d
"""
    ),
    19: _(
        """ Nombre de valeurs propres : %(i1)d
"""
    ),
    21: _(
        """La combinaison linéaire des vecteurs est impossible car les numérotations sont différentes."""
    ),
    22: _(
        """

       ===============================================
       =                                             =
       =          Opérateur CALC_ERC_DYN             =
       =                                             =
       =        Résolution d'un problème de          =
       =      minimisation d'une fonctionnelle       =
       =    d'erreur en relation de comportement     =
       =           en dynamique linéaire.            =
       =                                             =
       =        Formulation fréquentielle.           =
       =                                             =
       ===============================================

 Nombre de fréquences demandées    = %(i1)d

 Ordres de grandeur des matrices du problème:
 (informations combinant degrés de liberté physiques et de Lagrange)
 (informations des termes non nuls selon un stockage de la partie
  triangulaire supérieure des matrices)
 -------------------------------------------------------------------

 Dimension des matrices de masse (M) et raideur (K)  = %(i2)d x %(i2)d
 Termes non nuls des matrices de masse et raideur    = %(i3)d

 Dimension de la matrice d'observation (H)           = %(i4)d x %(i2)d
 Termes non nuls de la matrice d'observation         = %(i5)d
 Dimension de la matrice norme (G)                   = %(i4)d x %(i4)d

 Termes non nuls du sous-bloc H^T*G*H                = %(i6)d

 Dimension de la matrice du problème d'ERC           = %(i7)d x %(i7)d
 Termes non nuls de la matrice du problème d'ERC     = %(i8)d

 ===============================================

"""
    ),
}
