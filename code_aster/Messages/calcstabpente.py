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
from ..Utilities import _

cata_msg = {
    1: _(
        """
Les maillages dans le champ des matériaux et le modèle ne sont pas les mêmes. 
"""
    ),
    2: _(
        """
RESI_MAXI ne peut pas dépasser INCR_INIT !
Au cas contraire, RESI_MAXI est considéré égale à INCR_INIT, 
soit RESI_MAXI = %(r1)f. 
"""
    ),
    3: _(
        """
L'algorithme n'a pas convergé alors qu'on atteint le nombre maximum d'itération imposé.
Veuillez examiner les propriétés des matériaux, et augmenter le FS_INIT ou le ITER_MAXI.     
"""
    ),
    4: _(
        """
Le FS_INIT est probablement supérieure au vraie facteur de sécurité. 
Veuillez diminuer le FS_INIT.      

Il est également probable que le calcul non-linéaire a été mal configuré de sorte que le calcul diverge à 
la première itération. Il est donc favorable d'effectuer l'analyse statique non-linéaire via STAT_NON_LINE avant l'usage de la 
macro-commande. 
"""
    ),
    5: _(
        """
La loi de comportement %(k1)s n'est pas autorisée.  
"""
    ),
    6: _(
        """
Le mot clé MAILLE n'est pas autorisé, veuillez remplacer par GROUP_MA.     
"""
    ),
    7: _(
        """
Exactement un comportement dans la zone d'analyse autorisé !!
Non existence ou coexistence des comportements détectée.    
"""
    ),
    8: _(
        """
Le nom du group_ma %(k1)s est trop long, ce qui rend impossible le traitement du maillage dans la zone d'analyse.     
"""
    ),
    9: _(
        """
Le group_ma %(k1)s ne fait pas partie dans le maillage.    
"""
    ),
    10: _(
        """
La loi de comportement indiquée dans l'algorithme ne se trouve pas parmi les propriétés des matériaux affectant la zone d'analyse.     
"""
    ),
    11: _(
        """
Pour la loi MOHR_COULOMB, seule la loi d'écoulement plastique associée est autorisée, soit ANGDIL = PHI.
Systématiquement ANGDIL prend la valeur de PHI lorsqu'ils ne sont pas égaux.     
"""
    ),
    12: _(
        """
Un groupe_ma ne peut être affecté par qu'un seul matériau. 
"""
    ),
    13: _(
        """
Le FS initial doit être strictement positif.     
"""
    ),
    14: _(
        """
X1_MINI, X1_MAXI, X2_MINI, X2_MAXI doivent prendre les valeurs d'ordre croissante.     
"""
    ),
    15: _(
        """
Le champ PTOT doit affecter le modèle entier via TOUT = 'OUI' dans AFFE_VARC.     
"""
    ),
    16: _(
        """
Les zones des points d'extrémité sont en dehors du profil de la pente. 
"""
    ),
    17: _(
        """
La résolution du FS selon la procédure MORGENSTERN_PRICE diverge. 
"""
    ),
    18: _(
        """
Les opérandes Y_MINI et Y_MAXI doivent apparaître en paire.    
"""
    ),
    19: _(
        """
La valeur de Y_MINI est trop grande.   
"""
    ),
    20: _(
        """
La population des feux d'artifices doit être strictement positive.    
"""
    ),
    21: _(
        """
Impossible de déterminer les limites de l'ordonnée du point %(i1)d sur la surface de rupture. 
"""
    ),
}
