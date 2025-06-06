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
    1: _("""Liste des champs lus."""),
    2: _("""Champ %(k1)s."""),
    3: _("""Pour le numéro d'ordre %(i1)d, le paramètre d'accès de nom %(k1)s vaut %(r1)g."""),
    4: _("""Ce champ n'est pas autorisé dans le résultat."""),
    5: _(
        """le numéro d'archivage est inférieur au numéro précédent. Il doit être strictement croissant."""
    ),
    6: _(
        """Cette opération n'est pas possible si le modèle varie.

Risques & Conseils :
    Limiter le traitement sur des numéros d'ordre où le modèle ne change pas.
         """
    ),
    7: _(
        """Cette opération n'est pas possible si le champ de matériau varie.

Risques & Conseils :
    Limiter le traitement sur des numéros d'ordre où le champ de matériau ne change pas.
         """
    ),
    8: _(
        """Cette opération n'est pas possible si les caractéristiques élémentaires varient.

Risques & Conseils :
    Limiter le traitement sur des numéros d'ordre où les caractéristiques élémentaires
     ne changent pas.
         """
    ),
    9: _("""Cette structure de données n'est pas indexée par le temps."""),
    10: _(
        """L'index de la structure de données n'est pas disponible. Il faut allouer la structure de données avant d'y accéder."""
    ),
    14: _("""Le NUME_DDL a été déterminé à partir de la matrice de rigidité %(k1)s."""),
    15: _("""Les NUME_DDL associés aux matrices MATR_RIGI et MATR_MASS sont différents."""),
    16: _(
        """Vous essayez de stocker les chargements dans le résultat. Ce n'est pas possible pour un résultat de type %(k1)s, on ne stocke pas le chargement."""
    ),
    17: _(
        """Le modèle étant absent, on ne peut pas vérifier le comportement et la cohérence du nombre de variables internes.
 Il faut renseigner le mot-clé MODELE.
"""
    ),
    18: _(
        """Il n'est pas possible d'avoir plusieurs types de champ simultanément dans LIRE_RESU pour la structure de données des modes empiriques."""
    ),
    19: _(
        """Il n'est pas possible d'utiliser TOUT_ORDRE lorsqu'on enrichit un résultat dans LIRE_RESU. Il faut sélectionner les champs à lire."""
    ),
    20: _(
        """Le dernier instant lu était %(r1)g, il n'est donc pas possible de commencer à enrichir le résultat avec un champ à l'instant %(r2)g (car les instants doivent être strictement croissants)."""
    ),
    21: _(
        """La dernière fréquence lue était %(r1)g, il n'est donc pas possible de commencer à enrichir le résultat avec un champ à la fréquence %(r2)g (car les fréquences doivent être strictement croissantes)."""
    ),
    22: _(
        """Le dernier indice de rangement lu était %(i1)d, il n'est donc pas possible de commencer à enrichir le résultat avec un champ à l'indice de rangement %(i2)d (car ils doivent doivent être strictement croissants)."""
    ),
    24: _("""Le champ %(k2)s est incompatible avec le type de résultat %(k1)s."""),
    94: _("""Le champ %(k1)s n'est pas prévu dans LIRE_RESU. Vous pouvez demander l'évolution."""),
    95: _("""Le champ %(k1)s n'est pas prévu dans LIRE_CHAMP. Vous pouvez demander l'évolution."""),
    97: _(
        """On n'a pu lire aucun champ dans le fichier. La structure de données créée est vide.
Risques & Conseils :
  Si le fichier lu est au format IDEAS, et si la commande est LIRE_RESU,
  le problème vient peut-être d'une mauvaise utilisation (ou d'une absence d'utilisation)
  du mot clé FORMAT_IDEAS. Il faut examiner les "entêtes" des DATASET du fichier à lire.
"""
    ),
}
