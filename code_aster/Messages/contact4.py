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
    1: _(
        """
Contact LAC
Seuls les algorithmes en Newton sont utilisables (ALGO_RESO_GEOM et ALGO_RESO_CONT)
"""
    ),
    2: _(
        """
Contact LAC
    Le maillage %(k1)s ne contient pas les objets spécifiques à la méthode ALGO_CONT='LAC'.
Conseil:
    Il faut faire CREA_MAILLAGE/DECOUPE_LAC avant DEFI_CONTACT
"""
    ),
    4: _(
        """
Contact LAC
        ALGO_CONT='LAC' ne fonctionne pas avec le frottement.
"""
    ),
    5: _(
        """
Contact LAC
         On ne détecte pas le bon nombre de mailles esclaves.
         Conseil :
             Cette erreur est probablement dû au fait que vous avez inversé les rôles maîtres et esclaves.
             Vérifiez que votre GROUP_MA_ESCL est bien celui utilisé par DECOUPE_LAC de CREA_MAILLAGE.
"""
    ),
    6: _(
        """
Contact appariement de type LAC
         Une erreur non fatale est détectée lors de l'appariement de la maille esclave numéro %(i1)d et de la maille maître numéro %(i2)d .
         Ce couple est ignoré, il y a potentiellement interpénétration entre ces deux mailles.
         Conseil :
             Vérifiez la qualité de l'appariement visuellement au niveau des mailles concernées.
"""
    ),
}
