# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
 Des valeurs différentes de temps caractéristiques de viscosité ont été renseignées dans la définition des paramètres matériaux.
 On a retenu la plus grande de ces valeurs. Vérifier la mise en donnée.
"""
    ),
    2: _(
        """
 Aucune valeur de temps caractéristique de viscosité n'a été détectée dans la définition des paramètres matériaux.
 On utilisera TAU = 1.0
"""
    ),
    3: _(
        """
 Mot-clé OBSERVATION_VISC manquant dans ENDO_VISC.
 Il est obligatoire de définir une ou plusieurs quantités sur lesquelles faire porter le critère de stabilisation.
"""
    ),
    4: _(
        """
 CALC_ENDO ne sait pas traiter un chargement AFFE_CHAR_MECA_F fonction du temps.
 Il est recommandé d'utiliser FONC_MULT.
"""
    ),
    5: _(
        """
 CALC_ENDO n'a détecté aucun chargement ou variable de commande dépendant du temps.
 Vérifier la mise en données.
"""
    ),
    6: _(
        """
 Mot-clé OBSERVATION_VISC dans ENDO_VISC : les cas EVAL_CHAM = "VALE" et EVAL_ELGA = "VALE" ne sont pas prévus.
 Vérifier la mise en données.
"""
    ),
    7: _(
        """
 Mot-clé OBSERVATION_VISC dans ENDO_VISC : on ne peut définir qu'une composante NOM_CMP par mot-clé facteur.
 Il faut définir autant de mots-clés facteurs que de composantes souhaitées.
"""
    ),
    8: _(
        """
 Dans ETAT_INIT, les mots-clés NUME_ORDRE, INST_INIT et NUME_DIDI sont interdits.
 La poursuite d'un calcul se fait nécessairement à partir du dernier instant.
"""
    ),

}
