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

# person_in_charge: samuel.geniaut at edf.fr

from .cata_comportement import LoiComportement

loi = LoiComportement(
    nom="VISC_ENDO_LEMA",
    lc_type=("MECANIQUE",),
    doc="""Modèle viscoplastique couplé à l'endommagement isotrope de Lemaitre-Chaboche [R5.03.15].
   Ce modèle s'emploie avec les mots clés DEFORMATION = PETIT ou PETIT_REAC.""",
    num_lc=31,
    nb_vari=10,
    nom_vari=(
        "EPSPXX",
        "EPSPYY",
        "EPSPZZ",
        "EPSPXY",
        "EPSPXZ",
        "EPSPYZ",
        "EPSPEQ",
        "ECROISOT",
        "ENDO",
        "INDIPLAS",
    ),
    mc_mater=("ELAS", "VENDOCHAB"),
    modelisation=("3D", "AXIS", "D_PLAN"),
    deformation=("PETIT", "PETIT_REAC", "GROT_GDEP"),
    algo_inte=("SECANTE", "BRENT", "DEKKER"),
    type_matr_tang=("PERTURBATION", "VERIFICATION"),
    proprietes=None,
    syme_matr_tang=("Yes",),
    exte_vari=None,
    deform_ldc=("OLD",),
    regu_visc=("No",),
    post_incr=None,
)
