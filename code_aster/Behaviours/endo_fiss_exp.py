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


from .cata_comportement import LoiComportement

loi = LoiComportement(
    nom="ENDO_FISS_EXP",
    lc_type=("MECANIQUE",),
    doc="""Comportement élastique-fragile, à endommagement scalaire, seuil exponentiel et
   non local à gradient d'endommagement - R5.03.25""",
    num_lc=57,
    nb_vari=9,
    nom_vari=(
        "ENDO",
        "INDIENDO",
        "ENDORIGI",
        "EPSEXX",
        "EPSEYY",
        "EPSEZZ",
        "EPSEXY",
        "EPSEXZ",
        "EPSEYZ",
    ),
    mc_mater=("ELAS", "ENDO_SCALAIRE", "NON_LOCAL"),
    modelisation=("GRADVARI",),
    deformation=("PETIT",),
    algo_inte=("NEWTON",),
    type_matr_tang=("PERTURBATION", "VERIFICATION"),
    proprietes=("COMP_ELAS",),
    syme_matr_tang=("No",),
    exte_vari=None,
    deform_ldc=("OLD",),
    regu_visc=("REGU_VISC_ELAS",),
    post_incr=None,
)
