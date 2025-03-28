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
    nom="ELAS",
    lc_type=("MECANIQUE",),
    doc="""élasticité linéaire isotrope""",
    num_lc=1,
    nb_vari=1,
    nom_vari=("VIDE",),
    mc_mater=("ELAS", "ELAS_ORTH", "ELAS_ISTR"),
    modelisation=("3D", "AXIS", "C_PLAN", "D_PLAN", "1D"),
    deformation=("PETIT", "PETIT_REAC", "GROT_GDEP", "GDEF_LOG", "GREEN_LAGRANGE"),
    algo_inte=("ANALYTIQUE",),
    type_matr_tang=("PERTURBATION", "VERIFICATION", "IMPLEX"),
    proprietes=("COMP_ELAS",),
    syme_matr_tang=("Yes",),
    exte_vari=None,
    deform_ldc=("OLD",),
    regu_visc=("No",),
    post_incr=None,
)
