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

# person_in_charge: sofiane.hendili at edf.fr

from .cata_comportement import LoiComportement

loi = LoiComportement(
    nom="META_G_ISOT_PTRE",
    lc_type=("META_G_ISOT_PTRE",),
    doc="""Loi de comportement prenant en compte la métallurgie - Ecrouissage isotrope avec plasticité de transformation et resturation d'écrouissage""",
    num_lc=0,
    nb_vari=2,
    nom_vari=("EPSPEQ", "INDIPLAS"),
    mc_mater=None,
    modelisation=("3D", "AXIS", "D_PLAN"),
    deformation=("PETIT", "PETIT_REAC", "GROT_GDEP", "SIMO_MIEHE"),
    algo_inte=("SPECIFIQUE",),
    type_matr_tang=("PERTURBATION", "VERIFICATION"),
    proprietes=None,
    syme_matr_tang=("Yes",),
    exte_vari=None,
    deform_ldc=("OLD",),
    regu_visc=("No",),
    post_incr=None,
)
