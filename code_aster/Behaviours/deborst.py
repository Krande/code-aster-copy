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
    nom="DEBORST",
    lc_type=("UTILITAIRE",),
    doc="""Algo pour résolution en contraintes planes.""",
    num_lc=0,
    nb_vari=1,
    nom_vari=("EPZZ_CP",),
    mc_mater=None,
    modelisation=("C_PLAN",),
    deformation=("PETIT", "GROT_GDEP"),
    algo_inte=None,
    type_matr_tang=None,
    proprietes=None,
    syme_matr_tang=None,
    exte_vari=None,
    deform_ldc=("OLD",),
    regu_visc=("No",),
    post_incr=None,
)
