# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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

# person_in_charge: sylvie.michel-ponnelle at edf.fr

from .cata_comportement import LoiComportement

loi = LoiComportement(
    nom="SECH_RFT",
    lc_type=("SECHAGE",),
    doc="""Relation de comportement de thermique non lineaire pour modéliser le séchage du béton suivant le modèle de Richards Fick avec tempétarute (RFT)""",
    num_lc=0,
    nb_vari=0,
    nom_vari=None,
    mc_mater=None,
    modelisation=("3D", "AXIS", "PLAN", "3D_DIAG", "PLAN_DIAG", "AXIS_DIAG"),
    deformation=("PETIT", "PETIT_REAC", "GROT_GDEP"),
    algo_inte=("SANS_OBJET",),
    type_matr_tang=None,
    proprietes=None,
    syme_matr_tang=("Yes",),
    exte_vari=None,
    deform_ldc=("OLD",),
    regu_visc=("No",),
    post_incr=None,
)
