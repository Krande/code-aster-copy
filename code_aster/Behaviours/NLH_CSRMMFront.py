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

# person_in_charge: simon.raude at edf.fr

from .cata_comportement import LoiComportementMFront

loi = LoiComportementMFront(
    nom="NLH_CSRM",
    lc_type=("MECANIQUE",),
    doc="""To complete ...""",
    num_lc=58,
    nb_vari=0,
    nom_vari=None,
    mc_mater=None,
    modelisation=("3D", "AXIS", "D_PLAN"),
    deformation=("PETIT", "PETIT_REAC", "GDEF_LOG"),
    algo_inte=("NEWTON", "NEWTON_PERT"),
    type_matr_tang=("PERTURBATION", "VERIFICATION"),
    proprietes=None,
    syme_matr_tang=("No",),
    exte_vari=None,
    deform_ldc=("MECANIQUE",),
    regu_visc=("No",),
)
