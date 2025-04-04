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
    nom="BETON_GRANGER",
    lc_type=("MECANIQUE",),
    doc="""
        Comportement de fluage propre du beton, identique à BETON_GRANGER_V mais
        traitant uniquement un comportement isotherme. cf. R7.01.01
    """,
    num_lc=26,
    nb_vari=55,
    nom_vari=(
        "VG1",
        "VG2",
        "VG3",
        "VG4",
        "VG5",
        "VG6",
        "VG7",
        "VG8",
        "VG9",
        "VG10",
        "VG11",
        "VG12",
        "VG13",
        "VG14",
        "VG15",
        "VG16",
        "VG17",
        "VG18",
        "VG19",
        "VG20",
        "VG21",
        "VG22",
        "VG23",
        "VG24",
        "VG25",
        "VG26",
        "VG27",
        "VG28",
        "VG29",
        "VG30",
        "VG31",
        "VG32",
        "VG33",
        "VG34",
        "VG35",
        "VG36",
        "VG37",
        "VG38",
        "VG39",
        "VG40",
        "VG41",
        "VG42",
        "VG43",
        "VG44",
        "VG45",
        "VG46",
        "VG47",
        "VG48",
        "VG49",
        "VG50",
        "VG51",
        "VG52",
        "VG53",
        "VG54",
        "VG55",
    ),
    mc_mater=("ELAS", "BETON_GRANGER"),
    modelisation=("3D", "AXIS", "D_PLAN"),
    deformation=("PETIT", "PETIT_REAC", "GROT_GDEP"),
    algo_inte=("ANALYTIQUE",),
    type_matr_tang=("PERTURBATION", "VERIFICATION"),
    proprietes=None,
    syme_matr_tang=("Yes",),
    exte_vari=None,
    deform_ldc=("OLD",),
    regu_visc=("No",),
    post_incr=None,
)
