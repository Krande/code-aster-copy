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
    nom="KIT_THM",
    lc_type=("KIT_THM",),
    doc="""KIT associé au comportement des milieux poreux (modélisations thermo-hydro-mécanique).
   Pour plus de détails sur les modélisations thermo-hydro-mécaniques et les modèles de comportement,
   on pourra consulter les documents [R7.01.10] et [R7.01.11], ainsi que la notice d'utilisation [U2.04.05].
   Les relations KIT_XXXX permettent de résoudre simultanément de deux à quatre équations d'équilibre.
   Les équations considérées dépendent du suffixe XXXX avec la règle suivante :
   - M désigne l'équation d'équilibre mécanique,
   - T désigne l'équation d'équilibre thermique,
   - H désigne une équation d'équilibre hydraulique.
   - V désigne la présence d'une phase sous forme vapeur (en plus du liquide)
   Les problèmes thermo-hydro-mécaniques associés sont traités de facon totalement couplée.
   Une seule lettre H signifie que le milieu poreux est saturé (une seule variable de pression p),
   par exemple soit de gaz, soit de liquide, soit d'un mélange liquide/gaz (dont la pression du gaz est constante).
   Deux lettres H signifient que le milieu poreux est non saturé (deux variables de pression p), par exemple
   un mélange liquide/vapeur/gaz. La présence des deux lettres HV signifie que le milieu poreux est saturé par
   un composant (en pratique de l'eau), mais que ce composant peut être sous forme liquide ou vapeur.
   Il n'y a alors qu'une équation de conservation de ce composant, donc un seul degré de liberté pression,
   mais il y a un flux liquide et un flux vapeur.
   """,
    num_lc=0,
    nb_vari=0,
    nom_vari=None,
    mc_mater=None,
    modelisation=(
        "D_PLAN_THM",
        "D_PLAN_THMS",
        "D_PLAN_THMD",
        "AXIS_THM",
        "AXIS_THMS",
        "AXIS_THMD",
        "3D_THM",
        "3D_THMS",
        "3D_THMD",
    ),
    deformation=("PETIT",),
    algo_inte=("SANS_OBJET",),
    type_matr_tang=None,
    proprietes=None,
    syme_matr_tang=("Yes",),
    exte_vari=None,
    deform_ldc=("OLD",),
    regu_visc=("No",),
    post_incr=None,
)
