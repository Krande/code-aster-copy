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

# person_in_charge: nicolas.brie at edf.fr

from ...Cata.Syntax import _F
from ...Messages import UTMESS
from ...Objects import AssemblyMatrixDisplacementComplex, GeneralizedAssemblyMatrixComplex
from .mode_iter_inv import MODE_ITER_INV


def calc_modes_amelioration(self, modes, TYPE_RESU, INFO, **args):
    """
    Macro-command CALC_MODES, file for improving the quality of the eigenmodes
    """

    args = _F(args)
    SOLVEUR = args.get("SOLVEUR")
    VERI_MODE = args.get("VERI_MODE")
    TITRE = args.get("TITRE")

    # import the definitions of the commands to use in the macro-command
    # The name of the variable has to be the name of the command

    ##############################################################################
    # 1. CHECK IF THE COMPUTATION WITH MODE_ITER_INV CAN BE PERFORMED

    if TYPE_RESU == "DYNAMIQUE":
        type_vp = "FREQ"
        matr_A = "MATR_RIGI"
        matr_B = "MATR_MASS"
    elif TYPE_RESU == "MODE_FLAMB":
        type_vp = "CHAR_CRIT"
        matr_A = "MATR_RIGI"
        matr_B = "MATR_RIGI_GEOM"
    elif TYPE_RESU == "GENERAL":
        type_vp = "CHAR_CRIT"
        matr_A = "MATR_A"
        matr_B = "MATR_B"

    # 1.1 check if the input matrices are symmetric and real
    lsym = True
    lreel = True
    if not args[matr_A].isSymmetric():
        lsym = False
    if isinstance(args[matr_A], AssemblyMatrixDisplacementComplex) or isinstance(
        args[matr_A], GeneralizedAssemblyMatrixComplex
    ):
        lreel = False

    if not args[matr_B].isSymmetric():
        lsym = False
    if TYPE_RESU == "DYNAMIQUE":
        if args["MATR_AMOR"] is not None:
            if not args["MATR_AMOR"].isSymmetric():
                lsym = False
    if not lsym:
        UTMESS("I", "MODAL_15")
        return modes
    if not lreel:
        UTMESS("I", "MODAL_16", valk=matr_A)
        return modes

    # 1.2 detect too closed eigen values (gap < CALC_* / SEUIL_*)
    list_vp = modes.LIST_PARA()[type_vp]  # list of the eigen values
    seuil_vp = args["CALC_" + type_vp]["SEUIL_" + type_vp]
    # One reproduces here the detection performed in the routine vpgsmm.F90
    seuilr = 100.0 * seuil_vp
    seuilp = seuil_vp
    ltest = False
    for k in range(0, len(list_vp) - 1):
        if abs(list_vp[k]) < seuilr:
            ltest = abs(list_vp[k + 1]) < seuilr
        else:
            if abs(list_vp[k + 1]) >= seuilr:
                ltest = (
                    abs(2.0 * (list_vp[k] - list_vp[k + 1]) / (list_vp[k] + list_vp[k + 1]))
                    < seuilp
                )
    if ltest:
        UTMESS("I", "MODAL_17", valk=modes.getName())
        return modes

    ##############################################################################
    # 2. PERFORM THE IMPROVEMENT OF THE MODES WITH MODE_ITER_INV / OPTION='PROCHE'
    motcles = {}
    matrices = {}

    # read the input matrices
    if TYPE_RESU == "DYNAMIQUE":
        type_vp = "FREQ"
        matrices["MATR_RIGI"] = args["MATR_RIGI"]
        matrices["MATR_MASS"] = args["MATR_MASS"]
        if args["MATR_AMOR"] is not None:
            matrices["MATR_AMOR"] = args["MATR_AMOR"]

    elif TYPE_RESU == "MODE_FLAMB":
        type_vp = "CHAR_CRIT"
        matrices["MATR_RIGI"] = args["MATR_RIGI"]
        matrices["MATR_RIGI_GEOM"] = args["MATR_RIGI_GEOM"]

    elif TYPE_RESU == "GENERAL":
        type_vp = "CHAR_CRIT"
        matrices["MATR_A"] = args["MATR_A"]
        matrices["MATR_B"] = args["MATR_B"]

    motcles.update(matrices)

    #################################################################

    motcles_calc_vp = {}

    motcles_calc_vp[type_vp] = list_vp

    motcles["CALC_" + type_vp] = _F(OPTION="PROCHE", **motcles_calc_vp)

    #################################################################
    # read the keyword SOLVEUR (linear solver)
    solveur = SOLVEUR[0].cree_dict_valeurs(SOLVEUR[0].mc_liste)
    if "TYPE_RESU" in solveur:  # because TYPE_RESU is a keyword with a 'global' position
        solveur.pop("TYPE_RESU")
    if "OPTION" in solveur:  # because OPTION is a keyword with a 'global' position
        solveur.pop("OPTION")
    if "FREQ" in solveur:  # because FREQ can be a keyword with a 'global' position
        solveur.pop("FREQ")
    motcles["SOLVEUR"] = _F(**solveur)  # if this line is commented,
    # one will use the default keywords for SOLVEUR

    #################################################################
    # read the keyword VERI_MODE
    #
    sturm = VERI_MODE["STURM"]
    if sturm in ("GLOBAL", "LOCAL", "OUI"):
        # for MODE_ITER_INV, value for STURM can be only OUI or NON. Other
        # values are equivalent to OUI
        motveri = "OUI"
    elif sturm in ("NON"):
        # for keyword AMELIORATION
        motveri = "NON"
    else:
        assert False  # Pb parametrage STURM                              )
    motcles["VERI_MODE"] = _F(
        STOP_ERREUR=VERI_MODE["STOP_ERREUR"],
        SEUIL=VERI_MODE["SEUIL"],
        STURM=motveri,
        PREC_SHIFT=VERI_MODE["PREC_SHIFT"],
    )
    #################################################################
    if TITRE is not None:
        motcles["TITRE"] = TITRE

    #################################################################

    modes = MODE_ITER_INV(TYPE_RESU=TYPE_RESU, INFO=INFO, **motcles)

    return modes
