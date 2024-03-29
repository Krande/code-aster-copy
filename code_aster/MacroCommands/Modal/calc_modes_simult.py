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

# person_in_charge: nicolas.brie at edf.fr

from ...Cata.Syntax import _F
from .mode_iter_simult import MODE_ITER_SIMULT


def calc_modes_simult(self, stop_erreur, sturm, TYPE_RESU, OPTION, INFO, **args):
    """
    Macro-command CALC_MODES, case of the simultaneous iterations method
    """

    args = _F(args)
    SOLVEUR = args.get("SOLVEUR")
    SOLVEUR_MODAL = args.get("SOLVEUR_MODAL")
    VERI_MODE = args.get("VERI_MODE")
    TITRE = args.get("TITRE")
    CARA_ELEM = args.get("CARA_ELEM")
    CHAM_MATER = args.get("CHAM_MATER")

    motcles = {}
    matrices = {}

    if CARA_ELEM is not None:
        motcles["CARA_ELEM"] = CARA_ELEM
    if CHAM_MATER is not None:
        motcles["CHAM_MATER"] = CHAM_MATER

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

    #
    # read the keyword CALC_FREQ or CALC_CHAR_CRIT
    motcles_calc_vp = {}

    calc_vp = args["CALC_" + type_vp]
    nmax_vp = "NMAX_" + type_vp

    if OPTION in ("PLUS_PETITE", "PLUS_GRANDE"):
        motcles_calc_vp[nmax_vp] = calc_vp[nmax_vp]

    if OPTION == "CENTRE":
        motcles_calc_vp[type_vp] = calc_vp[type_vp]
        if type_vp == "FREQ":
            if calc_vp["AMOR_REDUIT"] is not None:
                motcles_calc_vp["AMOR_REDUIT"] = calc_vp["AMOR_REDUIT"]
        motcles_calc_vp[nmax_vp] = calc_vp[nmax_vp]

    if OPTION == "BANDE":
        motcles_calc_vp[type_vp] = calc_vp[type_vp]
        if calc_vp["TABLE_" + type_vp] is not None:
            motcles_calc_vp["TABLE_" + type_vp] = calc_vp["TABLE_" + type_vp]

    motcles_calc_vp["SEUIL_" + type_vp] = calc_vp["SEUIL_" + type_vp]

    if SOLVEUR_MODAL["DIM_SOUS_ESPACE"] is not None:
        motcles_calc_vp["DIM_SOUS_ESPACE"] = SOLVEUR_MODAL["DIM_SOUS_ESPACE"]
    if SOLVEUR_MODAL["COEF_DIM_ESPACE"] is not None:
        motcles_calc_vp["COEF_DIM_ESPACE"] = SOLVEUR_MODAL["COEF_DIM_ESPACE"]
    if SOLVEUR_MODAL["APPROCHE"] is not None:
        motcles_calc_vp["APPROCHE"] = SOLVEUR_MODAL["APPROCHE"]

    motcles["CALC_" + type_vp] = _F(
        OPTION=OPTION,
        NMAX_ITER_SHIFT=calc_vp["NMAX_ITER_SHIFT"],
        PREC_SHIFT=calc_vp["PREC_SHIFT"],
        **motcles_calc_vp
    )

    #
    # read the modal solver parameters
    motcles_solveur_modal = {}

    methode = SOLVEUR_MODAL["METHODE"]
    motcles_solveur_modal["METHODE"] = methode

    if methode == "TRI_DIAG":
        if SOLVEUR_MODAL["NMAX_ITER_ORTHO"] is not None:
            motcles_solveur_modal["NMAX_ITER_ORTHO"] = SOLVEUR_MODAL["NMAX_ITER_ORTHO"]
        if SOLVEUR_MODAL["PREC_ORTHO"] is not None:
            motcles_solveur_modal["PREC_ORTHO"] = SOLVEUR_MODAL["PREC_ORTHO"]
        if SOLVEUR_MODAL["PREC_LANCZOS"] is not None:
            motcles_solveur_modal["PREC_LANCZOS"] = SOLVEUR_MODAL["PREC_LANCZOS"]
        if SOLVEUR_MODAL["NMAX_ITER_QR"] is not None:
            motcles_solveur_modal["NMAX_ITER_QR"] = SOLVEUR_MODAL["NMAX_ITER_QR"]
        if SOLVEUR_MODAL["MODE_RIGIDE"] is not None:
            if SOLVEUR_MODAL["MODE_RIGIDE"] == "OUI":
                motcles["OPTION"] = "MODE_RIGIDE"
            else:
                motcles["OPTION"] = "SANS"
    elif methode == "JACOBI":
        if SOLVEUR_MODAL["NMAX_ITER_BATHE"] is not None:
            motcles_solveur_modal["NMAX_ITER_BATHE"] = SOLVEUR_MODAL["NMAX_ITER_BATHE"]
        if SOLVEUR_MODAL["PREC_BATHE"] is not None:
            motcles_solveur_modal["PREC_BATHE"] = SOLVEUR_MODAL["PREC_BATHE"]
        if SOLVEUR_MODAL["NMAX_ITER_JACOBI"] is not None:
            motcles_solveur_modal["NMAX_ITER_JACOBI"] = SOLVEUR_MODAL["NMAX_ITER_JACOBI"]
        if SOLVEUR_MODAL["PREC_JACOBI"] is not None:
            motcles_solveur_modal["PREC_JACOBI"] = SOLVEUR_MODAL["PREC_JACOBI"]
    elif methode == "SORENSEN":
        if SOLVEUR_MODAL["NMAX_ITER_SOREN"] is not None:
            motcles_solveur_modal["NMAX_ITER_SOREN"] = SOLVEUR_MODAL["NMAX_ITER_SOREN"]
        if SOLVEUR_MODAL["PARA_ORTHO_SOREN"] is not None:
            motcles_solveur_modal["PARA_ORTHO_SOREN"] = SOLVEUR_MODAL["PARA_ORTHO_SOREN"]
        if SOLVEUR_MODAL["PREC_SOREN"] is not None:
            motcles_solveur_modal["PREC_SOREN"] = SOLVEUR_MODAL["PREC_SOREN"]
    elif methode == "QZ":
        if SOLVEUR_MODAL["TYPE_QZ"] is not None:
            motcles_solveur_modal["TYPE_QZ"] = SOLVEUR_MODAL["TYPE_QZ"]

    motcles.update(motcles_solveur_modal)

    #
    # read the keyword SOLVEUR (linear solver)
    solveur = SOLVEUR[0].cree_dict_valeurs(SOLVEUR[0].mc_liste)
    if "TYPE_RESU" in solveur:  # because TYPE_RESU is a keyword with a 'global' position
        solveur.pop("TYPE_RESU")
    if "OPTION" in solveur:  # because OPTION is a keyword with a 'global' position
        solveur.pop("OPTION")
    if "FREQ" in solveur:  # because FREQ can be a keyword with a 'global' position
        solveur.pop("FREQ")
    motcles["SOLVEUR"] = _F(**solveur)

    #
    # read the keyword VERI_MODE
    if sturm in ("GLOBAL", "LOCAL", "OUI"):
        # for MODE_ITER_SIMULT, value for STURM can be only OUI or NON. Other
        # values are equivalent to OUI
        motveri = "OUI"
    elif sturm in ("NON"):
        # for keyword AMELIORATION
        motveri = "NON"
    else:
        assert False  # Pb parametrage STURM

    motcles["VERI_MODE"] = _F(
        STOP_ERREUR=stop_erreur,
        SEUIL=VERI_MODE["SEUIL"],
        STURM=motveri,
        PREC_SHIFT=VERI_MODE["PREC_SHIFT"],
    )

    #

    if args["STOP_BANDE_VIDE"] is not None:
        motcles["STOP_BANDE_VIDE"] = args["STOP_BANDE_VIDE"]

    if TITRE is not None:
        motcles["TITRE"] = TITRE

    modes = MODE_ITER_SIMULT(TYPE_RESU=TYPE_RESU, INFO=INFO, **motcles)

    materialField = None
    for matrice in matrices.values():
        try:
            if matrice.getNumberOfElementaryMatrix() != 0:
                modes.setMaterialField(matrice.getMaterialField())
                break
        except:
            pass

    return modes
