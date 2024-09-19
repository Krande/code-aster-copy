# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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
import libAsterGC
import aster
from ..CodeCommands import CREA_TABLE
from ..Messages import UTMESS


def trouver_marges_seuil(data, seuil):
    """Extracts the list of mesh groups where the margin is below a threshold.

    This function scans the list of MARGE values and identifies the corresponding
    mesh entities whose margin is strictly below the given threshold value.

    Arguments:
        data (tuple): A tuple where the first element is the list of MARGE values
            (floats), and the second element is a list containing the mesh group
            names associated to each margin.
        seuil (float): Threshold value below which margins are considered critical.

    Returns:
        list[str]: A list of cells numbers whose associated margin is
        strictly lower than the specified threshold.
    """
    marges = data[0]
    mailles = data[1][0]
    mailles_correspondantes = []
    for i, marge in enumerate(marges):
        if marge < seuil:
            mailles_correspondantes.append(mailles[i])

    return mailles_correspondantes


def trouver_NDIAG_marges_minimales(data, N_DIAG):
    """Finds the mesh groups with the N smallest margin values.

    This function sorts all mesh cells by their margin values in ascending order
    and extracts the names of those with the N_DIAG lowest margins.

    Arguments:
        data (tuple): A tuple where the first element is the list of MARGE values
            (floats), and the second element is a list containing the mesh group
            names associated to each margin.
        N_DIAG (int): The number of diagrams to consider. .

    Returns:
        list[str]: A list of cells numbers associated with the N_DIAG
        lowest margins.
    """
    marges = data[0]
    mailles = data[1][0]
    marges_mailles = list(zip(marges, mailles))

    # Trier la liste des tuples par marge
    marges_mailles.sort(key=lambda x: x[0])

    # Extraire les N_DIAG premières marges et leurs mailles correspondantes
    mailles_correspondantes = [maille for marge, maille in marges_mailles[:N_DIAG]]

    return mailles_correspondantes


def post_veri_ferraillage_ops(self, CHAM_VFER, **args):
    """
    macro POST_VERI_FERRAILLAGE
    """
    if (len(args)) > 1:
        UTMESS("F", "VERIFERRAILLAGE_15")

    Ma = CHAM_VFER.getMesh()
    lisgrma = Ma.getGroupsOfCells()

    if "GROUP_MA" in args:
        GROUP_MA = args.get("GROUP_MA")
        if not (all(element in lisgrma for element in GROUP_MA)):
            UTMESS("F", "VERIFERRAILLAGE_12")
        list_ma = CHAM_VFER.getValuesWithDescription("MARGE", GROUP_MA)[1][0]
    elif "SEUIL_MIN_MARGE" in args:
        marge_struct = CHAM_VFER.getValuesWithDescription("MARGE")
        SEUIL_MARGE = args["SEUIL_MIN_MARGE"]
        list_ma = trouver_marges_seuil(marge_struct, SEUIL_MARGE)
        if len(list_ma) == 0:
            UTMESS("F", "VERIFERRAILLAGE_13")
    elif "NB_DIAG" in args:
        marge_struct = CHAM_VFER.getValuesWithDescription("MARGE")
        NB_DIAG = args.get("NB_DIAG")
        list_ma = trouver_NDIAG_marges_minimales(marge_struct, NB_DIAG)
        if NB_DIAG > len(list_ma):
            UTMESS("F", "VERIFERRAILLAGE_14")
    else:
        raise KeyError("unexpected arguments")

    # lecture des paramètres CHAM_VFER
    cmps = CHAM_VFER.getComponents()
    dint_param = CHAM_VFER.getValuesWithDescription([], list_ma)[0]

    # Initialisation du dictionnaire de résultats
    resultats = {nom: [] for nom in cmps}
    n_composantes = len(cmps)
    # Remplissage du dictionnaire
    for i in range(len(list_ma)):
        for j, nom in enumerate(cmps):
            resultats[nom].append(dint_param[i * n_composantes + j])

    dnsinf_crit = resultats["DNSINF"]
    dnssup_crit = resultats["DNSSUP"]
    effn0_crit = resultats["N0"]
    effm0_crit = resultats["M0"]
    effn_crit = resultats["NED"]
    effm_crit = resultats["MED"]
    typcmb = resultats["TYPCMB"]
    typco = resultats["TYPCO"]
    cequi = resultats["CEQUI"]
    enrobi = resultats["ENROBI"]
    enrobs = resultats["ENROBS"]
    sigs = resultats["SIGS"]
    sigci = resultats["SIGCI"]
    sigcs = resultats["SIGCS"]
    alphacc = resultats["ALPHACC"]
    gammas = resultats["GAMMAS"]
    gammac = resultats["GAMMAC"]
    facier = resultats["FACIER"]
    eys = resultats["EYS"]
    typdiag = resultats["TYPDIAG"]
    fbeton = resultats["FBETON"]
    clacier = resultats["CLACIER"]
    uc = resultats["UC"]
    ht = resultats["HT"]
    bw = resultats["BW"]
    length = len(list_ma)

    # Appeler la fonction pour chaque maille
    # Listes pour stocker les résultats
    nom_results = []
    nrd_results = []
    mrd_results = []
    for i in range(length):
        if typcmb[i] == 0:
            typco[i] = int(typco[i])
            clacier[i] = int(clacier[i])
            typdiag[i] = int(typdiag[i])
            uc[i] = int(uc[i])
            nrd, mrd = libAsterGC.dintelu(
                typco[i],
                alphacc[i],
                ht[i],
                bw[i],
                enrobi[i],
                enrobs[i],
                facier[i],
                fbeton[i],
                gammas[i],
                gammac[i],
                clacier[i],
                eys[i],
                typdiag[i],
                uc[i],
                dnsinf_crit[i],
                dnssup_crit[i],
            )

        elif typcmb[i] == 1:
            uc[i] = int(uc[i])
            nrd, mrd = libAsterGC.dintels(
                cequi[i],
                ht[i],
                bw[i],
                enrobi[i],
                enrobs[i],
                sigci[i],
                sigcs[i],
                sigs[i],
                uc[i],
                dnsinf_crit[i],
                dnssup_crit[i],
            )
        # on ajoute 1 au numero de maille pour correspondre a la numerotation de salome
        nom_results.extend([list_ma[i] + 1] * (len(nrd) + 2))
        nrd_results.extend(nrd)
        nrd_results.append(effn0_crit[i])
        nrd_results.append(effn_crit[i])
        mrd_results.extend(mrd)
        mrd_results.append(effm0_crit[i])
        mrd_results.append(effm_crit[i])

    tabout = CREA_TABLE(
        LISTE=(
            _F(PARA="NUMERO_MAILLE", LISTE_I=nom_results),
            _F(PARA="N", LISTE_R=nrd_results),
            _F(PARA="M", LISTE_R=mrd_results),
        )
    )

    return tabout
