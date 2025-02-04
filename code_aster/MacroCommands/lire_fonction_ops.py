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

# person_in_charge: mathieu.courtois at edf.fr

import os.path as osp
from math import cos, sin

import numpy

from ..Messages import UTMESS

from ..Cata.Syntax import _F
from ..CodeCommands import DEFI_FONCTION, DEFI_NAPPE
from ..Helpers import LogicalUnitFile, ReservedUnitUsed
from ..Utilities.transpose import transpose


class LectureBlocError(Exception):
    pass


def lire_blocs(nomfich, SEPARATEUR, INFO=1):
    """Retourne la liste des blocs"""

    def info(ib, nlig, ncol):
        """Affiche les infos d'un bloc"""
        if INFO < 2:
            return
        print("   . Bloc %2d : %6d lignes, %6d colonnes" % (ib, nlig, ncol))

    print("  Lecture des blocs du fichier '%s'..." % nomfich)
    fich = open(nomfich, "r")
    if SEPARATEUR == "None":
        SEPARATEUR = None
    blocs = []
    lignes = []
    llen = 0
    il = 0
    for line in fich:
        il += 1
        try:
            line = line.strip()
            if line == "":
                raise ValueError
            splin = [i for i in line.split(SEPARATEUR) if i != ""]
            lignes.append(list(map(float, splin)))
            if llen == 0:
                llen = len(splin)
            elif len(splin) != llen:
                raise LectureBlocError(
                    "Ligne {0} : {1} champs au lieu de {2} " "attendus".format(il, len(splin), llen)
                )
        except ValueError:
            if lignes == []:
                pass  # dans ce cas, on a plusieurs lignes délimitant 2 fonctions
            else:
                blocs.append(numpy.array(lignes))
                info(len(blocs), len(lignes), llen)
                lignes = []
                llen = 0
    fich.close()
    if len(lignes) > 0:
        blocs.append(numpy.array(lignes))
        info(len(blocs), len(lignes), llen)
    return blocs


def liste_double(nomfich, INDIC_PARA, INDIC_RESU, SEPARATEUR, INFO=1):
    """Méthode de construction du VALE pour le format libre
    format LIBRE
    Les lignes contenant autre chose que des séquences de nombres
    réels et de séparateurs sont considérées comme délimitant deux
    fonctions différentes. Cette situation correspond à l exception
    ValueError levée par le map de float. Le deuxieme indice de
    INDIC_PARA et INDIC_RESU est l indice permettant de pointer sur la
    fonction voulue, au sens de ce découpage.
    """
    blocs = lire_blocs(nomfich, SEPARATEUR, INFO)

    # vérifications de cohérences lignes et colonnes
    nb_blocs = len(blocs)
    bloc_para = INDIC_PARA[0]
    col_para = INDIC_PARA[1]
    bloc_resu = INDIC_RESU[0]
    col_resu = INDIC_RESU[1]
    if bloc_para > nb_blocs:
        raise LectureBlocError(
            "Il y a {0} blocs or INDIC_PARA=({1}, .)".format(nb_blocs, bloc_para)
        )
    if bloc_resu > nb_blocs:
        raise LectureBlocError(
            "Il y a {0} blocs or INDIC_RESU=({1}, .)".format(nb_blocs, bloc_resu)
        )

    if col_para > len(blocs[bloc_para - 1][0]):
        raise LectureBlocError(
            "Le bloc {0} comporte {1} colonnes or "
            "INDIC_PARA=(., {2})".format(bloc_para, len(blocs[bloc_para - 1][0]), col_para)
        )
    if col_resu > len(blocs[bloc_resu - 1][0]):
        raise LectureBlocError(
            "Le bloc {0} comporte {1} colonnes or "
            "INDIC_RESU=(., {2})".format(bloc_resu, len(blocs[bloc_resu - 1][0]), col_resu)
        )

    # construction du VALE de la fonction par recherche des indices
    # de colonnes et de fonctions dans le tableau blocs
    vale_para = blocs[bloc_para - 1][:, col_para - 1]
    vale_resu = blocs[bloc_resu - 1][:, col_resu - 1]
    if len(vale_para) != len(vale_resu):
        print("VALE_PARA =", vale_para)
        print("VALE_RESU =", vale_resu)
        message = """Les deux colonnes extraites n'ont pas la meme longueur
         %d lignes pour PARA
         %d lignes pour RESU""" % (
            len(vale_para),
            len(vale_resu),
        )
        raise LectureBlocError(message)

    laux = transpose([vale_para, vale_resu])
    liste_vale = []
    for v in laux:
        liste_vale.extend(v)
    return liste_vale


def liste_simple(nomfich, INDIC_PARA, SEPARATEUR, INFO=1):
    """recherche d'une liste simple"""
    blocs = lire_blocs(nomfich, SEPARATEUR, INFO)

    # vérifications de cohérences lignes et colonnes
    nb_blocs = len(blocs)
    bloc_para = INDIC_PARA[0]
    col_para = INDIC_PARA[1]
    if bloc_para > nb_blocs:
        raise LectureBlocError(
            "Il y a {0} blocs or INDIC_PARA=({1}, .)".format(nb_blocs, bloc_para)
        )
    if col_para > len(blocs[bloc_para - 1][0]):
        raise LectureBlocError(
            "Le bloc {0} comporte {1} colonnes or "
            "INDIC_PARA=(., {2})".format(bloc_para, len(blocs[bloc_para - 1][0]), col_para)
        )

    # construction du VALE de la fonction par recherche des indices
    # de colonnes et de fonctions dans le tableau l_fonc
    vale_1 = blocs[bloc_para - 1][:, col_para - 1]
    return vale_1.tolist()


def column_values(fmt, filename, idbx, separ=" ", info=0):
    """Return the values of a column.

    Arguments:
        fmt (str): Format of the file to read.
        filename (str): File path.
        idbx ([int, int]): Indexes of the block to read and of the column in
            the block for the abscissas.
        separ (str): Text separator (if needed).
        info (int): Verbosity level.
    """
    if fmt == "LIBRE":
        try:
            values = liste_simple(filename, idbx, separ, info)
        except LectureBlocError as exc:
            UTMESS("F", "FONCT0_42", valk=exc.args)
    else:
        idx = idbx[1] - 1
        matrix = numpy.load(filename)
        values = matrix[:, idx]
    return values


def function_values(fmt, filename, idbx, idby, separ=" ", info=0):
    """Return the values to be passed to DEFI_FONCTION.

    Arguments:
        fmt (str): Format of the file to read.
        filename (str): File path.
        idbx ([int, int]): Indexes of the block to read and of the column in
            the block for the abscissas.
        idby ([int, int]): Some for the ordinates.
        separ (str): Text separator (if needed).
        info (int): Verbosity level.
    """
    if fmt == "LIBRE":
        try:
            values = liste_double(filename, idbx, idby, separ, info)
        except LectureBlocError as exc:
            UTMESS("F", "FONCT0_42", valk=exc.args)
    else:
        idx = idbx[1] - 1
        idy = idby[1] - 1
        matrix = numpy.load(filename)
        valx = matrix[:, idx]
        valy = matrix[:, idy]
        values = numpy.vstack((valx, valy)).transpose().ravel()
    return values


def complex_values(filename, idbx, idbr, idbi, module_phase=False):
    """Return the values to be passed to DEFI_FONCTION for a complex function.

    Arguments:
        filename (str): File path.
        idbx ([int, int]): Indexes of the block to read (unused) and of
            the column in the block for the abscissas.
        idbr ([int, int]): Some for the real part.
        idbi ([int, int]): Some for the imaginary part.
        module_phase (bool): Indicator that the columns contain module and
            phase values.
    """
    idx = idbx[1] - 1
    idr = idbr[1] - 1
    idi = idbi[1] - 1
    matrix = numpy.load(filename)
    valx = matrix[:, idx]
    valr = matrix[:, idr]
    vali = matrix[:, idi]
    if module_phase:
        module = valr
        phase = vali
        valr = module * numpy.cos(phase)
        vali = module * numpy.sin(phase)
    cols = numpy.vstack((valx, valr, vali)).transpose()
    return cols.ravel()


def lire_fonction_ops(
    self,
    UNITE,
    NOM_PARA,
    FORMAT=None,
    TYPE=None,
    SEPARATEUR=None,
    INDIC_PARA=None,
    NOM_RESU=None,
    INTERPOL=None,
    PROL_DROITE=None,
    PROL_GAUCHE=None,
    VERIF=None,
    INFO=None,
    TITRE=None,
    **args
):
    """Méthode corps de la macro"""

    # On recopie le mot cle defi_fonction pour le proteger
    if TYPE == "NAPPE":
        mc_DEFI_FONCTION = args["DEFI_FONCTION"]

    # On importe les definitions des commandes a utiliser dans la macro
    assert FORMAT in ("LIBRE", "NUMPY")

    # Lecture de la fonction dans un fichier d unité logique UNITE
    nomfich = LogicalUnitFile.filename_from_unit(UNITE)
    if not osp.isfile(nomfich):
        UTMESS("F", "FONCT0_41", valk=nomfich)

    if TYPE == "FONCTION":
        values = function_values(FORMAT, nomfich, INDIC_PARA, args["INDIC_RESU"], SEPARATEUR, INFO)
        # création de la fonction ASTER :
        ut_fonc = DEFI_FONCTION(
            NOM_PARA=NOM_PARA,
            NOM_RESU=NOM_RESU,
            PROL_DROITE=PROL_DROITE,
            PROL_GAUCHE=PROL_GAUCHE,
            INTERPOL=INTERPOL,
            INFO=INFO,
            TITRE=TITRE,
            VERIF=VERIF,
            VALE=values,
        )

    elif TYPE == "FONCTION_C":
        # mise en forme de la liste de valeurs suivant le format choisi :
        if FORMAT == "LIBRE":
            if "INDIC_REEL" in args:
                indic1 = args["INDIC_REEL"]
                indic2 = args["INDIC_IMAG"]
            if "INDIC_MODU" in args:
                indic1 = args["INDIC_MODU"]
                indic2 = args["INDIC_PHAS"]
            try:
                liste_vale_r = liste_double(nomfich, INDIC_PARA, indic1, SEPARATEUR, INFO)
            except LectureBlocError as exc:
                UTMESS("F", "FONCT0_42", valk=exc.args)

            try:
                liste_vale_i = liste_double(nomfich, INDIC_PARA, indic2, SEPARATEUR, INFO)
            except LectureBlocError as exc:
                UTMESS("F", "FONCT0_42", valk=exc.args)

            liste = []
            if "INDIC_REEL" in args:
                for i in range(len(liste_vale_r) // 2):
                    liste.extend(
                        [liste_vale_r[2 * i], liste_vale_r[2 * i + 1], liste_vale_i[2 * i + 1]]
                    )
            elif "INDIC_MODU" in args:
                for i in range(len(liste_vale_r) // 2):
                    module = liste_vale_r[2 * i + 1]
                    phase = liste_vale_i[2 * i + 1]
                    liste.extend([liste_vale_r[2 * i], module * cos(phase), module * sin(phase)])
        else:
            liste = complex_values(nomfich, INDIC_PARA, indic1, indic2, "INDIC_MODU" in args)

        # création de la fonction ASTER :
        ut_fonc = DEFI_FONCTION(
            NOM_PARA=NOM_PARA,
            NOM_RESU=NOM_RESU,
            PROL_DROITE=PROL_DROITE,
            PROL_GAUCHE=PROL_GAUCHE,
            INTERPOL=INTERPOL,
            INFO=INFO,
            TITRE=TITRE,
            VERIF=VERIF,
            VALE_C=liste,
        )

    elif TYPE == "NAPPE":
        # création de la nappe ASTER :
        motscles = {}
        motscles["DEFI_FONCTION"] = []
        for elem in mc_DEFI_FONCTION:
            values = function_values(
                FORMAT, nomfich, args["INDIC_ABSCISSE"], elem["INDIC_RESU"], SEPARATEUR, INFO
            )

            motscles["DEFI_FONCTION"].append(
                _F(
                    INTERPOL=args["INTERPOL_FONC"],
                    PROL_DROITE=args["PROL_DROITE_FONC"],
                    PROL_GAUCHE=args["PROL_GAUCHE_FONC"],
                    VALE=values,
                )
            )

        vale_para = column_values(FORMAT, nomfich, INDIC_PARA, SEPARATEUR, INFO)
        # création de la nappe
        ut_fonc = DEFI_NAPPE(
            PARA=vale_para,
            NOM_PARA=NOM_PARA,
            NOM_PARA_FONC=args["NOM_PARA_FONC"],
            NOM_RESU=NOM_RESU,
            PROL_DROITE=PROL_DROITE,
            PROL_GAUCHE=PROL_GAUCHE,
            INTERPOL=INTERPOL,
            INFO=INFO,
            TITRE=TITRE,
            VERIF=VERIF,
            **motscles
        )
    # remet UNITE dans son état initial
    ReservedUnitUsed(UNITE)
    return ut_fonc
