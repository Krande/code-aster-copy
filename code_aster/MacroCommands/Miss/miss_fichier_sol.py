# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

import os

from math import log

import aster
from .miss_utils import calc_param_auto, dict_format, l_coor_sort, verif_sol_homogene
from ...Messages import UTMESS

sfmt = dict_format["sR"]


def fichier_sol(tab, struct, param=None):
    """Retourne le contenu du fichier de sol construit à partir de
    la table des données de sol et éventuellement des paramètres du calcul.
    """
    # vérification de la table
    for p in (
        "NUME_COUCHE",
        "NUME_MATE",
        "E",
        "NU",
        "RHO",
        "EPAIS",
        "AMOR_HYST",
        "RECEPTEUR",
        "SUBSTRATUM",
    ):
        if not p in tab.para:
            UTMESS("F", "TABLE0_2", valk=(p, "de sol"))
    nb_couche = len(tab)
    if max(tab.NUME_COUCHE.values()) != nb_couche:
        UTMESS("F", "MISS0_5")
    tsubstr = tab.SUBSTRATUM == "OUI"
    if len(tsubstr) != 1:
        UTMESS("F", "MISS0_3")

    # complète la table
    tsol = tab.copy()
    # ... niveau récepteur

    def f_recep(v):
        res = ""
        if v.strip() == "OUI":
            res = "RECEP"
        return res

    tsol.fromfunction("s_RECEP", f_recep, "RECEPTEUR")
    # ... niveau source

    def f_force(num, v):
        res = 0
        if v.strip() == "OUI":
            res = num
        return res

    tsol.fromfunction("s_FORCE", f_force, ("NUME_COUCHE", "SOURCE"))
    # ... êta
    tsol.fromfunction("ETA", lambda x: 0.0, "NUME_COUCHE")

    content = []
    # titre de la table
    content.append("TITRE")
    content.append(tsol.titr)
    # materiaux
    tsol.sort(CLES=["NUME_MATE"])
    nb_mate = max(tsol.NUME_MATE.values())
    content.append("MATERIAU %8d" % nb_mate)
    content.append("RO           E            NU           BETA         ETA")
    format = "%%(RHO)%(R)s %%(E)%(R)s %%(NU)%(R)s %%(AMOR_HYST)%(R)s %%(ETA)%(R)s" % dict_format
    last_id_mate = 0
    for row in tsol:
        if row["NUME_MATE"] == last_id_mate:  # déjà vu, on saute
            continue
        last_id_mate = row["NUME_MATE"]
        content.append(format % row)
    # couches
    tsol.sort(CLES=["NUME_COUCHE"])
    content.append("COUCHE %8d" % (nb_couche - 1))
    format = "%%(EPAIS)%(R)s MATE %%(NUME_MATE)8d %%(s_RECEP)s" % dict_format
    for ic, row in enumerate(tsol):
        if ic == nb_couche - 1:
            continue
        content.append(format % row)
    # substratum
    substr = tsubstr.rows[0]
    content.append("SUBS   MATE %8d" % substr["NUME_MATE"])
    # sources
    nb_source = len(tsol.SOURCE == "OUI")
    content.append("SOURCE %8d 3D" % nb_source)
    # forces
    format = "FORCE HORIZ POSI %(s_FORCE)8d"
    if (
        param
        and param.get("TYPE_RESU") == "FICHIER_TEMPS"
        and param.get("FICHIER_SOL_INCI") is None
    ):
        # champ incident au substratum
        content.append("FORCE PLANE POSI   %8d" % nb_couche)
    else:
        for ic, row in enumerate(tsol):
            if row["s_FORCE"] != 0:
                content.append(format % row)
    # complément des paramètres du calcul
    spec_max_auto = False
    # Gestion automatique de certains parametres
    if param and param.get("AUTO") == "OUI" and param.get("_auto_first_LT") is None:
        l_coor_nodes = struct.noeud_coor
        l_coor_x, l_coor_y, l_coor_z = l_coor_sort(l_coor_nodes)
        surf = param["SURF"]
        coef_offset = param["COEF_OFFSET"]
        dref_auto, rfic_auto, offset_max_auto, offset_nb_auto = calc_param_auto(
            l_coor_x, l_coor_y, l_coor_z, surf, coef_offset
        )
        if param["OFFSET_MAX"] is None:
            param["OFFSET_MAX"] = offset_max_auto
        else:
            UTMESS("A", "MISS0_43", valk="OFFSET_MAX")
        if param["OFFSET_NB"] is None:
            param["OFFSET_NB"] = offset_nb_auto
        else:
            UTMESS("A", "MISS0_43", valk="OFFSET_NB")
        if param["OPTION_DREF"] == "OUI":
            if param["DREF"] is None:
                param["DREF"] = dref_auto
            else:
                UTMESS("A", "MISS0_43", valk="DREF")
        if param["OPTION_RFIC"] == "OUI":
            if param["RFIC"] is None:
                param["RFIC"] = rfic_auto
            else:
                UTMESS("A", "MISS0_43", valk="RFIC")
        else:
            param["RFIC"] = 0.0
        #
        print("Mode automatique :")
        print("-    OFFSET_MAX auto = ", param["OFFSET_MAX"])
        print("-    OFFSET_NB auto = ", param["OFFSET_NB"])
        print("-    DREF auto = ", param["DREF"])
        print("-    RFIC auto = ", param["RFIC"])
        # Pour Laplace-temps, on ne recalcule les parametre qu'une fois
        if param["_auto_first_LT"] is None:
            param["_auto_first_LT"] = False

    if param["SPEC_MAX"] is None and param["SPEC_NB"] is None:
        sol_homo, vs, vp, betamin = verif_sol_homogene(tab)
        if param["SURF"] == "OUI" and sol_homo:
            spec_max_auto = True
            spec_max = 15.0 / vs
            dp = betamin / vp
            k_nb = int(log(spec_max / dp + 1.0) / log(2.0))
            k_nb = k_nb + 1
            spec_nb = 2 ** k_nb
        if spec_max_auto:
            param["SPEC_MAX"] = spec_max
            param["SPEC_NB"] = spec_nb

    #
    if (param and param.get("OFFSET_MAX")) or (param and param.get("AUTO") == "OUI"):
        # ALGO
        if param and param.get("ALGO"):
            content.append("ALGO %s" % param["ALGO"])
        elif param["SURF"] == "OUI":
            content.append("ALGO DEPL")
        else:
            content.append("ALGO REGU")
        # DREF / SPEC / OFFSET
        if param["DREF"]:
            content.append(("DREF" + sfmt) % param["DREF"])
        if param["SPEC_MAX"] or param["SPEC_NB"]:
            content.append(("SPEC" + sfmt + " / %d") % (param["SPEC_MAX"], param["SPEC_NB"]))
        else:
            content.append("SPEC AUTO")
        content.append(("OFFSET" + sfmt + " / %d") % (param["OFFSET_MAX"], param["OFFSET_NB"]))

    content.append("FIND")

    # terminer le fichier par un retour chariot
    content.append("")
    return os.linesep.join(content)
