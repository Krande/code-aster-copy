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


def EXTR_MATR(matrrr, vall):
    from code_aster.SD.sd_stoc_morse import sd_stoc_morse
    from code_aster.SD.sd_nume_equa import sd_nume_equa

    # fonction permettant de recuperer les matrices assemblees au format numpy
    # attention a l'espace memoire
    import numpy as NP
    import aster

    # construction des vecteurs jeveux
    nom = matrrr.sdj.REFA.get()
    typm = nom[8]

    var = matrrr.sdj.VALM.get()
    nume = matrrr.getDOFNumbering().getName().strip()
    smos = sd_stoc_morse("{0:<14s}.SMOS".format(nume))
    neq = sd_nume_equa("{0:<14s}.NUME".format(nume))

    adia = smos.SMDI.get()
    numl = smos.SMHC.get()
    rtt = neq.DELG.get()
    rtt2 = neq.DEEQ.get()
    lili = neq.LILI.get()
    refn = neq.REFN.get()

    valr = var[1]
    if typm[0:2] == "MS":
        mult = 1
    if typm[0:2] == "MR":
        mult = -1
        valr2 = var[2]

    vc = len(rtt)
    ddltot = 0
    nnddl = [0]
    tryy = 0
    # numerotation des ddl de la matrice extraite

    nomail = refn[0]
    nomodel = lili[1]

    nomail = nomail[0:8]
    nomodel = nomodel[0:8]

    # calc nbre de ddl au total et vect contenant le nbre de ddl pr chaque noeud (ex : barre =3, poutre=6...)
    for jj in range(0, vc):
        if (rtt2[2 * jj + 1] == 1) and (tryy == 1):
            if tryy == 1:
                nnddl = nnddl + [0]
                tryy = 0
        if rtt2[2 * jj + 1] > 0:
            ddltot = ddltot + 1
            nnddl[len(nnddl) - 1] = nnddl[len(nnddl) - 1] + 1
            tryy = 1
    nnddl = nnddl + [0]
    ddlag = (vc - ddltot) // 2
    ddlphy = ddltot - ddlag
    indint = [None] * ddlphy

    # contruction du vect contenant les indices des ddls non lagr et non nuls
    gg = 0
    tt = 0
    pp = 0
    tt = 0
    limm = nnddl[pp]
    vecc = [1] * limm
    ste = limm
    vectddl = [""] * ddlphy
    refddl = ["DX", "DY", "DZ", "DRX", "DRY", "DRZ"]
    for ii in range(0, vc):
        if rtt2[2 * ii + 1] < 0:
            if rtt[ii] == -1:
                varr = -1 * rtt2[2 * ii + 1] - 1
                vecc[varr] = 0
            ste = limm
        if (rtt2[2 * ii + 1] > 0) and (ste == limm):
            for hh in range(0, limm):
                if vecc[hh] == 1:
                    indint[gg] = ii + hh
                    vectddl[gg] = refddl[hh]
                    gg = gg + 1
                vecc[hh] = 1
            ste = 0
            pp = pp + 1
            limm = nnddl[pp - 1]
            vecc = [1] * limm
        if (rtt2[2 * ii + 1] > 0) and (ste < limm):
            ste = ste + 1

    # constrcution de la matrice complete ddl lag compris (cette etape peut etre fusionnee avec la suivante mais ca devient complique)
    lon = len(adia)
    if vall == 1:
        valeu = NP.zeros([lon, lon])
    if vall == 2:
        valeu = NP.zeros([lon, lon], complex)
    opi = 0
    indc = 0
    indl = 0
    for ii in range(0, lon):
        if ii == 0:
            opi = adia[ii]
        else:
            opi = adia[ii] - adia[ii - 1]
        for jj in range(0, opi):
            if ii == 0:
                indc = jj
                indd = numl[indc] - 1
            else:
                indc = adia[ii - 1] + jj
                indd = numl[indc] - 1
            valeu[indd][ii] = valr[indc]
            if typm[0:2] == "MS":
                valeu[ii][indd] = valr[indc]
            if typm[0:2] == "MR":
                valeu[ii][indd] = valr2[indc]
    matrig = valeu

    # elimination des ddl de lagrange et des ddl nuls
    if vall == 1:
        rig11 = NP.zeros([ddlphy, ddlphy])
    if vall == 2:
        rig11 = NP.zeros([ddlphy, ddlphy], complex)
    for ii in range(0, ddlphy):
        for jj in range(0, ddlphy):
            rig11[ii][jj] = matrig[indint[ii]][indint[jj]]

    return [vectddl, rig11]
