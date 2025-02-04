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

"""
    Maquette calcul de fonctions de coherence

"""

from math import exp, log, pi, sqrt, tanh

import numpy as np


def calc_cohefromdata(acce1, acce2, dt, Mm):
    # IN : acce1,acce2: accelerograms at station1 and stationa 2 (matrices: nbsignal x nbstep),
    # dt: time step    Mm : points smoothing window
    # OUT : coherency function gam(wm)
    acce_size = np.array(acce1).shape
    nbval = acce_size[0]
    Nf = acce_size[1]
    TT = Nf * dt  #%duree du signal  (5s pour N=1000)
    OM = pi / dt  #% frequence de coupure
    dw = 2.0 * OM / Nf  # % pas de frequence
    #% discretisation freq
    w = np.arange(-OM + 0.5 * dw, OM, dw)
    #%  smoothening parameters and filter
    m = np.arange(-Mm, Mm + 1, 1)
    wm = w[Mm : Nf - Mm]
    N2 = len(wm)
    filt = 0.54 - 0.46 * np.cos(pi * (m + Mm) / Mm)
    dt2 = pi / (w[Nf - Mm - 1] + 0.5 * dw)
    t2 = np.arange(0.0, N2 * dt2, dt2)

    # FFT
    Xw1 = np.fft.fft(np.array(acce1), Nf, 1) * dt
    Xw2 = np.fft.fft(np.array(acce2), Nf, 1) * dt

    SX1 = []
    SX2 = []
    SX12 = []
    #%spectral estimate
    fact = 1.0 / (2.0 * pi) / TT
    for iii in range(nbval):
        SX12.append(fact * Xw1[iii] * np.conj(Xw2[iii]))
        SX1.append(np.real(fact * Xw1[iii] * np.conj(Xw1[iii])))
        SX2.append(np.real(fact * Xw2[iii] * np.conj(Xw2[iii])))

    fact = 1.0 / 1.08 / Mm
    Sxf1 = []
    for line in SX1:
        term = np.fft.fftshift(line)
        out = []
        for ii in range(Mm, Nf - Mm):
            nvale_m = m + ii
            out.append(fact * np.sum(term[nvale_m[0] : nvale_m[-1] + 1] * filt))
        Sxf1.append(out)

    Sxf2 = []
    for line in SX2:
        term = np.fft.fftshift(line)
        out = []
        for ii in range(Mm, Nf - Mm):
            nvale_m = m + ii
            out.append(fact * np.sum(term[nvale_m[0] : nvale_m[-1] + 1] * filt))
        Sxf2.append(out)

    Sxf12 = []
    for line in SX12:
        term = np.fft.fftshift(line)
        out = []
        for ii in range(Mm, Nf - Mm):
            nvale_m = m + ii
            out.append(fact * np.sum(term[nvale_m[0] : nvale_m[-1] + 1] * filt))
        Sxf12.append(out)

    # coherency
    gamw = []
    for i2 in range(nbval):
        gami2 = Sxf12[i2] / np.sqrt(np.array(Sxf1[i2]) * np.array(Sxf2[i2]))
        gamw.append(gami2)
    return wm / 2.0 / pi, (np.mean(gamw, 0))


# -------------------------------------------------------------------
# COHERENCY MATRIX
# --------------------------------------------------------------------
def CALC_COHE(freq, **kwargs):
    #    Frequency is in rad/s: freq= f*2*pi
    #    kwargs: VITE_ONDE, PARA_ALPHA, TYPE, MAILLAGE,
    model = kwargs["TYPE"]
    nom_mail = kwargs["MAILLAGE"]
    nom_group_inter = kwargs["GROUP_NO_INTERF"]
    if "NOEUDS_INTERF" in kwargs:
        noe_interf = kwargs["NOEUDS_INTERF"]
    else:
        liste_nom, noe_interf = get_group_nom_coord(nom_group_inter, nom_mail)
    if "DIST" in kwargs:
        DIST2 = kwargs["DIST"]
    else:
        DIST2 = calc_dist2(noe_interf)
    # # ----MITA & LUCO
    if model == "MITA_LUCO":
        # PARAMETRES fonction de coherence
        VITE_ONDE = kwargs["VITE_ONDE"]
        alpha = kwargs["PARA_ALPHA"]
        COHE = NP.exp(-(DIST2 * (alpha * freq / VITE_ONDE) ** 2.0))
    # ----ABRAHAMSON ROCK (EPRI)
    elif model == "ABRAHAMSON":
        p_a1 = 1.647
        p_a2 = 1.01
        p_a3 = 0.4
        p_n1 = 7.02
        nbno = len(noe_interf)
        freqk = freq / (2.0 * pi)
        COHE = NP.zeros((nbno, nbno))
        for no1 in range(nbno):
            for no2 in range(nbno):
                #                dist_xi = sqrt((XX[no1] - XX[no2])**2 + (YY[no1] - YY[no2])**2)
                dist_xi = sqrt(DIST2[no1, no2])
                p_n2 = 5.1 - 0.51 * log(dist_xi + 10.0)
                pfc = -1.886 + 2.221 * log(4000.0 / (dist_xi + 1.0) + 1.5)
                term1 = 1.0 + (freqk * tanh(p_a3 * dist_xi) / (p_a1 * pfc)) ** p_n1
                term2 = 1.0 + (freqk * tanh(p_a3 * dist_xi) / (p_a2 * pfc)) ** p_n2
                COHE[no1, no2] = 1.0 / sqrt(term1 * term2)
    elif model == "ABRA_ROCHER":
        p_a1 = 1.0
        p_a2 = 40.0
        p_a3 = 0.4
        p_n2 = 16.4
        nbno = len(noe_interf)
        freqk = freq / (2.0 * pi)
        COHE = NP.zeros((nbno, nbno))
        for no1 in range(nbno):
            for no2 in range(nbno):
                dist_xi = sqrt(DIST2[no1, no2])
                p_n1 = 3.8 - 0.04 * log(dist_xi + 1.0) + 0.0105 * (log(dist_xi + 1.0) - 3.6) ** 2
                pfc = 27.9 - 4.82 * log(dist_xi + 1.0) + 1.24 * (log(dist_xi + 1.0) - 3.6) ** 2
                term1 = 1.0 + (freqk * tanh(p_a3 * dist_xi) / (p_a1 * pfc)) ** p_n1
                term2 = 1.0 + (freqk * tanh(p_a3 * dist_xi) / (p_a2)) ** p_n2
                COHE[no1, no2] = 1.0 / sqrt(term1 * term2)
    elif model == "ABRA_SOLMOYEN":
        p_a1 = 1.0
        p_a3 = 0.4
        p_n1 = 3.0
        p_n2 = 15.0
        nbno = len(noe_interf)
        freqk = freq / (2.0 * pi)
        COHE = NP.zeros((nbno, nbno))
        for no1 in range(nbno):
            for no2 in range(nbno):
                dist_xi = sqrt(DIST2[no1, no2])
                p_a2 = 15.8 - 0.044 * dist_xi
                pfc = 14.3 - 2.35 * log(dist_xi + 1.0)
                term1 = 1.0 + (freqk * tanh(p_a3 * dist_xi) / (p_a1 * pfc)) ** p_n1
                term2 = 1.0 + (freqk * tanh(p_a3 * dist_xi) / (p_a2)) ** p_n2
                COHE[no1, no2] = 1.0 / sqrt(term1 * term2)
    return COHE
