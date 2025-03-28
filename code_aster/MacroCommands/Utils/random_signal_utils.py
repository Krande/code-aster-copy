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

# person_in_charge: irmela.zentner at edf.fr

# Routines for random signal generation
"""random_signal_utils.py

A collection of general-purpose routines using Numpy

DSP2ACCE1D        ---      generation of trajectories of a stationary Gaussian process
gene_traj_gauss_evol1D ---      generation of trajectories of a non stationary Gaussian process
calc_dsp_KT            ---      construct KT PSD
calc_dsp_FR            ---      construct rational PSD
acce_filtre_CP         ---       high pass filter for seismic signals
butterfilter             ---     Butterworth  lowpass filter for seismic signals
DSP2SRO                ---       identify PSD from given response spectrum
SRO2DSP                ---      determine response spectrum for a given PSD
Acce2SRO               ---      calculate response spectrum of a seismic signal
                                (accelerogram) by FFT
"""

from cmath import sqrt as csqrt
from math import ceil, cos, exp, log, pi, sqrt

import numpy as NP
import aster_fonctions
from ...Messages import UTMESS

from ...Objects.function_py import t_fonction
from .optimize import fmin


# -------------------------------------------------------------------
# ALGORITHME DE GENERATION DE SIGNAUX GAUSSIENS POUR LE CAS SCALAIRE
# --------------------------------------------------------------------
def DSP2ACCE1D(f_dsp, rv=None):
    # ----------------------------------
    # IN: f_dsp: dsp function for of list frequencies lw2 on (0, OM)
    #     rv: realisation du vecteur de variables aleatoires gaussiennes complexe
    # OUT: Xt: trajectoire du processus gaussien stationnaire normalise (m=0, ect=1)
    # ----------------------------------
    vale_dsp = f_dsp.vale_y
    lw2 = f_dsp.vale_x
    DW = lw2[1] - lw2[0]
    nbfreq2 = len(lw2)
    nbfreq = nbfreq2 * 2
    if rv is None:
        rv = NP.random.normal(0.0, 1.0, nbfreq) + 1j * NP.random.normal(0.0, 1.0, nbfreq)
    rv1 = rv[0:nbfreq2]
    rv2 = rv[nbfreq2:]
    CS = NP.array([0.0 + 0j] * nbfreq)
    for iifr in range(nbfreq2):
        vale_i = sqrt(vale_dsp[iifr])
        CS[nbfreq2 + iifr] = vale_i * rv1[iifr]
        CS[nbfreq2 - iifr - 1] = vale_i * rv2[iifr]
    SX = NP.fft.ifft(CS) * nbfreq
    ha = NP.exp(-1.0j * pi * NP.arange(nbfreq) * (1.0 - 1.0 / nbfreq))
    Xt = sqrt(DW) * (SX * ha).real
    return Xt


# --------------------------------------------------------------------------
#  ALGORITHME DE GENERATION DE SIGNAUX GAUSSIENS DSP evolutive non separable
# ---------------------------------------------------------------------------
def gene_traj_gauss_evol1D(self, rv=None, **kwargs):
    # ---------------------------------------------
    # IN: calc_dsp_KT: function for the definition of the PSD matrix (KT ou rational type)
    #      lw2: the list of frequencies corresponding to spec (0, OM)
    #      wg, wn: fond freq and evolution [rad/s],
    #      fcp [Hz]: corner frequency for Clough & Penzien filter
    # OUT: Xt: trajectoire du processus gaussien stationnaire normalise (m=0, ect=1)
    # ---------------------------------------------
    nbfreq2 = len(self.sampler.liste_w2)
    nbfreq = 2 * nbfreq2
    DW = self.sampler.DW
    Xt = []
    fg = kwargs["FREQ_FOND"]
    amo = kwargs["AMORT"]
    fp = kwargs["FREQ_PENTE"]
    TYPE = kwargs["TYPE_DSP"]
    if TYPE == "FR":
        R0 = kwargs["para_R0"]
        R2 = kwargs["para_R2"]
        l_FIT = kwargs["fonc_FIT"].vale_y
        assert len(l_FIT) == nbfreq2, "ERREUR listes frequences: emettre une fiche anomalie!"
        dsp_fr_refe = calc_dsp_FR(self.sampler.liste_w2, fg, amo, R0, R2, self.FREQ_CORNER)
        #   calcul de la variance (sigma^2) de normalisation mof
        if "ALEA_DSP" in kwargs:
            l_ALPHA = kwargs["ALEA_DSP"]
            mof = NP.trapz(dsp_fr_refe * l_FIT * l_ALPHA, self.sampler.liste_w2) * 2.0
            l_FIT = l_FIT * l_ALPHA
        else:
            mof = NP.trapz(dsp_fr_refe * l_FIT, self.sampler.liste_w2) * 2.0
    if rv is None:
        rv = NP.random.normal(0.0, 1.0, nbfreq) + 1j * NP.random.normal(0.0, 1.0, nbfreq)
    #      vecc1=(NP.random.normal(0.0,1.,nbfreq2)+1j*NP.random.normal(0.0,1.,nbfreq2))
    #      vecc2=(NP.random.normal(0.0,1.,nbfreq2)+1j*NP.random.normal(0.0,1.,nbfreq2))
    #   else:
    vecc1 = rv[0:nbfreq2]
    vecc2 = rv[nbfreq2:]
    t_mid = 0.5 * (self.modulator.T1 + self.modulator.T2)
    fg_fin = fg + fp * (self.modulator.T2 - t_mid)
    fg_ini = fg + fp * (self.modulator.T1 - t_mid)

    for tii in self.sampler.liste_temps:
        if tii < self.modulator.T1:
            fgt = fg_ini
        elif tii > self.modulator.T2:
            fgt = fg_fin
        else:
            fgt = fg + fp * (tii - t_mid)
        if fgt <= 0.0:
            UTMESS("F", "SEISME_35", valk=(str(tii)))
        # calcul du facteur de normalisation
        if TYPE == "KT":
            dsp = calc_dsp_KT(self, fgt, amo)
            # constante de normalisation pour que ecart_type=1 a pour tout t
            S_cst = 1.0 / NP.trapz(dsp, self.sampler.liste_w2) * 0.5
            MAT = calc_dsp_KT(self, fgt, amo, S_cst)
        elif TYPE == "FR":
            dsp = calc_dsp_FR(self.sampler.liste_w2, fgt, amo, R0, R2, self.FREQ_CORNER)
            # constante de normalisation pour que ecart_type=1 a pour tout t
            S_cst = mof / (NP.trapz(dsp * l_FIT, self.sampler.liste_w2) * 2.0)
            MAT = (
                calc_dsp_FR(self.sampler.liste_w2, fgt, amo, R0, R2, self.FREQ_CORNER, So=S_cst)
                * l_FIT
            )
        vale_xp = NP.sqrt(MAT) * vecc1 * NP.exp(1.0j * self.sampler.liste_w2 * tii)
        vale_xn = NP.sqrt(MAT) * vecc2 * NP.exp(-1.0j * self.sampler.liste_w2 * tii)
        vale_Xt = sum(vale_xp) + sum(vale_xn)
        Xt.append(vale_Xt.real * sqrt(DW))
    return Xt


# -----------------------------------------------------------------
#    filtre corner frequency wcp (modele Clough&Penzien)
# -----------------------------------------------------------------


def acce_filtre_CP(vale_acce, dt, fcorner, amoc=1.0):
    """Applies a high-pass filter to an accelerogram
    Args:
        vale_acce (ndarray or list): signal
        dt (float) : time step of signal
        fcorner (float) : the eigenfrequency (or corner frequency) of filter, Hz
        amoc (float): the damping ration of the filter
    Returns:
        ndarray : the filtered signal
    """
    # CP filter/corner frequency : wcp
    wcp = fcorner * 2.0 * pi
    N = len(vale_acce)
    ws = NP.fft.rfftfreq(N, d=dt) * 2 * pi

    acce_in = NP.fft.rfft(NP.array(vale_acce))

    hw2 = -(ws**2) * 1.0 / ((wcp**2 - ws**2) + 2.0 * amoc * 1j * wcp * ws)
    Yw = acce_in * hw2

    acce_out = NP.fft.irfft(Yw, n=N)
    #      f_out = t_fonction(vale_t, acce_out, para=f_in.para)
    return acce_out


# -----------------------------------------------------------------
#    filtre Butterworth passe bas
# -----------------------------------------------------------------
def butterfilter(cut_freq, dsp_input, order=10):
    # ---------------------------------------------------------
    # IN : acce: accelerogram m/s2 (optional)
    #      order : order of filter
    #      cut_freq: filter cut-off frequency if 'freq'
    #               need to give normalized cut_freq: cut_freq/fc if 'temp'
    # OUT: filtered time hist if domain = 'temp' or psd if "freq"
    #
    # ---------------------------------------------------------
    wc = 2.0 * pi * cut_freq
    w = dsp_input.vale_x
    f_h = 1.0 / NP.sqrt(1.0 + (w / wc) ** (2.0 * order))
    dsp_filt = dsp_input.vale_y * NP.abs(f_h) ** 2
    output = t_fonction(w, dsp_filt, para=dsp_input.para)
    return output


# ------------------------------------------------------------------------


def dsp_filtre_CP(f_in, fcorner, amoc=1.0):
    # ---------------------------------------------------------
    # IN : f_in: DSP (frequence rad/s),
    #         fcorner : corner frequency Hz,
    #         amoc: amortissement, l_freq: list of frequencies in Hz
    # OUT: f_out: DSP filtre (frequence rad/s),
    # attention: il faut de preference  2**N
    # ---------------------------------------------------------
    wcp = fcorner * 2.0 * pi
    vale_freq = f_in.vale_x
    vale_dsp = f_in.vale_y
    HW = 1.0 / ((wcp**2 - vale_freq**2) ** 2 + 4.0 * (amoc**2) * (wcp**2) * vale_freq**2)
    dsp_out = vale_freq**4 * vale_dsp * HW
    f_out = t_fonction(vale_freq, dsp_out, para=f_in.para)
    return f_out


# ------------------------------------------------------------------------


def calc_phase_delay(t, Xt, phase_data):
    # ---------------------------------------------------------
    # IN : 1D signal Xt, phase delay data dt =dx/Vs, coordinates of N points
    # OUT: ND signal with delay
    # ---------------------------------------------------------
    noe_interf = phase_data["NOEUDS_INTERF"]
    noe_refe = phase_data["COOR_REFE"]
    direction = phase_data["DIRECTION"]
    Vs = phase_data["VITE_ONDE"]

    #   if len(noe_interf[0])== 2:
    #        z = NP.zeros((len(noe_interf[:,0]),1))
    #        noe_interf = NP.append(noe_interf, z, axis=1)
    #        #noe_interf = interf3d
    dt = t[1] - t[0]
    DUREE = t[-1]
    X_out = []
    for line in noe_interf:
        Xtv = list(Xt)
        coord = NP.array(line) - noe_refe
        delay = 1.0 / Vs * NP.dot(NP.array(direction), coord)
        if delay < 0.0:
            UTMESS("F", "SEISME_78")
        if delay < dt:
            pass
        elif delay < DUREE:
            Nphase = int(delay / dt)
            for ii in range(Nphase):
                Xtv.insert(0, 0.0)
                Xtv.pop()
        else:
            Xtv = [0.0] * len(Xtv)
            print("ATTENTION : le délai de phase est supérieur à la durée du signal")

        X_out.append(Xtv)

    return X_out


# ------------------------------------------------------------------
#     PSD models
# -----------------------------------------------------------------
#
#     KANAI TAJIMI PSD
# -----------------------------------------------------------------


def calc_dsp_KT(self, freq_fond, amo, So=1.0):
    w0 = freq_fond * 2.0 * pi
    # KT model
    x11 = NP.array([4.0 * (amo**2) * (w0**2) * FREQ**2 for FREQ in self.sampler.liste_w2])
    xnum = x11 + w0**4
    denom = NP.array([(w0**2 - FREQ**2) ** 2 for FREQ in self.sampler.liste_w2])
    denom = denom + x11
    valkt = xnum / denom
    # CP filter
    amocp = 1.0
    valcp = NP.array(
        [
            FREQ**4
            / (
                4.0 * (amocp**2) * (2.0 * pi * self.FREQ_CORNER**2) * FREQ**2
                + (2.0 * pi * self.FREQ_CORNER**2 - FREQ**2) ** 2
            )
            for FREQ in self.sampler.liste_w2
        ]
    )
    dsp = valcp * valkt
    return dsp * So


#
#     FRACTION RATIONELLE
# -----------------------------------------------------------------
##
def calc_dsp_FR(lfreq, freq_fond, amor, R0, R1, FREQ_CORNER, So=1.0):
    # KT model parameters
    w0 = freq_fond * 2.0 * pi
    q0 = w0**2
    q1 = 2.0 * amor * w0
    valkt = NP.array(
        [
            (R0**2 + R1**2 * FREQ**2) / ((w0**2 - FREQ**2) ** 2 + q1**2 * FREQ**2)
            for FREQ in lfreq
        ]
    )
    # CP filter
    if FREQ_CORNER > 0.0:
        wcp = 2.0 * pi * FREQ_CORNER
        #      wcp=0.5*pi
        amocp = 1.0
        valcp = NP.array(
            [
                FREQ**4
                / (4.0 * (amocp**2) * (wcp**2) * FREQ**2 + (wcp**2 - FREQ**2) ** 2)
                for FREQ in lfreq
            ]
        )
        dsp = valcp * valkt * So
    else:
        dsp = valkt * So
    return dsp


# -----------------------------------------------------------------
#     ARIAS, duree phase forte TSM , T1 et T2 -----
# -----------------------------------------------------------------
def f_ARIAS(ta, acce, norme):
    acce2 = NP.array(acce) ** 2
    arias = NP.trapz(acce2, ta)  # energie
    arias = arias * pi / (2.0 * norme)  # indic Arias
    return arias


def f_ARIAS_TSM(ta, acce, norme):
    arias = f_ARIAS(ta, acce, norme)  # indic Arias
    ener = arias * (2.0 * norme) / pi
    acce2 = NP.array(acce) ** 2
    cumener = NP.array([NP.trapz(acce2[0 : ii + 1], ta[0 : ii + 1]) for ii in range(len(ta))])
    fract = cumener / ener
    n1 = NP.searchsorted(fract, 0.05)
    n2 = NP.searchsorted(fract, 0.95)
    #      n45= NP.searchsorted(fract,0.45)
    TSM = ta[n2] - ta[n1]
    T1 = ta[n1]
    T2 = ta[n2]
    return arias, TSM, T1, T2


def f_phase_forte(ta, acce, p1, p2):
    arias = f_ARIAS(ta, acce, 1.0)  # indic Arias
    ener = arias * (2.0 * 1.0) / pi
    acce2 = NP.array(acce) ** 2
    cumener = NP.array([NP.trapz(acce2[0 : ii + 1], ta[0 : ii + 1]) for ii in range(len(ta))])
    fract = cumener / ener
    n1 = NP.searchsorted(fract, p1)
    n2 = NP.searchsorted(fract, p2)
    return n1, n2


def f_ENER_qt(ta, acce, n1, n2):
    acce2 = acce**2
    ener = NP.trapz(acce2, ta)  # energie totale
    P1 = NP.trapz(acce2[0:n1], ta[0:n1]) / ener
    P2 = NP.trapz(acce2[0:n2], ta[0:n2]) / ener
    return ener, P1, P2


# -----------------------------------------------------------------
#     FONCTION DE MODULATION Gamma
# ----------------------------------------------------------------
# fonction de modulation gamma:  calcul pour liste de freq (normalisee si
# a1=1.0)


def fonctm_gam(ltemps, a1, a2, a3):
    qt = NP.array([a1 * tt ** (a2 - 1) * exp(-a3 * tt) for tt in ltemps])
    return qt


# fonction de modulation gamma: fonction cout pour identification des
# parametres


def f_opta(x0, ltemps, n1, n2):
    alpha = x0[0]
    beta = x0[1]
    if alpha <= 1.0:
        resu = 10.0
    elif beta < 0.0:
        resu = 1000.0
    else:
        qt = fonctm_gam(ltemps, 1.0, alpha, beta)
        ener, PINI, PFIN = f_ENER_qt(ltemps, qt, n1, n2)
        resu = sqrt((PINI - 0.05) ** 2 + (PFIN - 0.95) ** 2)
    return resu


# -----------------------------------------------------------------
#     FONCTION DE MODULATION Jennnings & Housner
# -----------------------------------------------------------------


def f_opt1(t1, ltemps, TS, a1, a2):
    T1 = t1[0]
    T2 = T1 + TS
    qt = fonctm_JetH(ltemps, T1, T2, a1, a2)
    n1 = NP.searchsorted(ltemps, T1)
    ener, PINI, PFIN = f_ENER_qt(ltemps, qt, n1, n1)
    residu = sqrt((PINI - 0.05) ** 2)
    return residu


def f_opt2(x0, ltemps, T1, TS):
    T2 = T1 + TS
    a1 = x0[0]
    a2 = x0[1]
    qt = fonctm_JetH(ltemps, T1, T2, a1, a2)
    n2 = NP.searchsorted(ltemps, T2)
    n1 = NP.searchsorted(ltemps, T1)
    ener, PINI, PFIN = f_ENER_qt(ltemps, qt, n1, n2)
    residu = sqrt((PFIN - 0.95) ** 2)
    return residu


# fonction de modulation Jennings & Housner normalisee


def fonctm_JetH(ltemps, T1, T2, a1, a2):
    qt = []
    for tt in ltemps:
        if tt < T1:
            qt.append((tt / T1) ** 2)
        elif tt < T2:
            qt.append(1.0)
        else:
            qt.append(exp(-a1 * (tt - T2) ** a2))
    return NP.array(qt)


# -----------------------------------------------------------------
#     FORMULES DE RICE
# -----------------------------------------------------------------


def Rice2(w2, DSP):
    #   MOMENTS
    m0 = NP.trapz(DSP, w2) * 2.0
    m1 = NP.trapz(DSP * abs(w2), w2) * 2.0
    m2 = NP.trapz(DSP * w2**2, w2) * 2.0
    #   FREQ_CENTRALE, BANDWIDTH
    vop = 1 / (2.0 * pi) * sqrt(m2 / m0)
    delta = sqrt(1.0 - m1**2 / (m0 * m2))
    return m0, m1, m2, vop, delta


# -----------------------------------------------------------------
#     FACTEUR DE PIC (VANMARCKE)
# -----------------------------------------------------------------
# calcul du facteur de peak par formule approche (Vanmarcke)


def peak(p, TSM, vop, amort):
    # ---------------------------------------------
    # IN: oscillator eigenfrequency  vop (Hz), reduced damping amort ,
    #     fractile p, duration TSM
    # OUT: peak factor
    # ---------------------------------------------
    omega0 = vop * 2.0 * pi
    deuxn = 2.0 * vop * TSM / (-log(p))
    if deuxn < 1.0:
        return 1.0
    else:
        xis = amort / (1.0 - exp(-2.0 * amort * omega0 * TSM))
        delta = sqrt(4.0 * xis / pi)
        sexp = -(delta**1.2) * sqrt(pi * log(deuxn))
        nup2 = 2.0 * log(deuxn * (1.0 - exp(sexp)))
        nup2 = max(1.0, nup2)
        return sqrt(nup2)


# calcul du facteur de peak par moments(formule Rice +Vanmarcke)


def peakm(p, TSM, w2, DSP):
    # ---------------------------------------------
    # IN   :  S(w) : DSP, w
    #          fractile p, duration TSM
    # OUT  :  peak factor
    # ---------------------------------------------
    m0 = NP.trapz(DSP, w2) * 2.0
    m1 = NP.trapz(DSP * abs(w2), w2) * 2.0
    m2 = NP.trapz(DSP * w2**2, w2) * 2.0
    vop = 1.0 / (2.0 * pi) * sqrt(m2 / m0)  # FREQ_CENTRALE
    delta = sqrt(1.0 - m1**2.0 / (m0 * m2))  # BANDWIDTH
    deuxn = 2.0 * vop * TSM / (-log(p))
    if deuxn < 1.0:
        return 1.0, m0
    else:
        sexp = -(delta**1.2) * sqrt(pi * log(deuxn))
        nup2 = 2.0 * log(deuxn * (1.0 - exp(sexp)))
        nup2 = max(1.0, nup2)
        return sqrt(nup2), m0


# -----------------------------------------------------------------
#     DSPSRO      SRO and DSP: functions of frequency (rad/s)
# -----------------------------------------------------------------
# conversion DSP en SRO par formule de Rice


def DSP2SRO(f_in, xig, TSM, liste_freq, ideb=2):
    # ---------------------------------------------
    # IN: f_in: DSP, function of frequency (rad/s),
    #     TSM: duree pase forte, xig: damping ratio
    #     liste_freq: list of freq SRO (Hz)
    # OUT: f_out: SRO, function of frequency (Hz), same norm. as DSP
    # ---------------------------------------------
    para_dsp = f_in.para
    vale_dsp_in = f_in.vale_y
    vale_sro = []
    vale_freq = f_in.vale_x
    vale_freq2 = vale_freq**2
    for f_0 in liste_freq:
        if f_0 == 0.0:
            vale_sro.append(0.0)
        else:
            w_0 = f_0 * 2.0 * pi
            vale_dsp_rep = vale_dsp_in / (
                (w_0**2 - vale_freq2) ** 2 + 4.0 * xig**2 * w_0**2 * vale_freq2
            )
            npeakm, m0i = peakm(0.5, TSM, vale_freq, vale_dsp_rep)
            vale_sro.append(w_0**ideb * npeakm * sqrt(m0i))
    f_out = t_fonction(liste_freq, vale_sro, para=para_dsp)
    return f_out


# -----------------------------------------------------------------
#     SRO2DSP     DSP: function of frequency (rad/s)
# -----------------------------------------------------------------
# iteration par formule de Rice pour mieux fitter le spectre cible


def iter_SRO(f_dsp, f_sro, amort, TS, Niter=10, nbliss=0):
    # ---------------------------------------------
    # IN  : f_in: DSP [rad/s], sro : spectre cible [Hz],
    #       amort: amortissement sro, TS: duree phase forte, meme disretisation
    # OUT : f_out: dsp apres iterations pour fitter au mieux le spectre sro
    # ---------------------------------------------
    para_dsp = f_dsp.para
    freq_dsp = f_dsp.vale_x
    vale_dsp = f_dsp.vale_y
    freq_sro = f_sro.vale_x
    vale_sro_ref = f_sro.vale_y
    nbvale = len(freq_dsp)
    ii = 0
    while ii < Niter:
        ii = ii + 1
        f_sroi = DSP2SRO(f_dsp, amort, TS, freq_sro)
        valesro = f_sroi.vale_y
        #  calcul de la correction des DSP
        nz = NP.nonzero(valesro)
        factm = NP.ones(nbvale)
        factm[nz] = vale_sro_ref[nz] / valesro[nz]
        vale_dspi = vale_dsp * factm**2

        if nbliss > 0:
            vale_dsp = smoothing(vale_dspi, nbliss)
        else:
            vale_dsp = vale_dspi

        f_dsp = t_fonction(freq_dsp, vale_dsp, para=para_dsp)
    f_out = f_dsp
    return f_out


def smoothing(yin, Mm):
    """Smoothes a function with a Hamming filter
    Args:
        yin (NP.ndarray):  the function to smooth
        Mm (int): the number of steps to apply the Hamming filter on
    Returns:
        NP.ndarray : the smoothed function
    """

    ysmoothed = NP.copy(yin)

    m = NP.arange(-Mm, Mm + 1, 1)

    filt = 0.54 - 0.46 * NP.cos(pi * (m + Mm) / Mm)
    fact = 1.0 / 1.08 / Mm

    for ii in range(Mm, len(yin) - Mm):
        nvale_m = m + ii
        ysmoothed[ii] = fact * NP.sum(yin[nvale_m[0] : nvale_m[-1] + 1] * filt)

    return ysmoothed


# iteration par simulation temporelle pour fitter le spectre cible sur une
# realisation (accelerogramme)


def itersim_SRO(self, FONC_DSP, NB_TIRAGE=1, **SRO_args):
    # (f_dsp, f_sro,nb_iter,f_modul, SRO_args ,dico_err,NB_TIRAGE=1 )
    # FONC_SPEC, AMORT, FMIN,  PAS=None, LIST_FREQ=None
    # ---------------------------------------------
    # IN  : FONC_DSP: DSP [rad/s], FONC_SPEC: spectre cible [Hz],
    #    amort: amortissement sro,  meme disretisation
    #    type_mod: type de fonction de modulation     niter: nombre d'iterations,
    #    FMIN: fequence min pour fit et filtrage ("corner frequency" Hz)
    # OUT : f_out: accelerogramme apres iterations pour fitter au mieux le spectre cible
    # ---------------------------------------------
    # parameters
    FMIN = SRO_args["FMIN"]
    amort = SRO_args["AMORT"]
    dico_err = SRO_args["DICO_ERR"]
    NB_ITER = SRO_args["NB_ITER"]
    # dsp initiale
    para_dsp = FONC_DSP.para
    freq_dsp = FONC_DSP.vale_x
    vale_dsp = FONC_DSP.vale_y
    nbfreq2 = len(freq_dsp)
    nbfreq = 2 * nbfreq2
    # sro cible
    freq_sro = freq_dsp / (2.0 * pi)
    vale_sro_ref = SRO_args["FONC_SPEC"].evalfonc(freq_sro).vale_y
    #  fonction de modulation
    hmod = self.modulator.fonc_modul.vale_y
    dt = self.sampler.DT

    FMINM = max(FMIN, 0.1)  #  FMIN pour le calcul de l'erreur relative
    FC = max(self.FREQ_FILTRE, FMINM)
    N1 = NP.searchsorted(freq_sro, FMINM) + 1
    FRED = freq_sro[N1:]
    ZPA = vale_sro_ref[-1]
    vpsum = sum([err_listes[0] for err_listes in list(dico_err.values())])
    coef_ZPA = dico_err["ERRE_ZPA"][0] / vpsum
    coef_MAX = dico_err["ERRE_MAX"][0] / vpsum
    coef_RMS = dico_err["ERRE_RMS"][0] / vpsum
    rv = NP.random.normal(0.0, 1.0, nbfreq) + 1j * NP.random.normal(0.0, 1.0, nbfreq)
    list_rv = [rv]
    if NB_TIRAGE > 1:
        ntir = 1
        while ntir < NB_TIRAGE:
            rv = NP.random.normal(0.0, 1.0, nbfreq) + 1j * NP.random.normal(0.0, 1.0, nbfreq)
            list_rv.append(rv)
            ntir = ntir + 1

    #  INITIALISATION
    errmult = []
    l_dsp = [FONC_DSP]

    if NB_TIRAGE == 1:
        acce = DSP2ACCE1D(FONC_DSP, rv) * hmod  # modulation
        if self.FREQ_FILTRE > 0.0:
            acce = acce_filtre_CP(acce, dt, self.FREQ_FILTRE)
        f_acce = t_fonction(self.sampler.liste_temps, acce, para=self.modulator.para_fonc_modul)
        l_acce = [f_acce]
        f_sroi = ACCE2SROM(self, f_acce, amort, freq_sro, 2, SRO_args["METHODE_SRO"])
        valesro = f_sroi.vale_y

    elif NB_TIRAGE > 1:
        liste_valesro = []
        for ntir in range(NB_TIRAGE):
            Xt = DSP2ACCE1D(FONC_DSP, list_rv[ntir])
            acce = Xt * hmod  # modulation
            if self.FREQ_FILTRE > 0.0:
                acce = acce_filtre_CP(acce, dt, self.FREQ_FILTRE)
            f_acce = t_fonction(self.sampler.liste_temps, acce, para=self.modulator.para_fonc_modul)
            f_sroi = ACCE2SROM(self, f_acce, amort, freq_sro, 2, SRO_args["METHODE_SRO"])
            liste_valesro.append(f_sroi.vale_y)
        if SRO_args["TYPE_ITER"] == "SPEC_MEDIANE":
            valesro = NP.median(NP.array(liste_valesro), axis=0)
        elif SRO_args["TYPE_ITER"] == "SPEC_MOYENNE":
            valesro = NP.mean(NP.array(liste_valesro), axis=0)

    l_sro = [valesro]
    err_zpa, err_max, err_min, err_rms, freq_err = erre_spectre(
        FRED, valesro[N1:], vale_sro_ref[N1:]
    )
    #  erreur multiobjectif
    err_ZPA = coef_ZPA * err_zpa
    err_MAX = coef_MAX * err_max
    err_RMS = coef_RMS * err_rms
    errmult.append(sqrt(1.0 / 3.0 * (err_ZPA**2 + err_MAX**2 + err_RMS**2)))
    if self.INFO == 2:
        UTMESS("I", "SEISME_43", valr=(err_zpa, err_max, err_rms, errmult[-1]))

    # ITERATIONS
    for kk in range(NB_ITER):
        #  CALCUL CORRECTION des DSP et mise a jour f_dsp
        nz = NP.nonzero(valesro)
        factm = NP.ones(nbfreq2)
        factm[nz] = vale_sro_ref[nz] / valesro[nz]
        vale_dspi = vale_dsp * factm**2
        vale_dsp = vale_dspi
        f_dsp = t_fonction(freq_dsp, vale_dsp, para=para_dsp)
        f_dsp = dsp_filtre_CP(f_dsp, FC)
        l_dsp.append(f_dsp)

        #  ITERATION DSP ACCE
        if NB_TIRAGE == 1:
            #  calcul accelerogramme et SRO
            Xt = DSP2ACCE1D(f_dsp, rv)
            acce = Xt * hmod  # modulation
            if self.FREQ_FILTRE > 0.0:
                acce = acce_filtre_CP(acce, dt, self.FREQ_FILTRE)
            f_acce = t_fonction(self.sampler.liste_temps, acce, para=self.modulator.para_fonc_modul)
            l_acce.append(f_acce)
            f_sroi = ACCE2SROM(self, f_acce, amort, freq_sro, 2, SRO_args["METHODE_SRO"])
            valesro = f_sroi.vale_y

        elif NB_TIRAGE > 1:
            liste_valesro = []
            for ntir in range(NB_TIRAGE):
                Xt = DSP2ACCE1D(f_dsp, list_rv[ntir])
                acce = Xt * hmod  # modulation
                if self.FREQ_FILTRE > 0.0:
                    acce = acce_filtre_CP(acce, dt, self.FREQ_FILTRE)
                f_acce = t_fonction(
                    self.sampler.liste_temps, acce, para=self.modulator.para_fonc_modul
                )
                f_sroi = ACCE2SROM(self, f_acce, amort, freq_sro, 2, SRO_args["METHODE_SRO"])
                liste_valesro.append(f_sroi.vale_y)
                if SRO_args["TYPE_ITER"] == "SPEC_MEDIANE":
                    valesro = NP.median(NP.array(liste_valesro), axis=0)
                elif SRO_args["TYPE_ITER"] == "SPEC_MOYENNE":
                    valesro = NP.mean(NP.array(liste_valesro), axis=0)

        #  CALCUL DES ERREURS
        l_sro.append(valesro)
        err_zpa, err_max, err_min, err_rms, freq_err = erre_spectre(
            FRED, valesro[N1:], vale_sro_ref[N1:]
        )
        #  erreur multiobjective
        err_ZPA = coef_ZPA * err_zpa
        err_MAX = coef_MAX * err_max
        err_RMS = coef_RMS * err_rms
        errmult.append(sqrt(1.0 / 3.0 * (err_ZPA**2 + err_MAX**2 + err_RMS**2)))
        if self.INFO == 2:
            UTMESS("I", "SEISME_42", vali=(kk + 1, NB_ITER), valr=errmult[-1])
    # OPTIMUM
    ind_opt = NP.argmin(NP.array(errmult))
    f_dsp_opt = l_dsp[ind_opt]
    valesro_opt = l_sro[ind_opt]
    err_zpa, err_max, err_min, err_rms, freq_err = erre_spectre(
        FRED, valesro_opt[N1:], vale_sro_ref[N1:]
    )
    dico_err["ERRE_ZPA"].append(err_zpa)
    dico_err["ERRE_MAX"].append(err_max)
    dico_err["ERRE_RMS"].append(err_rms)
    if self.INFO == 2:
        valargs = _F(
            vali=ind_opt,
            valr=(errmult[ind_opt], err_max, freq_err[0], err_min, freq_err[1], err_zpa, err_rms),
        )
        UTMESS("I", "SEISME_41", **valargs)
    for keys, listev in list(dico_err.items()):
        tole = listev[1] * 100.0
        erre = abs(listev[-1])
        if abs(erre) > tole:
            nbi = ind_opt
            UTMESS("A", "SEISME_36", vali=nbi, valk=keys, valr=(erre, tole))
    return f_dsp_opt, list_rv


# routines pour le calcul de l'erreur
# ---------------------------------------------
def erre_spectre(Freq, valesro, vale_sro_ref):
    errlin = (valesro - vale_sro_ref) / vale_sro_ref * 100.0
    errzpa = errlin[-1]
    errmax = max(abs(errlin))
    errmin = min(errlin)
    errms = sqrt(1.0 / len(Freq) * NP.sum(errlin**2))
    freqerr = [Freq[NP.argmax(abs(errlin))], Freq[NP.argmin((errlin))]]
    return errzpa, errmax, errmin, errms, freqerr


# conversion SRO en DSP equivalente par formule de Vanmarcke
# ---------------------------------------------
def SRO2DSP(
    FREQ_COUP,
    DUREE_PHASE_FORTE,
    FONC_SPEC,
    AMORT,
    FMIN,
    NORME,
    PAS=None,
    LIST_FREQ=None,
    NITER=10,
    FREQ_FILTRE_ZPA=0,
    NB_FREQ_LISS=0,
    **args
):
    # ---------------------------------------------
    #  f_in : SRO cible, frequency given in (Hz)
    #  f_out: DSP compatible avec SRO, frequency list lw in (rad/s)
    # ---------------------------------------------
    wmax = FREQ_COUP * 2.0 * pi
    f_in = FONC_SPEC
    fmin = max(FMIN, 0.01)
    wmin = fmin * 2.0 * pi
    #      wmin=1.001
    freq0 = 0.0
    freqi = freq0
    DSP = [0.0]
    lw = [freq0]
    lf = [freq0]
    lsro = [0.0]
    #      Sa_min=float(f_in.evalfonc([fmin]).vale_y*NORME)
    #      nupi=peak(0.5,  TSM, fmin ,  AMORT)
    #      DSP_min=Sa_min**2*2.*AMORT/(wmin*nupi**2)
    #      dsp_p= DSP_min/ wmin
    ii = 0
    while freqi < wmax:
        if PAS is not None:
            freqi = freqi + PAS * 2.0 * pi
        else:
            if ii < len(LIST_FREQ):
                freqi = LIST_FREQ[ii] * 2.0 * pi
                ii = ii + 1
            else:
                freqi = wmax
        if freqi <= wmin:
            assert freqi > 0.0
            fi = freqi / 2.0 / pi
            #            valsro=float(f_in.evalfonc([fi]).vale_y*NORME)
            #            lsro.append(valsro)
            lsro.append(0.0)
            #            valg = freqi*dsp_p
            DSP.append(0.0)
        else:
            fi = freqi / 2.0 / pi
            valsro = float(f_in.evalfonc([fi]).vale_y * NORME)
            lsro.append(valsro)
            nupi = peak(0.5, DUREE_PHASE_FORTE, fi, AMORT)
            nup2 = nupi**2
            v1 = 1.0 / (freqi * (pi / (2.0 * AMORT) - 2.0))
            v2 = (valsro**2) / nup2
            v3 = 2.0 * NP.trapz(NP.array(DSP), NP.array(lw))
            v4 = v1 * (v2 - v3)
            valg = max(v4, 0.0)
            DSP.append(valg)
        lw.append(freqi)
        lf.append(freqi / 2.0 / pi)

    para = dict(f_in.para)
    para.update({"INTERPOL": ["LIN", "LIN"]})
    f_out = t_fonction(lw, DSP, para=para)
    f_iter_sro_ref = t_fonction(lf, lsro, para=para)

    if FREQ_FILTRE_ZPA > 1.0:
        # filtrage passe_bas
        f_out = butterfilter(FREQ_FILTRE_ZPA, f_out)

    if NITER > 0:
        # iteration sans simulation: formule de rice
        # PSA for frequency list lw (rad/s), physical units (not g)!!
        f_dsp = iter_SRO(f_out, f_iter_sro_ref, AMORT, DUREE_PHASE_FORTE, NITER, NB_FREQ_LISS)
    else:
        f_dsp = f_out
    return f_dsp, f_iter_sro_ref


# -----------------------------------------------------------------
#     ACCE2SRO
# -----------------------------------------------------------------
# conversion ACCE en SRO par methode HARMO ou NIGAM
def ACCE2SROM(self, f_in, xig, l_freq, ideb, METHODE_SRO):
    if METHODE_SRO == "NIGAM":
        spectr = aster_fonctions.SPEC_OSCI(f_in.vale_x, f_in.vale_y, l_freq, [xig])
        vale_sro = spectr[0, ideb, :]
        f_out = t_fonction(l_freq, vale_sro, para=self.para_sro)
    elif METHODE_SRO == "HARMO":
        f_out = ACCE2SRO(f_in, xig, l_freq, ideb=2)
    else:
        print("ERROR METHODE SRO")
    return f_out


# conversion ACCE en SRO par fft et filtrage: METHODE_SRO=HARMO
def ACCE2SRO(f_in, xig, l_freq, ideb=2):
    """This function computes the response spectrum of an accelerogram
    Args:
        f_in(t_fonction):signal temporel (accelerogram)
        xig(float): damping ratio
        l_freq (list): list of frequencies for the response spectrum (Hz)
        ideb (int):
    Returns:
        t_fonction : the response spectrum function (as a function of freq in Hz)
    """
    #
    para_sro = {
        "INTERPOL": ["LIN", "LIN"],
        "NOM_PARA": "FREQ",
        "PROL_DROITE": "CONSTANT",
        "PROL_GAUCHE": "EXCLU",
        "NOM_RESU": "ACCE",
    }
    vale_t = f_in.vale_x
    vale_acce = f_in.vale_y
    N = len(vale_t)
    dt = vale_t[1] - vale_t[0]
    ws = NP.fft.rfftfreq(N, d=dt) * 2 * pi
    vale_sro = []
    acce_in = NP.fft.rfft(NP.array(vale_acce))
    for fi in l_freq:
        w_0 = fi * 2.0 * pi
        hw2 = 1.0 / ((w_0**2 - ws**2) + 2.0 * xig * 1j * w_0 * ws)
        Yw = acce_in * hw2
        acce_out = NP.fft.irfft(Yw, n=N)
        vale_sro.append(w_0**ideb * max(abs(acce_out)))
    f_out = t_fonction(l_freq, vale_sro, para=para_sro)
    return f_out


#
# -----------------------------------------------------------------
# DSP2FR
# -----------------------------------------------------------------
# Ajustement d'une DSP rationelle proche de KT


def DSP2FR(f_dsp_refe, FC):
    # ---------------------------------------------------------
    # IN : f_spec: SRO cible en fonction de la frequence en Hz
    #
    # OUT: f_out: DSP FR fonction de la frequence(rad/s)
    # ---------------------------------------------------------
    para_dsp = f_dsp_refe.para
    lfreq = f_dsp_refe.vale_x
    vale_dsp = f_dsp_refe.vale_y
    m0, m1, m2, vop, deltau = Rice2(lfreq, vale_dsp)
    # parametres initiales
    #      w0= vop*2.*pi
    xi0 = deltau ** (2.0 / 1.2) * pi / 4.0
    dsp_FR_ini = calc_dsp_FR(lfreq, vop, xi0, (vop * 2.0 * pi) ** 2, 4.0 * vop * pi * xi0, FC)
    const_ini = 2.0 * NP.trapz(dsp_FR_ini, lfreq)
    R0 = (vop * 2.0 * pi) ** 2 * sqrt(m0) / sqrt(const_ini)
    R2 = 4.0 * vop * pi * xi0 * sqrt(m0) / sqrt(const_ini)
    x0 = [R0, R2]
    para_opt = fmin(f_opt_FR1, x0, args=(f_dsp_refe, vop, xi0, FC))
    R0 = abs(para_opt[0])
    R2 = abs(para_opt[1])
    x0 = [vop, xi0]
    para_opt = fmin(f_opt_FR2, x0, args=(f_dsp_refe, R0, R2, FC))
    vop = para_opt[0]
    xi0 = para_opt[1]
    x0 = [R0, R2]
    para_opt = fmin(f_opt_FR1, x0, args=(f_dsp_refe, vop, xi0, FC))
    R0 = abs(para_opt[0])
    R2 = abs(para_opt[1])
    dsp_FR_fin = calc_dsp_FR(lfreq, vop, xi0, R0, R2, FC)
    FIT = NP.ones(len(lfreq))
    nz = NP.nonzero(dsp_FR_fin)
    FIT[nz] = vale_dsp[nz] / dsp_FR_fin[nz]
    f_fit = t_fonction(lfreq, FIT, para=para_dsp)
    return vop, xi0, R0, R2, f_fit


# ---------------------------------------------------------


def f_opt_FR1(para_ini, f_dsp_refe, vop, xi0, fcorner):
    R0 = para_ini[0]
    R2 = para_ini[1]
    lfreq = f_dsp_refe.vale_x
    sFR = calc_dsp_FR(lfreq, vop, xi0, R0, R2, fcorner, So=1.0)
    residu2 = NP.sum((sFR - f_dsp_refe.vale_y) ** 2)
    return sqrt(residu2)


# ---------------------------------------------------------


def f_opt_FR2(para_ini, f_dsp_refe, R0, R2, fcorner):
    vop = para_ini[0]
    xi0 = para_ini[1]
    lfreq = f_dsp_refe.vale_x
    sFR = calc_dsp_FR(lfreq, vop, xi0, R0, R2, fcorner, So=1.0)
    residu2 = NP.sum((sFR - f_dsp_refe.vale_y) ** 2)
    return sqrt(residu2)


# ---------------------------------------------------------
#
# -----------------------------------------------------------------
# RAND_DSP
# -----------------------------------------------------------------
# TIRAGE DSP ALEATOIRE AVEC LOI LOGNORMALE


def RAND_DSP(MAT_CHOL, Nbf, f_dsp):
    # ---------------------------------------------------------
    # IN : f_dsp: DSP mediane
    # MAT_CHOL  : chol(COV) pour la liste Periods
    # OUT: f_rand_dsp = f_dsp*rand_vec: realisation DSP aleatoire
    # ---------------------------------------------------------
    vale_dsp = f_dsp.vale_y
    freq_dsp = f_dsp.vale_x
    alpha2 = RAND_VEC(MAT_CHOL, Nbf, para=2.0)
    rand_dsp = vale_dsp * alpha2
    f_rand_dsp = t_fonction(freq_dsp, rand_dsp, para=f_dsp.para)
    return f_rand_dsp


#
#
def RAND_VEC(MAT_CHOL, Nbf, para=1.0):
    # ---------------------------------------------------------
    # IN : MAT_CHOL  : chol(COV) pour la liste Periods
    # OUT: alpha = vecteur aleatoire lognormal
    # ---------------------------------------------------------
    # on genere le vecteur Gaussien independ de moyenne 0 et COV=MAT_CHOL**2
    nbp = len(MAT_CHOL)
    rv = NP.random.normal(0.0, 1.0, nbp)
    rvec = NP.inner(MAT_CHOL, rv)
    if nbp < Nbf:  # il faut completer pour les tres basses frequences avec Period>10s
        nbv = Nbf - nbp
        vec0 = NP.ones(nbv) * rvec[0]
        rvec = NP.concatenate((vec0, rvec), axis=0)
    # on prend la variable lognormale de median 1 et sigma=beta,
    # on prend le carre car DSP: exp(rv)**2
    alpha = NP.exp(para * rvec)
    return alpha


# Coefficients de correlation (Baker)


def corrcoefmodel(Period, f_beta=None):
    # ---------------------------------------------------------
    # IN : Periods= liste des periodes 1/f  [s]
    #      optionnel: liste de beta (ecart-type)
    # OUT : mat_out= matrice de covariance pour periodes T
    #       >>>coef de correlation (ecart-type=1)  si beta=None
    #       >>>covariance si beta=tfonction
    #
    # REFERENCE     corrcoef selon le modele de Baker:
    #          Baker & Jayaram, Earthquake Spectra 24(1),299-317, 2008.
    #
    # ---------------------------------------------------------
    PMIN = min(Period)
    if PMIN < 0.01:
        UTMESS("F", "SEISME_37", valk=(str(1.0 / PMIN)))
    if max(Period) > 10.0:
        Periods = NP.extract(NP.array(Period) <= 10.0, NP.array(Period))
    else:
        Periods = Period
    nbT = len(Periods)
    Mat_Eps = NP.array([0.0] * nbT * nbT)
    Mat_Eps.resize(nbT, nbT)
    # Le modele de Baker est defini pour  max(Periods)<=10.

    if f_beta is not None:
        if min(f_beta.vale_x) > 0.1:
            UTMESS("F", "SEISME_82")
        else:
            f_beta = f_beta.evalfonc(1.0 / Periods)
            vale_beta = f_beta.vale_y

    for ii, Ti in enumerate(Periods):
        for jj, Tj in enumerate(Periods):
            Tmin = min(Ti, Tj)
            Tmax = max(Ti, Tj)
            C1 = 1.0 - cos(pi / 2.0 - 0.366 * log(Tmax / max(Tmin, 0.109)))
            C3 = C1
            if Tmax < 0.109:
                C2 = 1.0 - 0.105 * (1.0 - 1.0 / (1.0 + exp(100.0 * Tmax - 5.0))) * (
                    (Tmax - Tmin) / (Tmax - 0.0099)
                )
                Mat_Eps[ii, jj] = C2
            elif Tmin > 0.109:
                Mat_Eps[ii, jj] = C1
            elif Tmax < 0.2:
                C2 = 1.0 - 0.105 * (1 - 1 / (1 + exp(100 * Tmax - 5))) * (
                    (Tmax - Tmin) / (Tmax - 0.0099)
                )
                C4 = C1 + 0.5 * (sqrt(C3) - C3) * (1 + cos(pi * Tmin / 0.109))
                Mat_Eps[ii, jj] = min(C2, C4)
            else:
                C4 = C1 + 0.5 * (sqrt(C3) - C3) * (1.0 + cos(pi * Tmin / 0.109))
                Mat_Eps[ii, jj] = C4
            if f_beta is not None:
                Mat_Eps[ii, jj] = Mat_Eps[ii, jj] * vale_beta[ii] * vale_beta[jj]

    Mat_Gx = NP.linalg.cholesky(Mat_Eps)
    return Periods, Mat_Gx


#
# -----------------------------------------------------------------
# CORRECTION ZPA DES SIGNAUX
# -----------------------------------------------------------------
#
## Ces fonctions permettent de corriger les zpa des signaux acce


# create the Gaussian mask
def def_mask(signal: list, y00, epsilon: float):
    """Computes the Gaussian mask
    Args:
        signal : signal
        y00    : maximum requested
        epsilon : width of the Gaussian mask
    Returns:
        array: the Gaussian mask
    """
    t = signal[0]
    y = signal[1]
    t0_idx = NP.argmax(NP.abs(y))

    t0 = t[t0_idx]
    y0 = y[t0_idx] * NP.sign(y[t0_idx])
    mask = 1 - (1 - y00 / y0) * NP.exp(-0.5 * ((t - t0) / epsilon) ** 2)

    return mask


# correct the accelerogram to yield pga
def correct_signal(signal: list, pga: float, epsilon: float):
    """correct the signal by applying the Gaussian mask
    Args:
        signal : signal to modify
        pga    : target pga
        epsilon : width of the Gaussian mask
            (0.03 provides good results for seismic signals)
    Returns:
        list : the corrected signal
    """
    mask = def_mask(signal, pga, epsilon)

    sig = signal[1] * mask
    newsig = sig - NP.mean(sig)

    new_signal = (signal[0], newsig)

    return new_signal


# zpa match function
def zpa_match(signal: list, pga: float, epsilon: float = 0.03):
    """This function applies the zpa correction to an accelerogram,
        it corrects the signal to yield the target pga by applying a Gaussian mask.
        The correction function can be called several times
        when there are more than one exceedances.
    Args:
        signal : signal to modify
        pga    : target pga
        epsilon : width of the Gaussian mask
            (0.03 provides good results for seismic signals)
    Returns :
        list : the signal with corrected maximum
    """
    new_signal = correct_signal(signal, pga, epsilon)

    while NP.max(NP.abs(new_signal[1])) > pga * 1.001:
        new_signal = correct_signal(new_signal, pga, epsilon)

    return new_signal[1]
