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
import numpy as np


def fsolve(func, x0):
    epsfcn = np.finfo(np.dtype(float)).eps
    eps = epsfcn**0.5
    delta_x = 1.0
    x = x0
    while abs(delta_x) > eps:
        f_x = func(x)
        d_x = max(eps, eps * abs(x))
        df_x = (func(x + d_x) - f_x) / d_x
        delta_x = -f_x / df_x
        x += delta_x
    return x


class PostRocheAnalytic:
    """classe pour calcul analytique des valeurs de post roche"""

    # loi ramberg osgood
    _K = 0.01
    _n = 2
    _E = 2.0e11
    # Z : module d'inertie pour les poutres circulaires : Iy / R
    _Z = 7.36311e-09 / 0.01
    # B_2 = 1. pour les tronçons de partie droite
    _B_2 = 1.0
    # coefficients de contraintes
    _D_21 = 0.712
    _D_22 = 0.429 * (0.02 / 0.005) ** 0.16
    _D_23 = _D_22

    def __init__(
        self,
        mfy_poids_propre,
        mt_deplacement_impose,
        mfy_sismique_dyn,
        mfz_sismique_dyn,
        mfy_sismique_qs,
        mfz_sismique_qs,
        mfy_sismique,
        mfz_sismique,
        pression=0,
        coude=False,
        RCCM=False,
    ):
        ep = 0.005
        rayon = 0.01
        self._sigpres = 0.87 * pression * (rayon - ep) / ep

        if coude:
            rcourb = 0.01
            angle_rad = np.pi / 2
            rp02min = 410e6

            ff = ep * rcourb / (rayon - ep / 2) ** 2
            self._B_2 = max(1.0, 1.3 / ff ** (2.0 / 3))
            self._D_22 = max(
                1.07
                * pow(np.pi / angle_rad, -0.4)
                * pow(ff, -2 / 3)
                / (1 + 0.142 * pression * (rayon - ep) / ep / rp02min / pow(ff, 1.45)),
                1.02,
            )
            self._D_23 = max(
                0.809
                * pow(ff, -0.44)
                / (1 + 0.142 * pression * (rayon - ep) / ep / rp02min / pow(ff, 1.45)),
                1.02,
            )
            self._sigpres2 = 0.87 * pression * (rayon - ep) / ep
        else:
            self._sigpres2 = 0.0

        self.RCCM = RCCM
        self._mfy_poids_propre = mfy_poids_propre
        self._mt_deplacement_impose = mt_deplacement_impose
        self._mfy_sismique_dyn = mfy_sismique_dyn
        self._mfz_sismique_dyn = mfz_sismique_dyn
        self._mfy_sismique_qs = mfy_sismique_qs
        self._mfz_sismique_qs = mfz_sismique_qs
        self._mfy_sismique = mfy_sismique
        self._mfz_sismique = mfz_sismique

    def _epsi_p(self, sigma):
        """epsi plastique pour Ranberg osgood"""
        return self._K * (sigma / self._E) ** (1 / self._n)

    def calcul_ressort(self):
        """calcule self._sigma_deplacement_ref, self._sigma_sismique_ref, self._t, self._ts, self._T, self._T_s, self._r_m et self._r_s"""
        # moments en deplacements à abattres
        mt_deplacement = np.abs(self._mt_deplacement_impose)
        mf_sismique = np.sqrt(self._mfy_sismique**2 + self._mfz_sismique**2)
        # sigma ref
        self._sigma_deplacement_ref = 0.87 * mt_deplacement / self._Z
        self._sigma_sismique_ref = 0.79 * self._B_2 * mf_sismique / self._Z
        # reversibilites locales
        sig_over_E = self._sigma_deplacement_ref / self._E
        sig_over_E = (sig_over_E > 1e-6) * sig_over_E
        self._t = sig_over_E
        ii = 0
        for s_E, s in zip(sig_over_E, self._sigma_deplacement_ref):
            if s_E != 0.0:
                self._t[ii] = s_E / self._epsi_p(s)
            ii += 1
        sig_over_E = self._sigma_sismique_ref / self._E
        sig_over_E = (sig_over_E > 1e-6) * sig_over_E
        self._t_s = sig_over_E
        ii = 0
        for s_E, s in zip(sig_over_E, self._sigma_sismique_ref):
            if s_E != 0.0:
                self._t_s[ii] = s_E / self._epsi_p(s)
            ii += 1
        # reversibilité totales

        if sum(self._sigma_deplacement_ref[self._t != 0] ** 2) == 0:
            self._T = 0.0
        else:
            self._T = sum(self._sigma_deplacement_ref[self._t != 0] ** 2) / sum(
                self._sigma_deplacement_ref[self._t != 0] ** 2 / self._t[self._t != 0]
            )
        if sum(self._sigma_sismique_ref[self._t_s != 0] ** 2) == 0:
            self._T_s = 0.0
        else:
            self._T_s = sum(self._sigma_sismique_ref[self._t_s != 0] ** 2) / sum(
                self._sigma_sismique_ref[self._t_s != 0] ** 2 / self._t_s[self._t_s != 0]
            )
        # ressorts
        self._r_m = np.maximum(
            np.array([self._T / t - 1 if t != 0.0 else 0.0 for t in self._t]), np.zeros([4])
        )
        self._r_s = np.array([self._T_s / t - 1 if t != 0.0 else 0.0 for t in self._t_s])
        self._r_mmax = np.array([np.amax(self._r_m)] * 4)
        self._r_smax = np.array([np.amax(self._r_s)] * 4)

        return self._r_m, self._r_s, self._r_mmax, self._r_smax

    def calcul_abattement(self):
        """calcule self._sigma_deplacement, self._sigma_sismique, self._g, self._g_s"""

        def func(x):
            """equation de contrainte vrai"""
            if x < 0:
                x = 0.0
            return (
                (x / self._E + self._epsi_p(x) - (sigma_ref + self._sigpres) / self._E)
                - self._epsi_p(self._sigpres)
                + r * (x - sigma_ref - self._sigpres) / self._E
            )

        # resolution sigma_deplacement
        sigma_deplacement = []
        for sigma_ref, r in zip(self._sigma_deplacement_ref, self._r_m):
            if sigma_ref == 0.0:
                root = 0.0
            else:
                root = fsolve(func, sigma_ref)
            sigma_deplacement.append(root)
        self._sigma_deplacement = np.array(sigma_deplacement)

        # pour g +opt
        sigma_deplacement_opt = []
        for sigma_ref, r in zip(self._sigma_deplacement_ref, self._r_mmax):
            if sigma_ref == 0.0:
                root = 0.0
            else:
                root = fsolve(func, sigma_ref)
            sigma_deplacement_opt.append(root)
        self._sigma_deplacement_opt = np.array(sigma_deplacement_opt)

        # resolution sigma_sismique
        sigma_sismique = []
        for sigma_ref, r in zip(self._sigma_sismique_ref, self._r_s):
            if sigma_ref == 0.0:
                root = 0.0
            else:
                root = fsolve(func, sigma_ref)
            sigma_sismique.append(root)
        self._sigma_sismique = np.array(sigma_sismique)

        # pour gsopt
        sigma_sismique_opt = []
        for sigma_ref, r in zip(self._sigma_sismique_ref, self._r_smax):
            if sigma_ref == 0.0:
                root = 0.0
            else:
                root = fsolve(func, sigma_ref)
            sigma_sismique_opt.append(root)
        self._sigma_sismique_opt = np.array(sigma_sismique_opt)

        # abbatements g
        self._g = []
        for sigma_ref, sigma_vrai in zip(self._sigma_deplacement_ref, self._sigma_deplacement):
            if sigma_vrai < self._sigpres or sigma_vrai == 0:
                g = 1.0
            else:
                g = (sigma_vrai - self._sigpres) / sigma_ref
            self._g.append(g)
        self._g = np.array(self._g)

        # abbatements gopt
        self._gopt = []
        for sigma_ref, sigma_vrai in zip(self._sigma_deplacement_ref, self._sigma_deplacement_opt):
            if sigma_vrai < self._sigpres or sigma_vrai == 0:
                g = 1.0
            else:
                g = (sigma_vrai - self._sigpres) / sigma_ref
            self._gopt.append(g)
        self._gopt = np.array(self._gopt)

        # abbatements gs
        self._g_s = []
        for sigma_ref, sigma_vrai in zip(self._sigma_sismique_ref, self._sigma_sismique):
            if sigma_vrai < self._sigpres or sigma_vrai == 0:
                g = 1.0
            else:
                g = (sigma_vrai - self._sigpres) / sigma_ref
            self._g_s.append(g)
        self._g_s = np.array(self._g_s)

        # abbatements gsopt
        self._g_sopt = []
        for sigma_ref, sigma_vrai in zip(self._sigma_sismique_ref, self._sigma_sismique_opt):
            if sigma_vrai < self._sigpres or sigma_vrai == 0:
                g = 1.0
            else:
                g = (sigma_vrai - self._sigpres) / sigma_ref
            self._g_sopt.append(g)
        self._g_sopt = np.array(self._g_sopt)

        return self._g, self._g_s, self._gopt, self._g_sopt

    def donne_abattement(self, gopt, g_sopt):
        """fournir directement les coefficients d'abattement"""
        self._g = 0.0
        self._g_s = 0
        self._gopt = gopt
        self._g_sopt = g_sopt

    def calcul_sigma_eq(self):
        """calcule le sigma equivalent"""
        self._sigma_eq = (
            self._sigpres2**2
            + 1
            / self._Z**2
            * (
                self._D_21**2 * (self._g * self._mt_deplacement_impose) ** 2
                + self._D_22**2
                * (
                    abs(self._mfy_poids_propre)
                    + abs(self._mfy_sismique_qs)
                    + self._g_s * abs(self._mfy_sismique_dyn)
                )
                ** 2
                + self._D_23**2
                * (abs(self._mfz_sismique_qs) + self._g_s * abs(self._mfz_sismique_dyn)) ** 2
            )
        ) ** 0.5

        self._sigma_eq_opt = (
            self._sigpres2**2
            + 1
            / self._Z**2
            * (
                self._D_21**2 * (self._gopt * self._mt_deplacement_impose) ** 2
                + self._D_22**2
                * (
                    abs(self._mfy_poids_propre)
                    + abs(self._mfy_sismique_qs)
                    + self._g_sopt * abs(self._mfy_sismique_dyn)
                )
                ** 2
                + self._D_23**2
                * (abs(self._mfz_sismique_qs) + self._g_sopt * abs(self._mfz_sismique_dyn)) ** 2
            )
        ) ** 0.5

        return self._sigma_eq, self._sigma_eq_opt
