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


def crea_max_array(arr):
    _list = []
    for val in arr:
        _list.append(max(arr))
    _arr_resu = np.array(_list)
    return _arr_resu


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
        mfy_dilatation_thermique,
        mt_sismique_dyn,
        mt_sismique_qs,
        mt_sismique,
        mfy_sismique,
        mfz_sismique,
        pression=0.0,
    ):

        ep = 0.005
        rayon = 0.01
        self._sigpres = 0.87 * pression * (rayon - ep) / ep
        self._mfy_poids_propre = mfy_poids_propre
        self._mt_deplacement_impose = mt_deplacement_impose
        self._mfy_sismique_dyn = mfy_sismique_dyn
        self._mfz_sismique_dyn = mfz_sismique_dyn
        self._mfy_sismique_qs = mfy_sismique_qs
        self._mfz_sismique_qs = mfz_sismique_qs

        self._mfy_dilatation_thermique = mfy_dilatation_thermique
        self._mt_sismique_dyn = mt_sismique_dyn
        self._mt_sismique_qs = mt_sismique_qs

        self._mt_sismique = mt_sismique
        self._mfy_sismique = mfy_sismique
        self._mfz_sismique = mfz_sismique

    def _epsi_p(self, sigma):
        """epsi plastique pour Ranberg osgood"""
        return self._K * (sigma / self._E) ** (1 / self._n)

    def calcul_ressort(self):
        """calcule self._sigma_deplacement_ref, self._sigma_sismique_ref, self._t, self._ts, self._T, self._T_s, self._r_m et self._r_s"""
        # moments en deplacements à abattres
        mt_deplacement = np.abs(self._mt_deplacement_impose)
        mf_delpacement = np.abs(self._mfy_dilatation_thermique)
        mf_sismique = np.sqrt(self._mfy_sismique**2 + self._mfz_sismique**2)
        # sigma ref
        self._sigma_deplacement_ref = (
            (0.87 * mt_deplacement / self._Z) ** 2
            + (0.79 * self._B_2 * mf_delpacement / self._Z) ** 2
        ) ** 0.5
        self._sigma_sismique_ref = (
            (0.87 * self._mt_sismique / self._Z) ** 2
            + (0.79 * self._B_2 * mf_sismique / self._Z) ** 2
        ) ** 0.5
        # reversibilites locales
        sig_over_E = self._sigma_deplacement_ref / self._E
        sig_over_E = (sig_over_E > 1e-6) * sig_over_E
        if np.count_nonzero(sig_over_E) == 0:
            self._t = np.zeros([6])
        else:
            self._t = sig_over_E / self._epsi_p(self._sigma_deplacement_ref)
        sig_over_E = self._sigma_sismique_ref / self._E
        sig_over_E = (sig_over_E > 1e-6) * sig_over_E
        if np.count_nonzero(sig_over_E) == 0:
            self._t_s = np.zeros([6])
        else:
            self._t_s = sig_over_E / self._epsi_p(self._sigma_sismique_ref)

        if np.count_nonzero(self._t) == 0:
            # reversibilité totales
            T_1 = 0.0
            T_2 = 0.0
            self._T = np.array([T_1] * 4 + [T_2] * 2)
            # ressorts
            self._r_m = np.zeros([6])
            r_mmax1 = np.amax(self._r_m[:4])
            r_mmax2 = np.amax(self._r_m[4:])
            self._r_mmax = np.array([r_mmax1] * 4 + [r_mmax2] * 2)
        else:

            # reversibilité totales
            T_1 = sum(self._sigma_deplacement_ref[:4][self._t[:4] != 0] ** 2) / sum(
                self._sigma_deplacement_ref[:4][self._t[:4] != 0] ** 2
                / self._t[:4][self._t[:4] != 0]
            )
            T_2 = sum(self._sigma_deplacement_ref[4:][self._t[4:] != 0] ** 2) / max(
                sum(
                    self._sigma_deplacement_ref[4:][self._t[4:] != 0] ** 2
                    / self._t[4:][self._t[4:] != 0]
                ),
                1.0e-6,
            )
            self._T = np.array([T_1] * 4 + [T_2] * 2)
            self._r_m = np.maximum(
                np.array([T / t - 1 if t != 0.0 else 0.0 for t, T in zip(self._t, self._T)]),
                np.zeros([6]),
            )

            r_mmax1 = np.amax(self._r_m[:4])
            r_mmax2 = np.amax(self._r_m[4:])
            self._r_mmax = np.array([r_mmax1] * 4 + [r_mmax2] * 2)

        if np.count_nonzero(self._t_s) == 0:
            T_s_1 = 0.0
            T_s_2 = 0.0
            self._T_s = np.array([T_s_1] * 4 + [T_s_2] * 2)
            # ressorts
            self._r_s = np.zeros([6])

            r_smax1 = np.amax(self._r_s[:4])
            r_smax2 = np.amax(self._r_s[4:])
            self._r_smax = np.array([r_smax1] * 4 + [r_smax1] * 2)

        else:

            # reversibilité totales
            T_s_1 = sum(self._sigma_sismique_ref[:4][self._t_s[:4] != 0] ** 2) / sum(
                self._sigma_sismique_ref[:4][self._t_s[:4] != 0] ** 2
                / self._t_s[:4][self._t_s[:4] != 0]
            )
            T_s_2 = sum(self._sigma_sismique_ref[4:][self._t_s[4:] != 0] ** 2) / max(
                sum(
                    self._sigma_sismique_ref[4:][self._t_s[4:] != 0] ** 2
                    / self._t_s[4:][self._t_s[4:] != 0]
                ),
                1.0e-6,
            )
            self._T_s = np.array([T_s_1] * 4 + [T_s_2] * 2)
            # ressorts
            self._r_s = np.array(
                [T / t - 1 if t != 0.0 else 0.0 for t, T in zip(self._t_s, self._T_s)]
            )

            r_smax1 = np.amax(self._r_s[:4])
            r_smax2 = np.amax(self._r_s[4:])
            self._r_smax = np.array([r_smax1] * 4 + [r_smax1] * 2)

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

    def calcul_sigma_eq(self):
        """calcule le sigma equivalent"""

        self._sigma_eq = (
            1
            / self._Z
            * (
                self._D_21**2
                * (
                    abs(self._mt_sismique_qs)
                    + self._g * abs(self._mt_deplacement_impose)
                    + self._g_s * abs(self._mt_sismique_dyn)
                )
                ** 2
                + self._D_22**2
                * (
                    abs(self._mfy_poids_propre)
                    + abs(self._mfy_sismique_qs)
                    + self._g * abs(self._mfy_dilatation_thermique)
                    + self._g_s * abs(self._mfy_sismique_dyn)
                )
                ** 2
                + self._D_23**2
                * (abs(self._mfz_sismique_qs) + self._g_s * abs(self._mfz_sismique_dyn)) ** 2
            )
            ** 0.5
        )

        self._sigma_eq_opt = (
            1
            / self._Z
            * (
                self._D_21**2
                * (
                    abs(self._mt_sismique_qs)
                    + self._gopt * abs(self._mt_deplacement_impose)
                    + self._g_sopt * abs(self._mt_sismique_dyn)
                )
                ** 2
                + self._D_22**2
                * (
                    abs(self._mfy_poids_propre)
                    + abs(self._mfy_sismique_qs)
                    + self._gopt * abs(self._mfy_dilatation_thermique)
                    + self._g_sopt * abs(self._mfy_sismique_dyn)
                )
                ** 2
                + self._D_23**2
                * (abs(self._mfz_sismique_qs) + self._g_sopt * abs(self._mfz_sismique_dyn)) ** 2
            )
            ** 0.5
        )

        return self._sigma_eq, self._sigma_eq_opt
