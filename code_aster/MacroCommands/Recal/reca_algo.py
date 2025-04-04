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

# person_in_charge: mathieu.courtois@edf.fr

import copy
import os

import numpy as NP
import numpy.linalg as linalg

from ...Utilities import disable_fpe
from ...Messages import UTMESS

from ...Cata.Syntax import _F


def calcul_gradient(A, erreur):
    grad = NP.dot(NP.transpose(A), erreur)
    return grad


def calcul_norme2(V):
    a = NP.array(V)
    return NP.dot(a, NP.transpose(a)) ** 0.5


class Dimension:
    """
    Classe gérant l'adimensionnement et le dimensionnement
    """

    def __init__(self, val_initiales):
        """
        Le constructeur calcul la matrice D et son inverse
        """
        self.val_init = val_initiales
        dim = len(self.val_init)
        self.D = NP.zeros((dim, dim), float)
        for i in range(dim):
            self.D[i][i] = self.val_init[i]
        self.inv_D = linalg.inv(self.D)

    # ------------------------------------------------------------------------
    def adim_sensi(self, A):
        for i in range(A.shape[0]):
            for j in range(A.shape[1]):
                A[i, j] = A[i, j] * self.val_init[j]
        return A

    # ------------------------------------------------------------------------
    def redim_sensi(self, A):
        for i in range(A.shape[0]):
            for j in range(A.shape[1]):
                A[i, j] = A[i, j] / self.val_init[j]
        return A

    # ------------------------------------------------------------------------
    def adim(self, tab):
        tab_adim = NP.dot(self.inv_D, copy.copy(tab))
        return tab_adim

    # ------------------------------------------------------------------------
    def redim(self, tab_adim):
        tab = NP.dot(self.D, tab_adim)
        return tab


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------
def cond(matrix):
    e1 = linalg.eigvals(matrix)
    e = list(map(abs, e1))
    size = len(e)
    e = NP.sort(e)
    if NP.all(e[0] != 0):
        condi = e[size - 1] / e[0]
    else:
        condi = 0.0
    return condi, e[size - 1], e[0]


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------
def norm(matrix):
    e = linalg.eigvalsh(matrix)
    size = len(e)
    e = NP.sort(e)
    norm = e[size - 1]
    return norm


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------
def lambda_init(matrix):
    """
    Routine qui calcule la valeur initial du parametre de regularisation l.
    """
    condi, emax, emin = cond(matrix)
    id = NP.identity(matrix.shape[0])
    if condi == 0.0:
        l = 1.0e-3 * norm(matrix)
    elif condi <= 10000:
        l = 1.0e-16 * norm(matrix)
    elif condi > 10000:
        l = abs(10000.0 * emin - emax) / 10001.0
    return l


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------
def Levenberg_bornes(val, Dim, val_init, borne_inf, borne_sup, A, erreur, l, ul_out):
    """
    On resoud le système par contraintes actives:
       Q.dval + s + d =0
       soumis à :
       borne_inf < dval < borne_sup
               0 <  s
               s.(borne_inf - dval)=0
               s.(borne_sup - dval)=0
    """
    dim = len(val)
    id = NP.identity(dim)
    # Matrice du système
    Q = NP.dot(NP.transpose(A), A) + l * id
    # Second membre du système
    d = NP.dot(NP.transpose(A), erreur)
    # Ens. de liaisons actives
    Act = NP.array([], dtype=int)
    k = 0
    done = 0
    # Increment des parametres
    dval = NP.zeros(dim)
    while done < 1:
        k = k + 1
        I = NP.ones(dim, dtype=int)
        for i in Act:
            I[i] = 0
        I = NP.nonzero(NP.greater(I, 0))[0]
        s = NP.zeros(dim)
        for i in Act:
            # test sur les bornes (on stocke si on est en butée haute ou basse)
            if val[i] + dval[i] >= borne_sup[i]:
                dval[i] = borne_sup[i] - val[i]
                s[i] = 1.0
            if val[i] + dval[i] <= borne_inf[i]:
                dval[i] = borne_inf[i] - val[i]
                s[i] = -1.0
        if len(I) != 0:
            # xi=-Q(I)-1.(d(I)+Q(I,Act).dval(Act))
            t_QI = NP.take(Q, I, axis=0)
            t_tQI_Act = NP.take(t_QI, Act, axis=1)
            t_adim_Act = NP.take(Dim.adim(dval), Act)
            if NP.size(t_tQI_Act) > 0 and NP.size(t_adim_Act) > 0:
                smemb = NP.take(d, I) + NP.dot(t_tQI_Act, t_adim_Act)
            else:
                smemb = NP.take(d, I)
            xi = -linalg.solve(NP.take(t_QI, I, axis=1), smemb)
            for i in NP.arange(len(I)):
                dval[I[i]] = xi[i] * val_init[I[i]]
        if len(Act) != 0:
            # s(Av)=-d(Act)-Q(Act,:).dval
            sa = -NP.take(d, Act) - NP.dot(NP.take(Q, Act, axis=0), Dim.adim(dval))
            for i in range(len(Act)):
                if s[Act[i]] == -1.0:
                    s[Act[i]] = -sa[i]
                else:
                    s[Act[i]] = sa[i]
        # Nouvel ens. de liaisons actives
        Act = NP.concatenate(
            (
                NP.nonzero(NP.greater(dval, borne_sup - val))[0],
                NP.nonzero(NP.less(dval, borne_inf - val))[0],
                NP.nonzero(NP.greater(s, 0.0))[0],
            )
        ).astype(int)
        done = (
            (max(val + dval - borne_sup) <= 0)
            & (min(val + dval - borne_inf) >= 0)
            & (min(s) >= 0.0)
        )
        # Pour éviter le cyclage
        if k > 50:
            try:
                l = l * 2
                Q = NP.dot(NP.transpose(A), A) + l * id
                k = 0
            except:
                res = open(os.getcwd() + "/fort." + str(ul_out), "a")
                res.write("\n\nQ = \n" + NP.array2string(Q - l * id, array_output=1, separator=","))
                res.write("\n\nd = " + NP.array2string(d, array_output=1, separator=","))
                res.write("\n\nval = " + NP.array2string(val, array_output=1, separator=","))
                res.write(
                    "\n\nval_ini= " + NP.array2string(val_init, array_output=1, separator=",")
                )
                res.write(
                    "\n\nborne_inf= " + NP.array2string(borne_inf, array_output=1, separator=",")
                )
                res.write(
                    "\n\nborne_sup= " + NP.array2string(borne_sup, array_output=1, separator=",")
                )
                UTMESS("F", "RECAL0_18")
                return
    newval = copy.copy(val + dval)
    return newval, s, l, Act


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------
def actualise_lambda(l, val, new_val, A, erreur, new_J, old_J):
    dim = len(val)
    id = NP.identity(dim)
    # Matrice du système
    Q = NP.dot(NP.transpose(A), A) + l * id
    # Second membre du système
    d = NP.dot(NP.transpose(A), erreur)
    old_Q = old_J
    new_Q = (
        old_J
        + 0.5 * NP.dot(NP.transpose(new_val - val), NP.dot(Q, new_val - val))
        + NP.dot(NP.transpose(new_val - val), d)
    )
    # Ratio de la décroissance réelle et de l'approx. quad.
    if NP.all((old_Q - new_Q) != 0.0):
        R = (old_J - new_J) / (old_Q - new_Q)
        if R < 0.25:
            l = l * 10.0
        elif R > 0.75:
            l = l / 15.0
    else:
        if old_J > new_J:
            l = l * 10.0
        else:
            l = l / 10.0
    return l


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------
def test_convergence(gradient_init, erreur, A, s):
    """
    Renvoie le residu
    """
    gradient = calcul_gradient(A, erreur) + s
    try:
        epsilon = NP.dot(gradient, gradient) / NP.dot(gradient_init, gradient_init)
    except:
        UTMESS("F", "RECAL0_19")
        return
    epsilon = epsilon**0.5
    return epsilon


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------
def calcul_etat_final(para, A, iter, max_iter, prec, residu, Messg):
    """
    Fonction appelée quand la convergence est atteinte
    on calcule le Hessien et les valeurs propres et vecteurs
    propre associés au Hessien
    A    = sensibilite
    At*A = hessien
    """

    #   if ((iter < max_iter) or (residu < prec)):
    Hessien = NP.dot(NP.transpose(A), A)

    # Desactive temporairement les FPE qui pourraient etre generees (a
    # tord!) par blas
    with disable_fpe():
        valeurs_propres, vecteurs_propres = linalg.eig(Hessien)
        vecteurs_propres = NP.transpose(
            vecteurs_propres
        )  # numpy et Numeric n'ont pas la meme convention
        sensible = NP.nonzero(NP.greater(abs(valeurs_propres / max(abs(valeurs_propres))), 1.0e-1))[
            0
        ]
        insensible = NP.nonzero(NP.less(abs(valeurs_propres / max(abs(valeurs_propres))), 1.0e-2))[
            0
        ]

    Messg.affiche_calcul_etat_final(
        para, Hessien, valeurs_propres, vecteurs_propres, sensible, insensible
    )
