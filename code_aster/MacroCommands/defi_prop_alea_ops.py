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


"""Commande DEFI_PROP_ALEA"""

import sys
import traceback
from math import ceil, exp, log, pi, sqrt

import numpy as np

from ..CodeCommands import FORMULE
from ..Messages import UTMESS


def defi_prop_alea_ops(self, **kwargs):
    """Corps de la macro DEFI_PROP_ALEA"""
    # conteneur des paramètres du calcul
    params = Randomfield(**kwargs)
    np.random.seed(params.seed)
    # création de l'objet generator
    generator = Generator.factory(self, params)
    try:
        return generator.run()
    except Exception as err:
        trace = "".join(traceback.format_tb(sys.exc_info()[2]))
        UTMESS("F", "SUPERVIS2_5", valk=("DEFI_PROP_ALEA", trace, str(err)))


def evaluate_KL1D(X1, DIM, RANGE, XLISTE, Ux, beta, mediane, pseed):
    np.random.seed(pseed)
    nb = len(Ux[0])
    x1 = (X1 - RANGE[0][0]) / DIM[0]
    U1 = np.array([np.interp(x1, XLISTE[0], term) for term in Ux[0]])
    rand = np.random.normal(0.0, 1.0, nb)
    Ux_1 = mediane * np.exp(beta * np.sum(U1 * rand))
    return Ux_1


def evaluate_KL2D(X1, X2, DIM, RANGE, XLISTE, Ux, beta, mediane, pseed):
    np.random.seed(pseed)
    nb1 = len(Ux[0])
    nb2 = len(Ux[1])
    x1 = (X1 - RANGE[0][0]) / DIM[0]
    x2 = (X2 - RANGE[1][0]) / DIM[1]
    U1 = [np.interp(x1, XLISTE[0], term) for term in Ux[0]]
    U2 = [np.interp(x2, XLISTE[1], term) for term in Ux[1]]
    U1 = np.array(U1).reshape((nb1, 1))
    U2 = np.array(U2).reshape((1, nb2))
    KL_terms = (U1 * U2).ravel()

    rand = np.random.normal(0.0, 1.0, len(KL_terms))
    if Ux[2][0] != "All":
        KL_terms = np.array(KL_terms)[Ux[2]]
        rand = rand[Ux[2]]

    #        rand = np.random.normal(0., 1., len(KL_terms))
    Ux_12 = mediane * np.exp(beta * np.sum(KL_terms * rand))
    return Ux_12


def evaluate_KL3D(X1, X2, X3, DIM, RANGE, XLISTE, Ux, beta, mediane, pseed):
    np.random.seed(pseed)
    nb1 = len(Ux[0])
    nb2 = len(Ux[1])
    nb3 = len(Ux[2])
    x1 = (X1 - RANGE[0][0]) / DIM[0]
    x2 = (X2 - RANGE[1][0]) / DIM[1]
    x3 = (X3 - RANGE[2][0]) / DIM[2]
    U1 = [np.interp(x1, XLISTE[0], term) for term in Ux[0]]
    U2 = [np.interp(x2, XLISTE[1], term) for term in Ux[1]]
    U3 = [np.interp(x3, XLISTE[2], term) for term in Ux[2]]

    U1 = np.array(U1).reshape((nb1, 1, 1))
    U2 = np.array(U2).reshape((1, nb2, 1))
    U3 = np.array(U3).reshape((1, 1, nb3))
    KL_terms = (U1 * U2 * U3).ravel()

    if Ux[3][0] != "All":
        KL_terms = np.array(KL_terms)[Ux[3]]

    rand = np.random.normal(0.0, 1.0, len(KL_terms))
    Ux_123 = mediane * np.exp(beta * np.sum(KL_terms * rand))
    return Ux_123


class Randomfield:
    def __init__(self, **kwargs):
        """Enregistrement des valeurs des mots-clés dans un dictionnaire."""
        # GeneralKeys
        self.args = kwargs
        self.seed = kwargs.get("INIT_ALEA")
        self.mediane = kwargs.get("MEDIANE")
        self.cov = kwargs.get("COEF_VARI")
        self.beta = sqrt(log(1.0 + self.cov**2))
        self.cas = None
        self.nbtot = None
        self.precision = None

        cdict = {
            "RANGE": [],
            "NBTERMS": [],
            "DIM": [],
            "LONG_CORR": [],
            "COORD": None,
            "XLISTE": [],
        }
        # XYZKeys
        liste_coord = []
        if kwargs.get("LONG_CORR_X") is not None:
            self.dimx = kwargs.get("X_MAXI") - kwargs.get("X_MINI")
            self.Lcx = 0.5 * kwargs.get(
                "LONG_CORR_X"
            )  # the parameter of the Markov kernel is a=0.5*Lcorr
            cdict["RANGE"].append((kwargs.get("X_MINI"), kwargs.get("X_MAXI")))
            cdict["NBTERMS"].append(kwargs.get("NB_TERM_X"))
            cdict["DIM"].append(self.dimx)
            cdict["LONG_CORR"].append(self.Lcx)
            xliste = np.arange(0, 1.0 + self.Lcx / self.dimx / 100.0, self.Lcx / self.dimx / 100.0)
            cdict["XLISTE"].append(xliste)
            liste_coord.append("X")

        if kwargs.get("LONG_CORR_Y") is not None:
            self.dimy = kwargs.get("Y_MAXI") - kwargs.get("Y_MINI")
            self.Lcy = 0.5 * kwargs.get("LONG_CORR_Y")
            cdict["RANGE"].append((kwargs.get("Y_MINI"), kwargs.get("Y_MAXI")))
            cdict["NBTERMS"].append(kwargs.get("NB_TERM_Y"))
            cdict["DIM"].append(self.dimy)
            cdict["LONG_CORR"].append(self.Lcy)
            yliste = np.arange(0, 1.0 + self.Lcy / self.dimy / 100.0, self.Lcy / self.dimy / 100.0)
            cdict["XLISTE"].append(yliste)
            liste_coord.append("Y")

        if kwargs.get("LONG_CORR_Z") is not None:
            self.dimz = kwargs.get("Z_MAXI") - kwargs.get("Z_MINI")
            self.Lcz = 0.5 * kwargs.get("LONG_CORR_Z")
            cdict["RANGE"].append((kwargs.get("Z_MINI"), kwargs.get("Z_MAXI")))
            cdict["NBTERMS"].append(kwargs.get("NB_TERM_Z"))
            cdict["DIM"].append(self.dimz)
            cdict["LONG_CORR"].append(self.Lcz)
            zliste = np.arange(0, 1.0 + self.Lcz / self.dimz / 100.0, self.Lcz / self.dimz / 100.0)
            cdict["XLISTE"].append(zliste)
            liste_coord.append("Z")

        if len(liste_coord) > 1:
            if "PRECISION" in kwargs:
                self.precision = kwargs.get("PRECISION")
            if "NB_TERM" in kwargs:
                nbprod = int(np.prod(cdict["NBTERMS"]))
                print("NB_TERM kwargs", kwargs.get("NB_TERM"), nbprod)
                if int(kwargs.get("NB_TERM")) > nbprod:
                    dict_args = dict(
                        valk="NB_TERM must be smaller than the total number of terms computed"
                    )
                    UTMESS("F", "GENERIC_1", **dict_args)
                self.nbtot = kwargs.get("NB_TERM")
        elif len(liste_coord) == 1:
            if "PRECISION" or "NB_TERM" in kwargs:
                dict_args = dict(
                    valk="This is a 1D case. Keywords NB_TERM / PRECISION are used only for 2D or 3D fields"
                )
                UTMESS("A", "GENERIC_1", **dict_args)

        cdict["COORD"] = liste_coord
        self.coord = liste_coord
        self.cas = str(len(liste_coord)) + "D"
        self.data = cdict
        if len(liste_coord) > 3:
            raise ValueError("unknown configuration")


class Generator:
    """Base class Generator"""

    @staticmethod
    def factory(macro, params):
        """create an instance of the appropriated type of Generator"""
        if params.cas == "1D":
            return Generator1(macro, params)
        elif params.cas == "2D":
            return Generator2(macro, params)
        elif params.cas == "3D":
            return Generator3(macro, params)
        else:
            raise ValueError("unknown configuration")

    def __init__(self, macro, params):
        """Constructor Base class"""
        self.macro = macro
        self.mediane = params.mediane
        self.cas = params.cas
        self.data = params.data
        self.mediane = params.mediane
        self.beta = params.beta
        self.seed = params.seed
        self.coord = params.coord
        self.nbtot = params.nbtot
        self.precision = params.precision
        self.KL_terms = ["All"]

    def compute_KL(self):
        """specific to each method"""
        raise NotImplementedError("must be implemented in a subclass")

    def is_even(self, num):
        """Return whether the number num is even."""
        return num % 2 == 0

    def find_roots(self, x, vecf):
        roots = []
        signum = np.sign(vecf)
        for ii, vale in enumerate(signum):
            if ii == 0:
                pass
            elif vale != signum[ii - 1]:
                pos = (x[ii] - x[ii - 1]) * 0.5 + x[ii - 1]
                roots.append(pos)
        return roots

    def eval_eigfunc(self, x, Lc, nbmod):
        vmax = nbmod * 4.0
        v = np.arange(1.0, vmax, 0.01)
        veck = [
            (1 - Lc * vale * np.tan(0.5 * vale)) * (Lc * vale + np.tan(0.5 * vale)) for vale in v
        ]
        troots = self.find_roots(v, veck)
        while len(troots) < nbmod:
            vmax = vmax * 1.5
            v = np.arange(1.0, vmax, 0.01)
            veck = [
                (1 - Lc * vale * np.tan(0.5 * vale)) * (Lc * vale + np.tan(0.5 * vale))
                for vale in v
            ]
            troots = self.find_roots(v, veck)
        roots = troots[:nbmod]
        lamk = 2.0 * Lc * (1.0 + np.array(roots) ** 2 * Lc**2) ** (-1)
        print("NUMBER of ROOTS:", len(troots), "RETAINED EIGENVALUES:", nbmod)
        phik = []
        for ii, vk in enumerate(roots):
            if self.is_even(ii):  # %even
                phik.append(np.cos(vk * (x - 0.5)) / np.sqrt(0.5 * (1.0 + np.sin(vk) / vk)))
            else:  # %odd
                phik.append(np.sin(vk * (x - 0.5)) / np.sqrt(0.5 * (1.0 - np.sin(vk) / vk)))
        return list(zip(lamk, phik))


class Generator1(Generator):
    """1D class"""

    def compute_KL(self):
        Lcx1 = self.data["LONG_CORR"][0]
        dimx1 = self.data["DIM"][0]
        nbmod1 = int(self.data["NBTERMS"][0])
        KL_data1 = self.eval_eigfunc(self.data["XLISTE"][0], Lcx1 / dimx1, nbmod1)
        self.Ux1 = [np.sqrt(leig) * np.array(veig) for (leig, veig) in KL_data1]

    def run(self):
        self.compute_KL()

        if self.coord == ["X"]:
            formule_out = FORMULE(
                NOM_PARA=("X"),
                VALE="user_func(X,DIM,RANGE,XLISTE,Ux,beta,mediane, seed)",
                user_func=evaluate_KL1D,
                XLISTE=self.data["XLISTE"],
                DIM=self.data["DIM"],
                RANGE=self.data["RANGE"],
                Ux=(self.Ux1,),
                mediane=self.mediane,
                beta=self.beta,
                seed=self.seed,
            )
        elif self.coord == ["Y"]:
            formule_out = FORMULE(
                NOM_PARA=("Y"),
                VALE="user_func(Y,DIM,RANGE,XLISTE,Ux,beta,mediane, seed)",
                user_func=evaluate_KL1D,
                XLISTE=self.data["XLISTE"],
                DIM=self.data["DIM"],
                RANGE=self.data["RANGE"],
                Ux=(self.Ux1,),
                mediane=self.mediane,
                beta=self.beta,
                seed=self.seed,
            )
        elif self.coord == ["Z"]:
            formule_out = FORMULE(
                NOM_PARA=("Z"),
                VALE="user_func(Z,DIM,RANGE,XLISTE,Ux,beta,mediane, seed)",
                user_func=evaluate_KL1D,
                XLISTE=self.data["XLISTE"],
                DIM=self.data["DIM"],
                RANGE=self.data["RANGE"],
                Ux=(self.Ux1,),
                mediane=self.mediane,
                beta=self.beta,
                seed=self.seed,
            )
        else:
            raise ValueError("unknown configuration")
        return formule_out


class Generator2(Generator):
    """2D class"""

    def compute_KL(self):
        Lcx1 = self.data["LONG_CORR"][0]
        Lcx2 = self.data["LONG_CORR"][1]
        dimx1 = self.data["DIM"][0]
        dimx2 = self.data["DIM"][1]
        nbmod1 = int(self.data["NBTERMS"][0])
        nbmod2 = int(self.data["NBTERMS"][1])

        KL_data1 = self.eval_eigfunc(self.data["XLISTE"][0], Lcx1 / dimx1, nbmod1)
        self.Ux1 = [np.sqrt(leig) * np.array(veig) for (leig, veig) in KL_data1]
        eig1, vec = zip(*KL_data1)

        KL_data2 = self.eval_eigfunc(self.data["XLISTE"][1], Lcx2 / dimx2, nbmod2)
        self.Ux2 = [np.sqrt(leig) * np.array(veig) for (leig, veig) in KL_data2]
        eig2, vec = zip(*KL_data2)

        self.eigs = [eig1, eig2]

    def select_KL_terms(self):
        eig1 = np.array(self.eigs[0]).reshape((len(self.eigs[0]), 1))
        eig2 = np.array(self.eigs[1]).reshape((1, len(self.eigs[1])))
        eig12 = (eig1 * eig2).ravel()
        ind = np.flip(np.argsort(eig12), 0)
        eig12 = np.flip(np.sort(eig12), 0)

        #        print('squared sorted eigs')
        #        print(eig12)
        #        print('squared sorted summed eigs')
        #        print(np.cumsum(eig12))

        if self.precision:
            ind_cut = np.searchsorted(
                np.cumsum(eig12), self.precision * np.sum(eig12), side="right"
            )
            cutindlist = ind[:ind_cut]
        #            print('precision')
        #            print(ind, ind_cut, cutindlist)
        elif self.nbtot:
            cutindlist = ind[: self.nbtot]
        #            print('nbtot')
        #            print(cutindlist, self.nbtot)

        print("TOTAL NUMBER OF RETAINED EIGENVALUES:", len(cutindlist))
        self.KL_terms = cutindlist

    def run(self):
        self.compute_KL()
        if self.precision or self.nbtot:
            self.select_KL_terms()

        print("X,Y", self.coord)
        if self.coord == ["X", "Y"]:
            formule_out = FORMULE(
                NOM_PARA=("X", "Y"),
                VALE="user_func(X,Y,DIM,RANGE,XLISTE,Ux,beta,mediane, seed)",
                user_func=evaluate_KL2D,
                XLISTE=self.data["XLISTE"],
                DIM=self.data["DIM"],
                RANGE=self.data["RANGE"],
                Ux=(self.Ux1, self.Ux2, self.KL_terms),
                mediane=self.mediane,
                beta=self.beta,
                seed=self.seed,
            )
        elif self.coord == ["X", "Z"]:
            formule_out = FORMULE(
                NOM_PARA=("X", "Z"),
                VALE="user_func(X,Z,DIM,RANGE,XLISTE,Ux,beta,mediane, seed)",
                user_func=evaluate_KL2D,
                XLISTE=self.data["XLISTE"],
                DIM=self.data["DIM"],
                RANGE=self.data["RANGE"],
                Ux=(self.Ux1, self.Ux2, self.KL_terms),
                mediane=self.mediane,
                beta=self.beta,
                seed=self.seed,
            )
        elif self.coord == ["Y", "Z"]:
            formule_out = FORMULE(
                NOM_PARA=("Y", "Z"),
                VALE="user_func(Y,Z,DIM,RANGE,XLISTE,Ux,beta,mediane, seed)",
                user_func=evaluate_KL2D,
                XLISTE=self.data["XLISTE"],
                DIM=self.data["DIM"],
                RANGE=self.data["RANGE"],
                Ux=(self.Ux1, self.Ux2, self.KL_terms),
                mediane=self.mediane,
                beta=self.beta,
                seed=self.seed,
            )
        else:
            raise ValueError("unknown configuration")
        return formule_out


class Generator3(Generator):
    """3D class"""

    def compute_KL(self):
        Lcx1 = self.data["LONG_CORR"][0]
        Lcx2 = self.data["LONG_CORR"][1]
        Lcx3 = self.data["LONG_CORR"][2]
        dimx1 = self.data["DIM"][0]
        dimx2 = self.data["DIM"][1]
        dimx3 = self.data["DIM"][2]
        nbmod1 = int(self.data["NBTERMS"][0])
        nbmod2 = int(self.data["NBTERMS"][1])
        nbmod3 = int(self.data["NBTERMS"][2])

        KL_data1 = self.eval_eigfunc(self.data["XLISTE"][0], Lcx1 / dimx1, nbmod1)
        self.Ux1 = [np.sqrt(leig) * np.array(veig) for (leig, veig) in KL_data1]
        eig1, vec = zip(*KL_data1)

        KL_data2 = self.eval_eigfunc(self.data["XLISTE"][1], Lcx2 / dimx2, nbmod2)
        self.Ux2 = [np.sqrt(leig) * np.array(veig) for (leig, veig) in KL_data2]
        eig2, vec = zip(*KL_data2)

        KL_data3 = self.eval_eigfunc(self.data["XLISTE"][2], Lcx3 / dimx3, nbmod3)
        self.Ux3 = [np.sqrt(leig) * np.array(veig) for (leig, veig) in KL_data3]
        eig3, vec = zip(*KL_data3)

        self.eigs = [eig1, eig2, eig3]

    def select_KL_terms(self):
        eig1 = np.array(self.eigs[0]).reshape((len(self.eigs[0]), 1, 1))
        eig2 = np.array(self.eigs[1]).reshape((1, len(self.eigs[1]), 1))
        eig3 = np.array(self.eigs[2]).reshape((1, 1, len(self.eigs[2])))
        eig123 = (eig1 * eig2 * eig3).ravel()
        ind = np.flip(np.argsort(eig123), 0)
        eig123 = np.flip(np.sort(eig123), 0)

        if self.precision:
            ind_cut = np.searchsorted(
                np.cumsum(eig123), self.precision * np.sum(eig123), side="right"
            )
            cutindlist = ind[:ind_cut]
        elif self.nbtot:
            cutindlist = ind[: self.nbtot]

        print("TOTAL NUMBER OF RETAINED EIGENVALUES:", len(cutindlist))
        self.KL_terms = cutindlist

    def run(self):
        self.compute_KL()
        if self.precision or self.nbtot:
            self.select_KL_terms()

        formule_out = FORMULE(
            NOM_PARA=("X", "Y", "Z"),
            VALE="user_func(X, Y, Z, DIM, RANGE, XLISTE, Ux, beta, mediane, seed)",
            user_func=evaluate_KL3D,
            XLISTE=self.data["XLISTE"],
            DIM=self.data["DIM"],
            RANGE=self.data["RANGE"],
            Ux=(self.Ux1, self.Ux2, self.Ux3, self.KL_terms),
            mediane=self.mediane,
            beta=self.beta,
            seed=self.seed,
        )
        return formule_out
