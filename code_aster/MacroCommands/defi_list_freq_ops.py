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

import numpy as np

from ..Cata.Syntax import _F
from ..Cata.SyntaxUtils import remove_none
from ..CodeCommands import DEFI_LIST_REEL
from ..Messages import UTMESS
from ..Utilities import no_new_attributes


class Refinement:
    """Refinement procedure.

    Arguments:
        keywords (dict): Keywords passed to the RAFFINEMENT keyword.
    """

    epsilon = 1.0e-12
    crit = df_min = nbpts = freq = amor = None
    _disp = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, keywords) -> None:
        self.crit = keywords["CRITERE"]
        self.df_min = keywords["PAS_MINI"]
        self.nbpts = keywords["NB_POINTS"]
        self.freq = np.array(keywords["LIST_RAFFINE"])
        if self.crit in ("RELATIF", "ABSOLU"):
            amor = []
            self._disp = keywords["DISPERSION"]
        else:
            self._disp = None
            if keywords.get("AMOR_REDUIT") is not None:
                amor = list(keywords["AMOR_REDUIT"])
            else:
                amor = keywords["LIST_AMOR"].getValues()
        self.amor = np.array(amor)

    def dispersion(self, idx):
        """Return the dispersion at the given index (the intervalle width to be
        refined)."""
        if self.crit == "ABSOLU":
            return self._disp
        if self.crit == "RELATIF":
            return self._disp * self.freq[idx]
        if self.amor[idx] > self.epsilon:
            return 2.0 * self.amor[idx] * self.freq[idx]
        return 0.01 * self.freq[idx]

    def _extend_amor(self):
        """Extend the list of damping values with the last one if the size is
        less than the number of frequencies.
        """
        if self.amor.size == 0:
            return
        amor2 = np.ones(self.freq.size) * self.amor[-1]
        amor2[: self.amor.size] = self.amor
        self.amor = amor2

    @staticmethod
    def _remove_multiple(first, second=None, delta=0.0):
        """Remove multiple values from the first list and remove same indexes in
        the other lists.

        Arguments:
            first (list[float]): List to be filtered.
            second (list[float]): List to filter with the same indexes.
            delta (float): Precision.
        """
        idx = np.ones(first.size, dtype=bool)
        idx[1:] = np.diff(first) > delta
        first = first[idx]
        if second is not None and second.size > 0:
            second = second[idx]
        return first, second

    def _center_freq(self):
        """Shift the frequencies to be centered on the module resonance"""
        if self.amor.size == 0:
            return
        idx = self.amor > self.epsilon
        # f * sqrt(1 - 2 amor^2) if amor > eps else f * 1.
        coef = np.sqrt(1.0 - 2 * self.amor**2) * idx
        coef += np.ones(self.amor.size) * (1 - idx)
        self.freq *= coef

    def refine(self, list_init):
        """Refine the given list.

        Arguments:
            list_init (list[float]): Initial list of frequencies.
        """
        self._extend_amor()
        self.freq, self.amor = self._remove_multiple(self.freq, self.amor, delta=self.df_min)

        freq0 = self.freq.copy()
        self._center_freq()

        # add refinement around initial frequencies
        cum = np.array([])
        for i, frq in enumerate(self.freq):
            dfr = self.dispersion(i)
            # check for overlap
            if i > 0 and cum.max() > frq - dfr / 2.0:
                UTMESS("I", "DYNAMIQUE_26", valr=(freq0[i - 1], freq0[i]))
            steps = (np.arange(self.nbpts, dtype=float) / (self.nbpts - 1) - 0.5) * dfr
            cum = np.hstack((cum, frq + steps))

        if self.nbpts % 2 == 0:
            # add a point at the center
            cum = np.hstack((cum, self.freq, self.freq / np.sqrt(1 - 2 * self.amor**2)))
        cum = np.hstack((cum, list_init))
        cum.sort()

        self.freq = cum
        # should we keep double frequencies if present in list_init or LIST_RAFFINE?
        result = self._remove_multiple(self.freq, delta=self.df_min)[0]
        return result


def defi_list_freq_ops(self, **args):
    """Definition of a list of frequencies."""
    args = _F(args)
    remove_none(args)
    refine_args = args.pop("RAFFINEMENT", [None])[0]

    listr = DEFI_LIST_REEL(**args)
    list_init = listr.getValues()
    builder = Refinement(refine_args)
    refined_list = builder.refine(list_init)

    return DEFI_LIST_REEL(VALE=refined_list)
