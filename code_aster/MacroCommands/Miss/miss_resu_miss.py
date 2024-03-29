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

# person_in_charge: mathieu.courtois at edf.fr

"""Module permettant de lire un fichier produit par Miss"""

import re

import numpy as NP

from .miss_utils import double, lire_nb_valeurs
from ...Messages import UTMESS


class MissCsolReader:

    """Lit un fichier csol"""

    def __init__(self, nbpc, freq_nb):
        """Initialisation"""
        self.fobj = None
        self.ln = 0
        self.nbpc = nbpc
        self.freq_nb = freq_nb
        self.listfreq = None
        self.values = [ResultatPC() for i in range(nbpc)]

    def read(self, fname):
        """Read the file line per line."""
        try:
            self.fobj = open(fname, "r")
            self._read_all()
            self.fobj.close()
        except (ValueError, IOError, AssertionError) as err:
            UTMESS("F", "MISS0_7", vali=self.ln, valk=str(err))
        return self.listfreq, self.values

    def _read_all(self):
        """Read the file line per line."""
        re_lab = re.compile("POINT *([0-9]+) *CHAMP NO *([0-9]+) *" "DDL NO *([0-9]+)")
        nbf = self.freq_nb
        for ich in range(3):
            for iddl in range(3):
                if ich != iddl:
                    for ifrq in range(self.nbpc * (nbf + 2)):
                        self.fobj.readline()
                        self.ln += 1
                    continue
                for ipc in range(self.nbpc):
                    self.ln += 2
                    self.fobj.readline()
                    mat = re_lab.search(self.fobj.readline())
                    assert mat is not None, "unexpected line"
                    nums = list(map(int, mat.groups()))
                    assert nums == [ipc + 1, ich + 1, iddl + 1], "(%d, %d, %d) expected" % (
                        ipc + 1,
                        ich + 1,
                        iddl + 1,
                    )
                    val = []
                    self.ln += lire_nb_valeurs(self.fobj, 4 * nbf, val, double)
                    array = NP.array(val).reshape((nbf, 4))
                    self.listfreq = array[:, 0]
                    real = array[:, 1]
                    imag = array[:, 2]
                    self.values[ipc].set(iddl, real, imag)


class ResultatPC:

    """Simple conteneur des valeurs relus en un point de contrôle"""

    def __init__(self):
        """Initialisation"""
        self.comp = [None] * 3

    def set(self, component, real, imag):
        """register the value for a component"""
        self.comp[component] = (real, imag)
