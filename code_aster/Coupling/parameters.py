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

"""
Definition of object that stores the parameters of a coupled simulation
between code_saturne and code_aster.
"""

import sys


class SchemeParams:
    """Object thats holds the values of the parameters of the scheme."""

    __slots__ = ("nb_steps", "init_time", "delta_t", "final_time", "epsilon", "nb_iter", "_")

    def __init__(self):
        self.nb_steps = None
        self.init_time = None
        self.delta_t = None
        self.final_time = None
        self.epsilon = 1.0e-6
        self.nb_iter = 1
        self._ = 0

    def update(self, kwargs):
        """Set values from keyword arguments.

        Arguments:
            kwargs (dict): Dict of parameters values.
        """
        for key, value in kwargs.items():
            try:
                setattr(self, key, value)
            except AttributeError:
                print(f"unknown parameter: {key!r}, ignored", file=sys.stderr)

    def check(self):
        """Check consistency of time steps values."""
        values = [self.nb_steps, self.init_time, self.delta_t, self.final_time]
        if values.count(None) > 1:
            raise ValueError(f"missing values to define time steps: {values}")
        if self.nb_steps is None:
            self.nb_steps = (self.final_time - self.init_time) / self.delta_t
        if self.init_time is None:
            self.init_time = self.final_time - self.delta_t * self.nb_steps
        if self.delta_t is None:
            self.delta_t = (self.final_time - self.init_time) / self.nb_steps
        if self.final_time is None:
            self.final_time = self.init_time + self.delta_t * self.nb_steps
        if (
            abs(self.final_time - (self.init_time + self.delta_t * self.nb_steps))
            > 1.0e6 * self.final_time
        ):
            raise ValueError("inconsistent definition of time steps!")
