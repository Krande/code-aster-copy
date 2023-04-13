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

# person_in_charge: mathieu.courtois@edf.fr
"""
:py:class:`AssemblyMatrixDisplacementReal` --- Assembly matrix
****************************************************
"""

import numpy as NP
from libaster import (
    AssemblyMatrixDisplacementComplex,
    AssemblyMatrixDisplacementReal,
    AssemblyMatrixEliminatedReal,
    AssemblyMatrixPressureComplex,
    AssemblyMatrixPressureReal,
    AssemblyMatrixTemperatureComplex,
    AssemblyMatrixTemperatureReal,
    assemblyMatrixToPetsc,
)

from ..SD.sd_stoc_morse import sd_stoc_morse
from ..Utilities import injector
from ..Objects.Serialization import InternalStateBuilder


class AssemblyMatrixStateBuilder(InternalStateBuilder):
    """Class that returns the internal state of a *AssemblyMatrix* to be pickled."""

    def save(self, matrix):
        """Return the internal state of a *AssemblyMatrix* to be pickled.

        Arguments:
            matrix (*AssemblyMatrix*): The *AssemblyMatrix* object to be pickled.

        Returns:
            *InternalStateBuilder*: The internal state itself.
        """
        super().save(matrix)
        self._st["numbering"] = matrix.getDOFNumbering()
        return self

    def restore(self, matrix):
        """Restore the *AssemblyMatrix* content from the previously saved internal
        state.

        Arguments:
            matrix (*AssemblyMatrix*): The *DataStructure* object to be pickled.
        """
        super().restore(matrix)
        matrix.setDOFNumbering(self._st["numbering"])


class BaseAssemblyMatrix:
    """Base object for AssemblyMatrix."""

    cata_sdj = "SD.sd_matr_asse.sd_matr_asse"
    internalStateBuilder = AssemblyMatrixStateBuilder

    def toPetsc(self):
        """Convert the matrix to a PETSc matrix object.

        Returns:
            PetscMat: PETSc matrix.
        """

        if isinstance(self, (AssemblyMatrixDisplacementReal, AssemblyMatrixTemperatureReal)):
            return assemblyMatrixToPetsc(self)
        else:
            raise NotImplementedError("Type not supported by Petsc")


class BaseAssemblyMatrixReal(BaseAssemblyMatrix):
    """Base object for real AssemblyMatrix."""

    def EXTR_MATR(self, sparse=False, epsilon=None):
        """Returns the matrix values as `numpy.array`.

        Arguments:
            sparse (bool): By default, the returned matrix is dense. If *True*
                the returned matrix is sparse.
            epsilon (float): Terms less than this value is considered null.
                By default, no filtering is done.
                Only used if *sparse=True*.

        Returns:
            misc: A single `numpy.array` of the dense matrix if *sparse=False*.
            Or if *sparse=True* a tuple `(data, rows, cols, dim)`. `data`
            contains the values, `rows` the rows indices, `cols` the columns
            indices and `dim` the number of terms.
        """
        refa = NP.array(self.sdj.REFA.get())
        ma = refa[0]
        nu = refa[1]
        smos = sd_stoc_morse(nu[:14] + ".SMOS")
        valm = self.sdj.VALM.get()
        smhc = smos.SMHC.get()
        smdi = smos.SMDI.get()
        ccid = self.sdj.cine.CCID.get()
        ccid = ccid[:-1] if ccid else None  # suppress last entry containing number of bc
        sym = refa[8].strip() == "MS"
        dim = len(smdi)
        nnz = len(smhc)
        triang_sup = NP.array(valm[1])
        if sym:
            triang_inf = triang_sup
        else:
            triang_inf = NP.array(valm[2])
        if type(valm[1][0]) == complex:
            dtype = complex
        else:
            dtype = float

        if sparse:
            rows = NP.array(smhc) - 1
            diag = NP.array(smdi) - 1

            # generate the columns indices
            cols = NP.zeros(nnz, dtype=int)
            jcol = 0
            for i in range(1, nnz):
                if i > diag[jcol]:
                    jcol += 1
                cols[i] = jcol

            # diag is where "row == col"
            helper = rows - cols
            diag_indices = NP.where(helper == 0)[0]

            # transpose indices in 'inf' part and remove diagonal terms
            cols_inf = NP.delete(rows, diag_indices)
            rows_inf = NP.delete(cols, diag_indices)
            triang_inf = NP.delete(triang_inf, diag_indices)

            # join 'sup' and 'inf' parts
            rows = NP.concatenate((rows, rows_inf))
            cols = NP.concatenate((cols, cols_inf))
            data = NP.concatenate((triang_sup, triang_inf))

            # filter terms
            if epsilon is not None:
                nulls = NP.where(abs(data) < epsilon)
                rows = NP.delete(rows, nulls)
                cols = NP.delete(cols, nulls)
                data = NP.delete(data, nulls)

            # apply kinematic boundary conditions
            if ccid:
                elim = NP.where(NP.array(ccid) == 1)[0]
                keep = NP.isin(rows, elim, invert=True) & NP.isin(cols, elim, invert=True)
                data = data[keep]
                rows = rows[keep]
                cols = cols[keep]
                rows = NP.concatenate((elim, rows))
                cols = NP.concatenate((elim, cols))
                data = NP.concatenate((NP.ones(len(elim)), data))
            return data, rows, cols, dim
        else:
            data = NP.zeros([dim, dim], dtype=dtype)
            jcol = 1
            for kterm in range(1, nnz + 1):
                ilig = smhc[kterm - 1]
                if smdi[jcol - 1] < kterm:
                    jcol += 1
                data[jcol - 1, ilig - 1] = triang_inf[kterm - 1]
                data[ilig - 1, jcol - 1] = triang_sup[kterm - 1]
            # apply kinematic boundary conditions
            if ccid:
                elim = NP.where(NP.array(ccid) == 1)
                data[elim, :] = 0.0
                data[:, elim] = 0.0
                data[elim, elim] = 1.0
            return data


class BaseAssemblyMatrixComplex(BaseAssemblyMatrix):
    """Base object for complex AssemblyMatrix."""

    def EXTR_MATR(self, sparse=False, epsilon=None):
        """Returns the matrix values as `numpy.array`.

        Arguments:
            sparse (bool): By default, the returned matrix is dense. If *True*
                the returned matrix is sparse.
            epsilon (float): Terms less than this value is considered null.
                By default, no filtering is done.
                Only used if *sparse=True*.

        Returns:
            misc: A single `numpy.array` of the dense matrix if *sparse=False*.
            Or if *sparse=True* a tuple `(data, rows, cols, dim)`. `data`
            contains the values, `rows` the rows indices, `cols` the columns
            indices and `dim` the number of terms.
        """
        refa = NP.array(self.sdj.REFA.get())
        ma = refa[0]
        nu = refa[1]
        smos = sd_stoc_morse(nu[:14] + ".SMOS")
        valm = self.sdj.VALM.get()
        smhc = smos.SMHC.get()
        smdi = smos.SMDI.get()
        ccid = self.sdj.cine.CCID.get()
        ccid = ccid[:-1] if ccid else None  # suppress last entry containing number of bc
        sym = refa[8].strip() == "MS"
        dim = len(smdi)
        nnz = len(smhc)
        triang_sup = NP.array(valm[1])
        if sym:
            triang_inf = triang_sup
        else:
            triang_inf = NP.array(valm[2])
        if type(valm[1][0]) == complex:
            dtype = complex
        else:
            dtype = float

        if sparse:
            rows = NP.array(smhc) - 1
            diag = NP.array(smdi) - 1

            # generate the columns indices
            cols = NP.zeros(nnz, dtype=int)
            jcol = 0
            for i in range(1, nnz):
                if i > diag[jcol]:
                    jcol += 1
                cols[i] = jcol

            # diag is where "row == col"
            helper = rows - cols
            diag_indices = NP.where(helper == 0)[0]

            # transpose indices in 'inf' part and remove diagonal terms
            cols_inf = NP.delete(rows, diag_indices)
            rows_inf = NP.delete(cols, diag_indices)
            triang_inf = NP.delete(triang_inf, diag_indices)

            # join 'sup' and 'inf' parts
            rows = NP.concatenate((rows, rows_inf))
            cols = NP.concatenate((cols, cols_inf))
            data = NP.concatenate((triang_sup, triang_inf))

            # filter terms
            if epsilon is not None:
                nulls = NP.where(abs(data) < epsilon)
                rows = NP.delete(rows, nulls)
                cols = NP.delete(cols, nulls)
                data = NP.delete(data, nulls)

            # apply kinematic boundary conditions
            if ccid:
                elim = NP.where(NP.array(ccid) == 1)[0]
                keep = NP.isin(rows, elim, invert=True) & NP.isin(cols, elim, invert=True)
                data = data[keep]
                rows = rows[keep]
                cols = cols[keep]
                rows = NP.concatenate((elim, rows))
                cols = NP.concatenate((elim, cols))
                data = NP.concatenate((NP.ones(len(elim)), data))
            return data, rows, cols, dim
        else:
            data = NP.zeros([dim, dim], dtype=dtype)
            jcol = 1
            for kterm in range(1, nnz + 1):
                ilig = smhc[kterm - 1]
                if smdi[jcol - 1] < kterm:
                    jcol += 1
                data[jcol - 1, ilig - 1] = triang_inf[kterm - 1]
                data[ilig - 1, jcol - 1] = triang_sup[kterm - 1]
            # apply kinematic boundary conditions
            if ccid:
                elim = NP.where(NP.array(ccid) == 1)
                data[elim, :] = 0.0
                data[:, elim] = 0.0
                data[elim, elim] = 1.0
            return data


_orig_DisplReal_getType = AssemblyMatrixDisplacementReal.getType


@injector(AssemblyMatrixDisplacementReal)
class ExtendedAssemblyMatrixDisplacementReal(BaseAssemblyMatrixReal):
    pass


@injector(AssemblyMatrixEliminatedReal)
class ExtendedAssemblyMatrixEliminatedReal(BaseAssemblyMatrixReal):
    pass


@injector(AssemblyMatrixDisplacementComplex)
class ExtendedAssemblyMatrixDisplacementComplex(BaseAssemblyMatrixComplex):
    pass


@injector(AssemblyMatrixTemperatureReal)
class ExtendedAssemblyMatrixTemperatureReal(BaseAssemblyMatrixReal):
    pass


@injector(AssemblyMatrixTemperatureComplex)
class ExtendedAssemblyMatrixTemperatureComplex(BaseAssemblyMatrixComplex):
    pass


@injector(AssemblyMatrixPressureReal)
class ExtendedAssemblyMatrixPressureReal(BaseAssemblyMatrixReal):
    pass


@injector(AssemblyMatrixPressureComplex)
class ExtendedAssemblyMatrixPressureComplex(BaseAssemblyMatrixComplex):
    pass
