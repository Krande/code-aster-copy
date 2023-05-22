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
:py:class:`GeneralizedAssemblyMatrixReal` --- Generalized Assembly matrix
****************************************************
"""

import numpy as NP

from libaster import GeneralizedAssemblyMatrixComplex, GeneralizedAssemblyMatrixReal

from ..Utilities import injector, deprecated
from ..Objects.Serialization import InternalStateBuilder


class GeneralizedAssemblyMatrixStateBuilder(InternalStateBuilder):
    """Class that returns the internal state of a *GeneralizedAssemblyMatrix* to be pickled."""

    def save(self, matrix):
        """Return the internal state of a *GeneralizedAssemblyMatrix* to be pickled.

        Arguments:
            matrix (*GeneralizedAssemblyMatrix*): The *GeneralizedAssemblyMatrix* object to be pickled.

        Returns:
            *InternalStateBuilder*: The internal state itself.
        """
        super().save(matrix)
        self._st["numbering"] = matrix.getGeneralizedDOFNumbering()
        self._st["base"] = matrix.getModalBasis()
        return self

    def restore(self, matrix):
        """Restore the *GeneralizedAssemblyMatrix* content from the previously saved internal
        state.

        Arguments:
            matrix (*GeneralizedAssemblyMatrix*): The *DataStructure* object to be pickled.
        """
        super().restore(matrix)
        matrix.setGeneralizedDOFNumbering(self._st["numbering"])
        matrix.setModalBasis(self._st["base"])


class BaseGeneralizedAssemblyMatrix:
    """Base object for AssemblyMatrix."""

    cata_sdj = "SD.sd_matr_asse_gene.sd_matr_asse_gene"
    internalStateBuilder = GeneralizedAssemblyMatrixStateBuilder

    def toNumpy(self):
        """Returns the matrix values as `numpy.array`.

        Returns:
            numpy.array: A simple `numpy.array` of the dense matrix.
        """

        # On teste si la matrix existe
        if not self.exists():
            raise AsException("L'objet matrix {0!r} n'existe pas".format(self.getName()))

        if isinstance(self, (GeneralizedAssemblyMatrixReal,)):
            dtype = float
        else:
            dtype = complex

        dim = self.size()
        valeur = NP.zeros([dim, dim], dtype=dtype)

        # Si le stockage est plein
        if self.isDense():
            triang_sup = NP.array(self.getUpperValues())
            assert dim * (dim + 1) // 2 == len(
                triang_sup
            ), "matrix non pleine : %d*(%d+1)/2 != %d" % (dim, dim, len(triang_sup))

            if self.isSymmetric():
                triang_inf = triang_sup
            else:
                triang_inf = NP.array(self.getLowerValues())

            for i in range(dim):
                for j in range(i + 1):
                    k = i * (i + 1) // 2 + j
                    valeur[i, j] = triang_inf[k]
                    valeur[j, i] = triang_sup[k]

        # Si le stockage est diagonal
        elif self.isDiagonal():
            diag = NP.array(self.getUpperValues())
            assert dim == len(diag), "Dimension incorrecte : %d != %d" % (dim, len(diag))
            for i in range(dim):
                valeur[i, i] = diag[i]

        # Sinon on arrete tout
        else:
            raise KeyError

        return valeur

    def fromNumpy(self, matrix):
        """Replace inplace the matrix values by the given `numpy.array`.
           The matrix has to exist and be allocated before.

        Arguments:
            matrix [numpy.array]: A simple `numpy.array` matrix.
        """

        # On teste si le DESC de la matrix existe
        if not self.exists():
            raise AsException("L'objet matrix {0!r} n'existe pas".format(self.getName()))

        NP.asarray(matrix)

        if isinstance(self, (GeneralizedAssemblyMatrixReal,)):
            dtype = float
        else:
            dtype = complex

        # Symétrique ou non
        sym = self.isSymmetric()
        dim = self.size()

        # On teste si la dimension de la matrix python est 2
        if len(NP.shape(matrix)) != 2:
            raise AsException("La dimension de la matrix est incorrecte ")

        # On teste si les tailles des matrixs jeveux et python sont identiques
        if tuple([dim, dim]) != NP.shape(matrix):
            raise AsException("La taille de la matrix est incorrecte ")

        # Si le stockage est plein
        if self.isDense():
            taille = int(dim * dim / 2.0 + dim / 2.0)
            # Triangulaire supérieure
            tmp = NP.zeros(taille, dtype=dtype)
            for j in range(dim):
                for i in range(j + 1):
                    k = j * (j + 1) // 2 + i
                    tmp[k] = matrix[i, j]
            self.setUpperValues(tmp)

            # Cas non-symétrique
            if not sym:
                # Triangulaire inférieure
                tmp = NP.zeros(taille, dtype=dtype)
                for j in range(dim):
                    for i in range(j + 1):
                        k = j * (j + 1) // 2 + i
                        tmp[k] = matrix[j, i]
                self.setLowerValues(tmp)

        # Si le stockage est diagonal
        elif self.isDiagonal():
            tmp = NP.zeros(dim, dtype=dtype)
            for j in range(dim):
                tmp[j] = matrix[j, j]
            self.setUpperValues(tmp)
        # Sinon on arrete tout
        else:
            raise KeyError

        self.build()

    @deprecated(case=1, help="Use 'toNumpy() instead.")
    def EXTR_MATR(self):
        """Returns the matrix values as `numpy.array`."""

        return self.toNumpy()

    @deprecated(case=4, help="Use 'toNumpy() instead.")
    def EXTR_MATR_GENE(self):
        """Returns the matrix values as `numpy.array`."""

        raise RuntimeError("EXTR_MATR_GENE() is replaced by toNumpy()")

        return None

    @deprecated(case=1, help="Use 'fromNumpy() instead.")
    def RECU_MATR(self, matrix):
        """Returns the matrix values as `numpy.array`."""

        self.fromNumpy(matrix)

    @deprecated(case=4, help="Use 'fromNumpy() instead.")
    def RECU_MATR_GENE(self):
        """Returns the matrix values as `numpy.array`."""

        raise RuntimeError("RECU_MATR_GENE() is replaced by fromNumpy()")

        return None


@injector(GeneralizedAssemblyMatrixComplex)
class ExtendedGeneralizedAssemblyMatrixComplex(BaseGeneralizedAssemblyMatrix):
    pass


@injector(GeneralizedAssemblyMatrixReal)
class ExtendedGeneralizedAssemblyMatrixReal(BaseGeneralizedAssemblyMatrix):
    pass
