# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

import aster
from libaster import (GeneralizedAssemblyMatrixComplex,
                      GeneralizedAssemblyMatrixReal)

from ..Utilities import injector
from ..Objects.Serialization import InternalStateBuilder


def VALM_triang2array(dict_VALM, dim, dtype=None):
    # stockage symetrique ou non (triang inf+sup)
    sym = len(dict_VALM) == 1
    triang_sup = NP.array(dict_VALM[1])
    assert dim*(dim+1) // 2 == len(triang_sup), \
            'Matrice non pleine : %d*(%d+1)/2 != %d' % (dim, dim, len(triang_sup))
    if sym:
        triang_inf = triang_sup
    else:
        triang_inf = NP.array(dict_VALM[2])
    valeur=NP.zeros([dim, dim], dtype=dtype)
    for i in range(1, dim+1):
        for j in range(1, i+1):
            k = i*(i-1) // 2 + j
            valeur[i-1, j-1]=triang_inf[k-1]
            valeur[j-1, i-1]=triang_sup[k-1]
    return valeur

def VALM_diag2array(dict_VALM, dim, dtype=None):
    diag = NP.array(dict_VALM[1])
    assert dim == len(diag), 'Dimension incorrecte : %d != %d' % (dim, len(diag))
    valeur=NP.zeros([dim, dim], dtype=dtype)
    for i in range(dim):
        valeur[i,i] =  diag[i]
    return valeur


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


@injector(GeneralizedAssemblyMatrixComplex)
class ExtendedGeneralizedAssemblyMatrixComplex:
    cata_sdj = "SD.sd_matr_asse_gene.sd_matr_asse_gene"
    internalStateBuilder = GeneralizedAssemblyMatrixStateBuilder

    def EXTR_MATR_GENE(self) :
        desc = self.sdj.DESC.get()
        # On teste si le DESC de la matrice existe
        if not desc:
            raise AsException("L'objet matrice {0!r} n'existe pas"
                              .format(self.sdj.DESC.nomj()))
        desc = NP.array(desc)
        # Si le stockage est plein
        if desc[2] == 2 :
            valeur = VALM_triang2array(self.sdj.VALM.get(), desc[1], complex)

        # Si le stockage est diagonal
        elif desc[2]==1 :
            valeur = VALM_diag2array(self.sdj.VALM.get(), desc[1], complex)

        # Sinon on arrete tout
        else:
            raise KeyError
        return valeur

    def RECU_MATR_GENE(self,matrice) :
        NP.asarray(matrice)
        ncham=self.getName()
        desc = self.sdj.DESC.get()
        # On teste si le DESC de la matrice existe
        if not desc:
            raise AsException("L'objet matrice {0!r} n'existe pas"
                            .format(self.sdj.DESC.nomj()))
        desc = NP.array(desc)
        NP.asarray(matrice)

        # On teste si la dimension de la matrice python est 2
        if (len(NP.shape(matrice)) != 2) :
            raise AsException("La dimension de la matrice est incorrecte ")

        # On teste si la taille de la matrice jeveux et python est identique
        if (tuple([desc[1],desc[1]]) != NP.shape(matrice)) :
            raise AsException("La taille de la matrice est incorrecte ")

        # Si le stockage est plein
        if desc[2]==2 :
            taille=desc[1]*desc[1]/2.0+desc[1]/2.0
            tmpr=NP.zeros([int(taille)])
            tmpc=NP.zeros([int(taille)])
            for j in range(desc[1]+1):
                for i in range(j):
                    k=j*(j-1) // 2+i
                    tmpr[k]=matrice[j-1,i].real
                    tmpc[k]=matrice[j-1,i].imag
            aster.putvectjev('%-19s.VALM' % ncham, len(tmpr), tuple((\
                            list(range(1,len(tmpr)+1)))),tuple(tmpr),tuple(tmpc),1)
        # Si le stockage est diagonal
        elif desc[2]==1 :
            tmpr=NP.zeros(desc[1])
            tmpc=NP.zeros(desc[1])
            for j in range(desc[1]):
                tmpr[j]=matrice[j,j].real
                tmpc[j]=matrice[j,j].imag
            aster.putvectjev('%-19s.VALM' % ncham,len(tmpr),tuple((\
                            list(range(1,len(tmpr)+1)))),tuple(tmpr),tuple(tmpc),1)
        # Sinon on arrete tout
        else:
            raise KeyError
        return

@injector(GeneralizedAssemblyMatrixReal)
class ExtendedGeneralizedAssemblyMatrixReal:
    cata_sdj = "SD.sd_matr_asse_gene.sd_matr_asse_gene"
    internalStateBuilder = GeneralizedAssemblyMatrixStateBuilder

    def EXTR_MATR(self, sparse=False):

        refa = NP.array(self.sdj.REFA.get())

        stock = "diag" if self.sdj.DESC.get()[2]==1 else "full"
        valm = self.sdj.VALM.get()
        sym = len(valm) == 1
        if not sym:
            raise Accas.AsException(
                "Not implemented for non symetric matrix")
        dim = len(valm[1]) if stock=="diag" else int((-1+NP.sqrt(1+8*len(valm[1])))/2.)
        if stock=="diag":
            return NP.diag(valm[1])
        else:
            def make_sym_matrix(n,vals,ntype):
                m = NP.zeros([n,n], dtype=ntype)
                xs,ys = NP.triu_indices(n)
                m[xs,ys] = vals
                m[ys,xs] = vals
                return m
            triang_sup = NP.array(valm[1])
            if type(valm[1][0]) == complex:
                ntype = complex
            else:
                ntype = float
            return make_sym_matrix(dim,triang_sup,ntype)

    def EXTR_MATR_GENE(self):
        desc = self.sdj.DESC.get()
        # On teste si le DESC du vecteur existe
        if not desc:
            raise AsException("L'objet vecteur {0!r} n'existe pas"
                              .format(self.sdj.DESC.nomj()))
        desc = NP.array(desc)

        # Si le stockage est plein
        if desc[2]==2:
            valeur = VALM_triang2array(self.sdj.VALM.get(), desc[1])

        # Si le stockage est diagonal
        elif desc[2]==1:
            valeur = VALM_diag2array(self.sdj.VALM.get(), desc[1])

        # Sinon on arrete tout
        else:
            raise KeyError
        return valeur

    def RECU_MATR_GENE(self,matrice) :
        ncham=self.getName()

        desc = self.sdj.DESC.get()
        # On teste si le DESC de la matrice existe
        if not desc:
            raise AsException("L'objet matrice {0!r} n'existe pas"
                            .format(self.sdj.DESC.nomj()))
        desc = NP.array(desc)

        NP.asarray(matrice)

        # On teste si la dimension de la matrice python est 2
        if (len(NP.shape(matrice)) != 2) :
            raise AsException("La dimension de la matrice est incorrecte ")

        # On teste si les tailles des matrices jeveux et python sont identiques
        if (tuple([desc[1],desc[1]]) != NP.shape(matrice)) :
            raise AsException("La taille de la matrice est incorrecte ")

        # Si le stockage est plein
        if desc[2]==2 :
            taille=desc[1]*desc[1]/2.0+desc[1]/2.0
            tmp=NP.zeros([int(taille)])
            for j in range(desc[1]+1):
                for i in range(j):
                    k=j*(j-1) // 2+i
                    tmp[k]=matrice[j-1,i]
            aster.putcolljev('%-19s.VALM' % ncham,len(tmp),tuple((\
            list(range(1,len(tmp)+1)))),tuple(tmp),tuple(tmp),1)
        # Si le stockage est diagonal
        elif desc[2]==1 :
            tmp=NP.zeros(desc[1])
            for j in range(desc[1]):
                tmp[j]=matrice[j,j]
            aster.putcolljev('%-19s.VALM' % ncham,len(tmp),tuple((\
            list(range(1,len(tmp)+1)))),tuple(tmp),tuple(tmp),1)
        # Sinon on arrete tout
        else:
            raise KeyError
        return
