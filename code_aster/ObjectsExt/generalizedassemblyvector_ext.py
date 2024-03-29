# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
:py:class:`GeneralizedAssemblyVectorOnMesh` --- Assignment of material properties on mesh
************************************************************************
"""

import numpy

import aster
from libaster import GeneralizedAssemblyVectorComplex, GeneralizedAssemblyVectorReal

from ..Utilities import injector


@injector(GeneralizedAssemblyVectorComplex)
class ExtendedGeneralizedAssemblyVectorComplex:
    cata_sdj = "SD.sd_cham_gene.sd_cham_gene"

    def EXTR_VECT_GENE_C(self):

        valeur = numpy.array(self.sdj.VALE.get(), complex)

        return valeur

    def RECU_VECT_GENE_C(self, vecteur):
        numpy.asarray(vecteur)
        ncham = self.getName()
        ncham = ncham + (8 - len(ncham)) * " "
        desc = self.sdj.DESC.get()
        # On teste si le DESC de la matrice existe
        if not desc:
            raise AsException("L'objet vecteur {0!r} n'existe pas".format(self.sdj.DESC.nomj()))
        desc = numpy.array(desc)
        # On teste si la taille de la matrice jeveux et python est identique
        if desc[1] != numpy.shape(vecteur)[0]:
            raise AsException("La taille du vecteur python est incorrecte")
        tmpr = vecteur.real
        tmpc = vecteur.imag
        aster.putvectjev(
            ncham + (19 - len(ncham)) * " " + ".VALE",
            len(tmpr),
            tuple(range(1, len(tmpr) + 1)),
            tuple(tmpr),
            tuple(tmpc),
            1,
        )
        return


@injector(GeneralizedAssemblyVectorReal)
class ExtendedGeneralizedAssemblyVectorReal:
    cata_sdj = "SD.sd_cham_gene.sd_cham_gene"

    def EXTR_VECT_GENE_R(self):
        valeur = numpy.array(self.sdj.VALE.get())
        return valeur

    def EXTR_VECT_GENE_C(self):
        valeur = numpy.array(self.sdj.VALE.get(), complex)

        return valeur

    def RECU_VECT_GENE_C(self, vecteur):
        numpy.asarray(vecteur)
        ncham = self.getName()
        ncham = ncham + (8 - len(ncham)) * " "
        desc = self.sdj.DESC.get()
        # On teste si le DESC de la matrice existe
        if not desc:
            raise AsException("L'objet vecteur {0!r} n'existe pas".format(self.sdj.DESC.nomj()))
        desc = numpy.array(desc)
        # On teste si la taille de la matrice jeveux et python est identique
        if desc[1] != numpy.shape(vecteur)[0]:
            raise AsException("La taille du vecteur python est incorrecte")
        tmpr = vecteur.real
        tmpc = vecteur.imag
        aster.putvectjev(
            ncham + (19 - len(ncham)) * " " + ".VALE",
            len(tmpr),
            tuple(range(1, len(tmpr) + 1)),
            tuple(tmpr),
            tuple(tmpc),
            1,
        )
        return
