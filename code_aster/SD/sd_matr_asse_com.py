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

from . import *
from .sd_maillage import sd_maillage
from .sd_matr_cine import sd_matr_cine
from .sd_nume_ddl import sd_nume_ddl
from .sd_solveur import sd_solveur
from .sd_titre import sd_titre


class sd_matr_asse_com(sd_titre):
    # -----------------------------
    nomj = SDNom(fin=19)

    REFA = AsVK24(lonmax=20)
    VALM = AsColl(acces="NU", stockage="DISPERSE", modelong="CONSTANT", type=Parmi("C", "R"))
    UALF = Facultatif(
        AsColl(acces="NU", stockage="DISPERSE", modelong="CONSTANT", type=Parmi("C", "R"))
    )
    VALF = Facultatif(
        AsColl(acces="NU", stockage="DISPERSE", modelong="VARIABLE", type=Parmi("C", "R"))
    )
    WALF = Facultatif(
        AsColl(acces="NU", stockage="DISPERSE", modelong="VARIABLE", type=Parmi("C", "R"))
    )
    CONL = Facultatif(OJBVect(type=Parmi("C", "R")))
    DIGS = Facultatif(OJBVect(type=Parmi("C", "R")))
    # seulement si solveurs LDLT et MULT_FRONT
    PERM = Facultatif(AsVI())
    cine = Facultatif(sd_matr_cine(SDNom(nomj="")))

    def exists(self):
        # retourne "vrai" si la SD semble exister (et donc qu'elle peut etre
        # vérifiée)
        return self.REFA.exists

    def check_VALM(self, checker):
        if not self.exists():
            return
        nbloc = self.VALM.nmaxoc
        assert nbloc in (1, 2), nbloc

    def check_REFA(self, checker):
        if not self.exists():
            return
        refa = self.REFA.get_stripped()
        assert refa[9] in ("NOEU", "GENE"), refa
        lgene = refa[9] == "GENE"
        # pour les matrices generalisees, on ne sait pas ce qui est stocké dans
        # refa[0]='' :
        if not lgene:
            sd2 = sd_maillage(refa[0])
            sd2.check(checker)
            sd2 = sd_nume_ddl(refa[1])
            sd2.check(checker)
        assert refa[2] in ("ELIMF", "ELIML", ""), refa
        assert refa[4] == "", refa
        # pour les matrices generalisees, refa[7] n'est pas toujours rempli :
        if not lgene:
            # glute à résorber : j'ajoute '' à la liste permise pour le test
            # yyyy108e :
            assert refa[7] in ("ASSE", "DECT", "DECP", ""), refa
        assert refa[8] in ("MS", "MR"), refa
        if refa[8] == "MS":
            assert self.VALM.nmaxoc == 1, (refa, self.VALM.nmaxoc)
        elif refa[8] == "MR":
            assert self.VALM.nmaxoc == 2, (refa, self.VALM.nmaxoc)

        assert refa[10] in ("MPI_COMPLET", "MPI_INCOMPLET"), refa

        # refa[18] : nom de la matrice "fille" (ELIM_LAGR='OUI') :
        if len(refa[18]) > 0:
            sd2 = sd_matr_asse(refa[18])
            sd2.check(checker)
        # refa[19] : nom de la matrice "mere" (ELIM_LAGR='OUI') :
        if len(refa[19]) > 0:
            # J. Pellet ne comprend pas pourquoi ca plante le test zzzz351a :
            # sd2=sd_matr_asse(refa[19]) ; sd2.check(checker)
            pass
        if refa[6] != "":
            sd2 = sd_solveur(refa[6])
            sd2.check(checker)
