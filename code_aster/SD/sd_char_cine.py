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

from . import *
from .sd_fonction import sd_fonction
from .sd_modele import sd_modele


class sd_char_cine(AsBase):
    # ===========================
    nomj = SDNom(fin=19)

    AFCK = AsVK8(lonmax=3)
    AFCI = AsVI()
    AFCV = Facultatif(OJBVect(type=Parmi("C", "R", "K")))

    def exists(self):
        # retourne "vrai" si la SD semble exister (et donc qu'elle peut etre
        # vérifiée)
        return self.AFCK.exists

    def u_veri1(self):  # retourne (CIME/CITH/CIAC, RE/CX/FT)
        # ---------------------------------------------------------------
        if not self.exists():
            return
        afck = self.AFCK.get()
        l1 = afck[0].strip().split("_")
        assert len(l1) == 2, afck
        phen, tsca = l1[0], l1[1]
        assert phen in ("CIME", "CITH", "CIAC"), afck
        assert tsca in ("RE", "CX", "FT"), tsca
        return phen, tsca

    def check_AFCK(self, checker):
        # ---------------------------------------------
        if not self.exists():
            return
        phen, tsca = self.u_veri1()
        afck = self.AFCK.get()
        nomo = afck[1].strip()
        sd2 = sd_modele(nomo)
        sd2.check(checker)
        if afck[2].strip() != "":
            assert (phen == "CIME" or phen == "CITH") and tsca == "FT", afck

    def check_AFCI(self, checker):
        # ---------------------------------------------
        if not self.exists():
            return
        phen, tsca = self.u_veri1()
        afci = self.AFCI.get()
        nbloc = afci[0]
        assert len(afci) >= 3 * nbloc + 1, afci
        for k in range(nbloc):
            nuno = afci[3 * k + 1]
            nucmp = afci[3 * k + 2]
            assert afci[3 * k + 3] == 0, (k, afci)
            assert nuno > 0, (k, afci)
            assert nucmp > 0, (k, afci)

    def check_AFCV(self, checker):
        # -------------------------------------------------
        if not self.exists():
            return
        phen, tsca = self.u_veri1()
        afci = self.AFCI.get()
        nbloc = afci[0]
        if not self.AFCV.exists:
            assert tsca == "FT", tsca
            afck = self.AFCK.get()
            assert afck[2].strip() != "", afck
        else:
            tsca2 = self.AFCV.type.strip()
            assert self.AFCV.lonmax >= nbloc, (nbloc, self.AFCV.lonmax)

            if tsca == "RE":
                assert tsca2 == "R", tsca2  # champ de réels
            if tsca == "FT":
                assert tsca2 == "K", tsca2  # champ de fonctions
            if tsca == "CX":
                assert tsca2 == "C", tsca2  # champ de complexes

            # vérification des fonctions :
            if tsca == "FT":
                afcv = self.AFCV.get()
                for fonc in afcv[:nbloc]:
                    sd2 = sd_fonction(fonc)
                    sd2.check(checker)
