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
from .sd_interf_dyna_clas import sd_interf_dyna_clas
from .sd_macr_elem_dyna import sd_macr_elem_dyna

# from .sd_base_modale import sd_base_modale
from .sd_mode_meca import sd_mode_meca
from .sd_util import *


class sd_modele_gene(AsBase):
    # -----------------------------
    nomj = SDNom(fin=14)
    MODG_LIPR = AsVI(SDNom(nomj=".MODG.LIPR"))
    MODG_LIDF = AsColl(
        SDNom(nomj=".MODG.LIDF"),
        acces="NU",
        stockage="DISPERSE",
        modelong="CONSTANT",
        type="K",
        ltyp=8,
    )
    MODG_SSTR = AsColl(
        SDNom(nomj=".MODG.SSTR"),
        acces="NU",
        stockage="CONTIG",
        modelong="CONSTANT",
        type="R",
        ltyp=8,
    )
    MODG_SSOR = AsColl(
        SDNom(nomj=".MODG.SSOR"),
        acces="NU",
        stockage="CONTIG",
        modelong="CONSTANT",
        type="R",
        ltyp=8,
    )
    MODG_SSNO = AsPn(SDNom(nomj=".MODG.SSNO"), ltyp=8)
    MODG_SSME = AsColl(
        SDNom(nomj=".MODG.SSME"),
        acces="NU",
        stockage="CONTIG",
        modelong="CONSTANT",
        type="K",
        ltyp=8,
    )
    MODG_DESC = AsVI(SDNom(nomj=".MODG.DESC"), lonmax=3)
    MODG_LIMA = AsColl(
        SDNom(nomj=".MODG.LIMA"),
        acces="NU",
        stockage="DISPERSE",
        modelong="VARIABLE",
        type="R",
        ltyp=8,
    )

    def check_dimensions(self, checker):
        nb_struc = self.MODG_SSME.nmaxoc
        nb_liaison = self.MODG_LIDF.nmaxoc

        assert self.MODG_LIPR.lonmax == 9 * nb_liaison
        assert self.MODG_LIMA.nmaxoc == 3 * nb_liaison
        assert self.MODG_LIMA.nutioc == 3 * nb_liaison

        assert self.MODG_SSNO.nomuti == nb_struc
        assert self.MODG_SSNO.nommax == nb_struc
        assert self.MODG_SSOR.nmaxoc == nb_struc
        assert self.MODG_SSOR.nutioc == nb_struc
        assert self.MODG_SSTR.nmaxoc == nb_struc
        assert self.MODG_SSTR.nutioc == nb_struc

    def check_SSME(self, checker):
        nb_struc = self.MODG_SSME.nmaxoc
        ssme = self.MODG_SSME.get()
        for k in range(nb_struc):
            sd2 = sd_macr_elem_dyna(ssme[k + 1][0].strip())
            sd2.check

    def check_DESC(self, checker):
        desc = self.MODG_DESC.get()
        nomgd = sdu_nom_gd(desc[2])
        assert nomgd == "DEPL_R", (nomgd, desc)
        assert desc[0] > 2 and desc[0] < 15, desc
        assert desc[1] > 2 * 30 and desc[1] < 15 * 30, desc

    def check_SSOR(self, checker):
        nb_struc = self.MODG_SSME.nmaxoc
        ssor = self.MODG_SSOR.get()
        for k in range(nb_struc):
            assert len(ssor[k + 1]) == 3, ssor

    def check_SSTR(self, checker):
        nb_struc = self.MODG_SSME.nmaxoc
        sstr = self.MODG_SSTR.get()
        for k in range(nb_struc):
            assert len(sstr[k + 1]) == 3, sstr

    def check_LIDF(self, checker):
        lidf = self.MODG_LIDF.get()
        nb_liaison = self.MODG_LIDF.nmaxoc
        for k in range(nb_liaison):
            assert len(lidf[k + 1]) == 5, lidf
            assert lidf[k + 1][4].strip() in ("OUI", "NON"), lidf

    def check_LIPR_LIMA(self, checker):
        lipr = self.MODG_LIPR.get()
        lima = self.MODG_LIMA.get()
        nb_liaison = self.MODG_LIDF.nmaxoc
        for k in range(nb_liaison):
            mat1_nlig = lipr[9 * k + 0]
            assert mat1_nlig > 0
            mat1_ncol = lipr[9 * k + 1]
            assert mat1_ncol > 0
            mat1_nume = lipr[9 * k + 2]
            assert mat1_nume == 3 * k + 1, (mat1_nume, k)
            assert len(lima[3 * k + 1]) == mat1_nlig * mat1_ncol

            mat2_nlig = lipr[9 * k + 3]
            assert mat2_nlig > 0
            mat2_ncol = lipr[9 * k + 4]
            assert mat2_ncol > 0
            mat2_nume = lipr[9 * k + 5]
            assert mat2_nume == 3 * k + 2, (mat2_nume, k)
            assert len(lima[3 * k + 2]) == mat2_nlig * mat2_ncol

            mat3_nlig = lipr[9 * k + 6]
            assert mat3_nlig > 0
            mat3_ncol = lipr[9 * k + 7]
            assert mat3_ncol > 0
            mat3_nume = lipr[9 * k + 8]
            assert mat3_nume == 3 * k + 3, (mat3_nume, k)
            assert len(lima[3 * k + 3]) == mat3_nlig * mat3_ncol
