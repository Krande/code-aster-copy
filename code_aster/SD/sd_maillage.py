# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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
from .sd_carte import sd_carte
from .sd_cham_geom import sd_cham_geom
from .sd_l_table import sd_l_table
from .sd_titre import sd_titre


class sd_maillage(sd_titre):
    # -------------------------------
    nomj = SDNom(fin=8)

    DIME = AsVI(lonmax=6)

    # un sd_maillage a toujours des noeuds :
    COORDO = sd_cham_geom()
    NOMNOE = Facultatif(AsPn(ltyp=8))

    # normalement, un sd_maillage a toujours une "sd_l_table" contenant des
    # caractéristiques géométriques :
    lt = sd_l_table(SDNom(nomj=""))

    # si le sd_maillage a des groupes :
    GROUPENO = Facultatif(AsColl(stockage="DISPERSE", modelong="VARIABLE", type="I"))
    GROUPEMA = Facultatif(AsColl(stockage="DISPERSE", modelong="VARIABLE", type="I"))
    PTRNOMNOE = Facultatif(AsPn(type="K", ltyp=24))
    PTRNOMMAI = Facultatif(AsPn(type="K", ltyp=24))

    # si le sd_maillage a des mailles :
    CONNEX = Facultatif(AsColl(acces="NU", stockage="CONTIG", modelong="VARIABLE", type="I"))
    TYPMAIL = Facultatif(AsVI())
    NOMMAI = Facultatif(AsPn(ltyp=8))

    # si le sd_maillage a des patchs:
    PATCH = Facultatif(AsColl(acces="NU", stockage="CONTIG", modelong="VARIABLE", type="I"))
    COMAPA = Facultatif(AsVI())
    CONOPA = Facultatif(AsVI())
    PTRNOMPAT = Facultatif(AsVK24())

    # si le sd_maillage a des super-mailles :
    NOMACR = Facultatif(AsVK8())
    SUPMAIL = Facultatif(AsColl(acces="NO", stockage="DISPERSE", modelong="VARIABLE", type="I"))
    PARA_R = Facultatif(AsVR())
    TYPL = Facultatif(AsVI())

    # si le sd_maillage est linéique (tube_GV) :
    absc_curv = Facultatif(sd_carte(SDNom(nomj=".ABSC_CURV")))

    ADAPTATION = Facultatif(AsVI(lonmax=1))

    # Ces objets sont nécessaires pour CREA_MAILLAGE RESTREINT
    CRNO = Facultatif(AsVI())
    CRMA = Facultatif(AsVI())
    MAOR = Facultatif(AsVK8())

    def u_dime(self):
        dime = self.DIME.get()
        nb_no = dime[0]
        nb_nl = dime[1]
        nb_ma = dime[2]
        nb_sm = dime[3]
        nb_sm_mx = dime[4]
        dim_coor = dime[5]
        return nb_no, nb_nl, nb_ma, nb_sm, nb_sm_mx, dim_coor

    # remarque :  la sd_maillage pouvant etre "volumineuse", on s'interdit (pour des raisons de temps CPU)
    #             de vérifier le contenu des gros objets.

    def check_DIME(self, checker):
        nb_no, nb_nl, nb_ma, nb_sm, nb_sm_mx, dim_coor = self.u_dime()
        assert nb_sm <= nb_sm_mx, (nb_sm, nb_sm_mx)
        if nb_nl > 0:
            assert nb_sm > 0
        assert nb_no > 0, nb_no
        assert dim_coor in (2, 3), dim_coor

    def check_NOEUDS(self, checker):
        nb_no, nb_nl, nb_ma, nb_sm, nb_sm_mx, dim_coor = self.u_dime()
        assert self.COORDO.VALE.lonmax == 3 * nb_no, nb_no

    def check_MAILLES(self, checker):
        nb_no, nb_nl, nb_ma, nb_sm, nb_sm_mx, dim_coor = self.u_dime()
        if nb_ma == 0:
            return
        assert self.TYPMAIL.lonmax == nb_ma, nb_ma
        assert self.CONNEX.nmaxoc == nb_ma, nb_ma

    def check_SSS(self, checker):
        nb_no, nb_nl, nb_ma, nb_sm, nb_sm_mx, dim_coor = self.u_dime()
        if nb_sm == 0:
            return
        assert self.NOMACR.lonmax == nb_sm, nb_sm
        assert self.PARA_R.lonmax == 14 * nb_sm, nb_sm
        assert self.SUPMAIL.nmaxoc == nb_sm, nb_sm

    def check_TYPL(self, checker):
        nb_no, nb_nl, nb_ma, nb_sm, nb_sm_mx, dim_coor = self.u_dime()
        if nb_nl == 0:
            return
        assert self.TYPL.lonmax == nb_nl, nb_nl
        typl = self.TYPL.get()
        for k in typl:
            assert k in (-1, -2), typl


class sd_connection_mesh(AsBase):
    nomj = SDNom(fin=8)
