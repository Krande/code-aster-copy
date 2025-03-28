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
from .sd_matr_asse_gd import sd_matr_asse_gd
from .sd_proj_mesu import sd_proj_mesu
from .sd_stoc_lciel import sd_stoc_lciel


class sd_macr_elem_stat(AsBase):
    # ----------------------------------------------
    nomj = SDNom(fin=8)

    # description géométrique et topolique :
    DESM = AsVI(lonmax=10)
    REFM = AsVK8()
    LINO = AsVI()
    VARM = AsVR(lonmax=2)
    CONX = Facultatif(AsVI())
    # l'objet devient obligatoire dès l'étape de condensation
    # de la rigidité

    # rigidité condensée :
    rigimeca = Facultatif(sd_matr_asse_gd(SDNom(nomj=".RIGIMECA", fin=19)))
    #   MAEL_RAID_VALE = Facultatif(AsVR())
    MAEL_RAID_VALE = Facultatif(
        AsColl(acces="NU", stockage="DISPERSE", modelong="CONSTANT", type="R", ltyp=8)
    )
    PHI_IE = Facultatif(
        AsColl(acces="NU", stockage="DISPERSE", modelong="CONSTANT", type="R", ltyp=8)
    )

    # masse condensée :
    massmeca = Facultatif(sd_matr_asse_gd(SDNom(nomj=".MASSMECA", fin=19)))
    #   MAEL_MASS_VALE = Facultatif(AsVR())
    MAEL_MASS_VALE = Facultatif(
        AsColl(acces="NU", stockage="DISPERSE", modelong="CONSTANT", type="R", ltyp=8)
    )

    # amortissement condensé :
    #   MAEL_AMOR_VALE = Facultatif(AsVR())
    MAEL_AMOR_VALE = Facultatif(
        AsColl(acces="NU", stockage="DISPERSE", modelong="CONSTANT", type="R", ltyp=8)
    )

    # chargements condensés :
    LICA = Facultatif(
        AsColl(acces="NO", stockage="DISPERSE", modelong="CONSTANT", type="R", ltyp=8)
    )
    LICH = Facultatif(AsColl(acces="NO", stockage="CONTIG", modelong="CONSTANT", type="K", ltyp=8))

    # si utilisation de PROJ_MESU_MODAL :
    PROJM = Facultatif(sd_proj_mesu())

    def check_longueurs(self, checker):
        # ------------------------------------
        # vérifs existence, longueurs, ...

        desm = self.DESM.get()
        refm = self.REFM.get()
        assert desm[0] == 0, desm
        nbnoe, nbnoi, nddle, nddli, nbchar, nbcas, nlage, nlagl, nlagi = desm[1:10]
        assert nbnoe > 0, desm
        assert nbchar >= 0, desm

        # si on n'a pas encore condensé la rigidité, certaines valeurs ne sont
        # pas encore calculées :
        if self.MAEL_RAID_VALE.exists:
            assert nbnoi > 0, desm
            assert nddle > 1, desm
            assert nddli > 0, desm
            assert nbcas >= 0, desm
            assert nlage >= 0, desm
            assert nlagl >= 0, desm
            assert nlagi >= 0, desm
            assert self.CONX.lonmax == 3 * (nbnoe + nlage + nlagl), (desm, self.CONX.get())
            assert refm[5] == "OUI_RIGI"

        assert self.REFM.lonmax == 9 + nbchar, (desm, self.REFM.get())
        assert self.LINO.lonmax == nbnoe, (desm, self.LINO.get())

        # rigidité condensée :
        if self.MAEL_RAID_VALE.exists:
            #           assert self.MAEL_RAID_VALE.lonmax == (nddle * (nddle + 1)) // 2
            assert len(self.MAEL_RAID_VALE.get()[1]) == (nddle * (nddle + 1)) // 2
            assert self.PHI_IE.exists
            phi_ie = self.PHI_IE.get()
            # on ne sait pas faire autrement que
            # ramener l'objet entier ...
            nlblph = len(phi_ie[1]) // nddli  # nombre de lignes de phi_ie par bloc
            assert self.PHI_IE.nmaxoc == (nddle - 1) // nlblph + 1, (nddle, self.PHI_IE.nmaxoc)
            for ke in list(phi_ie.keys()):
                assert len(phi_ie[ke]) == nddli * nlblph, (
                    nddli,
                    nlblph,
                    nddle,
                    len(phi_ie[ke]),
                    ke,
                )

        # masse condensée :
        if self.MAEL_MASS_VALE.exists:
            assert len(self.MAEL_MASS_VALE.get()[1]) == (nddle * (nddle + 1)) // 2
            #           assert self.MAEL_MASS_VALE.lonmax == (nddle * (nddle + 1)) // 2
            assert refm[6] == "OUI_MASS"

        # amortissement condensé : (JP : je ne sais pas si ca peut exister ?)
        if self.MAEL_AMOR_VALE.exists:
            assert len(self.MAEL_AMOR_VALE.get()[1]) == (nddle * (nddle + 1)) // 2
            #           assert self.MAEL_AMOR_VALE.lonmax == (nddle * (nddle + 1)) // 2

            assert refm[7] == "OUI_AMOR"

        # chargements condensés :
        if nbcas > 0:
            assert self.LICA.exists
            assert self.LICA.nmaxoc >= nbcas
            lica = self.LICA.get()
            for k in list(lica.keys()):
                assert len(lica[k]) == 2 * (nddli + nddle)

            assert self.LICH.exists
            assert self.LICH.nmaxoc == self.LICA.nmaxoc
            assert self.LICH.nutioc == self.LICA.nutioc
            lich = self.LICH.get()
            for k in list(lich.keys()):
                assert len(lich[k]) >= nbchar + 1
