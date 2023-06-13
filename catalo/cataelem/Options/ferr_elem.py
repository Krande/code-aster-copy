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


from cataelem.Tools.base_objects import InputParameter, OutputParameter, Option, CondCalcul
import cataelem.Commons.physical_quantities as PHY
import cataelem.Commons.parameters as SP
import cataelem.Commons.attributes as AT


PEFFORR = InputParameter(phys=PHY.SIEF_R, comment=""" PEFFORR : EFFORTS GENERALISES (EFGE_ELNO) """)


FERR_ELEM = Option(
    para_in=(
        SP.PCACOQU,  # codeaster-civilengine%  Caracteristiques des coques
        PEFFORR,  # Etat de contrainte (ou d'effort interne)
        SP.PFERRA1,  # DONNEES UTILISATEUR POUR LE CALCUL DE FERRAILLAGE
        SP.PCAGEPO,  # Caracteristiques des poutre
    ),
    para_out=(SP.PFERRA2,),  # Component fields DNSXI, DNSXS in FERRAILLAGE
    condition=(
        CondCalcul(
            "+",
            (
                (AT.PHENO, "ME"),
                (AT.BORD, "0"),
                (AT.COQUE, "OUI"),
                (AT.EFGE, "OUI"),
                (AT.DIM_TOPO_MODELI, "2"),
            ),
        ),
        CondCalcul(
            "+",
            (
                (AT.PHENO, "ME"),
                (AT.BORD, "0"),
                (AT.POUTRE, "OUI"),
                (AT.EFGE, "OUI"),
                (AT.DIM_TOPO_MODELI, "1"),
            ),
        ),
    ),
    comment="""  FERRAILLAGE : CALCUL DES TAUX DE FERRAILLAGE
           A PARTIR DES EFFORTS GENERALISES (COQUE&POUTRE).
           LICITE EN LINEAIRE SEULEMENT. """,
)
