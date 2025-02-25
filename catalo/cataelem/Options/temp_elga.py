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

#

from cataelem.Tools.base_objects import InputParameter, OutputParameter, Option, CondCalcul
import cataelem.Commons.physical_quantities as PHY
import cataelem.Commons.parameters as SP
import cataelem.Commons.attributes as AT


# parametres specifiques X-FEM

PPINTTO = InputParameter(phys=PHY.N132_R)

PCNSETO = InputParameter(phys=PHY.N1280I)

PHEAVTO = InputParameter(phys=PHY.N512_I)

PHEA_NO = InputParameter(phys=PHY.N120_I)

PLONCHA = InputParameter(phys=PHY.N120_I)

PBASLOR = InputParameter(phys=PHY.NEUT_R)

PLSN = InputParameter(phys=PHY.NEUT_R)

PLST = InputParameter(phys=PHY.NEUT_R)

# parametres specifiques aux éléments à sous-points

PNBSP_I = InputParameter(
    phys=PHY.NBSP_I, container="CARA!.CANBSP", comment="""  PNBSP_I : NOMBRE DE SOUS_POINTS  """
)

PVARCPR = InputParameter(
    phys=PHY.VARI_R,
    container="VOLA!&&CCPARA.VARI_INT_N",
    comment="""  VARIABLES DE COMMANDES POUR T+  """,
)

PCHHOBS = InputParameter(phys=PHY.N3600R, comment=""" HHO - coefficient base locale""")


TEMP_ELGA = Option(
    para_in=(
        PBASLOR,
        PCNSETO,
        SP.PGEOMER,
        PHEAVTO,
        PHEA_NO,
        PLONCHA,
        PLSN,
        PLST,
        PPINTTO,
        SP.PTEMPER,
        SP.PCACOQU,
        PVARCPR,
        PNBSP_I,
        PCHHOBS,
    ),
    para_out=(SP.PTEMP_R,),
    condition=(
        CondCalcul("+", ((AT.PHENO, "TH"), (AT.BORD, "0"), (AT.LXFEM, "OUI"))),
        CondCalcul("+", ((AT.PHENO, "TH"), (AT.BORD, "0"), (AT.HHO, "OUI"))),
        CondCalcul("+", ((AT.PHENO, "ME"), (AT.BORD, "0"))),
        CondCalcul("-", ((AT.PHENO, "ME"), (AT.ABSO, "OUI"))),
        CondCalcul("-", ((AT.PHENO, "ME"), (AT.FLUIDE, "OUI"))),
        CondCalcul("-", ((AT.PHENO, "ME"), (AT.INTERFACE, "OUI"))),
    ),
    comment="""
        POUR LES ELEMENTS X-FEM
            CALCUL DE LA TEMPERATURE ET DE SON GRADIENT AUX POINTS DE GAUSS
        POUR LES ELEMENTS A SOUS-POINTS
             CALCUL DE LA TEMPERATURE AUX SOUS-POINTS
        POUR LES ELEMENTS HHO
            CALCUL DE LA TEMPERATURE AUX POINTS DE GAUSS""",
)
