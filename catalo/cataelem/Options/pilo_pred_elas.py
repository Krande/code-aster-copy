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

# person_in_charge: kyrylo.kazymyrenko at edf.fr


from cataelem.Tools.base_objects import InputParameter, OutputParameter, Option, CondCalcul
import cataelem.Commons.physical_quantities as PHY
import cataelem.Commons.parameters as SP
import cataelem.Commons.attributes as AT


PCOMPOR = InputParameter(phys=PHY.COMPOR)


PCONTMR = InputParameter(phys=PHY.SIEF_R)


PVARIMR = InputParameter(phys=PHY.VARI_R)


PLST = InputParameter(phys=PHY.NEUT_R)


PLSN = InputParameter(phys=PHY.NEUT_R)


PPINTER = InputParameter(phys=PHY.N816_R)


PAINTER = InputParameter(phys=PHY.N1360R)


PCFACE = InputParameter(phys=PHY.N720_I)


PLONGCO = InputParameter(
    phys=PHY.N120_I,
    container="MODL!.TOPOSE.LON",
    comment="""  XFEM - NBRE DE TETRAEDRES ET DE SOUS-ELEMENTS  """,
)


PBASECO = InputParameter(phys=PHY.N2448R)


PCOPILO = OutputParameter(phys=PHY.PILO_R, type="ELGA")


PILO_PRED_ELAS = Option(
    para_in=(
        PAINTER,
        PBASECO,
        SP.PBORNPI,
        SP.PCAMASS,
        SP.PCDTAU,
        PCFACE,
        SP.PCOHES,
        PCOMPOR,
        SP.PCARCRI,
        PCONTMR,
        SP.PDDEPLR,
        SP.PDEPL0R,
        SP.PDEPL1R,
        SP.PDEPLMR,
        SP.PDONCO,
        SP.PGEOMER,
        SP.PINDCOI,
        PLONGCO,
        PLSN,
        PLST,
        SP.PMATERC,
        PPINTER,
        SP.PTYPEPI,
        PVARIMR,
    ),
    para_out=(PCOPILO,),
    condition=(
        CondCalcul("+", ((AT.PHENO, "ME"), (AT.BORD, "0"))),
        #  -- Pour les elements XFEM, cette option ne concerne que les elements 'XHC':
        #     Mais il est difficile pour l'utilisateur de ne designer QUE les elements XHC,
        #     c'est pour cela que l'on retire les autres elements XFEM qui les bordent :
        CondCalcul("-", ((AT.PHENO, "ME"), (AT.LXFEM, "OUI"))),
        CondCalcul("+", ((AT.PHENO, "ME"), (AT.XFEM, "XHC"))),
    ),
)
