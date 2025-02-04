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

# person_in_charge: sylvie.granet at edf.fr


from cataelem.Tools.base_objects import InputParameter, OutputParameter, Option, CondCalcul
import cataelem.Commons.physical_quantities as PHY
import cataelem.Commons.parameters as SP
import cataelem.Commons.attributes as AT


PPINTTO = InputParameter(phys=PHY.N132_R)


PCNSETO = InputParameter(
    phys=PHY.N1280I,
    container="MODL!.TOPOSE.CNS",
    comment="""  XFEM - CONNECTIVITE DES SOUS-ELEMENTS  """,
)


PHEAVTO = InputParameter(phys=PHY.N512_I)


PLONCHA = InputParameter(
    phys=PHY.N120_I,
    container="MODL!.TOPOSE.LON",
    comment="""  XFEM - NBRE DE TETRAEDRES ET DE SOUS-ELEMENTS  """,
)


PLSN = InputParameter(phys=PHY.NEUT_R)


PLST = InputParameter(phys=PHY.NEUT_R)

PSTANO = InputParameter(phys=PHY.N120_I)


PPMILTO = InputParameter(phys=PHY.N792_R)


PFISNO = InputParameter(phys=PHY.NEUT_I)


PHEA_NO = InputParameter(phys=PHY.N120_I)


PHEA_SE = InputParameter(phys=PHY.N512_I)


CHAR_MECA_FLUX_F = Option(
    para_in=(
        PCNSETO,
        PFISNO,
        SP.PFLUXF,
        SP.PGEOMER,
        PHEAVTO,
        PHEA_NO,
        PHEA_SE,
        PLONCHA,
        PLSN,
        PLST,
        PPINTTO,
        PPMILTO,
        PSTANO,
        SP.PINSTR,
        SP.PHEAVNO,
    ),
    para_out=(SP.PVECTUR,),
    condition=(
        CondCalcul("+", ((AT.PHENO, "ME"), (AT.TYPMOD2, "THM"), (AT.BORD, "-1"))),
        CondCalcul("+", ((AT.PHENO, "ME"), (AT.TYPMOD2, "JHMS"), (AT.BORD, "-1"))),
        CondCalcul("+", ((AT.PHENO, "ME"), (AT.TYPMOD2, "XFEM_HM"), (AT.BORD, "-1"))),
        CondCalcul("+", ((AT.PHENO, "ME"), (AT.TYPMOD2, "XFEM_HM"), (AT.CONTACT, "OUI"))),
    ),
    comment=""" CHAR_MECA_FLUX_F (MOT-CLE FLUX_THM_REP) : CALCUL DU SECOND MEMBRE
           CORRESPONDANT A UN FLUX DE CHALEUR ET/OU UN APPORT DE MASSE FLUIDE
           (FLUX HYDRAULIQUE) APPLIQUER A UN DOMAINE DE MILIEU CONTINU 2D OU 3D.
           LES FLUX SONT DONNES SOUS FORME DE FONCTION """,
)
