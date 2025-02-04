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

# person_in_charge: jacques.pellet at edf.fr


from cataelem.Tools.base_objects import InputParameter, OutputParameter, Option, CondCalcul
import cataelem.Commons.physical_quantities as PHY
import cataelem.Commons.parameters as SP
import cataelem.Commons.attributes as AT


PVARCPR = InputParameter(
    phys=PHY.VARI_R,
    container="VOLA!&&CCPARA.VARI_INT_N",
    comment="""  PVARCPR : VARIABLES DE COMMANDE  """,
)


PCAORIE = InputParameter(
    phys=PHY.CAORIE_R,
    container="CARA!.CARORIEN",
    comment="""  PCAORIE : ORIENTATION LOCALE D'UN ELEMENT DE POUTRE OU DE TUYAU  """,
)


PNBSP_I = InputParameter(
    phys=PHY.NBSP_I, container="CARA!.CANBSP", comment="""  PNBSP_I :  NOMBRE DE SOUS_POINTS  """
)


PCOMPOR = InputParameter(
    phys=PHY.COMPOR,
    container="CHMA!.COMPOR",
    comment="""  PCOMPOR :  DESCRIPTION DU COMPORTEMENT DE CHAQUE GROUPE DE FIBRES
           NECESSITE DE FOURNIR LE CONCEPT PRODUIT PAR AFFE_MATERIAU  """,
)


PCNSETO = InputParameter(
    phys=PHY.N1280I,
    container="MODL!.TOPOSE.CNS",
    comment="""  XFEM - CONNECTIVITE DES SOUS-ELEMENTS  """,
)


PLONCHA = InputParameter(
    phys=PHY.N120_I,
    container="MODL!.TOPOSE.LON",
    comment="""  XFEM - NBRE DE TETRAEDRES ET DE SOUS-ELEMENTS  """,
)


PPINTTO = InputParameter(
    phys=PHY.N132_R,
    container="MODL!.TOPOSE.PIN",
    comment=""" XFEM - COORD. POINTS SOMMETS DES SOUS-ELEMENTS """,
)


PHEAVTO = InputParameter(
    phys=PHY.N512_I,
    container="MODL!.TOPOSE.HEA",
    comment=""" XFEM - VALEUR FONCTION HEAVISIDE SUR LES SOUS-ELEMENTS """,
)


PBASLOR = InputParameter(
    phys=PHY.NEUT_R, container="MODL!.BASLOC", comment=""" XFEM - BASE LOCALE AU FOND DE FISSURE """
)


PLSN = InputParameter(
    phys=PHY.NEUT_R, container="MODL!.LNNO", comment=""" XFEM - VALEURS DE LA LEVEL SET NORMALE """
)


PLST = InputParameter(
    phys=PHY.NEUT_R, container="MODL!.LTNO", comment=""" XFEM - VALEURS DE LA LEVEL SET TANGENTE """
)


PSTANO = InputParameter(
    phys=PHY.N120_I,
    container="MODL!.STNO",
    comment=""" XFEM - STATUT DES NOEUDS (ENRICHISSEMENT) """,
)


PPMILTO = InputParameter(phys=PHY.N792_R, container="MODL!.TOPOSE.PMI")


PHEA_NO = InputParameter(phys=PHY.N120_I, container="MODL!.TOPONO.HNO")


PFISNO = InputParameter(
    phys=PHY.NEUT_I,
    container="MODL!.FISSNO",
    comment=""" PFISNO : CONNECTIVITE DES FISSURES ET DES DDL HEAVISIDE """,
)


PSTRXRR = InputParameter(
    phys=PHY.STRX_R,
    container="RESU!STRX_ELGA!N",
    comment="""  PSTRXRR : CHAMPS ELEMENTS DE STRUCTURE INSTANT ACTUEL """,
)


PCONTRR = OutputParameter(phys=PHY.SIEF_R, type="ELGA")


SIEF_ELGA = Option(
    para_in=(
        PBASLOR,
        SP.PCAARPO,
        SP.PCACOQU,
        SP.PCADISK,
        SP.PCAGEPO,
        SP.PCAGNBA,
        SP.PCAGNPO,
        SP.PCAMASS,
        PCAORIE,
        SP.PCINFDI,
        PCNSETO,
        PCOMPOR,
        SP.PDEPLAR,
        SP.PFIBRES,
        PFISNO,
        SP.PGEOMER,
        SP.PHARMON,
        PHEAVTO,
        PHEA_NO,
        PLONCHA,
        PLSN,
        PLST,
        SP.PMATERC,
        PNBSP_I,
        PPINTTO,
        PPMILTO,
        PSTANO,
        PSTRXRR,
        SP.PINSTR,
        PVARCPR,
        SP.PVARCRR,
    ),
    para_out=(SP.PCONTRC, PCONTRR),
    condition=(
        CondCalcul("+", ((AT.PHENO, "ME"), (AT.BORD, "0"))),
        CondCalcul("-", ((AT.PHENO, "ME"), (AT.ABSO, "OUI"))),
        CondCalcul("-", ((AT.PHENO, "ME"), (AT.FLUIDE, "OUI"))),
        CondCalcul("-", ((AT.PHENO, "ME"), (AT.INTERFACE, "OUI"), (AT.TYPMOD2, "INTERFAC"))),
    ),
    comment=""" CALCUL DES CONTRAINTES ET/OU EFFORTS GENERALISES AUX POINTS DE GAUSS
   A PARTIR DES DEPLACEMENTS. LICITE EN LINEAIRE SEULEMENT. """,
)
