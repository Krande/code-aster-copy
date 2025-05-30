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

# person_in_charge: nicolas.greffet at edf.fr


from cataelem.Tools.base_objects import InputParameter, OutputParameter, Option, CondCalcul
import cataelem.Commons.physical_quantities as PHY
import cataelem.Commons.parameters as SP
import cataelem.Commons.attributes as AT


PVARCPR = InputParameter(phys=PHY.VARI_R)


PCAORIE = InputParameter(
    phys=PHY.CAORIE_R,
    container="CARA!.CARORIEN",
    comment="""  PCAORIE : ORIENTATION LOCALE D'UN ELEMENT DE POUTRE OU DE TUYAU  """,
)


PNBSP_I = InputParameter(
    phys=PHY.NBSP_I, container="CARA!.CANBSP", comment="""  PNBSP_I :  NOMBRE DE SOUS_POINTS  """
)


PCOMPOR = InputParameter(phys=PHY.COMPOR)


MASS_MECA_DIAG = Option(
    para_in=(
        SP.PCACOQU,
        SP.PCADISM,
        SP.PCAGNBA,
        SP.PCAGNPO,
        PCAORIE,
        SP.PCINFDI,
        PCOMPOR,
        SP.PFIBRES,
        SP.PGEOMER,
        SP.PMATERC,
        PNBSP_I,
        PVARCPR,
    ),
    para_out=(SP.PMATUNS, SP.PMATUUR),
    condition=(
        CondCalcul("+", ((AT.PHENO, "ME"), (AT.BORD, "0"))),
        CondCalcul("-", ((AT.PHENO, "ME"), (AT.ABSO, "OUI"))),
    ),
    comment=""" matrice de masse "diagonale".
   Cette matrice diagonale permet (en theorie) d'accelerer les calculs
   lorsqu'il faut resoudre des systemes lineaires avec la matrice de masse
   (puisque la resolution est alors une simple multiplication).

   Aujourd'hui (Novembre 2006), ce gain theorique n'est pas mis en oeuvre
   car toutes les matrices (rigi_meca, mass_meca,mass_meca_diag, ...)
   ont la meme topologie (non-diagonale). La matrice mass_meca_diag contient
   seulement beucoup plus de 0. que les autres.
""",
)
