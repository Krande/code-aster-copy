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

# person_in_charge: xavier.desroches at edf.fr


from cataelem.Tools.base_objects import InputParameter, OutputParameter, Option, CondCalcul
import cataelem.Commons.physical_quantities as PHY
import cataelem.Commons.parameters as SP
import cataelem.Commons.attributes as AT


PCHHOBS = InputParameter(phys=PHY.N3600R, comment=""" HHO - coefficient base locale""")


CHAR_MECA_FRSU23 = Option(
    para_in=(SP.PDEPLMR, SP.PDEPLPR, SP.PFR2D3D, SP.PGEOMER, SP.PMATERC, PCHHOBS),
    para_out=(SP.PVECTUR,),
    condition=(CondCalcul("+", ((AT.PHENO, "ME"), (AT.BORD, "-1"), (AT.DIM_TOPO_MODELI, "3"))),),
    comment=""" CHAR_MECA_FRSU23 (MOT-CLE : FORCE_FACE): CALCUL DU SECOND MEMBRE
           ELEMENTAIRE CORRESPONDANT A DES FORCES SUIVEUSES SURFACIQUES
           APPLIQUEES SUR UNE FACE D'UN DOMAINE 3D""",
)
