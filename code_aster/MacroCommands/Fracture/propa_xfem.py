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

# person_in_charge: samuel.geniaut at edf.fr

from ...Cata.Commons import *
from ...Cata.DataStructure import *
from ...Cata.Syntax import *
from ...Objects import XfemCrack
from ...Supervis import ExecuteCommand

PROPA_XFEM_CATA = OPER(
    nom="PROPA_XFEM",
    op=10,
    sd_prod=fiss_xfem,
    reentrant="n",
    fr=tr("Propagation de fissure avec X-FEM"),
    METHODE=SIMP(
        statut="f", typ="TXM", into=("SIMPLEXE", "GEOMETRIQUE", "UPWIND"), defaut="GEOMETRIQUE"
    ),
    OPERATION=SIMP(
        statut="f", typ="TXM", into=("RIEN", "DETECT_COHESIF", "PROPA_COHESIF"), defaut="RIEN"
    ),
    MODELE=SIMP(statut="o", typ=modele_sdaster),
    TEST_MAIL=SIMP(statut="f", typ="TXM", into=("NON", "OUI"), defaut="NON"),
    b_detec=BLOC(condition="OPERATION != 'DETECT_COHESIF' ", DA_MAX=SIMP(statut="o", typ="R")),
    FISS_PROP=SIMP(statut="o", typ=fiss_xfem),
    ZONE_MAJ=SIMP(statut="f", typ="TXM", into=("TOUT", "TORE"), defaut="TORE"),
    RAYON_TORE=SIMP(statut="f", typ="R"),
    LISTE_FISS=SIMP(statut="o", typ=fiss_xfem, min=1, max="**"),
    ANGLE_BETA=SIMP(statut="f", typ="R", max="**"),
    ANGLE_GAMMA=SIMP(statut="f", typ="R", max="**"),
    NOM_PARA_ANGLE=SIMP(statut="f", typ="TXM", into=("BETA", "BETA_GAMMA"), defaut="BETA"),
    VITESSE=SIMP(statut="f", typ="R", max="**"),
    b_pas_cohe=BLOC(
        condition="(OPERATION!= 'PROPA_COHESIF') and (OPERATION != 'DETECT_COHESIF')",
        DA_FISS=SIMP(statut="f", typ="R"),
        NB_CYCLES=SIMP(statut="f", typ="R"),
        RAYON=SIMP(statut="o", typ="R"),
    ),
    b_test_mail_const=BLOC(
        condition="TEST_MAIL == 'OUI' ",
        FISS_INITIALE=SIMP(statut="o", typ=fiss_xfem, max=1),
        DISTANCE=SIMP(statut="o", typ="R", max=1),
        TOLERANCE=SIMP(statut="o", typ="R", max=1),
    ),
    INFO=SIMP(statut="f", typ="I", defaut=0, into=(0, 1, 2)),
)


class CrackPropagation(ExecuteCommand):
    """Command that defines :class:`~code_aster.Objects.XfemCrack`."""

    command_name = "PROPA_XFEM"
    command_cata = PROPA_XFEM_CATA

    def create_result(self, keywords):
        """Initialize the result.

        Arguments:
            keywords (dict): Keywords arguments of user's keywords.
        """
        self._result = XfemCrack(keywords["MODELE"].getMesh())

    def post_exec(self, keywords):
        self._result.setDiscontinuityType(keywords["FISS_PROP"].getDiscontinuityType())
        self._result.updateInternalState()


PROPA_XFEM = CrackPropagation.run
