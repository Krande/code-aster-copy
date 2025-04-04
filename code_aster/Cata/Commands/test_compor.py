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

# person_in_charge: david.haboussa at edf.fr

from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *
from ..Commons.c_comportement import compat_syntax

TEST_COMPOR = MACRO(
    nom="TEST_COMPOR",
    op=OPS("code_aster.MacroCommands.test_compor_ops.test_compor_ops"),
    compat_syntax=compat_syntax,
    sd_prod=table_sdaster,
    docu="",
    reentrant="n",
    fr=tr("macro de test des comportements incrementaux dependant de la temperature"),
    OPTION=SIMP(statut="f", typ="TXM", into=("THER", "MECA"), defaut="THER"),
    COMPORTEMENT=C_COMPORTEMENT("SIMU_POINT_MAT"),
    NEWTON=C_NEWTON(),
    CONVERGENCE=C_CONVERGENCE("SIMU_POINT_MAT"),
    b_ther=BLOC(
        condition="""equal_to("OPTION", 'THER')""",
        regles=(EXCLUS("C_PRAG", "D_SIGM_EPSI"),),
        MATER=SIMP(
            statut="o", typ=mater_sdaster, max=1, fr=tr("materiau dependant de la temperature")
        ),
        ALPHA=SIMP(
            statut="o",
            typ=fonction_sdaster,
            fr=tr("coefficient de dilatation fonction de la temperature"),
        ),
        YOUNG=SIMP(
            statut="o", typ=fonction_sdaster, fr=tr("module d'Young fonction de la temperature")
        ),
        LIST_MATER=SIMP(
            statut="o",
            typ=mater_sdaster,
            max="**",
            fr=tr("liste des materiaux constants interpolés à chaque température"),
        ),
        TEMP_INIT=SIMP(statut="o", typ="R", fr=tr("temperature initiale et de reference")),
        TEMP_FIN=SIMP(statut="o", typ="R", fr=tr("temperature finale")),
        INST_FIN=SIMP(statut="f", typ="R", defaut=1.0, fr=tr("instant final")),
        SUPPORT=SIMP(statut="f", typ="TXM", max=1, into=("POINT", "ELEMENT"), defaut=("POINT")),
        NB_VARI=SIMP(statut="o", typ="I", fr=tr("nombre de variables internes - 0 en elasticité")),
        VARI_TEST=SIMP(
            statut="f",
            typ="TXM",
            max="**",
            fr=tr("liste de variables internes à tester - par defaut, toutes"),
        ),
        #           special ecrouissage cinematique
        D_SIGM_EPSI=SIMP(
            statut="f",
            typ=fonction_sdaster,
            fr=tr("module tangent fonction de la temperature- VMIS_CINE_LINE"),
        ),
        C_PRAG=SIMP(
            statut="f",
            typ=fonction_sdaster,
            fr=tr("constante de Prager fonction de la temperature- VMIS_ECMI_*"),
        ),
    ),
    b_meca=BLOC(
        condition="""equal_to("OPTION", 'MECA')""",
        LIST_MATER=SIMP(
            statut="o",
            typ=mater_sdaster,
            max=2,
            min=2,
            fr=tr("liste des materiaux en Pa puis MPa "),
        ),
        YOUNG=SIMP(statut="o", typ="R", fr=tr("module d'Young")),
        POISSON=SIMP(statut="o", typ="R", fr=tr("coef de Poisson")),
        LIST_NPAS=SIMP(
            statut="f",
            typ="I",
            max="**",
            fr=tr("nombre de pas de temps pour chaque discretisation"),
        ),
        LIST_TOLE=SIMP(statut="f", typ="R", max="**"),
        PREC_ZERO=SIMP(statut="f", typ="R", max="**"),
        VARI_TEST=SIMP(
            statut="f",
            typ="TXM",
            max="**",
            defaut=("V1", "VMIS", "TRACE"),
            fr=tr("liste des CMP à tester "),
        ),
        SUPPORT=SIMP(statut="f", typ="TXM", max=1, into=("POINT", "ELEMENT")),
        MODELISATION=SIMP(statut="f", typ="TXM", max=1, into=("3D", "C_PLAN"), defaut="3D"),
        ANGLE=SIMP(
            statut="f",
            typ="R",
            max=1,
            defaut=0.0,
            fr=tr(
                "Rotation de ANGLE autour de Z uniquement, et seulement pour les déformations imposées"
            ),
        ),
        MASSIF=FACT(
            statut="f",
            fr=tr("orientation du materiau (monocristal, orthotropie)"),
            regles=(UN_PARMI("ANGL_REP", "ANGL_EULER"),),
            ANGL_REP=SIMP(statut="f", typ="R", min=1, max=3),
            ANGL_EULER=SIMP(statut="f", typ="R", min=1, max=3),
        ),
        TEST_TANGENTE=SIMP(statut="f", typ="TXM", max=1, into=("OUI", "NON"), defaut="OUI"),
        VERI_MATR_OPTION=FACT(
            statut="f",
            max=1,
            fr=tr("options pour le test de la matrice tangente"),
            VALE_PERT_RELA=SIMP(statut="f", typ="R", defaut=1.0e-5),
            PRECISION=SIMP(statut="f", typ="R", defaut=1.0e-4),
            PREC_ZERO=SIMP(statut="f", typ="R", defaut=1.0e-12),
        ),
    ),
    INFO=SIMP(statut="f", typ="I", defaut=1, into=(1, 2)),
)
