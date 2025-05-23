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

# person_in_charge: david.haboussa at edf.fr

from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *

DEFI_COMPOR = OPER(
    nom="DEFI_COMPOR",
    op=59,
    sd_prod=compor_sdaster,
    fr=tr("Définir le comportement d'un monocristal, d'un polycristal ou de groupes de fibres"),
    reentrant="n",
    #   Exclusion de MULTIFBRE MONOCRISTAL POLYCRISTAL car la structure de données n'est pas organisée pareil pour ces cas
    regles=(
        UN_PARMI("MONOCRISTAL", "POLYCRISTAL", "MULTIFIBRE"),
        PRESENT_PRESENT("MULTIFIBRE", "GEOM_FIBRE", "MATER_SECT"),
    ),
    #
    # ===========================================================================================================================
    MONOCRISTAL=FACT(
        statut="f",
        max=5,
        MATER=SIMP(statut="o", typ=mater_sdaster, max=1),
        ELAS=SIMP(
            statut="f",
            typ="TXM",
            max=1,
            fr=tr(
                "Donner le nom du mot-clé facteur de DEFI_MATERIAU précisant le comportement élastique (un et un seul)"
            ),
        ),
        ECOULEMENT=SIMP(
            statut="o",
            typ="TXM",
            max=1,
            into=(
                "MONO_VISC1",
                "MONO_VISC2",
                "MONO_DD_KR",
                "MONO_DD_CFC",
                "MONO_DD_CFC_IRRA",
                "MONO_DD_CC",
                "MONO_DD_CC_IRRA",
                "MONO_DD_FAT",
            ),
            fr=tr(
                "Donner le nom du mot-clé facteur de DEFI_MATERIAU précisant le type d'écoulement viscoplastique"
            ),
        ),
        #
        b_non_dd=BLOC(
            condition="""equal_to("ECOULEMENT", 'MONO_VISC1') or equal_to("ECOULEMENT", 'MONO_VISC2')""",
            ECRO_ISOT=SIMP(
                statut="o",
                typ="TXM",
                max=1,
                into=("MONO_ISOT1", "MONO_ISOT2"),
                fr=tr(
                    "Donner le nom du mot-clé facteur de DEFI_MATERIAU précisant le type d'écrouissage isotrope"
                ),
            ),
            ECRO_CINE=SIMP(
                statut="o",
                typ="TXM",
                max=1,
                into=("MONO_CINE1", "MONO_CINE2"),
                fr=tr(
                    "Donner le nom du mot-clé facteur de DEFI_MATERIAU précisant le type d'écrouissage cinématique"
                ),
            ),
            FAMI_SYST_GLIS=SIMP(
                statut="f",
                typ="TXM",
                max=1,
                into=(
                    "OCTAEDRIQUE",
                    "BCC24",
                    "CUBIQUE1",
                    "CUBIQUE2",
                    "ZIRCONIUM",
                    "UNIAXIAL",
                    "UTILISATEUR",
                ),
            ),
            b_util=BLOC(
                condition="""equal_to("FAMI_SYST_GLIS", 'UTILISATEUR')""",
                TABL_SYST_GLIS=SIMP(statut="f", typ=table_sdaster, max=1),
            ),
        ),
        #
        b_dd_kr=BLOC(
            condition="""equal_to("ECOULEMENT", 'MONO_DD_KR')""",
            FAMI_SYST_GLIS=SIMP(
                statut="f", typ="TXM", max=1, into=("BCC24", "UTILISATEUR"), defaut="BCC24"
            ),
            b_util=BLOC(
                condition="""equal_to("FAMI_SYST_GLIS", 'UTILISATEUR') """,
                TABL_SYST_GLIS=SIMP(statut="f", typ=table_sdaster, max=1),
            ),
        ),
        #
        b_ecp_cfc=BLOC(
            condition="""equal_to("ECOULEMENT", 'MONO_DD_FAT')""",
            FAMI_SYST_GLIS=SIMP(
                statut="f", typ="TXM", max=1, into=("OCTAEDRIQUE",), defaut="OCTAEDRIQUE"
            ),
        ),
        #
        b_dd_cfc=BLOC(
            condition="""equal_to("ECOULEMENT", 'MONO_DD_CFC') or equal_to("ECOULEMENT", 'MONO_DD_CFC_IRRA')""",
            FAMI_SYST_GLIS=SIMP(
                statut="f",
                typ="TXM",
                max=1,
                into=("OCTAEDRIQUE", "UTILISATEUR"),
                defaut="OCTAEDRIQUE",
            ),
            b_util=BLOC(
                condition="""equal_to("FAMI_SYST_GLIS", 'UTILISATEUR')""",
                TABL_SYST_GLIS=SIMP(statut="f", typ=table_sdaster, max=1),
            ),
        ),
        #
        b_dd_cc=BLOC(
            condition="""equal_to("ECOULEMENT", 'MONO_DD_CC') or equal_to("ECOULEMENT", 'MONO_DD_CC_IRRA') """,
            FAMI_SYST_GLIS=SIMP(
                statut="f", typ="TXM", max=1, into=("CUBIQUE1", "UTILISATEUR"), defaut="CUBIQUE1"
            ),
            b_util=BLOC(
                condition="""equal_to("FAMI_SYST_GLIS", 'UTILISATEUR')""",
                TABL_SYST_GLIS=SIMP(statut="f", typ=table_sdaster, max=1),
            ),
        ),
    ),
    b_mono=BLOC(
        condition="""exists("MONOCRISTAL")""",
        MATR_INTER=SIMP(statut="f", typ=table_sdaster, max=1),
        ROTA_RESEAU=SIMP(
            statut="f",
            typ="TXM",
            max=1,
            into=("NON", "POST", "CALC"),
            defaut="NON",
            fr=tr("rotation de réseau : NON, POST, CALC"),
        ),
    ),
    #
    # ===========================================================================================================================
    POLYCRISTAL=FACT(
        statut="f",
        max="**",
        regles=(UN_PARMI("ANGL_REP", "ANGL_EULER"),),
        MONOCRISTAL=SIMP(statut="o", typ=compor_sdaster, max=1),
        FRAC_VOL=SIMP(
            statut="o",
            typ="R",
            max=1,
            fr=tr("fraction volumique de la phase correspondant au monocristal"),
        ),
        ANGL_REP=SIMP(
            statut="f",
            typ="R",
            max=3,
            fr=tr("orientation du monocristal : 3 angles nautiques en degrés"),
        ),
        ANGL_EULER=SIMP(
            statut="f",
            typ="R",
            max=3,
            fr=tr("orientation du monocristal : 3 angles d'Euler   en degrés"),
        ),
    ),
    b_poly=BLOC(
        condition="""exists("POLYCRISTAL")""",
        MU_LOCA=SIMP(statut="o", typ="R", max=1),
        LOCALISATION=SIMP(
            statut="f",
            typ="TXM",
            max=1,
            into=("BZ", "BETA"),
            fr=tr("Donner le nom de la règle de localisation"),
        ),
        b_beta=BLOC(
            condition="""equal_to("LOCALISATION", 'BETA')""",
            DL=SIMP(statut="o", typ="R", max=1),
            DA=SIMP(statut="o", typ="R", max=1),
        ),
    ),
    #
    # ===========================================================================================================================
    GEOM_FIBRE=SIMP(
        statut="f",
        max=1,
        typ=gfibre_sdaster,
        fr=tr(
            "Donner le nom du concept regroupant tous les groupes de fibres (issu de DEFI_GEOM_FIBRE)"
        ),
    ),
    MATER_SECT=SIMP(
        statut="f",
        max=1,
        typ=mater_sdaster,
        fr=tr("Donner le nom du matériau pour les caractéristiques homogénéisées sur la section"),
    ),
    MULTIFIBRE=FACT(
        statut="f",
        max="**",
        GROUP_FIBRE=SIMP(statut="o", typ="TXM", max="**"),
        MATER=SIMP(
            statut="o",
            typ=mater_sdaster,
            max=1,
            fr=tr("Donner le nom du matériau pour le groupe de fibres"),
        ),
        RELATION=SIMP(
            statut="f",
            typ="TXM",
            max=1,
            defaut="ELAS",
            into=C_RELATION("DEFI_COMPOR"),
            fr=tr("Donner le nom de la relation incrémentale pour le groupe de fibres"),
        ),
    ),
)
