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

# person_in_charge: mathieu.courtois at edf.fr


from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *

DEFI_SOL_MISS = MACRO(
    nom="DEFI_SOL_MISS",
    op=OPS("code_aster.MacroCommands.defi_sol_miss_ops.defi_sol_miss_ops"),
    sd_prod=table_sdaster,
    fr=tr("Définition des données de sol pour Miss"),
    reentrant="n",
    regles=(UN_PARMI("COUCHE", "COUCHE_AUTO"),),
    MATERIAU=FACT(
        statut="o",
        max="**",
        fr=tr("Définition des matériaux"),
        E=SIMP(statut="o", typ="R", fr=tr("Module d'Young")),
        NU=SIMP(statut="o", typ="R", fr=tr("Coefficient de Poisson")),
        RHO=SIMP(statut="o", typ="R", fr=tr("Masse volumique")),
        AMOR_HYST=SIMP(statut="o", typ="R", fr=tr("Coefficient d'amortissement")),
    ),
    COUCHE=FACT(
        statut="f",
        max="**",
        fr=tr("Définition des couches"),
        regles=(AU_MOINS_UN("EPAIS", "SUBSTRATUM"),),
        SUBSTRATUM=SIMP(statut="f", typ="TXM", into=("OUI", "NON")),
        EPAIS=SIMP(statut="f", typ="R", fr=tr("Epaisseur de la couche")),
        RECEPTEUR=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
        SOURCE=SIMP(statut="f", typ="TXM", defaut="NON", into=("OUI", "NON")),
        NUME_MATE=SIMP(statut="o", typ="I", fr=tr("Numéro du matériau")),
    ),
    COUCHE_AUTO=FACT(
        statut="f",
        max=1,
        fr=tr("Définition automatique des couches"),
        Z0=SIMP(statut="f", typ="R", max=1, fr=tr("Position de la surface libre")),
        HOMOGENE=SIMP(statut="o", typ="TXM", into=("OUI", "NON")),
        SURF=SIMP(statut="f", typ="TXM", into=("OUI", "NON"), defaut="NON"),
        EPAIS_PHYS=SIMP(statut="o", typ="R", max="**", fr=tr("Epaisseur des couches")),
        b_stratifie=BLOC(
            condition="""equal_to("HOMOGENE", 'NON')""",
            NUME_MATE_SUBSTRATUM=SIMP(statut="o", typ="I", fr="Numéro du matériau du substratum"),
            NUME_MATE=SIMP(statut="o", typ="I", max="**", fr=tr("Numéro du matériau")),
        ),
        b_surf=BLOC(
            condition="""equal_to("SURF", 'OUI')""",
            regles=(PRESENT_PRESENT("GROUP_MA_CONTROL", "MAILLAGE"),),
            GROUP_MA_CONTROL=SIMP(
                statut="f", typ=grma, max=1, fr=tr("Groupe de mailles des points de contrôle")
            ),
            MAILLAGE=SIMP(statut="f", typ=maillage_sdaster),
        ),
        b_enfonce=BLOC(
            condition="""equal_to("SURF", 'NON')""",
            regles=(UN_PARMI("GROUP_MA", "GROUP_NO"),),
            GROUP_MA=SIMP(
                statut="f", typ=grma, max=1, fr=tr("Groupe de mailles donnant les cotes verticales")
            ),
            GROUP_NO=SIMP(
                statut="f", typ=grno, max=1, fr=tr("Groupe de noeuds donnant les cotes verticales")
            ),
            NOMBRE_RECEPTEUR=SIMP(
                statut="f", typ="I", defaut=4, max=1, fr="Nombre de récépteurs par element"
            ),
            GROUP_MA_INTERF=SIMP(
                statut="f", typ=grma, max=1, fr="Groupe de mailles de l'interface"
            ),
            GROUP_MA_CONTROL=SIMP(
                statut="f", typ=grma, max=1, fr=tr("Groupe de mailles des points de contrôle")
            ),
            MAILLAGE=SIMP(statut="o", typ=maillage_sdaster),
            TOLERANCE=SIMP(statut="f", typ="R", defaut=1.0e-5),
            DECALAGE_AUTO=SIMP(statut="f", typ="TXM", into=("OUI", "NON"), defaut="OUI"),
        ),
    ),
    TITRE=SIMP(statut="f", typ="TXM", fr=tr("Titre de la table produite")),
    INFO=SIMP(statut="f", typ="I", defaut=1, into=(1, 2)),
)
