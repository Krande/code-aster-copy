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


def calc_miss_sdprod(self, TYPE_RESU, **kwargs):
    """Typage des structures de données produites"""
    if kwargs.get("__all__"):
        return (None, table_sdaster, harm_gene, tran_gene, char_meca)

    if TYPE_RESU in ("TABLE", "TABLE_CONTROL"):
        return table_sdaster
    elif TYPE_RESU == "HARM_GENE":
        return harm_gene
    elif TYPE_RESU == "TRAN_GENE":
        return tran_gene
    elif TYPE_RESU == "CHARGE":
        return char_meca
    else:
        return None


CALC_MISS = MACRO(
    nom="CALC_MISS",
    op=OPS("code_aster.MacroCommands.calc_miss_ops.calc_miss_ops"),
    sd_prod=calc_miss_sdprod,
    fr=tr("Préparation des données, exécution du logiciel Miss3D, et post-traitement"),
    regles=(EXCLUS("TABLE_SOL", "MATER_SOL"),),
    TYPE_RESU=SIMP(
        statut="o",
        typ="TXM",
        into=(
            "HARM_GENE",
            "TRAN_GENE",
            "TABLE",
            "TABLE_CONTROL",
            "FICHIER",
            "FICHIER_TEMPS",
            "CHARGE",
        ),
        fr=tr("Type de résultat produit en post-traitement. FICHIER : pas de post-traitement"),
    ),
    PROJET=SIMP(statut="f", typ="TXM", defaut="MODELE", fr=tr("Nom de l'étude Miss")),
    REPERTOIRE=SIMP(statut="f", typ="TXM", fr=tr("Répertoire de travail de Miss")),
    VERSION=SIMP(
        statut="f",
        typ="TXM",
        into=("V6.7", "V6.6", "V6.5"),
        defaut="V6.7",
        fr=tr("Version de Miss utilisée"),
    ),
    TABLE_SOL=SIMP(statut="f", typ=table_sdaster, fr=tr("Table des propriétés du sol stratifié")),
    MATER_SOL=FACT(
        statut="f",
        fr=tr("Propriétés du sol homogène"),
        E=SIMP(statut="o", typ="R", val_min=0.0),
        NU=SIMP(statut="o", typ="R", val_min=-1.0, val_max=0.5),
        RHO=SIMP(statut="o", typ="R", val_min=0.0),
        AMOR_HYST=SIMP(statut="f", typ="R", val_min=0.0, val_max=1.0),
    ),
    MATER_FLUIDE=FACT(
        statut="f",
        fr=tr("Propriétés du fluide (requis si ISSF='OUI')"),
        RHO=SIMP(statut="o", typ="R", val_min=0.0),
        CELE=SIMP(statut="o", typ="R", val_min=0.0),
        AMOR_BETA=SIMP(statut="f", typ="R", val_min=0.0, val_max=1.0),
        DEMI_ESPACE=SIMP(
            statut="f",
            typ="TXM",
            defaut="OUI",
            into=("OUI", "NON"),
            fr=tr("Demi-espace de fluide avec surface libre ou non"),
        ),
    ),
    # pas de post-traitement
    b_basic=BLOC(
        condition="""is_in("TYPE_RESU", ('FICHIER', 'TABLE_CONTROL'))""",
        regles=(
            UN_PARMI("MACR_ELEM_DYNA", "BASE_MODALE"),
            ENSEMBLE("GROUP_MA_FLU_STR", "GROUP_MA_FLU_SOL"),
            EXCLUS("SOURCE_SOL", "SOURCE_FLUIDE"),
        ),
        MACR_ELEM_DYNA=SIMP(
            statut="f", typ=macr_elem_dyna, fr=tr("Macro élément produit en amont")
        ),
        BASE_MODALE=SIMP(statut="f", typ=mode_meca, fr=tr("Base de modes")),
        b_base_modale=BLOC(
            condition="""exists("BASE_MODALE")""",
            MATR_RIGI=SIMP(statut="f", typ=(matr_asse_depl_r, matr_asse_depl_c)),
            MATR_MASS=SIMP(statut="f", typ=matr_asse_depl_r),
        ),
        AMOR_REDUIT=SIMP(statut="f", typ="R", max="**"),
        GROUP_MA_INTERF=SIMP(
            statut="o", typ=grma, max="**", fr=tr("Groupe de mailles de l'interface")
        ),
        GROUP_MA_FLU_STR=SIMP(
            statut="f", typ=grma, max="**", fr=tr("Groupe de mailles fluide-structure")
        ),
        GROUP_MA_FLU_SOL=SIMP(
            statut="f", typ=grma, max="**", fr=tr("Groupe de mailles fluide-sol")
        ),
        GROUP_MA_SOL_SOL=SIMP(statut="f", typ=grma, max="**", fr=tr("Groupe de mailles sol-sol")),
        UNITE_IMPR_ASTER=SIMP(
            statut="f",
            typ=UnitType(),
            inout="out",
            fr=tr("Unité des résultats transmis par code_aster à Miss"),
        ),
        UNITE_RESU_IMPE=SIMP(
            statut="f",
            typ=UnitType(),
            inout="out",
            fr=tr("Unité logique des impédances écrites par Miss"),
        ),
        UNITE_RESU_FORC=SIMP(
            statut="f",
            typ=UnitType(),
            inout="out",
            fr=tr("Unité logique des forces sismiques écrites par Miss"),
        ),
        SOURCE_SOL=FACT(
            statut="f",
            max="**",
            fr=tr("Source ponctuelle dans le sol"),
            POINT=SIMP(statut="o", typ="R", min=3, max=3, fr=tr("Position de la source")),
            DIRECTION=SIMP(statut="o", typ="R", min=3, max=3, fr=tr("Direction de la source")),
        ),
        SOURCE_FLUIDE=FACT(
            statut="f",
            max="**",
            fr=tr("Source ponctuelle dans le fluide"),
            POINT=SIMP(statut="o", typ="R", min=3, max=3, fr=tr("Position de la source")),
        ),
    ),
    # post-traitement : passage du domaine de Laplace au domaine temporel
    b_fichier_temps=BLOC(
        condition="""equal_to("TYPE_RESU", 'FICHIER_TEMPS')""",
        regles=(
            UN_PARMI("MACR_ELEM_DYNA", "BASE_MODALE"),
            ENSEMBLE("GROUP_MA_FLU_STR", "GROUP_MA_FLU_SOL"),
            AU_MOINS_UN("UNITE_RESU_RIGI", "UNITE_RESU_AMOR", "UNITE_RESU_MASS"),
            PRESENT_PRESENT("UNITE_RESU_AMOR", "MATR_GENE"),
            PRESENT_PRESENT("UNITE_RESU_MASS", "MATR_GENE"),
        ),
        MACR_ELEM_DYNA=SIMP(
            statut="f", typ=macr_elem_dyna, fr=tr("Macro élément produit en amont")
        ),
        BASE_MODALE=SIMP(statut="f", typ=mode_meca, fr=tr("Base de modes")),
        b_base_modale=BLOC(
            condition="""exists("BASE_MODALE")""",
            MATR_RIGI=SIMP(statut="f", typ=(matr_asse_depl_r, matr_asse_depl_c)),
            MATR_MASS=SIMP(statut="f", typ=matr_asse_depl_r),
        ),
        AMOR_REDUIT=SIMP(statut="f", typ="R", max="**"),
        GROUP_MA_INTERF=SIMP(
            statut="o", typ=grma, max="**", fr=tr("Groupe de mailles de l'interface")
        ),
        GROUP_MA_FLU_STR=SIMP(
            statut="f", typ=grma, max="**", fr=tr("Groupe de mailles fluide-structure")
        ),
        GROUP_MA_FLU_SOL=SIMP(
            statut="f", typ=grma, max="**", fr=tr("Groupe de mailles fluide-sol")
        ),
        GROUP_MA_SOL_SOL=SIMP(statut="f", typ=grma, max="**", fr=tr("Groupe de mailles sol-sol")),
        UNITE_IMPR_ASTER=SIMP(
            statut="f",
            typ=UnitType(),
            inout="out",
            fr=tr("Unité des résultats transmis par code_aster à Miss"),
        ),
        UNITE_RESU_RIGI=SIMP(statut="f", typ=UnitType(), inout="out"),
        UNITE_RESU_AMOR=SIMP(statut="f", typ=UnitType(), inout="out"),
        UNITE_RESU_MASS=SIMP(statut="f", typ=UnitType(), inout="out"),
        INST_FIN=SIMP(statut="f", typ="R", fr=tr("Instant final du calcul")),
        PAS_INST=SIMP(statut="f", typ="R", fr=tr("Pas de temps du calcul")),
        INST_ECRI_FIN=SIMP(statut="f", typ="R", fr=tr("Instant final d'écriture des impédances")),
        FACTEUR_INTERPOL=SIMP(
            statut="f",
            typ="I",
            fr=tr(
                "Valeur du pas d'échantillonnage et \
                                du facteur de réduction du temps de calcul"
            ),
            defaut=1,
            val_min=1,
        ),
        PCENT_FREQ_CALCUL=SIMP(
            statut="f",
            typ="R",
            fr=tr(
                "Valeur correspondante au ratio 100*Ns/Nt, où \
                 Ns est le nombre d'échantillons sans interpoler et Nt le nombre total d'échantillons"
            ),
            defaut=0,
            val_min=0,
            val_max=100,
        ),
        PRECISION=SIMP(statut="f", typ="R", defaut=0.000001),
        COEF_SURECH=SIMP(statut="f", typ="R", defaut=1.35),
        COEF_MULT=SIMP(statut="f", typ="R", defaut=1.0),
        TYPE_FICHIER_TEMPS=SIMP(statut="f", typ="TXM", into=("BINAIRE", "ASCII"), defaut="ASCII"),
        MATR_GENE=FACT(
            statut="f",
            max=1,
            DECOMP_IMPE=SIMP(
                statut="f", typ="TXM", defaut="PRODUIT", into=("PRODUIT", "SANS_PRODUIT")
            ),
            AMOR_HYST=SIMP(
                statut="o",
                typ="TXM",
                into=("DANS_IMPEDANCE", "DANS_MATR_AMOR"),
                fr=tr("Indique comment l'amortissement hysteretique est pris en compte"),
            ),
            b_amor_nohyst=BLOC(
                condition="""equal_to("AMOR_HYST", 'DANS_MATR_AMOR')""",
                MATR_MASS=SIMP(statut="f", typ=(matr_asse_gene_r, matr_asse_depl_r)),
                MATR_RIGI=SIMP(
                    statut="f", typ=(matr_asse_gene_r, matr_asse_gene_c, matr_asse_depl_r)
                ),
                MATR_AMOR=SIMP(
                    statut="o", typ=(matr_asse_gene_r, matr_asse_gene_c, matr_asse_depl_r)
                ),
            ),
            b_amor_hyst=BLOC(
                condition="""equal_to("AMOR_HYST", 'DANS_IMPEDANCE')""",
                regles=(AU_MOINS_UN("MATR_MASS", "MATR_RIGI", "MATR_AMOR"),),
                MATR_MASS=SIMP(statut="f", typ=(matr_asse_gene_r, matr_asse_depl_r)),
                MATR_RIGI=SIMP(
                    statut="f", typ=(matr_asse_gene_r, matr_asse_gene_c, matr_asse_depl_r)
                ),
                MATR_AMOR=SIMP(
                    statut="f", typ=(matr_asse_gene_r, matr_asse_gene_c, matr_asse_depl_r)
                ),
            ),
        ),
        EXCIT_SOL=FACT(
            statut="f",
            max=1,
            regles=(AU_MOINS_UN("CHAM_X", "CHAM_Y", "CHAM_Z"),),
            UNITE_RESU_FORC=SIMP(
                statut="o",
                typ=UnitType(),
                inout="out",
                fr=tr("Unité logique des forces sismiques écrites par Miss"),
            ),
            NOM_CHAM=SIMP(statut="f", typ="TXM", defaut="DEPL", into=("ACCE", "VITE", "DEPL")),
            CHAM_X=SIMP(statut="f", typ=fonction_sdaster),
            CHAM_Y=SIMP(statut="f", typ=fonction_sdaster),
            CHAM_Z=SIMP(statut="f", typ=fonction_sdaster),
        ),
    ),
    # si post-traitement
    b_donnees=BLOC(
        condition="""not is_in("TYPE_RESU", ('FICHIER', 'FICHIER_TEMPS', 'TABLE_CONTROL', 'CHARGE'))""",
        regles=(
            ENSEMBLE("GROUP_MA_FLU_STR", "GROUP_MA_FLU_SOL"),
            UN_PARMI("MATR_AMOR", "AMOR_REDUIT"),
        ),
        MACR_ELEM_DYNA=SIMP(
            statut="f", typ=macr_elem_dyna, fr=tr("Macro élément produit en amont")
        ),
        BASE_MODALE=SIMP(statut="o", typ=mode_meca, fr=tr("Base de modes")),
        MATR_RIGI=SIMP(statut="o", typ=(matr_asse_depl_r, matr_asse_depl_c)),
        MATR_MASS=SIMP(statut="o", typ=matr_asse_depl_r),
        MATR_AMOR=SIMP(statut="f", typ=matr_asse_depl_r),
        AMOR_REDUIT=SIMP(statut="f", typ="R", max="**"),
        GROUP_MA_INTERF=SIMP(
            statut="o", typ=grma, max="**", fr=tr("Groupe de mailles de l'interface")
        ),
        GROUP_MA_FLU_STR=SIMP(
            statut="f", typ=grma, max="**", fr=tr("Groupe de mailles fluide-structure")
        ),
        GROUP_MA_FLU_SOL=SIMP(
            statut="f", typ=grma, max="**", fr=tr("Groupe de mailles fluide-sol")
        ),
        GROUP_MA_SOL_SOL=SIMP(statut="f", typ=grma, max="**", fr=tr("Groupe de mailles sol-sol")),
        UNITE_IMPR_ASTER=SIMP(
            statut="f",
            typ=UnitType(),
            inout="out",
            fr=tr("Unité des résultats transmis par code_aster à Miss"),
        ),
        UNITE_RESU_IMPE=SIMP(
            statut="f", typ=UnitType(), inout="in", fr=tr("Unité logique des impédances à relire.")
        ),
        UNITE_RESU_FORC=SIMP(
            statut="f",
            typ=UnitType(),
            inout="in",
            fr=tr("Unité logique des forces sismiques à relire"),
        ),
    ),
    # Paramètres du calcul Miss
    PARAMETRE=FACT(
        statut="f",
        regles=(
            PRESENT_PRESENT("OFFSET_MAX", "OFFSET_NB"),
            PRESENT_PRESENT("FREQ_MIN", "FREQ_MAX", "FREQ_PAS"),
            UN_PARMI("FREQ_MIN", "LIST_FREQ", "FREQ_IMAG"),
        ),
        FREQ_MIN=SIMP(statut="f", typ="R"),
        FREQ_MAX=SIMP(statut="f", typ="R"),
        FREQ_PAS=SIMP(statut="f", typ="R"),
        LIST_FREQ=SIMP(statut="f", typ="R", max="**"),
        FREQ_IMAG=SIMP(statut="f", typ="R"),
        Z0=SIMP(statut="f", typ="R", defaut=0.0),
        TYPE=SIMP(statut="f", typ="TXM", into=("BINAIRE", "ASCII"), defaut="ASCII"),
        ISSF=SIMP(statut="f", typ="TXM", into=("OUI", "NON"), defaut="NON"),
        ALLU=SIMP(statut="f", typ="R", defaut=0.0),
        SURF=SIMP(statut="f", typ="TXM", into=("OUI", "NON"), defaut="NON"),
        DREF=SIMP(statut="f", typ="R"),
        OFFSET_MAX=SIMP(statut="f", typ="R"),
        OFFSET_NB=SIMP(statut="f", typ="I"),
        AUTO=SIMP(statut="f", typ="TXM", into=("OUI", "NON"), defaut="NON"),
        b_auto=BLOC(
            condition="""equal_to("AUTO", 'OUI')""",
            regles=(ENSEMBLE("SPEC_MAX", "SPEC_NB"),),
            OPTION_DREF=SIMP(statut="f", typ="TXM", into=("OUI", "NON"), defaut="NON"),
            OPTION_RFIC=SIMP(statut="f", typ="TXM", into=("OUI", "NON"), defaut="NON"),
            RFIC=SIMP(statut="f", typ="R"),
            SPEC_MAX=SIMP(statut="f", typ="R"),
            SPEC_NB=SIMP(statut="f", typ="I"),
            COEF_OFFSET=SIMP(statut="f", typ="I", defaut=12),
        ),
        b_noauto=BLOC(
            condition="""equal_to("AUTO", 'NON')""",
            regles=(ENSEMBLE("SPEC_MAX", "SPEC_NB"),),
            ALGO=SIMP(statut="f", typ="TXM", into=("DEPL", "REGU")),
            RFIC=SIMP(statut="f", typ="R", defaut=0.0),
            SPEC_MAX=SIMP(statut="f", typ="R"),
            SPEC_NB=SIMP(statut="f", typ="I"),
        ),
    ),
    # Post-traitement type 1 - tran_gene
    b_post_tran_gene=BLOC(
        condition="""equal_to("TYPE_RESU", 'TRAN_GENE')""",
        regles=(ENSEMBLE("INST_FIN", "PAS_INST"),),
        MODELE=SIMP(statut="o", typ=(modele_sdaster)),
        GROUP_NO=SIMP(statut="f", typ=grno, max="**"),
        INST_FIN=SIMP(statut="f", typ="R", fr=tr("Instant final du calcul")),
        PAS_INST=SIMP(statut="f", typ="R", fr=tr("Pas de temps du calcul")),
        b_post_tran_gene_temp=BLOC(
            condition="""exists("INST_FIN")""",
            regles=(
                AU_MOINS_UN("ACCE_X", "ACCE_Y", "ACCE_Z", "DEPL_X", "DEPL_Y", "DEPL_Z"),
                PRESENT_ABSENT("ACCE_X", "DEPL_X", "DEPL_Y", "DEPL_Z"),
                PRESENT_ABSENT("ACCE_Y", "DEPL_X", "DEPL_Y", "DEPL_Z"),
                PRESENT_ABSENT("ACCE_Z", "DEPL_X", "DEPL_Y", "DEPL_Z"),
            ),
            ACCE_X=SIMP(statut="f", typ=fonction_sdaster),
            ACCE_Y=SIMP(statut="f", typ=fonction_sdaster),
            ACCE_Z=SIMP(statut="f", typ=fonction_sdaster),
            DEPL_X=SIMP(statut="f", typ=fonction_sdaster),
            DEPL_Y=SIMP(statut="f", typ=fonction_sdaster),
            DEPL_Z=SIMP(statut="f", typ=fonction_sdaster),
        ),
        b_post_tran_gene_frreq=BLOC(
            condition="""not exists("INST_FIN")""",
            regles=(
                AU_MOINS_UN("ACCE_X", "ACCE_Y", "ACCE_Z", "DEPL_X", "DEPL_Y", "DEPL_Z"),
                PRESENT_ABSENT("ACCE_X", "DEPL_X", "DEPL_Y", "DEPL_Z"),
                PRESENT_ABSENT("ACCE_Y", "DEPL_X", "DEPL_Y", "DEPL_Z"),
                PRESENT_ABSENT("ACCE_Z", "DEPL_X", "DEPL_Y", "DEPL_Z"),
            ),
            ACCE_X=SIMP(statut="f", typ=fonction_c),
            ACCE_Y=SIMP(statut="f", typ=fonction_c),
            ACCE_Z=SIMP(statut="f", typ=fonction_c),
            DEPL_X=SIMP(statut="f", typ=fonction_c),
            DEPL_Y=SIMP(statut="f", typ=fonction_c),
            DEPL_Z=SIMP(statut="f", typ=fonction_c),
        ),
    ),
    # Post-traitement type 1 - harm_gene
    b_post_harm_gene=BLOC(
        condition="""equal_to("TYPE_RESU", 'HARM_GENE')""",
        regles=(PRESENT_ABSENT("EXCIT_HARMO", "INST_FIN"), ENSEMBLE("INST_FIN", "PAS_INST")),
        MODELE=SIMP(statut="o", typ=(modele_sdaster)),
        GROUP_NO=SIMP(statut="f", typ=grno, max="**"),
        INST_FIN=SIMP(statut="f", typ="R", fr=tr("Instant final du calcul")),
        PAS_INST=SIMP(statut="f", typ="R", fr=tr("Pas de temps du calcul")),
        # identique à EXCIT de DYNA_LINE_HARM au type attendu pour VECT_ASSE près
        EXCIT_HARMO=FACT(
            statut="f",
            max="**",
            regles=(
                UN_PARMI("VECT_ASSE", "CHARGE"),
                UN_PARMI("FONC_MULT", "FONC_MULT_C", "COEF_MULT", "COEF_MULT_C"),
            ),
            VECT_ASSE=SIMP(statut="f", typ=cham_no_sdaster),
            CHARGE=SIMP(statut="f", typ=char_meca),
            FONC_MULT_C=SIMP(statut="f", typ=(fonction_c, formule_c)),
            COEF_MULT_C=SIMP(statut="f", typ="C"),
            FONC_MULT=SIMP(statut="f", typ=(fonction_sdaster, nappe_sdaster, formule)),
            COEF_MULT=SIMP(statut="f", typ="R"),
            PHAS_DEG=SIMP(statut="f", typ="R", defaut=0.0),
            PUIS_PULS=SIMP(statut="f", typ="I", defaut=0),
        ),
        b_post_harm_gene_temp=BLOC(
            condition="""not exists("EXCIT_HARMO") and exists("INST_FIN")""",
            regles=(
                AU_MOINS_UN("ACCE_X", "ACCE_Y", "ACCE_Z", "DEPL_X", "DEPL_Y", "DEPL_Z"),
                PRESENT_ABSENT("ACCE_X", "DEPL_X", "DEPL_Y", "DEPL_Z"),
                PRESENT_ABSENT("ACCE_Y", "DEPL_X", "DEPL_Y", "DEPL_Z"),
                PRESENT_ABSENT("ACCE_Z", "DEPL_X", "DEPL_Y", "DEPL_Z"),
            ),
            ACCE_X=SIMP(statut="f", typ=fonction_sdaster),
            ACCE_Y=SIMP(statut="f", typ=fonction_sdaster),
            ACCE_Z=SIMP(statut="f", typ=fonction_sdaster),
            DEPL_X=SIMP(statut="f", typ=fonction_sdaster),
            DEPL_Y=SIMP(statut="f", typ=fonction_sdaster),
            DEPL_Z=SIMP(statut="f", typ=fonction_sdaster),
        ),
        b_post_harm_gene_freq=BLOC(
            condition="""not exists("EXCIT_HARMO") and not exists("INST_FIN")""",
            regles=(
                AU_MOINS_UN("ACCE_X", "ACCE_Y", "ACCE_Z", "DEPL_X", "DEPL_Y", "DEPL_Z"),
                PRESENT_ABSENT("ACCE_X", "DEPL_X", "DEPL_Y", "DEPL_Z"),
                PRESENT_ABSENT("ACCE_Y", "DEPL_X", "DEPL_Y", "DEPL_Z"),
                PRESENT_ABSENT("ACCE_Z", "DEPL_X", "DEPL_Y", "DEPL_Z"),
            ),
            ACCE_X=SIMP(statut="f", typ=fonction_c),
            ACCE_Y=SIMP(statut="f", typ=fonction_c),
            ACCE_Z=SIMP(statut="f", typ=fonction_c),
            DEPL_X=SIMP(statut="f", typ=fonction_c),
            DEPL_Y=SIMP(statut="f", typ=fonction_c),
            DEPL_Z=SIMP(statut="f", typ=fonction_c),
        ),
    ),
    # Post-traitement type 2
    b_post_table=BLOC(
        condition="""equal_to("TYPE_RESU", 'TABLE')""",
        regles=(ENSEMBLE("INST_FIN", "PAS_INST"),),
        MODELE=SIMP(statut="o", typ=(modele_sdaster)),
        GROUP_NO=SIMP(
            statut="o", typ=grno, max="**", fr=tr("Liste des groupes de noeud de post-traitement")
        ),
        INST_FIN=SIMP(statut="f", typ="R", fr=tr("Instant final du calcul")),
        PAS_INST=SIMP(statut="f", typ="R", fr=tr("Pas de temps du calcul")),
        NORME=SIMP(statut="o", typ="R", fr=tr("Valeur de la norme du spectre d'oscillateur")),
        AMOR_SPEC_OSCI=SIMP(
            statut="o", typ="R", max="**", fr=tr("Amortissement du spectre d'oscillateur")
        ),
        LIST_FREQ_SPEC_OSCI=SIMP(
            statut="f",
            typ=listr8_sdaster,
            fr=tr("Fréquences utilisées pour le calcul du spectre d'oscillateur"),
        ),
        b_post_table_temp=BLOC(
            condition="""exists("INST_FIN")""",
            regles=(AU_MOINS_UN("ACCE_X", "ACCE_Y", "ACCE_Z"),),
            ACCE_X=SIMP(statut="f", typ=fonction_sdaster),
            ACCE_Y=SIMP(statut="f", typ=fonction_sdaster),
            ACCE_Z=SIMP(statut="f", typ=fonction_sdaster),
        ),
        b_post_table_freq=BLOC(
            condition="""not exists("INST_FIN")""",
            regles=(AU_MOINS_UN("ACCE_X", "ACCE_Y", "ACCE_Z"),),
            ACCE_X=SIMP(statut="f", typ=fonction_c),
            ACCE_Y=SIMP(statut="f", typ=fonction_c),
            ACCE_Z=SIMP(statut="f", typ=fonction_c),
        ),
    ),
    # Post-traitement type 3 - points de controle
    b_post_control=BLOC(
        condition="""equal_to("TYPE_RESU", 'TABLE_CONTROL')""",
        regles=(ENSEMBLE("INST_FIN", "PAS_INST"), ENSEMBLE("NORME", "AMOR_SPEC_OSCI")),
        GROUP_MA_CONTROL=SIMP(
            statut="f", typ=grma, max="**", fr=tr("Groupe de mailles des points de contrôle")
        ),
        INST_FIN=SIMP(statut="f", typ="R", fr=tr("Instant final du calcul")),
        PAS_INST=SIMP(statut="f", typ="R", fr=tr("Pas de temps du calcul")),
        NORME=SIMP(statut="f", typ="R", fr=tr("Valeur de la norme du spectre d'oscillateur")),
        AMOR_SPEC_OSCI=SIMP(
            statut="f", typ="R", max="**", fr=tr("Amortissement du spectre d'oscillateur")
        ),
        LIST_FREQ_SPEC_OSCI=SIMP(
            statut="f",
            typ=listr8_sdaster,
            fr=tr("Fréquences utilisées pour le calcul du spectre d'oscillateur"),
        ),
        TOUT_CHAM=SIMP(statut="f", typ="TXM", into=("OUI",)),
        b_post_controle_temp=BLOC(
            condition="""exists("INST_FIN")""",
            ACCE_X=SIMP(statut="f", typ=fonction_sdaster),
            ACCE_Y=SIMP(statut="f", typ=fonction_sdaster),
            ACCE_Z=SIMP(statut="f", typ=fonction_sdaster),
        ),
        b_post_controle_freq=BLOC(
            condition="""not exists("INST_FIN")""",
            ACCE_X=SIMP(statut="f", typ=fonction_c),
            ACCE_Y=SIMP(statut="f", typ=fonction_c),
            ACCE_Z=SIMP(statut="f", typ=fonction_c),
        ),
    ),
    # post-traitement : creation d'une charge sismique temporelle
    b_charge=BLOC(
        condition="""equal_to("TYPE_RESU", 'CHARGE')""",
        regles=(EXCLUS("NOEUD_AFFE", "GROUP_NO_AFFE")),
        MODELE=SIMP(statut="o", typ=(modele_sdaster)),
        FONC_SIGNAL=SIMP(statut="o", typ=(fonction_sdaster)),
        NOM_CMP=SIMP(statut="o", typ="TXM", into=("DX", "DY", "DZ")),
        ISSF=SIMP(statut="f", typ="TXM", into=("OUI", "NON"), defaut="NON"),
        NOEUD_AFFE=SIMP(statut="c", typ=no, validators=NoRepeat(), max="**"),
        GROUP_NO_AFFE=SIMP(statut="f", typ=grno, validators=NoRepeat(), max="**"),
        UNITE_RESU_FORC=SIMP(statut="o", typ=UnitType(), inout="in"),
        FREQ_MAX=SIMP(statut="f", typ="R"),
        VARI=SIMP(statut="f", typ="TXM", into=("OUI", "NON"), defaut="NON"),
        MATR_GENE=FACT(
            statut="f",
            BASE=SIMP(statut="o", typ=mode_meca),
            NUME_DDL_GENE=SIMP(statut="o", typ=nume_ddl_gene),
        ),
        PRECISION=SIMP(statut="f", typ="R", defaut=0.999),
        INTERF=FACT(
            statut="f",
            GROUP_NO_INTERF=SIMP(statut="o", typ=grno),
            MODE_INTERF=SIMP(statut="o", typ="TXM", into=("CORP_RIGI", "TOUT", "QUELCONQUE")),
        ),
        MATR_COHE=FACT(
            statut="f",
            TYPE=SIMP(statut="o", typ="TXM", into=("MITA_LUCO", "ABRAHAMSON")),
            b_type_coh=BLOC(
                condition="""equal_to("TYPE", 'MITA_LUCO') """,
                VITE_ONDE=SIMP(statut="f", typ="R", defaut=600.0),
                PARA_ALPHA=SIMP(statut="f", typ="R", defaut=0.5),
            ),
        ),
        UNITE_RESU_IMPE=SIMP(statut="f", typ=UnitType(), inout="in"),
        TYPE=SIMP(statut="f", typ="TXM", into=("BINAIRE", "ASCII"), defaut="ASCII"),
    ),
    INFO=SIMP(statut="f", typ="I", defaut=1, into=(1, 2)),
)
