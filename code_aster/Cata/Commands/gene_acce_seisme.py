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

# person_in_charge: irmela.zentner at edf.fr

from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *

GENE_ACCE_SEISME = MACRO(
    nom="GENE_ACCE_SEISME",
    op=OPS("code_aster.MacroCommands.gene_acce_seisme_ops.gene_acce_seisme_ops"),
    sd_prod=table_fonction,
    fr=tr("Generation d'accelerogrammes sismiques "),
    reentrant="n",
    regles=(
        UN_PARMI("DSP", "SPEC_MEDIANE", "SPEC_MOYENNE", "SPEC_UNIQUE", "SPEC_FRACTILE"),
        EXCLUS("MATR_COHE", "COEF_CORR", "PHASE"),
    ),
    #                EXCLUS('DSP', 'PHASE')),
    INIT_ALEA=SIMP(statut="f", typ="I"),
    TITRE=SIMP(statut="f", typ="TXM"),
    INFO=SIMP(statut="f", typ="I", defaut=1, into=(1, 2)),
    PAS_INST=SIMP(statut="o", typ="R"),
    NB_POIN=SIMP(statut="f", typ="I", fr=tr("nombre de points")),
    PESANTEUR=SIMP(
        statut="o",
        typ="R",
        fr=tr("constante de normalisation de ACCE_MAX, ECART_TYPE et INTE_ARIAS (g) ou le spectre"),
    ),
    DUREE_PHASE_FORTE=SIMP(statut="o", typ="R", fr=tr("durée phase forte du signal")),
    NB_TIRAGE=SIMP(statut="f", typ="I", defaut=1, val_min=1, fr=tr("nombre accelerogrammes")),
    FREQ_CORNER=SIMP(
        statut="f", typ="R", val_min=0.0, fr=tr("frequence du filtre frequentiel: corner frequency")
    ),
    FREQ_FILTRE=SIMP(
        statut="f",
        typ="R",
        defaut=0.0,
        val_min=0.0,
        fr=tr("frequence du filtre temporel passe-haut"),
    ),
    FREQ_FILTRE_ZPA=SIMP(
        statut="f",
        typ="R",
        defaut=0.0,
        val_min=0.0,
        fr=tr("frequence du filtre Butterworth passe-bas"),
    ),
    FREQ_PENTE=SIMP(statut="f", typ="R", fr=tr("pente pour l'evolution de la frequence centrale")),
    COEF_CORR=SIMP(
        statut="f",
        typ="R",
        val_min=-1.0,
        val_max=1.0,
        fr=tr("Correlation coefficient for horizontal motion"),
    ),
    a_ratio=BLOC(
        condition="""exists('COEF_CORR') and exists('SPEC_FRACTILE')""",
        RATIO_HV=SIMP(
            statut="f", typ="R", val_min=0.0, fr=tr("Ratio H/V pour le spectre vertical")
        ),
    ),
    MATR_COHE=FACT(
        statut="f",
        MAILLAGE=SIMP(statut="o", typ=maillage_sdaster),
        GROUP_NO_INTERF=SIMP(statut="o", typ=grno),
        TYPE=SIMP(
            statut="o", typ="TXM", into=("MITA_LUCO", "ABRAHAMSON", "ABRA_ROCHER", "ABRA_SOLMOYEN")
        ),
        b_type_coh=BLOC(
            condition="""equal_to('TYPE', 'MITA_LUCO')""",
            VITE_ONDE=SIMP(statut="o", typ="R", val_min=0.0),
            PARA_ALPHA=SIMP(statut="f", typ="R", defaut=0.1),
        ),
    ),
    PHASE=FACT(
        statut="f",
        max=1,
        MAILLAGE=SIMP(statut="o", typ=maillage_sdaster),
        GROUP_NO_INTERF=SIMP(statut="o", typ=grno, fr=tr("groupe des noeud des appuis")),
        DIRECTION=SIMP(statut="o", typ="R", min=3, max=3, fr=tr("direction de propagation")),
        VITE_ONDE=SIMP(
            statut="o", typ="R", val_min=0.0, fr=tr("vitesse de propagation dans la direction")
        ),
        COOR_REFE=SIMP(
            statut="f", typ="R", min=3, max=3, fr=tr("coord de reference pour les phases")
        ),
    ),
    DSP=FACT(
        statut="f",
        max=1,
        AMOR_REDUIT=SIMP(statut="o", typ="R"),
        FREQ_FOND=SIMP(statut="o", typ="R", fr=tr("frequence centrale")),
        #           FREQ_PENTE    =SIMP(statut='f',typ='R',  fr=tr("pente pour l'evolution de la frequence centrale")),
    ),
    SPEC_MEDIANE=FACT(
        statut="f",
        max=1,
        regles=EXCLUS("FREQ_PAS", "LIST_FREQ"),
        SPEC_OSCI=SIMP(statut="o", typ=(fonction_sdaster)),
        AMOR_REDUIT=SIMP(statut="o", typ="R", val_min=0.00001, val_max=1.0),
        FREQ_PAS=SIMP(statut="f", typ="R", fr=tr("pas")),
        LIST_FREQ=SIMP(statut="f", typ=listr8_sdaster),
        NB_ITER=SIMP(
            statut="f", typ="I", val_min=0, fr=tr("nombre d'iterations pour fitter le spectre")
        ),
        ERRE_ZPA=SIMP(
            statut="f", typ="R", defaut=(1.0, 0.2), min=1, max=2, fr=tr("coef et erreur maxi ZPA")
        ),
        ERRE_MAX=SIMP(
            statut="f",
            typ="R",
            defaut=(0.5, 0.2),
            min=1,
            max=2,
            fr=tr("coef et erreur maxi global"),
        ),
        ERRE_RMS=SIMP(
            statut="f", typ="R", defaut=(0.5, 0.2), min=1, max=2, fr=tr("coef et erreur maxi rms")
        ),
        METHODE=SIMP(statut="f", typ="TXM", defaut="HARMO", into=("NIGAM", "HARMO")),
    ),
    SPEC_MOYENNE=FACT(
        statut="f",
        max=1,
        regles=EXCLUS("FREQ_PAS", "LIST_FREQ"),
        SPEC_OSCI=SIMP(statut="o", typ=(fonction_sdaster)),
        AMOR_REDUIT=SIMP(statut="o", typ="R", val_min=0.00001, val_max=1.0),
        FREQ_PAS=SIMP(statut="f", typ="R", fr=tr("pas")),
        LIST_FREQ=SIMP(statut="f", typ=listr8_sdaster),
        NB_ITER=SIMP(
            statut="f", typ="I", val_min=0, fr=tr("nombre d'iterations pour fitter le spectre")
        ),
        ERRE_ZPA=SIMP(
            statut="f", typ="R", defaut=(1.0, 0.2), min=1, max=2, fr=tr("coef et erreur maxi ZPA")
        ),
        ERRE_MAX=SIMP(
            statut="f",
            typ="R",
            defaut=(0.5, 0.2),
            min=1,
            max=2,
            fr=tr("coef et erreur maxi global"),
        ),
        ERRE_RMS=SIMP(
            statut="f", typ="R", defaut=(0.5, 0.2), min=1, max=2, fr=tr("coef et erreur maxi rms")
        ),
        METHODE=SIMP(statut="f", typ="TXM", defaut="HARMO", into=("NIGAM", "HARMO")),
    ),
    SPEC_UNIQUE=FACT(
        statut="f",
        max=1,
        regles=EXCLUS("FREQ_PAS", "LIST_FREQ"),
        ERRE_ZPA=SIMP(
            statut="f", typ="R", defaut=(1.0, 0.2), min=1, max=2, fr=tr("coef et erreur maxi ZPA")
        ),
        ERRE_MAX=SIMP(
            statut="f",
            typ="R",
            defaut=(0.5, 0.2),
            min=1,
            max=2,
            fr=tr("coef et erreur maxi global"),
        ),
        ERRE_RMS=SIMP(
            statut="f", typ="R", defaut=(0.5, 0.2), min=1, max=2, fr=tr("coef et erreur maxi rms")
        ),
        SPEC_OSCI=SIMP(statut="o", typ=(fonction_sdaster)),
        CORR_ZPA=SIMP(statut="f", typ="TXM", defaut="NON", into=("NON", "OUI")),
        AMOR_REDUIT=SIMP(statut="o", typ="R", val_min=0.00001, val_max=1.0),
        FREQ_PAS=SIMP(statut="f", typ="R", fr=tr("pas")),
        LIST_FREQ=SIMP(statut="f", typ=listr8_sdaster),
        NB_ITER=SIMP(
            statut="f", typ="I", val_min=0, fr=tr("nombre d'iterations pour fitter le spectre")
        ),
        METHODE=SIMP(statut="f", typ="TXM", defaut="HARMO", into=("NIGAM", "HARMO")),
    ),
    #
    SPEC_FRACTILE=FACT(
        statut="f",
        max=1,
        regles=(EXCLUS("FREQ_PAS", "LIST_FREQ"),),
        SPEC_OSCI=SIMP(statut="o", typ=(fonction_sdaster)),
        SPEC_1_SIGMA=SIMP(statut="o", typ=(fonction_sdaster)),
        AMOR_REDUIT=SIMP(statut="o", typ="R", val_min=0.00001, val_max=1.0),
        FREQ_PAS=SIMP(statut="f", typ="R", fr=tr("pas")),
        LIST_FREQ=SIMP(statut="f", typ=listr8_sdaster),
    ),
    a_type_dsp=BLOC(
        condition="""exists('DSP')""",
        MODULATION=FACT(
            statut="o",
            max=1,
            regles=(EXCLUS("ACCE_MAX", "INTE_ARIAS", "ECART_TYPE"),),
            TYPE=SIMP(statut="o", typ="TXM", into=("GAMMA", "JENNINGS_HOUSNER", "CONSTANT")),
            ACCE_MAX=SIMP(statut="f", typ="R", fr=tr("PGA: acceleration max au sol (g)")),
            ECART_TYPE=SIMP(statut="f", typ="R", fr=tr("ecart-type")),
            INTE_ARIAS=SIMP(statut="f", typ="R", fr=tr("intensite d'Arias")),
            a_type_mod=BLOC(
                condition="""equal_to('TYPE', 'GAMMA')""",
                INST_INI=SIMP(statut="o", typ="R", fr=tr("instant debut phase forte")),
            ),
        ),
    ),
    b_type_spec=BLOC(
        condition="""not exists('DSP')""",
        MODULATION=FACT(
            statut="o",
            max=1,
            TYPE=SIMP(statut="o", typ="TXM", into=("GAMMA", "JENNINGS_HOUSNER", "CONSTANT")),
            b_type_mod=BLOC(
                condition="""equal_to('TYPE', 'GAMMA')""",
                INST_INI=SIMP(statut="o", typ="R", fr=tr("instant debut phase forte")),
            ),
        ),
    ),
)
