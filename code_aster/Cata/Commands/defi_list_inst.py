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

# person_in_charge: mickael.abbas at edf.fr

from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *

# Bloc pour decoupe automatique
bloc_auto = BLOC(
    fr=tr("Subdivision de type automatique"),
    condition="""equal_to("SUBD_METHODE", 'AUTO')""",
    SUBD_PAS_MINI=SIMP(
        fr=tr("Pas de temps en dessous duquel on ne subdivise plus"),
        statut="f",
        typ="R",
        val_min=0.0,
        max=1,
        defaut=0.0,
    ),
)

# Bloc pour decoupe manuel
bloc_manu = BLOC(
    fr=tr("Subdivision de type manuel"),
    condition="""equal_to("SUBD_METHODE", 'MANUEL')""",
    SUBD_PAS=SIMP(
        fr=tr("Nombre de subdivision d'un pas de temps"),
        statut="f",
        typ="I",
        val_min=2,
        max=1,
        defaut=4,
    ),
    SUBD_NIVEAU=SIMP(
        fr=tr("Nombre maximum de niveau de subdivision d'un pas de temps"),
        statut="f",
        typ="I",
        val_min=0,
        max=1,
        defaut=3,
    ),
    SUBD_PAS_MINI=SIMP(
        fr=tr("Pas de temps en dessous duquel on ne subdivise plus"),
        statut="f",
        typ="R",
        val_min=0.0,
        max=1,
        defaut=0.0,
    ),
)

# Bloc pour decoupe automatique - Cas de la collision
bloc_auto2 = BLOC(
    fr=tr("Subdivision de type automatique"),
    condition="""equal_to("SUBD_METHODE", 'AUTO')""",
    SUBD_INST=SIMP(
        fr=tr("Parametre de decoupe fine du pas de temps"), statut="o", typ="R", val_min=0.0, max=1
    ),
    SUBD_DUREE=SIMP(
        fr=tr("Duree de decoupe apres collision"), statut="o", typ="R", val_min=0.0, max=1
    ),
)

# Bloc pour decoupe du pas de temps
bloc_deco = BLOC(
    fr=tr("Action de decoupe du pas temps"),
    condition="""equal_to("ACTION", 'DECOUPE')""",
    SUBD_METHODE=SIMP(
        fr=tr("Méthode de subdivision des pas de temps en cas de divergence"),
        statut="f",
        typ="TXM",
        max=1,
        into=("MANUEL", "AUTO"),
        defaut="MANUEL",
    ),
    b_deco_manu=bloc_manu,
    b_deco_auto=bloc_auto,
)


# Bloc pour decoupe du pas de temps - special pour collision
bloc_deco2 = BLOC(
    fr=tr("Action de decoupe du pas temps"),
    condition="""equal_to("ACTION", 'DECOUPE')""",
    SUBD_METHODE=SIMP(
        fr=tr("Méthode de subdivision des pas de temps en cas de collision"),
        statut="f",
        typ="TXM",
        max=1,
        into=("MANUEL", "AUTO"),
        defaut="AUTO",
    ),
    b_deco_manu=bloc_manu,
    b_deco_auto=bloc_auto2,
)

# Bloc pour extrapolation du nombre d'iterations de Newton
bloc_supp = BLOC(
    fr=tr("Action d'extrapolation du nombre d'iterations de Newton"),
    condition="""equal_to("ACTION", 'ITER_SUPPL')""",
    PCENT_ITER_PLUS=SIMP(
        fr=tr("Pourcentage d'itérations autorisées en plus"),
        statut="f",
        typ="I",
        val_min=20,
        max=1,
        defaut=50,
    ),
    SUBD_METHODE=SIMP(
        fr=tr("Méthode de subdivision des pas de temps en cas de divergence"),
        statut="f",
        typ="TXM",
        max=1,
        into=("MANUEL", "AUTO"),
        defaut="MANUEL",
    ),
    b_deco_manu=bloc_manu,
    b_deco_auto=bloc_auto,
)

# Bloc pour adaptation du coefficient de penalisation
bloc_pene = BLOC(
    fr=tr("Action d' adaptation du coefficient de penalisation"),
    condition="""equal_to("ACTION", 'ADAPT_COEF_PENA')""",
    COEF_MAXI=SIMP(
        fr=tr("Coefficient multiplicateur maximum du coefficient de penalisation"),
        statut="f",
        typ="R",
        val_min=1.0,
        max=1,
        defaut=1e12,
    ),
)


DEFI_LIST_INST = OPER(
    nom="DEFI_LIST_INST",
    op=28,
    sd_prod=list_inst,
    reentrant="n",
    fr=tr("Définition de la gestion de la liste d'instants"),
    MODELE=SIMP(statut="f", typ=(modele_sdaster)),
    METHODE=SIMP(
        fr=tr("Methode de definition de la liste d'instants"),
        statut="f",
        typ="TXM",
        max=1,
        into=("MANUEL", "AUTO"),
        defaut="MANUEL",
    ),
    # ----------------------------------------------------------------------------------------------------------------------------------
    # mot-cle pour la definition a priori de la liste d'instant
    # ----------------------------------------------------------------------------------------------------------------------------------
    b_manuel=BLOC(
        fr=tr("Liste d'instants donnée par l'utilisateur"),
        condition="""equal_to("METHODE", 'MANUEL') """,
        DEFI_LIST=FACT(
            fr=tr("Definition a priori de la liste d'instants"),
            statut="o",
            max=1,
            regles=(
                UN_PARMI("LIST_INST", "VALE", "RESULTAT"),
                PRESENT_PRESENT("RESULTAT", "SUBD_PAS"),
            ),
            VALE=SIMP(statut="f", typ="R", max="**"),
            LIST_INST=SIMP(statut="f", typ=listr8_sdaster),
            RESULTAT=SIMP(statut="f", typ=resultat_sdaster),
            SUBD_PAS=SIMP(statut="f", typ="I", max=1, val_min=1),
        ),  # end fkw_defi_list
    ),  # end b_manuel
    b_auto=BLOC(
        fr=tr("Gestion automatique de la liste d'instants"),
        condition="""(equal_to("METHODE", 'AUTO')) """,
        DEFI_LIST=FACT(
            fr=tr("Definition a priori de la liste d'instants"),
            statut="o",
            max=1,
            regles=(UN_PARMI("LIST_INST", "VALE"),),
            VALE=SIMP(statut="f", typ="R", max="**"),
            LIST_INST=SIMP(statut="f", typ=listr8_sdaster),
            PAS_MINI=SIMP(statut="f", typ="R", max=1, val_min=1.0e-12),
            PAS_MAXI=SIMP(statut="f", typ="R", max=1),
            NB_PAS_MAXI=SIMP(statut="f", typ="I", max=1, val_max=1000000, defaut=1000000),
        ),  # end fkw_defi_list
    ),  # end b_auto
    # ----------------------------------------------------------------------------------------------------------------------------------
    # mot-cle pour le comportement en cas d'echec (on doit recommencer le meme instant)
    # ----------------------------------------------------------------------------------------------------------------------------------
    ECHEC=FACT(
        fr=tr("Comportement en cas d'echec"),
        statut="d",
        max="**",
        EVENEMENT=SIMP(
            fr=tr("Type de l'evenement"),
            statut="f",
            typ="TXM",
            max=1,
            into=(
                "ERREUR",
                "DELTA_GRANDEUR",
                "COLLISION",
                "RESI_MAXI",
                "INTERPENETRATION",
                "DIVE_RESI",
                "INSTABILITE",
            ),
            defaut="ERREUR",
        ),
        b_erreur=BLOC(
            fr=tr("Event: erreur ou iter_maxi"),
            condition="""equal_to("EVENEMENT", 'ERREUR') """,
            ACTION=SIMP(
                fr=tr("Actions possibles"),
                statut="f",
                max=1,
                typ="TXM",
                into=("ARRET", "DECOUPE", "ITER_SUPPL"),
                defaut="DECOUPE",
            ),
            b_deco=bloc_deco,
            b_supp=bloc_supp,
        ),
        b_edelta=BLOC(
            fr=tr("Event: l'increment d'une composante d'un champ depasse le seuil"),
            condition="""equal_to("EVENEMENT", 'DELTA_GRANDEUR') """,
            regles=(EXCLUS("GROUP_MA", "GROUP_NO"),),
            GROUP_MA=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
            GROUP_NO=SIMP(statut="f", typ=grno, validators=NoRepeat(), max="**"),
            VALE_REF=SIMP(fr=tr("Valeur du seuil"), statut="o", typ="R", max=1),
            NOM_CHAM=SIMP(
                fr=tr("Nom du champ"),
                statut="o",
                typ="TXM",
                max=1,
                into=("DEPL", "VARI_ELGA", "SIEF_ELGA"),
            ),
            NOM_CMP=SIMP(fr=tr("Nom de la composante"), statut="o", typ="TXM", max=1),
            ACTION=SIMP(
                fr=tr("Actions possibles"),
                statut="f",
                max=1,
                typ="TXM",
                into=("ARRET", "DECOUPE"),
                defaut="DECOUPE",
            ),
            b_deco=bloc_deco,
        ),
        b_colli=BLOC(
            fr=tr("Event: collision"),
            condition="""equal_to("EVENEMENT", 'COLLISION') """,
            ACTION=SIMP(
                fr=tr("Actions possibles"),
                statut="f",
                max=1,
                typ="TXM",
                into=("ARRET", "DECOUPE"),
                defaut="DECOUPE",
            ),
            b_deco2=bloc_deco2,
        ),
        b_penetration=BLOC(
            fr=tr("Event: interpenetration des surfaces en contact"),
            condition="""equal_to("EVENEMENT", 'INTERPENETRATION') """,
            PENE_MAXI=SIMP(fr=tr("Valeur maxi de l'interpenetration"), statut="o", typ="R", max=1),
            ACTION=SIMP(
                fr=tr("Actions possibles"),
                statut="f",
                max=1,
                typ="TXM",
                into=("ARRET", "ADAPT_COEF_PENA"),
                defaut="ADAPT_COEF_PENA",
            ),
            b_pene=bloc_pene,
        ),
        b_dive_resi=BLOC(
            fr=tr("Event: divergence du residu"),
            condition="""equal_to("EVENEMENT", 'DIVE_RESI') """,
            ACTION=SIMP(
                fr=tr("Actions possibles"),
                statut="f",
                max=1,
                typ="TXM",
                into=("DECOUPE",),
                defaut="DECOUPE",
            ),
            b_deco=bloc_deco,
        ),
        b_resi_maxi=BLOC(
            fr=tr("Event: residu troup grande"),
            condition="""equal_to("EVENEMENT", 'RESI_MAXI') """,
            RESI_GLOB_MAXI=SIMP(fr=tr("Valeur du seuil"), statut="o", typ="R", max=1),
            ACTION=SIMP(
                fr=tr("Actions possibles"),
                statut="f",
                max=1,
                typ="TXM",
                into=("DECOUPE",),
                defaut="DECOUPE",
            ),
            b_deco=bloc_deco,
        ),
        b_instabilite=BLOC(
            fr=tr("Event: instabilite"),
            condition="""equal_to("EVENEMENT", 'INSTABILITE') """,
            ACTION=SIMP(
                fr=tr("Actions possibles"),
                statut="f",
                max=1,
                typ="TXM",
                into=("ARRET", "CONTINUE"),
                defaut="CONTINUE",
            ),
        ),
    ),
    # ----------------------------------------------------------------------------------------------------------------------------------
    # Mot-cle pour le comportement en cas de succes (on a bien converge)
    # ----------------------------------------------------------------------------------------------------------------------------------
    b_adap=BLOC(
        condition="""equal_to("METHODE", 'AUTO')""",
        ADAPTATION=FACT(
            fr=tr("Parametres de l'evenement declencheur de l'adaptation du pas de temps"),
            statut="d",
            max="**",
            EVENEMENT=SIMP(
                fr=tr("Nom de l'evenement declencheur de l'adaptation"),
                statut="f",
                max=1,
                typ="TXM",
                into=("SEUIL", "TOUT_INST", "AUCUN"),
                defaut="SEUIL",
            ),
            b_adap_seuil=BLOC(
                fr=tr("Seuil d'adaptation"),
                condition="""equal_to("EVENEMENT", 'SEUIL') """,
                NB_INCR_SEUIL=SIMP(statut="f", typ="I", defaut=2),
                NOM_PARA=SIMP(
                    statut="f", typ="TXM", into=("NB_ITER_NEWTON",), defaut="NB_ITER_NEWTON"
                ),
                CRIT_COMP=SIMP(statut="f", typ="TXM", into=("LT", "GT", "LE", "GE"), defaut="LE"),
                b_vale_I=BLOC(
                    fr=tr("Valeur entiere"),
                    condition="""equal_to("NOM_PARA", 'NB_ITER_NEWTON') """,
                    VALE_I=SIMP(statut="f", typ="I"),
                ),
            ),
            #
            #  Parametres du mode de calcul de dt+
            #      dans le cas FIXE            :(deltaT+) = (deltaT-)x(1+PCENT_AUGM/100)
            #      dans le cas DELTA_GRANDEUR : (deltaT+) = (deltaT-)x(VALREF/deltaVAL) : l'acceleration est inversement proportionnelle
            #                                                                             a la variation de la grandeur
            #      dans le cas ITER_NEWTON    : (deltaT+) = (deltaT-) x sqrt(VALREF/N)  : l'acceleration est inversement proportionnelle
            #                                                                             au nombre d'iterations de Newton precedent
            MODE_CALCUL_TPLUS=SIMP(
                fr=tr("Parametres du mode de calcul de dt+"),
                statut="f",
                max=1,
                typ="TXM",
                into=("FIXE", "DELTA_GRANDEUR", "ITER_NEWTON", "IMPLEX"),
                defaut="FIXE",
            ),
            b_mfixe=BLOC(
                fr=tr("Mode de calcul de dt+: fixe"),
                condition="""equal_to("MODE_CALCUL_TPLUS", 'FIXE') """,
                PCENT_AUGM=SIMP(statut="f", max=1, typ="R", defaut=100.0, val_min=-100.0),
            ),
            b_mdelta=BLOC(
                fr=tr("Mode de calcul de dt+: increment d'une grandeur"),
                condition="""equal_to("MODE_CALCUL_TPLUS", 'DELTA_GRANDEUR') """,
                regles=(EXCLUS("GROUP_MA", "GROUP_NO"),),
                GROUP_MA=SIMP(statut="f", typ=grma, validators=NoRepeat(), max="**"),
                GROUP_NO=SIMP(statut="f", typ=grno, validators=NoRepeat(), max="**"),
                VALE_REF=SIMP(statut="o", max=1, typ="R"),
                NOM_CHAM=SIMP(
                    statut="o", max=1, typ="TXM", into=("DEPL", "VARI_ELGA", "SIEF_ELGA")
                ),
                NOM_CMP=SIMP(statut="o", max=1, typ="TXM"),
            ),
            b_mitnew=BLOC(
                fr=tr("Mode de calcul de dt+: nb iterations de Newton"),
                condition="""equal_to("MODE_CALCUL_TPLUS", 'ITER_NEWTON') """,
                NB_ITER_NEWTON_REF=SIMP(statut="o", max=1, typ="I"),
            ),
            # les schemas pre-definis :
            #  abaqus :
            #      EVENEMENT       ='SEUIL'
            #      NB_INCR_SEUIL     = 2
            #      NOM_PARA          ='NB_ITER_NEWTON'
            #      CRIT_COMP         ='LE'
            #      VALE_I            = 5
            #      MODE_CALCUL_TPLUS ='FIXE'
            #      PCENT_AUGM        = 50.
            #  Zebulon 1 :
            #      EVENEMENT       ='TOUT_INST'
            #      MODE_CALCUL_TPLUS ='DELTA_GRANDEUR'
            #      VALE_REF          = valref
            #      NOM_CHAM          ='VARI_ELGA'
            #      NOM_CMP           ='V1'
            #  Zebulon 2 :
            #      EVENEMENT       ='TOUT_INST'
            #      MODE_CALCUL_TPLUS ='ITER_NEWTON'
            #      NB_ITER_NEWTON_REF= nc
            #  Tough2 :
            #      EVENEMENT       ='SEUIL'
            #      NB_INCR_SEUIL     = 1
            #      NOM_PARA          ='NB_ITER_NEWTON'
            #      CRIT_COMP         ='LE'
            #      VALE_I            = n
            #      MODE_CALCUL_TPLUS ='FIXE'
            #      PCENT_AUGM        = 100.
            #  Oliver :
            #      EVENEMENT       ='TOUT_INST'
            #      MODE_CALCUL_TPLUS ='FORMULE'
            #      NOM_SCHEMA        ='OLIVER'
        ),
    ),
    # ----------------------------------------------------------------------------------------------------------------------------------
    INFO=SIMP(statut="f", typ="I", defaut=1, into=(1, 2)),
)
