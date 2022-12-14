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

from ..Commons import *
from ..Language.DataStructure import *
from ..Language.Syntax import *


def macr_adap_mail_prod(self, MAJ_CHAM, ADD_CHAM, ADAPTATION, **args):
    if args.get("__all__"):
        return (
            [None],
            [None, maillage_sdaster],
            [None, maillage_sdaster],
            [None, cham_no_sdaster],
            [None, cham_elem],
            [None, carte_sdaster],
        )
    #
    # 0. Typage des structures produites
    #
    # print args
    if "MAILLAGE_NP1" in args:
        if args["MAILLAGE_NP1"] is not None:
            maillage_np1 = args["MAILLAGE_NP1"]
            self.type_sdprod(maillage_np1, maillage_sdaster)
    #
    if "MAILLAGE_NP1_ANNEXE" in args:
        if args["MAILLAGE_NP1_ANNEXE"] is not None:
            maillage_np1_annexe = args["MAILLAGE_NP1_ANNEXE"]
            self.type_sdprod(maillage_np1_annexe, maillage_sdaster)
    #
    # print "MAJ_CHAM =", MAJ_CHAM
    if MAJ_CHAM is not None:
        # Remarque : la liste qui suit doit etre conforme à son homologue de LIRE_CHAMP
        for ch in MAJ_CHAM:
            t = ch["TYPE_CHAM"]
            location, quantity, typ = t.split("_")

            if location == "NOEU":
                self.type_sdprod(ch["CHAM_MAJ"], cham_no_sdaster)
            else:
                self.type_sdprod(ch["CHAM_MAJ"], cham_elem)
    #
    # print "ADD_CHAM =", ADD_CHAM
    if ADD_CHAM is not None:
        for ch in ADD_CHAM:
            self.type_sdprod(ch["CHAM_GD"], carte_sdaster)
    #
    return None


MACR_ADAP_MAIL = MACRO(
    nom="MACR_ADAP_MAIL",
    op=OPS("code_aster.MacroCommands.macr_adap_mail_ops.macr_adap_mail_ops"),
    sd_prod=macr_adap_mail_prod,
    fr=tr("Adapter un maillage avec le logiciel HOMARD."),
    #
    # 1. Le niveau d'information
    #
    INFO=SIMP(statut="f", typ="I", defaut=1, into=(1, 2, 3, 4)),
    #
    # 2. Les maillages
    # 2.1. Quel que soit le type de traitement, il faut donner le concept du maillage initial (entree)
    #
    MAILLAGE_N=SIMP(statut="o", typ=maillage_sdaster, fr=tr("Maillage avant adaptation")),
    #
    # 2.2. Si ce n'est pas une simple lecture :
    #
    b_maillage=BLOC(
        condition=""" (not equal_to("ADAPTATION", 'LECTURE')) """,
        fr=tr("Lectures de champs aux points de Gauss ou aux noeuds par element."),
        IDENTIFICATION=SIMP(
            statut="f",
            typ="TXM",
            fr=tr(
                "Nom du maillage dont on souhaite reprendre "
                "l'historique de raffinement dans le cas où "
                "MAILLAGE_N en est une légère modification."
            ),
        ),
        #
        # 2.2.1. Le concept du maillage final (sortie)
        #
        MAILLAGE_NP1=SIMP(statut="o", typ=CO, fr=tr("Maillage après adaptation")),
        #
        # 2.2.2. Eventuellement, on peut produire un maillage annexe
        #      Actuellement, c'est le maillage n+1, mais de degré différent.
        #
        MAILLAGE_NP1_ANNEXE=SIMP(statut="f", typ=CO, fr=tr("Maillage annexe après adaptation")),
        #
    ),
    #
    # 3. Le pilotage de l'adaptation, avec les variantes suivantes :
    #  . Raffinement et deraffinement, selon un champ
    #  . Raffinement seul, selon un champ
    #  . Deraffinement seul, selon un champ
    #  . Raffinement et deraffinement, selon des zones geometriques
    #  . Raffinement uniforme : toutes les mailles sont divisées
    #  . Deraffinement uniforme : toutes les mailles sont regroupées
    #  . Modification : le maillage subit des transformations specifiques
    #  . Rien : le maillage est le meme a la sortie et a l'entree
    #
    ADAPTATION=SIMP(
        statut="o",
        typ="TXM",
        into=(
            "RAFF_DERA",
            "RAFFINEMENT",
            "DERAFFINEMENT",
            "RAFF_DERA_ZONE",
            "RAFFINEMENT_UNIFORME",
            "DERAFFINEMENT_UNIFORME",
            "MODIFICATION",
            "LECTURE",
            "RIEN",
        ),
        fr=tr("Pilotage de l'adaptation : selon un champ ou uniforme."),
    ),
    #
    # 4. Pour de l'adaptation libre, il faut un champ d'indicateur
    #
    b_champ=BLOC(
        condition="""is_in("ADAPTATION", ('RAFF_DERA','RAFFINEMENT','DERAFFINEMENT') )""",
        fr=tr(
            "Pour une adaptation libre, choix du champ définissant la région à raffiner/déraffiner"
        ),
        #
        regles=(UN_PARMI("CHAM_GD", "RESULTAT_N")),
        #
        # 4.1. Reperage de la région a raffiner a l'aide d'un champ
        #
        # 4.1.1. Sous forme de champ de grandeur
        #
        CHAM_GD=SIMP(
            statut="f",
            typ=cham_gd_sdaster,
            fr=tr("Champ de grandeur Code_Aster pilotant l'adaptation"),
        ),
        #
        # 4.1.2. Sous forme de concept resultat_sdaster
        #
        RESULTAT_N=SIMP(
            statut="f",
            typ=(evol_elas, evol_noli, evol_ther),
            fr=tr("Concept résultat Code_Aster contenant le champ"),
        ),
        #
        b_champ_adaptation=BLOC(
            condition="""(exists("RESULTAT_N"))""",
            NOM_CHAM=SIMP(
                statut="o",
                typ="TXM",
                validators=NoRepeat(),
                into=C_NOM_CHAM_INTO(),
                fr=tr("Champ dans le résultat"),
            ),
        ),
        #
        # 4.1.4. La ou les composantes retenues
        #
        NOM_CMP=SIMP(
            statut="f",
            typ="TXM",
            validators=NoRepeat(),
            max="**",
            fr=tr("Liste des composante(s) retenue(s) pour le champ."),
        ),
        #
        # 4.1.5. Le paramètre temporel pour le champ
        #
        b_parametre_temporel=BLOC(
            condition="""(exists("RESULTAT_N"))""",
            fr=tr("Choix éventuel du paramètre temporel pour le champ"),
            #
            regles=(EXCLUS("NUME_ORDRE", "INST"),),
            #
            # 4.1.5.1. Soit le numero d'ordre
            #
            NUME_ORDRE=SIMP(statut="f", typ="I", val_min=0, fr=tr("Numéro d ordre")),
            #
            # 4.1.5.2. Soit l'instant
            # 4.1.5.2.1. Sa valeur
            #
            INST=SIMP(statut="f", typ="R", fr=tr("Instant associé")),
            #
            # 4.1.5.2.2. La précision du choix de l'instant
            #
            b_precision=BLOC(
                condition="""(exists("INST"))""",
                fr=tr("Choix de la précision du choix de l'instant"),
                CRITERE=SIMP(
                    statut="f",
                    typ="TXM",
                    defaut="RELATIF",
                    into=("RELATIF", "ABSOLU"),
                    fr=tr("Critère de précision sur le choix de l'instant associé"),
                ),
                b_prec_rela=BLOC(
                    condition="""(equal_to("CRITERE", 'RELATIF'))""",
                    PRECISION=SIMP(
                        statut="f",
                        typ="R",
                        defaut=1.0e-6,
                        fr=tr("Précision relative sur le choix de l'instant associé"),
                    ),
                ),
                b_prec_abso=BLOC(
                    condition="""(equal_to("CRITERE", 'ABSOLU'))""",
                    PRECISION=SIMP(
                        statut="o",
                        typ="R",
                        fr=tr("Précision absolue sur le choix de l'instant associé"),
                    ),
                ),
            ),
            #
        ),
        #
        # 4.1.6. Usage des composantes : maximum, maximum de la valeur absolue, ou de la norme L2, ou de la norme infinie
        #
        USAGE_CMP=SIMP(
            statut="f",
            typ="TXM",
            defaut="NORME_L2",
            into=("ABSOLU", "NORME_L2", "NORME_INFINIE", "RELATIF"),
            fr=tr(
                "Valeur absolue de la composante, ou norme du champ, ou valeur relative de la composante"
            ),
        ),
        #
        # 4.1.7. Usage du champ : la valeur par maille ou le max du saut entre mailles
        #
        USAGE_CHAMP=SIMP(
            statut="f",
            typ="TXM",
            defaut="MAILLE",
            into=("MAILLE", "SAUT"),
            fr=tr("Usage du champ : la valeur par maille ou le saut entre mailles voisines"),
        ),
        #
        # 4.1.8. Initialisation de l'adaptation : raffinement ou déraffinement
        #
        ADAP_INIT=SIMP(
            statut="f",
            typ="TXM",
            defaut="GARDER",
            into=("GARDER", "RAFFINER", "DERAFFINER"),
            fr=tr(
                "Initialisation de l'adaptation dans les régions sans indicateur : garder, raffiner ou déraffiner"
            ),
        ),
        #
    ),
    #
    # 5. Les criteres pour de l'adaptation libre avec un champ :
    #        absolu, relatif, en proportion d'entite
    # 5.1. Pour le raffinement :
    #
    b_critere_de_raffinement=BLOC(
        condition="""is_in("ADAPTATION", ('RAFF_DERA','RAFFINEMENT') )""",
        fr=tr("Critère de raffinement."),
        #
        regles=(UN_PARMI("CRIT_RAFF_ABS", "CRIT_RAFF_REL", "CRIT_RAFF_PE", "CRIT_RAFF_MS"),),
        #
        CRIT_RAFF_ABS=SIMP(statut="f", typ="R", fr=tr("Critère absolu")),
        CRIT_RAFF_REL=SIMP(
            statut="f",
            typ="R",
            val_min=0.0,
            val_max=1.0,
            fr=tr("Critère relatif : fraction réelle entre 0. et 1."),
        ),
        CRIT_RAFF_PE=SIMP(
            statut="f",
            typ="R",
            val_min=0.0,
            val_max=1.0,
            fr=tr("Pourcentage de mailles : fraction réelle entre 0. et 1."),
        ),
        CRIT_RAFF_MS=SIMP(
            statut="f", typ="R", fr=tr("Critère absolu valant moyenne + n*sigma, n étant > 0")
        ),
    ),
    #
    # 5.2. Pour le deraffinement :
    #
    b_critere_de_deraffinement=BLOC(
        condition="""is_in("ADAPTATION", ('RAFF_DERA','DERAFFINEMENT'))""",
        fr=tr("Critère de déraffinement."),
        #
        regles=(UN_PARMI("CRIT_DERA_ABS", "CRIT_DERA_REL", "CRIT_DERA_PE", "CRIT_DERA_MS"),),
        #
        CRIT_DERA_ABS=SIMP(statut="f", typ="R", fr=tr("Critère absolu")),
        CRIT_DERA_REL=SIMP(
            statut="f",
            typ="R",
            val_min=0.0,
            val_max=1.0,
            fr=tr("Critère relatif : fraction réelle entre 0. et 1."),
        ),
        CRIT_DERA_PE=SIMP(
            statut="f",
            typ="R",
            val_min=0.0,
            val_max=1.0,
            fr=tr("Pourcentage de mailles : fraction réelle entre 0. et 1."),
        ),
        CRIT_DERA_MS=SIMP(
            statut="f", typ="R", fr=tr("Critère absolu valant moyenne - n*sigma, n étant > 0")
        ),
    ),
    #
    # 6. Pour de l'adaptation par zone, définitions des zones
    #
    b_zone=BLOC(
        condition=""" (equal_to("ADAPTATION", 'RAFF_DERA_ZONE')) """,
        fr=tr("Pour une adaptation selon une zone à raffiner"),
        #
        ZONE=FACT(
            statut="o",
            min=1,
            max="**",
            fr=tr("Définition de zones à raffiner ou déraffiner."),
            #
            # 6.1. Type de la zone
            #
            TYPE=SIMP(
                statut="o",
                typ="TXM",
                into=(
                    "RECTANGLE",
                    "BOITE",
                    "DISQUE",
                    "SPHERE",
                    "CYLINDRE",
                    "DISQUE_PERCE",
                    "TUYAU",
                ),
                fr=tr("Type de la zone"),
            ),
            #
            # 6.2. Usage de la zone
            #
            USAGE=SIMP(
                statut="f",
                typ="TXM",
                into=("RAFFINEMENT", "DERAFFINEMENT"),
                defaut="RAFFINEMENT",
                fr=tr("Zone pour raffiner ou déraffiner"),
            ),
            #
            # Ne sachant pas exploiter les blocs, je mets des regles
            #
            regles=(
                AU_MOINS_UN("X_MINI", "X_CENTRE", "HAUTEUR"),
                EXCLUS("X_MINI", "X_CENTRE", "HAUTEUR"),
                EXCLUS("Z_MINI", "X_CENTRE", "HAUTEUR"),
                EXCLUS("X_MINI", "Z_CENTRE", "HAUTEUR"),
                EXCLUS("Z_MINI", "Z_CENTRE", "HAUTEUR"),
                EXCLUS("X_MINI", "RAYON"),
                EXCLUS("Z_MINI", "RAYON"),
                EXCLUS("RAYON", "RAYON_INT"),
            ),
            #
            # 6.3. Une boite rectangulaire ou parallelepipedique
            # 6.3.1. Incontournables
            #
            X_MINI=SIMP(statut="f", typ="R", fr=tr("Abscisse minimum de la boite")),
            X_MAXI=SIMP(statut="f", typ="R", fr=tr("Abscisse maximum de la boite")),
            Y_MINI=SIMP(statut="f", typ="R", fr=tr("Ordonnée minimum de la boite")),
            Y_MAXI=SIMP(statut="f", typ="R", fr=tr("Abscisse maximum de la boite")),
            #
            # 6.3.2. Complement pour une boite parallelepipedique
            #
            Z_MINI=SIMP(statut="f", typ="R", fr=tr("Cote minimum de la boite")),
            Z_MAXI=SIMP(statut="f", typ="R", fr=tr("Cote maximum de la boite")),
            #
            # 6.4. Rayon pour un disque, une sphere ou un cylindre
            #
            RAYON=SIMP(statut="f", typ="R", val_min=0.0, fr=tr("Rayon")),
            #
            # 6.5. Pour un disque plein ou perce, une sphere
            # 6.5.1. Incontournables
            #
            X_CENTRE=SIMP(
                statut="f", typ="R", fr=tr("Abscisse du centre du disque ou de la sphère")
            ),
            Y_CENTRE=SIMP(
                statut="f", typ="R", fr=tr("Ordonnée du centre du disque ou de la sphère")
            ),
            #
            # 6.5.2. Complement pour une sphere
            #
            Z_CENTRE=SIMP(statut="f", typ="R", fr=tr("Cote du centre de la sphère")),
            #
            # 6.6. Rayons interieur et exterieur pour un disque perce ou un tuyau
            #
            RAYON_INT=SIMP(statut="f", typ="R", val_min=0.0, fr=tr("Rayon intérieur")),
            RAYON_EXT=SIMP(statut="f", typ="R", val_min=0.0, fr=tr("Rayon extérieur")),
            #
            # 6.7. Un cylindre ou un tuyau
            #
            X_AXE=SIMP(statut="f", typ="R", fr=tr("Abscisse du vecteur directeur de l'axe")),
            Y_AXE=SIMP(statut="f", typ="R", fr=tr("Ordonnée du vecteur directeur de l'axe")),
            Z_AXE=SIMP(statut="f", typ="R", fr=tr("Cote du vecteur directeur de l'axe")),
            X_BASE=SIMP(statut="f", typ="R", fr=tr("Abscisse d'un point de la base, sur l'axe")),
            Y_BASE=SIMP(statut="f", typ="R", fr=tr("Ordonnée d'un point de la base, sur l'axe")),
            Z_BASE=SIMP(statut="f", typ="R", fr=tr("Cote d'un point de la base, sur l'axe")),
            HAUTEUR=SIMP(statut="f", typ="R", val_min=0.0, fr=tr("Hauteur")),
            #
        ),
        #
    ),
    #
    # 7. Les niveaux extremes pour le maillage adapte
    # 7.1. Pour le raffinement :
    #
    b_niveau_maximum=BLOC(
        condition="""is_in("ADAPTATION", ('RAFF_DERA','RAFF_DERA_ZONE','RAFFINEMENT','RAFFINEMENT_UNIFORME') )""",
        fr=tr("Profondeur maximale de raffinement"),
        #
        NIVE_MAX=SIMP(
            statut="f", typ="I", val_min=1, fr=tr("Niveau maximum de profondeur de raffinement")
        ),
        #
        DIAM_MIN=SIMP(statut="f", typ="R", val_min=0.0, fr=tr("Diamètre minimal des mailles")),
        #
    ),
    #
    # 7.2. Pour le deraffinement :
    #
    b_niveau_minimum=BLOC(
        condition="""is_in("ADAPTATION", ('RAFF_DERA','RAFF_DERA_ZONE','DERAFFINEMENT','DERAFFINEMENT_UNIFORME') )""",
        fr=tr("Niveau minimum de profondeur de déraffinement"),
        NIVE_MIN=SIMP(
            statut="f", typ="I", val_min=0, fr=tr("Niveau minimum de profondeur de déraffinement")
        ),
    ),
    #
    # 8. Filtrage de l'adaptation par des groupes
    #
    b_filtrage_par_des_groupes=BLOC(
        condition="""is_in("ADAPTATION", ('RAFF_DERA','RAFFINEMENT','RAFFINEMENT_UNIFORME','RAFF_DERA_ZONE','DERAFFINEMENT','DERAFFINEMENT_UNIFORME') )""",
        fr=tr("Filtrage de l'adaptation par des groupes."),
        #
        GROUP_MA=SIMP(
            statut="f",
            typ=grma,
            validators=NoRepeat(),
            max="**",
            fr=tr("Liste des groupes de mailles pour le filtrage de l'adaptation."),
        ),
        #
        GROUP_NO=SIMP(
            statut="f",
            typ=grno,
            validators=NoRepeat(),
            max="**",
            fr=tr("Liste des groupes de noeuds pour le filtrage de l'adaptation."),
        ),
    ),
    #
    # 9. Suivi d'une frontière
    #
    # 9.1. Definition d'une frontière par un maillage (valable seulement pour des frontières 1D)
    #
    MAILLAGE_FRONTIERE=SIMP(
        statut="f", typ=maillage_sdaster, fr=tr("Maillage de la frontière discrète à suivre")
    ),
    #
    b_FRONTIERE=BLOC(
        condition=""" exists("MAILLAGE_FRONTIERE") """,
        fr=tr("Information complémentaire sur la frontière discrète"),
        #
        GROUP_MA_FRONT=SIMP(
            statut="f",
            typ=grma,
            validators=NoRepeat(),
            max="**",
            fr=tr("Liste des groupes de mailles définissant la frontière discrète"),
        ),
        #
    ),
    #
    # 9.2. Definition analytique d'une frontière
    #
    FRONTIERE_ANALYTIQUE=FACT(
        statut="f",
        max="**",
        fr=tr("Definition analytique de frontières à suivre."),
        #
        # 9.2.1. Nom de la frontière
        #
        NOM=SIMP(statut="o", typ="TXM", max=1, fr=tr("Nom de la frontière analytique")),
        #
        # 9.2.2. Type de la frontière
        #
        TYPE=SIMP(
            statut="o",
            typ="TXM",
            into=("SPHERE", "CYLINDRE", "CONE_A", "CONE_R", "TORE"),
            fr=tr("Type de la frontière analytique"),
        ),
        #
        # 9.2.3. Pour une sphere, un cylindre, un cone defini par ses rayons ou un tore : rayon et centre
        #        Pour un tore, c'est le rayon de revolution
        #
        b_fr_rayon=BLOC(
            condition=""" (equal_to("TYPE", 'SPHERE')) or (equal_to("TYPE", 'CYLINDRE')) or (equal_to("TYPE", 'CONE_R')) or (equal_to("TYPE", 'TORE')) """,
            fr=tr("Le rayon et le centre d'une sphère, d'un cylindre, d'un cone ou d'un tore."),
            RAYON=SIMP(statut="o", typ="R", val_min=0.0, fr=tr("Rayon")),
            X_CENTRE=SIMP(statut="o", typ="R", fr=tr("Abscisse du centre")),
            Y_CENTRE=SIMP(statut="o", typ="R", fr=tr("Ordonnée du centre")),
            Z_CENTRE=SIMP(statut="o", typ="R", fr=tr("Cote du centre")),
        ),
        #
        # 9.2.4. Pour un cylindre, un cone defini par axe et angle ou un tore : axe
        #
        b_fr_cylindre=BLOC(
            condition=""" (equal_to("TYPE", 'CYLINDRE')) or (equal_to("TYPE", 'CONE_A')) or (equal_to("TYPE", 'TORE')) """,
            fr=tr("Pour un cylindre, un cone ou un tore."),
            X_AXE=SIMP(statut="o", typ="R", fr=tr("Abscisse du vecteur directeur de l'axe")),
            Y_AXE=SIMP(statut="o", typ="R", fr=tr("Ordonnée du vecteur directeur de l'axe")),
            Z_AXE=SIMP(statut="o", typ="R", fr=tr("Cote du vecteur directeur de l'axe")),
        ),
        #
        # 9.2.5. Pour un cone defini par ses rayons : second couple (rayon et centre)
        #
        b_fr_rayon2=BLOC(
            condition=""" (equal_to("TYPE", 'CONE_R')) """,
            fr=tr("Le second rayon et centre d'un cone."),
            RAYON2=SIMP(statut="o", typ="R", val_min=0.0, fr=tr("Rayon numéro 2")),
            X_CENTRE2=SIMP(statut="o", typ="R", fr=tr("Abscisse du centre numéro 2")),
            Y_CENTRE2=SIMP(statut="o", typ="R", fr=tr("Ordonnée du centre numéro 2")),
            Z_CENTRE2=SIMP(statut="o", typ="R", fr=tr("Cote du centre numéro 2")),
        ),
        #
        # 9.2.6. Pour un cone defini par axe et angle : centre et angle
        #
        b_fr_angle=BLOC(
            condition=""" (equal_to("TYPE", 'CONE_A')) """,
            fr=tr("L'angle et le centre d'un cone."),
            ANGLE=SIMP(statut="o", typ="R", val_min=0.0, fr=tr("Angle en degré")),
            X_CENTRE=SIMP(statut="o", typ="R", fr=tr("Abscisse du centre")),
            Y_CENTRE=SIMP(statut="o", typ="R", fr=tr("Ordonnée du centre")),
            Z_CENTRE=SIMP(statut="o", typ="R", fr=tr("Cote du centre")),
        ),
        #
        # 9.2.7. Pour un tore : rayon primaire
        #
        b_fr_rayon3=BLOC(
            condition=""" (equal_to("TYPE", 'TORE')) """,
            fr=tr("Le rayon primaire du tore."),
            RAYON2=SIMP(statut="o", typ="R", val_min=0.0, fr=tr("Rayon primaire")),
        ),
        # 9.2.8. Groupe(s) lie(s) a la frontière
        #
        GROUP_MA=SIMP(
            statut="o",
            typ=grma,
            validators=NoRepeat(),
            max="**",
            fr=tr("Liste des groupes de mailles placées sur la frontière"),
        ),
        #
    ),
    #
    # 10. mise à jour de champs sur le nouveau maillage
    #
    MAJ_CHAM=FACT(
        statut="f",
        max="**",
        fr=tr("Mise à jour de champs sur le nouveau maillage."),
        #
        # 10.1. Le nom du champ de grandeur qui contiendra le resultat de la mise a jour
        #
        CHAM_MAJ=SIMP(
            statut="o", typ=CO, fr=tr("Nom du champ de grandeur qui contiendra le champ mis à jour")
        ),
        #
        # 10.2. Le type du champ qui contiendra le resultat de la mise a jour
        #
        TYPE_CHAM=SIMP(
            statut="o",
            typ="TXM",
            into=C_TYPE_CHAM_INTO(("NOEU", "ELEM", "ELNO", "ELGA")),
            fr=tr("Type du champ qui contiendra le champ mis à jour"),
        ),
        #
        # 10.3. Le champ a interpoler
        #
        regles=(UN_PARMI("CHAM_GD", "RESULTAT")),
        #
        # 10.3.1. Sous forme de champ de grandeur
        #
        CHAM_GD=SIMP(
            statut="f",
            typ=cham_gd_sdaster,
            fr=tr("Champ de grandeur Code_Aster contenant le champ à mettre à jour"),
        ),
        #
        # 10.3.2. Sous forme de champ dans un resultat
        #
        RESULTAT=SIMP(
            statut="f",
            typ=(evol_elas, evol_noli, evol_ther),
            fr=tr("Résultat contenant le champ à mettre à jour"),
        ),
        #
        b_nom_du_champ=BLOC(
            condition="""(exists("RESULTAT"))""",
            fr=tr("Choix éventuel du nom du champ à interpoler"),
            #
            NOM_CHAM=SIMP(
                statut="o",
                typ="TXM",
                validators=NoRepeat(),
                into=C_NOM_CHAM_INTO(),
                fr=tr("Nom du champ à mettre à jour"),
            ),
            #
        ),
        #
        # 10.4. Les composantes
        #
        NOM_CMP=SIMP(
            statut="f",
            typ="TXM",
            validators=NoRepeat(),
            max="**",
            fr=tr("Liste des composante(s) retenue(s) pour le champ."),
        ),
        #
        # 10.5. Le paramètre temporel pour le champ a interpoler
        #
        b_parametre_temporel=BLOC(
            condition="""(exists("RESULTAT"))""",
            fr=tr("Choix éventuel du paramètre temporel pour le champ à interpoler"),
            #
            regles=(EXCLUS("NUME_ORDRE", "INST"),),
            #
            # 10.5.1. Soit le numero d'ordre
            #
            NUME_ORDRE=SIMP(
                statut="f", typ="I", val_min=0, fr=tr("Numéro d ordre du champ à mettre à jour")
            ),
            #
            # 10.5.2. Soit l'instant
            # 10.5.2.1. Sa valeur
            #
            INST=SIMP(statut="f", typ="R", fr=tr("Instant associé")),
            #
            # 10.5.2.2. La précision du choix de l'instant
            #
            b_precision=BLOC(
                condition="""(exists("INST"))""",
                fr=tr("Choix de la précision du choix de l'instant"),
                #
                CRITERE=SIMP(
                    statut="f",
                    typ="TXM",
                    defaut="RELATIF",
                    into=("RELATIF", "ABSOLU"),
                    fr=tr("Critère de précision sur le choix de l'instant associé"),
                ),
                b_prec_rela=BLOC(
                    condition="""(equal_to("CRITERE", 'RELATIF'))""",
                    PRECISION=SIMP(
                        statut="f",
                        typ="R",
                        defaut=1.0e-6,
                        fr=tr("Précision relative sur le choix de l'instant associé"),
                    ),
                ),
                b_prec_abso=BLOC(
                    condition="""(equal_to("CRITERE", 'ABSOLU'))""",
                    PRECISION=SIMP(
                        statut="o",
                        typ="R",
                        fr=tr("Précision absolue sur le choix de l'instant associé"),
                    ),
                ),
                #
            ),
            #
        ),
        #
        # 10.6. Type d'interpolation
        #
        TYPE_MAJ=SIMP(
            statut="f",
            typ="TXM",
            defaut="AUTO",
            into=("AUTO", "ISOP2"),
            fr=tr("Type de mise à jour : automatique ou iso-P2"),
        ),
    ),
    #
    # 11. Les modifications
    #
    b_modifications=BLOC(
        condition=""" (equal_to("ADAPTATION", 'MODIFICATION')) """,
        fr=tr("Modification de maillage."),
        #
        # regles=(AU_MOINS_UN('DEGRE','JOINT'),),
        #
        # 11.1. Changement de degre
        #
        DEGRE=SIMP(statut="o", typ="TXM", into=("OUI",), fr=tr("Changement de degré du maillage")),
        #
        # 11.2. Création des joints
        #
        # JOINT         = SIMP(statut='f',typ='TXM',defaut="NON",into=("OUI", "NON"),
        # fr=tr("Creations des joints"),
        # ),
        #
    ),
    #
    # 12. Le modele pour les lectures de champs aux points de Gauss ou aux noeuds par element
    #
    b_lectures=BLOC(
        condition=""" (equal_to("ADAPTATION", 'LECTURE')) """,
        fr=tr("Lectures de champs aux points de Gauss."),
        #
        MODELE=SIMP(statut="o", typ=modele_sdaster, fr=tr("Modèle")),
        #
    ),
    #
    # 13. Les options d'analyse de maillage ; par defaut, on ne fait que les nombres
    #
    b_analyses=BLOC(
        condition=""" (not equal_to("ADAPTATION", 'LECTURE')) """,
        fr=tr("Analyse du maillage."),
        #
        # 13.1. Nombre de noeuds et mailles
        #
        NOMBRE=SIMP(
            statut="f",
            typ="TXM",
            defaut="OUI",
            into=("OUI", "NON"),
            fr=tr("Nombre de noeuds et de mailles du maillage"),
        ),
        #
        # 13.2. Determination de la qualité des mailles du maillage
        #
        QUALITE=SIMP(
            statut="f",
            typ="TXM",
            defaut="NON",
            into=("OUI", "NON"),
            fr=tr("Qualité des mailles du maillage."),
        ),
        #
        # 13.3. Connexite du maillage
        #
        CONNEXITE=SIMP(
            statut="f",
            typ="TXM",
            defaut="NON",
            into=("OUI", "NON"),
            fr=tr("Connexité du maillage."),
        ),
        #
        # 13.4. Taille des sous-domaines du maillage
        #
        TAILLE=SIMP(
            statut="f",
            typ="TXM",
            defaut="NON",
            into=("OUI", "NON"),
            fr=tr("Tailles des sous-domaines du maillage."),
        ),
        #
        # 13.5. Controle de la non-interpenetration des mailles
        #
        INTERPENETRATION=SIMP(
            statut="f",
            typ="TXM",
            defaut="NON",
            into=("OUI", "NON"),
            fr=tr("Controle de la non interpénétration des mailles."),
        ),
        #
        # 13.6. Propriétés du maillage de calcul
        #
        PROP_CALCUL=SIMP(
            statut="f",
            typ="TXM",
            defaut="NON",
            into=("OUI", "NON"),
            fr=tr("Propriétés du maillage de calcul."),
        ),
        #
        # 13.7. Determination des diametres des mailles du maillage
        #
        DIAMETRE=SIMP(
            statut="f",
            typ="TXM",
            defaut="NON",
            into=("OUI", "NON"),
            fr=tr("Diamètre des mailles du maillage."),
        ),
        #
    ),
    #
    # 14. Champs supplémentaires sur le nouveau maillage
    #
    ADD_CHAM=FACT(
        statut="f",
        max="**",
        fr=tr("Champs supplémentaires sur le nouveau maillage."),
        #
        # 14.1. Le nom du champ de grandeur qui contiendra le nouveau champ
        #
        CHAM_GD=SIMP(
            statut="o",
            typ=CO,
            fr=tr("Nom du champ de grandeur qui contiendra le champ supplémentaire"),
        ),
        #
        # 14.2. La categorie du champ supplementaire
        #
        CHAM_CAT=SIMP(
            statut="o",
            typ="TXM",
            into=("NIVEAU", "QUALITE", "DIAMETRE"),
            fr=tr("Catégorie du champ supplémentaire : niveau, qualité ou diamètre"),
        ),
        #
    ),
    #
    # 16. Les options avancées
    # 16.1. Langue des messages issus de HOMARD
    #
    LANGUE=SIMP(
        statut="f",
        typ="TXM",
        defaut="FRANCAIS",
        into=("FRANCAIS", "FRENCH", "ANGLAIS", "ENGLISH"),
        fr=tr("Langue des messages issus de HOMARD."),
    ),
    #
    # 16.2. Gestion des mailles acceptees dans le maillage initial
    #       "HOMARD" : exclusivement les mailles pouvant etre decoupees (defaut)
    #       "IGNORE_PYRA" : elles sont ignorées
    #
    b_autres_mailles=BLOC(
        condition=""" (not equal_to("ADAPTATION", 'LECTURE')) """,
        fr=tr("Gestion des pyramides."),
        #
        ELEMENTS_ACCEPTES=SIMP(
            statut="f",
            typ="TXM",
            defaut="HOMARD",
            into=("HOMARD", "IGNORE_PYRA"),
            fr=tr("Acceptation des mailles dans le maillage initial"),
        ),
        #
    ),
    #
    # 16.3. Version de HOMARD
    #
    VERSION_HOMARD=SIMP(
        statut="f",
        typ="TXM",
        defaut="V11_10",
        into=("V11_10", "V11_N", "V11_N_PERSO"),
        fr=tr("Version de HOMARD"),
    ),
    #
    # 17. Les options de débogage
    # 17.1. Exécutable pilotant HOMARD
    #
    LOGICIEL=SIMP(statut="f", typ="TXM", fr=tr("Logiciel pilotant HOMARD")),
    #
    # 17.2. Unite logique d'un fichier à ajouter à HOMARD.Configuration
    #
    b_unite=BLOC(
        condition="""is_in("VERSION_HOMARD", ('V11_N','V11_N_PERSO') )""",
        fr=tr("Fichier supplémentaire."),
        #
        UNITE=SIMP(
            statut="f",
            typ=UnitType(),
            inout="in",
            fr=tr("Unite logique à ajouter à HOMARD.Configuration"),
        ),
        #
    ),
    #
)
