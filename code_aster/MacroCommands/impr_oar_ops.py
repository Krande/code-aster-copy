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

from ..Objects.table_py import Table
from ..Messages import UTMESS
import numpy as np
import time
from ..Utilities.version import get_version


class OAR_EF:
    def __init__(self, **args):

        """donnes d entrees : tables issues de la commande MACR_LIGN_COU"""

        self.charg = {
            "PRESSION": None,
            "TORSION": None,
            "FLEXION_HP": None,
            "FLEXION_P": None,
            "TEMP": None,
            "CONTRAINTE": None,
        }
        self.etiquette = {
            "PRESSION": "SIGM_P",
            "TORSION": "SIGM_T",
            "FLEXION_HP": "SIGM_HP",
            "FLEXION_P": "SIGM_FP",
            "TEMP": "TEMP",
            "CONTRAINTE": "SIGM",
        }

        self.nbAzim = 0
        self.typefic = "meca"  # 2 possibilités : sortie 'meca' ou 'temp'
        self.listeinstant_temp = None
        self.ficsortie = None
        self.ordre = None

        self.titre = None

        self.type_unit = None
        self.type_unit_etiquette = {
            "SI": "Modele en S.I., Impression en (m, MPa)",
            "MM-MPA": "Modele en (mm,MPa) , Impression en (m, MPa)",
        }

    def Tablenonvide(self):
        """renvoie une table non vide du dictionnaire des chargements"""
        for key in self.charg.keys():
            if self.charg[key] is not None:
                table = self.charg[key]
                break
        return table

    def ajoutINST(self):
        """dans le cas d'un MACRO_LIGN_COUPE après MACRO_ELAS_MULT
         il n'y a pas d instants  on ajoute une liste d instants
        factice pour eviter le bug"""
        for key in self.charg.keys():
            if self.charg[key] is not None:
                try:
                    a = self.charg[key].values()["INST"][0]
                except:
                    self.charg[key]["INST"] = np.zeros(len(self.charg[key].values()["ABSC_CURV"]))

    def ListeAzim(self):
        """renvoie la liste des azimuts à imprimer sous forme de liste
        d'intitulés"""
        table = self.Tablenonvide()
        tablereduit = table.INST == table.values()["INST"][0]
        intituleliste = tablereduit.values()["INTITULE"]
        liste = [intituleliste[0]]
        i = 0
        for item in intituleliste:
            if item != liste[i]:
                i = i + 1
                liste.append(item)
        return liste

    def AbscisseAzim(self, macroc, intitule):
        """Renvoie les abscisses en m d'un intitule/azimut"""
        if self.type_unit == "SI":
            conversion = 1.0
        if self.type_unit == "MM-MPA":
            conversion = 1e-3
        macroreduit = macroc.INST == macroc.values()["INST"][0]
        macrocr = macroreduit.INTITULE == intitule
        abscisse = conversion * np.array(macrocr.values()["ABSC_CURV"])
        repere = range(1, len(abscisse) + 1)
        join = [repere, abscisse]
        join = np.transpose(join)
        return join

    def CoupeAzim(self, macroc, intitule, instant):
        """renvoie le tableau des contraintes
        pour un azimut/intitule donné"""
        macroreduit = macroc.INST == instant
        macrocr = macroreduit.INTITULE == intitule
        if self.type_unit == "SI":
            conversion = 1e-6
        if self.type_unit == "MM-MPA":
            conversion = 1.0
        cles_sig = ["SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ"]
        tableau = []
        for cles in cles_sig:
            try:
                tableau.append(conversion * np.array(macrocr.values()[cles]))
            except:
                # cas axisymétriques : certaines composantes peuvent manquer.
                # on les remplace alors par des listes nulles.
                tableau.append(np.zeros(len(macrocr.values()["ABSC_CURV"])))

        repere = range(1, len(tableau[0]) + 1)
        tableau.insert(0, repere)
        tableau = np.transpose(tableau)
        return tableau

    def CoupeAzim_temp(self, macroc, intitule, instant):
        """renvoie le tableau des contraintes
        pour un azimut/intitule donné"""
        macroreduit = macroc.INST == instant
        macrocr = macroreduit.INTITULE == intitule
        TEMP = np.array(macrocr.values()["TEMP"])
        repere = range(1, len(TEMP) + 1)
        tableau = [repere, TEMP]
        tableau = np.transpose(tableau)
        return tableau

    def ecrire_preambule(self, nomfichier):
        """écriture des lignes d'introduction du fichier de sortie"""
        with open(nomfichier, "a") as f:
            f.write("; CODE_ASTER version " + get_version() + "\n")
            f.write("; IMPRESSION AU FORMAT OAR \n")
            f.write("; CREATION " + (time.strftime("LE %m/%d/%Y A %H:%M:%S")) + "\n \n")
            if self.titre is not None:
                f.write("; {} \n \n".format(self.titre))
            f.write("; UNITES : " + self.type_unit_etiquette[self.type_unit])
            f.write("\n \n")

    def ecrire_abs(self, nomfichier, tableau):
        """écriture des abscisses"""
        with open(nomfichier, "a") as f:
            f.write("ABSC \n ")
            for i in range(len(tableau)):
                f.write("  %.6e \n " % (tableau[i][1]))
            f.write("\n \n")

    def ecrire_contr(self, nomfichier, tableau, etiquette):
        """écriture des contraintes - cas mécanique"""
        with open(nomfichier, "a") as f:
            f.write(etiquette + "\n")
            f.write(
                ";       SIXX            SIYY"
                "            SIZZ            SIXY            SIXZ"
                "            SIYZ \n"
            )
            for row in tableau:
                f.write(
                    f"{row[1]:15.4E} {row[2]:15.4E}"
                    f" {row[3]:15.4E} {row[4]:15.4E} {row[5]:15.4E} "
                    f"{row[6]:15.4E} \n"
                )
            f.write("\n")

    def ecrire_temp(self, nomfichier, tabl_temp, tabl_tempsig):
        """écriture des températures et contraintes - cas thermique"""
        with open(nomfichier, "a") as f:
            f.write(
                "; TEMPERATURE             SIXX            SIYY"
                "            SIZZ            SIXY            SIXZ"
                "            SIYZ \n"
            )
            for row1, row2 in zip(tabl_temp, tabl_tempsig):
                f.write(
                    f"    {row1[1]:5.1f}         {row2[1]:15.4E} "
                    f"{row2[2]:15.4E} {row2[3]:15.4E} {row2[4]:15.4E} "
                    f"{row2[5]:15.4E} {row2[6]:15.4E} \n"
                )
            f.write("\n")

    def Impression(self):
        """Impression proprement dite du rapport de sortie"""
        self.ajoutINST()
        listeAzim = self.ListeAzim()
        nbAzim = len(listeAzim)
        tablenonvide = self.Tablenonvide()

        if self.typefic == "meca":
            self.ecrire_preambule(self.ficsortie)
            for intitule in listeAzim:
                with open(self.ficsortie, "a") as f:
                    f.write("; " + intitule + "\n \n")
                tableau_abs = self.AbscisseAzim(tablenonvide, intitule)
                self.ecrire_abs(self.ficsortie, tableau_abs)

                for key in self.ordre:
                    if self.charg[key] is not None:
                        instant = self.charg[key].values()["INST"][0]
                        tableau = self.CoupeAzim(self.charg[key], intitule, instant)
                        self.ecrire_contr(self.ficsortie, tableau, self.etiquette[key])

        if self.typefic == "temp":
            # determination liste des instants
            premiereabsc = self.charg["TEMP"].values()["ABSC_CURV"][0]
            tableintro = self.charg["TEMP"].ABSC_CURV == premiereabsc
            tableintror = tableintro.INTITULE == tableintro.values()["INTITULE"][0]
            NUME_ORDRE = np.array(tableintror.values()["NUME_ORDRE"])
            INST = np.array(tableintror.values()["INST"])
            self.listeinstant_temp = INST

            # impression du preambule et de la liste des instants
            self.ecrire_preambule(self.ficsortie)
            with open(self.ficsortie, "a") as f:
                f.write("; LISTE DES INSTANTS \n \n")
                f.write("INST \n \n")
                for instant in self.listeinstant_temp:
                    f.write("{} \n".format(instant))
                f.write("\n \n")
                f.write("; LISTE DES TEMPERATURES ET CONTRAINTES PAR COUPE")
                f.write("\n \n")

            # impression
            for intitule in listeAzim:
                tableau_abs = self.AbscisseAzim(self.charg["TEMP"], intitule)
                with open(self.ficsortie, "a") as f:
                    f.write("; " + intitule + "\n \n")
                self.ecrire_abs(self.ficsortie, tableau_abs)
                with open(self.ficsortie, "a") as f:
                    f.write("TEMP SIGM \n \n")
                for instant in self.listeinstant_temp:
                    with open(self.ficsortie, "a") as f:
                        f.write("; INSTANT {} \n \n".format(instant))
                    tableau_temp = self.CoupeAzim_temp(self.charg["TEMP"], intitule, instant)
                    tableau_tempsig = self.CoupeAzim(self.charg["CONTRAINTE"], intitule, instant)
                    self.ecrire_temp(self.ficsortie, tableau_temp, tableau_tempsig)
                    with open(self.ficsortie, "a") as f:
                        f.write("\n \n")


def impr_oar_ops(self, **args):
    """
    Macro IMPR_OAR

    """
    resultat = OAR_EF()

    TABL_MECA = args.get("TABL_MECA")
    TABL_THER = args.get("TABL_THER")
    resultat.titre = args.get("TITRE")
    resultat.type_unit = "SI" if args.get("TYPE_UNIT") is None else args.get("TYPE_UNIT")

    try:
        unite = args["UNITE"]
    except:
        UTMESS("F", "OAR0_1")

    name = "fort." + str(unite)
    resultat.ficsortie = name
    open(resultat.ficsortie, "w").close()  # on vide le fichier

    if (TABL_MECA is not None) and (TABL_THER is not None):
        UTMESS("F", "OAR0_2")

    if resultat.type_unit not in ["SI", "MM-MPA"]:
        UTMESS("F", "OAR0_3")

    if TABL_MECA is not None:
        resultat.ordre = TABL_MECA[0].keys()
        for key in resultat.ordre:
            resultat.charg[key] = TABL_MECA[0][key].EXTR_TABLE()
        resultat.typefic = "meca"
        resultat.Impression()

    if TABL_THER is not None:
        resultat.charg["TEMP"] = TABL_THER[0]["TEMP"].EXTR_TABLE()
        resultat.charg["CONTRAINTE"] = TABL_THER[0]["CONTRAINTE"].EXTR_TABLE()
        resultat.typefic = "temp"
        resultat.Impression()
