# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

# person_in_charge: albert.alarcon at edf.fr

from ..Cata.SyntaxUtils import remove_none
from ..Commands import (
    ASSE_MATRICE,
    ASSE_VECTEUR,
    CALC_MATR_ELEM,
    CALC_VECT_ELEM,
    NUME_DDL,
)
from ..Messages import UTMESS


def create_nume(self, numeddl_status, option, numeddl, matr_rigi, CHARGE, INFO, MODELE):
    if numeddl_status == "To_Read":
        num = numeddl
    elif numeddl_status == "To_Create":
        if option in ("RIGI_MECA", "RIGI_THER", "RIGI_ACOU", "RIGI_FLUI_STRU"):
            num = NUME_DDL(MATR_RIGI=matr_rigi, INFO=INFO)
        else:
            if CHARGE != None:
                num = NUME_DDL(MODELE=MODELE, CHARGE=CHARGE, INFO=INFO)
            else:
                num = NUME_DDL(MODELE=MODELE, INFO=INFO)
    return num


def assemblage_ops(self, MODELE, NUME_DDL, INFO, **args):
    """
    Ecriture de la macro MACRO_ASSE
    """
    CHAM_MATER = args.get("CHAM_MATER")
    CARA_ELEM = args.get("CARA_ELEM")
    CHARGE = args.get("CHARGE")
    CHAR_CINE = args.get("CHAR_CINE")
    INST = args.get("INST")
    MATR_ASSE = args.get("MATR_ASSE")
    VECT_ASSE = args.get("VECT_ASSE")

    num = None

    # On met le mot-cl?? NUME_DDL dans une variable locale pour le prot??ger
    numeddl = NUME_DDL
    info = INFO
    # Le nom de la variable doit ??tre obligatoirement le nom de la commande

    if numeddl in self.sdprods:
        # Si le concept numeddl est dans self.sdprods,
        # il doit ??tre produit par la macro,
        # il faudra donc appeler la commande NUME_DDL.
        if MATR_ASSE is None:
            UTMESS("F", "MATRICE0_5")

        numeddl_status = "To_Create"
    else:
        numeddl_status = "To_Read"
        # Dans le cas o?? on assemble des vecteurs, le mot-cl?? charge utilis??
        # pr??alablement pour la construction des conditions de Dirichlet est
        # indispensable... sauf si initialement les matrices ??taient construites
        # sans aucune charge (peut-??tre le cas en dynamique)
        if VECT_ASSE and not CHARGE:
            UTMESS("A", "MATRICE0_6")

    lrigel = 0
    lmasel = 0

    # Assemblage des matrices
    if MATR_ASSE:
        # D??calage ??ventuel en premi??re position dans la liste de l'occurence
        # de MATR_ASSE contenant l'option de rigidit??
        MATR_ASSE = list(MATR_ASSE)
        for m in MATR_ASSE[:]:
            option = m["OPTION"]
            if option in (
                "RIGI_MECA",
                "RIGI_MECA_HYST",
                "RIGI_THER",
                "RIGI_ACOU",
                "RIGI_FLUI_STRU",
            ):
                MATR_ASSE.remove(m)
                MATR_ASSE.insert(0, m)
                break

        iocc = 0
        for m in MATR_ASSE:
            iocc = iocc + 1
            option = m["OPTION"]

            motscles = {"OPTION": option}
            if option == "RIGI_MECA_HYST":
                if not lrigel:
                    UTMESS("F", "MATRICE0_10")
                motscles["RIGI_MECA"] = rigel
            if option == "AMOR_MECA":
                if not lrigel:
                    UTMESS("F", "MATRICE0_11")
                if CHAM_MATER:
                    motscles["RIGI_MECA"] = rigel
                    if lmasel:
                        motscles["MASS_MECA"] = masel
            if CHARGE:
                if option[0:9] != "RIGI_GEOM":
                    motscles["CHARGE"] = CHARGE
            if CHAM_MATER:
                motscles["CHAM_MATER"] = CHAM_MATER
            if CARA_ELEM:
                motscles["CARA_ELEM"] = CARA_ELEM
            if INST:
                motscles["INST"] = INST
            if "SIEF_ELGA" in m:
                motscles["SIEF_ELGA"] = m["SIEF_ELGA"]
            if "MODE_FOURIER" in m:
                motscles["MODE_FOURIER"] = m["MODE_FOURIER"]
            if "GROUP_MA" in m:
                motscles["GROUP_MA"] = m["GROUP_MA"]

            remove_none(motscles)
            _a = CALC_MATR_ELEM(MODELE=MODELE, **motscles)

            if option == "RIGI_MECA":
                rigel = _a
                lrigel = 1
            if option == "MASS_MECA":
                masel = _a
                lmasel = 1

            # Create NUME_DDL
            if numeddl_status != "OK":
                num = create_nume(
                    self, numeddl_status, option, numeddl, _a, CHARGE, info, MODELE
                )
                if numeddl_status == "To_Create":
                    self.register_result(num, numeddl)
                numeddl_status = "OK"

            motscles = {"OPTION": option}
            if CHAR_CINE:
                mm = ASSE_MATRICE(MATR_ELEM=_a, NUME_DDL=num, CHAR_CINE=CHAR_CINE)
            else:
                mm = ASSE_MATRICE(MATR_ELEM=_a, NUME_DDL=num)
            self.register_result(mm, m["MATRICE"])

    # ASSEMBLAGE DES VECTEURS
    if VECT_ASSE:
        for v in VECT_ASSE:
            option = v["OPTION"]
            motscles = {"OPTION": option}
            liste_toutes_charges = []
            if CHARGE:
                liste_toutes_charges = list(CHARGE)

            # On ajoute ?? la liste de charges de Dirichlet les charges de Newmann
            for char_vv in v["CHARGE"] or []:
                if char_vv not in liste_toutes_charges:
                    liste_toutes_charges.append(char_vv)
            # On ne peut pas construire de vecteur si on a aucune charge
            if not liste_toutes_charges:
                UTMESS("F", "MATRICE0_8")

            motscles["CHARGE"] = tuple(liste_toutes_charges)
            # remplissage du mot-cle CHARGE pour toutes les options
            if option == "CHAR_MECA":
                if CARA_ELEM:
                    motscles["CARA_ELEM"] = CARA_ELEM
                if CHAM_MATER:
                    motscles["CHAM_MATER"] = CHAM_MATER
                if INST:
                    motscles["INST"] = INST
                if "MODE_FOURIER" in v:
                    motscles["MODE_FOURIER"] = v["MODE_FOURIER"]

            elif option == "CHAR_THER":
                if CARA_ELEM:
                    motscles["CARA_ELEM"] = CARA_ELEM
                if INST:
                    motscles["INST"] = INST

            elif option == 'CHAR_ACOU':
                if CHAM_MATER is None:
                    UTMESS("F", "MATRICE0_7")
                else:
                    motscles["CHAM_MATER"] = CHAM_MATER
            else:
                assert False, f"unexpected option {option}"

            _b = CALC_VECT_ELEM(**motscles)

            # Create NUME_DDL
            if numeddl_status != "OK":
                num = create_nume(
                    self, numeddl_status, option, numeddl, None, CHARGE, info, MODELE
                )
                if numeddl_status == "To_Create":
                    self.register_result(num, numeddl)
                numeddl_status = "OK"

            # Les vecteurs assembl??s sont des concepts sortants.
            # Assemblage des vecteurs
            vv = ASSE_VECTEUR(VECT_ELEM=_b, NUME_DDL=num)
            self.register_result(vv, v["VECTEUR"])
    return
