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

# person_in_charge: jean-luc.flejou at edf.fr

import math
import numpy as NP

from ..CodeCommands import DEFI_MATERIAU
from ..Messages import UTMESS


def FaitMessage(Dico):
    message = ""
    cpt = 1
    for xk, xv in Dico.items():
        message += "%s = %15.8E; " % (xk, xv)
        if len(message) > 80 * cpt:
            message += "\n"
            cpt += 1
    return message


def BetonEC2(Classe):
    Dico = {}
    # Décodage de la classe
    sfck, sfckc = Classe[1:].split("/")
    #
    Dico["fck"] = float(sfck)
    Dico["fckc"] = float(sfckc)
    # Contrainte en MPa et déformation en "pour mille"
    Dico["fcm"] = Dico["fck"] + 8.0
    Dico["nu"] = 0.20
    #
    if Dico["fck"] <= 50.0:
        Dico["fctm"] = 0.30 * NP.power(Dico["fck"], 2.0 / 3.0)
    else:
        Dico["fctm"] = 2.12 * NP.log(1.0 + Dico["fcm"] / 10.0)
    #
    Dico["ecm"] = 22.0e03 * NP.power(Dico["fcm"] / 10.0, 0.3)
    #
    Dico["epsi_c1"] = 0.7 * NP.power(Dico["fcm"], 0.31)
    if Dico["epsi_c1"] >= 2.8:
        Dico["epsi_c1"] = 2.8
    #
    if Dico["fck"] > 50.0:
        Dico["epsi_cu1"] = 2.80 + 27.000 * NP.power((98.0 - Dico["fcm"]) / 100.0, 4.0)
        Dico["epsi_c2"] = 2.00 + 0.085 * NP.power(Dico["fck"] - 50.0, 0.53)
        Dico["epsi_cu2"] = 2.60 + 35.000 * NP.power((90.0 - Dico["fck"]) / 100.0, 4.0)
        Dico["n"] = 1.40 + 23.400 * NP.power((90.0 - Dico["fck"]) / 100.0, 4.0)
        Dico["epsi_c3"] = 1.75 + 0.550 * (Dico["fck"] - 50.0) / 40.0
        Dico["epsi_cu3"] = 2.60 + 35.000 * NP.power((90.0 - Dico["fck"]) / 100.0, 4.0)
    else:
        Dico["epsi_cu1"] = 3.50
        Dico["epsi_c2"] = 2.00
        Dico["epsi_cu2"] = 3.50
        Dico["n"] = 2.00
        Dico["epsi_c3"] = 1.75
        Dico["epsi_cu3"] = 3.50
    #
    return Dico


def BetonBAEL91(fcj):
    Dico = {}
    # Contrainte en MPa
    Dico["fcj"] = fcj
    #
    Dico["eij"] = 11000.0 * NP.power(Dico["fcj"], 0.333333)
    Dico["ftj"] = 0.60 + 0.06 * Dico["fcj"]
    Dico["epsi_c"] = 0.620e-3 * NP.power(Dico["fcj"], 0.333333)
    Dico["nu"] = 0.20
    #
    return Dico


def Mazars_Unil(DMATER, args):
    """
    MAZARS_UNIL = Paramètres de la loi de comportement
        UNITE_LONGUEUR = unité du problème [M|MM]
        FCJ    [Unite] = Contrainte au pic en compression
        EIJ    [Unite] = Module d'young
        EPSI_C         = Déformation au pic en compression
        FTJ    [Unite] = Contrainte au pic en traction
        NU             = Coefficient de poisson
        EPSD0          = Déformation, seuil d'endommagement
        K              = Paramètre de décroissance post-pic en cisaillement
        AC             = Paramètre de décroissance post-pic en compression
        BC             = 1/(Déformation au pic en compression)
        AT             = Paramètre de décroissance post-pic en traction
        BT             = 1/(Déformation au pic en traction)
        SIGM_LIM       = Contrainte limite pour post-traitement
        EPSI_LIM       = Déformation limite pour post-traitement

    Masse volumique, dilatation, amortissements
        RHO            = Masse volumique
        ALPHA          = Coefficient de dilatation
        AMOR_ALPHA     =
        AMOR_BETA      =
        AMOR_HYST      =
    """
    #
    MATER = DMATER.cree_dict_valeurs(DMATER.mc_liste)
    #
    # Obligatoire : Règlement de codification
    Regle = MATER["CODIFICATION"]
    # Liste des paramètres matériaux facultatifs mais nécessaires pour calculer
    # les valeurs des paramètres de MAZARS
    listepara = [
        "NU",
        "EPSD0",
        "EPST0",
        "EPSC0",
        "BT",
        "K",
        "AT",
        "BC",
        "AC",
        "SIGM_LIM",
        "EPSI_LIM",
    ]
    #
    # Les paramètres
    #  Obligatoire   : FCJ , EIJ, FTJ, EPSI_C
    #  Donnée ou pas avec des valeurs par défaut : NU , EPSI_LIM
    NU = None
    EPSI_LIM = None
    if Regle == "BAEL91":
        # Obligatoire : FCJ UNITE_CONTRAINTE
        if MATER["UNITE_CONTRAINTE"] == "MPa":
            coeff = 1.0
        elif MATER["UNITE_CONTRAINTE"] == "Pa":
            coeff = 1.0e06
        beton = BetonBAEL91(MATER["FCJ"] / coeff)
        #
        FCJ = beton["fcj"] * coeff
        EIJ = beton["eij"] * coeff
        FTJ = beton["ftj"] * coeff
        EPSI_C = beton["epsi_c"]
        NU = beton["nu"]
        #
        for xx in listepara:
            MATER[xx] = None
        #
    elif Regle == "EC2":
        # Obligatoire CLASSE UNITE_CONTRAINTE
        if MATER["UNITE_CONTRAINTE"] == "MPa":
            coeff = 1.0
        elif MATER["UNITE_CONTRAINTE"] == "Pa":
            coeff = 1.0e06
        beton = BetonEC2(MATER["CLASSE"])
        #
        FCJ = beton["fcm"] * coeff
        EIJ = beton["ecm"] * coeff
        FTJ = beton["fctm"] * coeff
        EPSI_C = beton["epsi_c1"] / 1000.0
        NU = beton["nu"]
        EPSI_LIM = beton["epsi_cu1"] / 1000.0
        #
        for xx in listepara:
            MATER[xx] = None
        #
    elif Regle == "ESSAI":
        # Obligatoire FCJ , EIJ, FTJ, EPSI_C
        FCJ = MATER["FCJ"]
        EIJ = MATER["EIJ"]
        FTJ = MATER["FTJ"]
        EPSI_C = MATER["EPSI_C"]
        MATER["UNITE_CONTRAINTE"] = ""

    # L'ordre dans la liste est important à cause des dépendances des relations
    # Les coefficients FCJ , EIJ, FTJ, EPSI_C doivent déjà être définis
    # Optional keywords
    def mater_value(name, default, value=None):
        if MATER.get(name) is not None:
            return MATER[name]
        else:
            if value is not None:
                return value
        return default

    #
    NU = mater_value("NU", 0.2, NU)
    NUB = NU * (2.0**0.5)
    EPSD0 = mater_value("EPSD0", FTJ / EIJ)
    EPST0 = mater_value("EPST0", FTJ / EIJ)
    EPSC0 = mater_value("EPSC0", EPST0 / NUB)
    K = mater_value("K", 0.70)
    BT = mater_value("BT", 1.00 / EPST0)
    AT = mater_value("AT", 0.90)
    BC = mater_value("BC", 1.00 / EPSI_C)
    AC = mater_value("AC", (FCJ / EIJ - EPSC0) / (EPSI_C * NP.exp(BC * (EPSC0 - EPSI_C)) - EPSC0))
    SIGM_LIM = mater_value("SIGM_LIM", 0.60 * FCJ)
    EPSI_LIM = mater_value("EPSI_LIM", 3.5 / 1000.0, EPSI_LIM)
    #
    # Mot clef MATER
    mclef = elastic_properties(EIJ, NU, args)
    mclef["MAZARS"] = {
        "K": K,
        "EPSC0": EPSC0,
        "EPST0": EPST0,
        "AC": AC,
        "AT": AT,
        "BC": BC,
        "BT": BT,
        "SIGM_LIM": SIGM_LIM,
        "EPSI_LIM": EPSI_LIM,
    }
    #
    # On affiche dans tous les cas
    if len(MATER["UNITE_CONTRAINTE"]) > 0:
        message0 = "MAZARS [%s]" % MATER["UNITE_CONTRAINTE"]
    else:
        message0 = "MAZARS"
    #
    message1 = FaitMessage(mclef["ELAS"])
    message2 = FaitMessage(mclef["MAZARS"])
    Dico = {"FCJ": FCJ, "FTJ": FTJ, "EPSI_C": EPSI_C, "EPSD0": EPSD0}
    message3 = FaitMessage(Dico)
    #
    UTMESS("I", "COMPOR1_75", valk=(message0, message1, message2, message3))
    #
    return mclef


def Beton_GLRC(DMATER, args):
    """
    BETON_GLRC = Paramètres de la loi de comportement
        UNITE_LONGUEUR = unité du problème [M|MM]
        FCJ    [Unite] = Contrainte au pic en compression
        EIJ    [Unite] = Module d'young
        EPSI_C         = Déformation au pic en compression
        FTJ    [Unite] = Contrainte au pic en traction
        NU             = Coefficient de poisson

    """
    #
    MATER = DMATER.cree_dict_valeurs(DMATER.mc_liste)
    #
    # Obligatoire : Règlement de codification
    Regle = MATER["CODIFICATION"]

    if Regle == "EC2":
        # Obligatoire CLASSE UNITE_CONTRAINTE
        if MATER["UNITE_CONTRAINTE"] == "MPa":
            coeff = 1.0
        elif MATER["UNITE_CONTRAINTE"] == "Pa":
            coeff = 1.0e06
        beton = BetonEC2(MATER["CLASSE"])
        #
        FCJ = beton["fcm"] * coeff
        EIJ = beton["ecm"] * coeff
        FTJ = beton["fctm"] * coeff
        EPSI_C = beton["epsi_c1"] / 1000.0
        NU = beton["nu"]

    elif Regle == "ESSAI":
        # Obligatoire FCJ , EIJ, FTJ, EPSI_C
        FCJ = MATER["FCJ"]
        EIJ = MATER["EIJ"]
        FTJ = MATER["FTJ"]
        EPSI_C = MATER["EPSI_C"]
        NU = MATER["NU"]  # 0.2 par defaut dans le catalogue
        MATER["UNITE_CONTRAINTE"] = ""

    # Mot clef MATER
    # voir si ELAS est necessaire
    mclef = elastic_properties(EIJ, NU, args)
    mclef["BETON_GLRC"] = {"FCJ": FCJ, "FTJ": FTJ, "EPSI_C": EPSI_C}
    #
    # On affiche dans tous les cas
    if len(MATER["UNITE_CONTRAINTE"]) > 0:
        message0 = "BETON_GLRC [%s]" % MATER["UNITE_CONTRAINTE"]
    else:
        message0 = "BETON_GLRC"
    #
    message1 = FaitMessage(mclef["ELAS"])
    message2 = FaitMessage(mclef["BETON_GLRC"])
    # Dico = {'FCJ': FCJ, 'FTJ': FTJ, 'EPSI_C': EPSI_C}
    # message3 = FaitMessage(Dico)
    #
    UTMESS("I", "COMPOR1_75", valk=(message0, message1, message2))
    #
    return mclef


def Acier_Cine_Line(DMATER, args):
    """
    ACIER = Paramètres matériaux de l'acier
        E              = Module d'Young
        D_SIGM_EPSI    = Module plastique
        SY             = Limite élastique
        SIGM_LIM       = Contrainte limite pour post-traitement
        EPSI_LIM       = Déformation limite pour post-traitement

    Masse volumique, dilatation, amortissements
        RHO            = Masse volumique
        ALPHA          = Coefficient de dilatation
        AMOR_ALPHA     =
        AMOR_BETA      =
        AMOR_HYST      =
    """
    #
    MATER = DMATER.cree_dict_valeurs(DMATER.mc_liste)
    # Obligatoire E
    E = MATER["E"]
    # Obligatoire SY
    SY = MATER["SY"]

    # Optional keywords
    def mater_value(name, default):
        return MATER[name] if MATER.get(name) is not None else default

    NU = mater_value("NU", 0.3)
    D_SIGM_EPSI = mater_value("D_SIGM_EPSI", E / 1.0e04)
    SIGM_LIM = mater_value("SIGM_LIM", SY / 1.1)
    EPSI_LIM = mater_value("EPSI_LIM", 10.0 / 1000.0)

    # Mot clef MATER
    mclef = elastic_properties(E, NU, args)
    mclef["ECRO_LINE"] = {
        "D_SIGM_EPSI": D_SIGM_EPSI,
        "SY": SY,
        "SIGM_LIM": SIGM_LIM,
        "EPSI_LIM": EPSI_LIM,
    }
    # On affiche dans tous les cas
    message1 = FaitMessage(mclef["ELAS"])
    message2 = FaitMessage(mclef["ECRO_LINE"])
    Dico = {"EPSI_ELAS": SY / E}
    message3 = FaitMessage(Dico)
    #
    UTMESS("I", "COMPOR1_75", valk=("ECRO_LINE", message1, message2, message3))
    #
    return mclef


def Ident_Endo_Fiss_Exp(ft, fc, beta, prec=1e-10, itemax=100):
    # Estimation initiale
    A = (2.0 / 3.0 + 3 * beta**2) ** 0.5
    r = fc / ft
    C = 3**0.5
    L = A * (r - 1)
    p0 = 2 * (1 - C)
    pp = 1 - L
    delta = pp**2 - p0
    x = -pp + delta**0.5

    # Resolution de l'equation par methode de Newton
    for i in range(itemax):
        f = L * x + (2 + NP.exp(-2 * r * x)) ** 0.5 - (2 + NP.exp(2 * x)) ** 0.5
        if abs(f) < prec:
            break
        df = (
            L
            - r * NP.exp(-2 * r * x) / (2 + NP.exp(-2 * r * x)) ** 0.5
            - NP.exp(2 * x) / (2 + NP.exp(2 * x)) ** 0.5
        )
        x = x - f / df
    else:
        UTMESS("F", "COMPOR1_87")
    #
    tau = A * x + (2 + NP.exp(2 * x)) ** 0.5
    sig0 = ft / x
    return (sig0, tau)


def ConfinedTension(nu, sig0, tau, beta, prec=1e-10, itemax=100):
    # Initialisation
    s = NP.array((1 - nu, nu, nu))
    L = (2.0 / 3.0 * (1 - 2 * nu) ** 2 + 3 * beta**2 * (1 + nu) ** 2) ** 0.5
    # Estimation initiale
    xe = NP.log(tau**2 - 2) / (2 * s[0])
    xl = tau / L
    x = min(xe, xl)
    # Résolution de l'équation par méthode de Newton
    for i in range(itemax):
        ep = NP.exp(x * s)
        epr = NP.dot(ep, ep) ** 0.5
        f = L * x + epr - tau
        if abs(f) < prec:
            break
        df = L + NP.add.reduce(ep * ep * s) / epr
        x = x - f / df
    else:
        UTMESS("F", "COMPOR1_87")
    #
    sig1 = x * sig0 * (1 - nu)
    return sig1


def Endo_Fiss_Exp(DMATER, args):
    """
    ENDO_FISS_EXP = Paramètres utilisateurs de la loi ENDO_FISS_EXP
      E              = Module de Young
      NU             = Coefficient de Poisson
      FT             = Limite en traction simple
      FC             = Limite en compression simple
      GF             = Energie de fissuration
      P              = Parametre dominant de la loi cohésive asymptotique
      G_INIT         = Energie de fissuration initiale (via la pente initiale de la loi cohésive)
      Q              = Parametre secondaire de la loi cohesive asymptotique
      Q_REL          = Parametre Q exprime de maniere relative par rapport a Qmax(P)
      LARG_BANDE     = Largeur de bande d'endommagement (2*D)
      REST_RIGI_FC   = Restauration de rigidité pour eps=fc/E (0=sans)
    """
    #
    MATER = DMATER.cree_dict_valeurs(DMATER.mc_liste)
    # Lecture et interprétation des paramètres utilisateurs
    E = float(MATER["E"])
    NU = float(MATER["NU"])
    GF = float(MATER["GF"])
    FT = float(MATER["FT"])
    FC = float(MATER["FC"])
    CRM = float(MATER["COEF_RIGI_MINI"])
    D = float(MATER["LARG_BANDE"] / 2.0)
    rrc = float(MATER["REST_RIGI_FC"])
    # Valeur par défaut
    beta = 0.1
    # Paramètres de la fonction seuil
    if FC / FT < 5.83:
        UTMESS("F", "COMPOR1_86", valr=(float(FC) / float(FT),))
    (sig0, tau) = Ident_Endo_Fiss_Exp(FT, FC, beta)
    sigc = ConfinedTension(NU, sig0, tau, beta)
    # Paramètres de la fonction d'adoucissement
    if MATER["P"] is not None:
        P = float(MATER["P"])
    else:
        G1 = float(MATER["G_INIT"])
        P = ((3 * NP.pi * GF) / (4 * G1)) ** (2.0 / 3.0) - 2
        if P < 1:
            UTMESS("F", "COMPOR1_93")
    #
    if MATER["Q"] is not None:
        Q = float(MATER["Q"])
    elif MATER["Q_REL"] is not None:
        qmax = (1.11375 + 0.565239 * P - 0.003322 * P**2) * (1 - NP.exp(-1.98935 * P)) - 0.01
        Q = qmax * float(MATER["Q_REL"])
    else:
        Q = 0.0
    #
    # Paramètres internes au modèle
    rig = E * (1 - NU) / ((1 + NU) * (1 - 2 * NU))
    K = 0.75 * GF / D
    C = 0.375 * GF * D
    M = 1.5 * rig * GF / (D * sigc**2)
    #
    if M < P + 2:
        UTMESS("F", "COMPOR1_94", valr=(float(M), float(P)))
    # Restauration de rigidité
    if rrc == 0.0:
        gamma = 0
    else:
        gamma = -1.0 / (FC / E * NP.log(rrc))

    # Paramètres pour DEFI_MATERIAU
    mclef = elastic_properties(E, NU, args)
    mclef.update(
        {
            "ENDO_FISS_EXP": {
                "M": M,
                "P": P,
                "Q": Q,
                "K": K,
                "TAU": tau,
                "SIG0": sig0,
                "BETA": beta,
                "COEF_RIGI_MINI": CRM,
                "REST_RIGIDITE": gamma,
            },
            "NON_LOCAL": {"C_GRAD_VARI": C, "PENA_LAGR": 1.0e3 * K},
        }
    )

    return mclef


def Endo_Loca_Exp(DMATER, args):
    """
    ENDO_LOCA_EXP = Paramètres utilisateurs de la loi ENDO_LOCA_EXP
      E              = Module de Young
      NU             = Coefficient de Poisson
      FT             = Limite en traction simple
      FC             = Limite en compression simple
      GF             = Energie de fissuration
      P              = Parametre dominant de la loi cohésive asymptotique
      DIST_FISSURE   = Distance moyenne inter-fissure
      REST_RIGI_FC   = Restauration de rigidité pour eps=fc/E (0=sans)
    """

    MATER = DMATER.cree_dict_valeurs(DMATER.mc_liste)

    # Lecture et interprétation des paramètres utilisateurs
    E = float(MATER["E"])
    NU = float(MATER["NU"])
    GF = float(MATER["GF"])
    FT = float(MATER["FT"])
    FC = float(MATER["FC"])
    P = float(MATER["P"])
    D = float(MATER["DIST_FISSURE"])
    rrc = float(MATER["REST_RIGI_FC"])

    # Valeur par défaut
    beta = 0.1

    # Paramètres de la fonction seuil
    if FC / FT < 5.83:
        UTMESS("F", "COMPOR1_86", valr=(float(FC) / float(FT),))
    (sig0, tau) = Ident_Endo_Fiss_Exp(FT, FC, beta)
    sigc = ConfinedTension(NU, sig0, tau, beta)

    # Paramètres internes au modèle
    Ec = E * (1 - NU) / ((1 + NU) * (1 - 2 * NU))
    wc = sigc**2 / (2 * Ec)
    kappa = GF / (D * wc)
    if rrc == 0.0:
        gamma = 0
    else:
        gamma = E / (FC * 2.0 * (1 - rrc))

    # Controle de la distance inter-fissure
    Dc = 0.75 * NP.pi * GF / wc * (P + 2.0) ** (-1.5)
    if D > Dc:
        UTMESS("F", "COMPOR1_97", valr=(Dc,))

    # Paramètres pour DEFI_MATERIAU
    mclef = elastic_properties(E, NU, args)
    mclef.update(
        {
            "ENDO_LOCA_EXP": {
                "SIGC": sigc,
                "SIG0": sig0,
                "KAPPA": kappa,
                "P": P,
                "BETA0": beta,
                "REST_RIGIDITE": gamma,
            }
        }
    )

    return mclef


def Endo_Loca_TC(DMATER, args):
    """
    ENDO_LOCA_TC = Paramètres utilisateurs de la loi ENDO_LOCA_TC
      E              = Module de Young
      NU             = Coefficient de Poisson
      FT             = Limite en traction simple
      FC             = Limite en compression simple
      SIG0           = Limite de linearité compression simple
      GF             = Energie de fissuration
      P              = Parametre dominant de la loi cohésive asymptotique
      DIST_FISSURE   = Distance moyenne inter-fissure
      REST_RIGI_FC   = Restauration de rigidité pour eps=fc/E (0=sans)
      REGU_REDU_SEUIL= Facteur de reduction du seuil par regularisation
      TAU_REGU_VISC  = Temps caractéristique de la régularisation visqueuse
    """

    MATER = DMATER.cree_dict_valeurs(DMATER.mc_liste)

    # Lecture et interprétation des paramètres utilisateurs
    young = float(MATER["E"])
    nu = float(MATER["NU"])
    gf = float(MATER["GF"])
    ft = float(MATER["FT"])
    fc = float(MATER["FC"])
    p = float(MATER["P"])
    lf = float(MATER["DIST_FISSURE"])
    rrc = float(MATER["REST_RIGI_FC"])
    sig0 = float(MATER["SIG0"])
    eta = float(MATER["COEF_REDU_SEUIL"])
    tauv = float(MATER["TAU_REGU_VISC"])

    # Paramètres internes au modèle
    lbd = young * nu / ((1 + nu) * (1 - 2 * nu))
    dmu = young / (1 + nu)
    ec = lbd + dmu
    wc = ft**2 / (2 * ec)
    kappa = gf / (lf * wc)

    # Controle de la distance inter-fissure
    dc = 0.75 * math.pi * gf / wc * (p + 2.0) ** (-1.5)
    if lf > dc:
        UTMESS("F", "COMPOR1_97", valr=(dc,))

    # Paramètres pour DEFI_MATERIAU

    prms = dict(FT=ft, KAPPA=kappa, P=p, SIG0=sig0, FC=fc)

    gamma = young / (fc * 2.0 * (1 - rrc))
    prms["REST_RIGIDITE"] = gamma

    regu_p = math.log(3.0) / math.log(1.0 / eta + (1.0 / eta - 1) * 3 * nu / (1 - 2 * nu))
    prms["CRIT_REGU"] = regu_p

    prms["TAU_REGU_VISC"] = tauv

    mclef = elastic_properties(young, nu, args)
    mclef["ENDO_LOCA_TC"] = prms

    return mclef


def elastic_properties(E, NU, args):
    """Returns the properties for elasticity.

    Also add generic keywords (INFO).

    Arguments:
        args (dict): User arguments.

    Returns:
        dict: Updated `mclef` dict.
    """
    mclef = {"INFO": args.get("INFO", 1), "ELAS": {"E": E, "NU": NU}}
    list_para = ["RHO", "ALPHA", "AMOR_ALPHA", "AMOR_BETA", "AMOR_HYST"]
    for para in list_para:
        if args.get(para) is not None:
            mclef["ELAS"][para] = args[para]
    return mclef


def defi_mater_gc_ops(
    self,
    MAZARS=None,
    ACIER=None,
    ENDO_FISS_EXP=None,
    ENDO_LOCA_EXP=None,
    ENDO_LOCA_TC=None,
    BETON_GLRC=None,
    **args
):
    """
    C'est : un parmi : ACIER  MAZARS  ENDO_FISS_EXP, ENDO_LOCA_EXP, ENDO_LOCA_TC, BETON_GLRC
    """

    #
    if MAZARS is not None:
        mclef = Mazars_Unil(MAZARS[0], args)
    if ACIER is not None:
        mclef = Acier_Cine_Line(ACIER[0], args)
    if ENDO_FISS_EXP is not None:
        mclef = Endo_Fiss_Exp(ENDO_FISS_EXP[0], args)
    if ENDO_LOCA_EXP is not None:
        mclef = Endo_Loca_Exp(ENDO_LOCA_EXP[0], args)
    if ENDO_LOCA_TC is not None:
        mclef = Endo_Loca_TC(ENDO_LOCA_TC[0], args)
    if BETON_GLRC is not None:
        mclef = Beton_GLRC(BETON_GLRC[0], args)

    Materiau = DEFI_MATERIAU(**mclef)
    return Materiau
