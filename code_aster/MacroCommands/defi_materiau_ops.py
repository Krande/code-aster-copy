# coding: utf-8

# Copyright (C) 1991 - 2022  EDF R&D                www.code-aster.org
#
# This file is part of Code_Aster.
#
# Code_Aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Code_Aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Code_Aster.  If not, see <http://www.gnu.org/licenses/>.

import numpy
from libaster import createEnthalpy, resetFortranLoggingLevel, setFortranLoggingLevel

from ..Cata.Syntax import _F
from ..Messages import UTMESS
from ..Objects import Function, Material


def defi_materiau_ops(self, **args):
    """Execute the command.

    Arguments:
        **args (dict): User's keywords.
    """

    args = _F(args)

    setFortranLoggingLevel(args["INFO"])

    # reuse, copy from or new object
    if args.get("reuse"):
        assert args["reuse"] == args["MATER"]
        mater = args["MATER"]
    elif args.get("MATER"):
        mater = Material(args["MATER"])
    else:
        mater = Material()

    # In this function, we can check the value of keywords and add some properties
    check_keywords(args)

    visit = mater.Builder(mater)
    self._cata.accept(visit, args)

    resetFortranLoggingLevel()
    check_young_consistency(mater)
    if args["INFO"] == 2:
        mater.debugPrint()

    return mater


def check_young_consistency(mater):
    """Check the consistency between E provided as a constant under 'ELAS'
    and the first point of the 'TRACTION' function."""
    matNames = mater.getMaterialNames()
    if "ELAS" not in matNames or "TRACTION" not in matNames:
        return
    trac = mater.getFunction("TRACTION", "SIGM")
    typ = trac.getProperties()[0]
    if typ == "NAPPE":
        UTMESS("I", "MATERIAL1_6")
    if typ != "FONCTION":
        # formula not supported here
        return
    moduleE = mater.getValueReal("ELAS", "E")
    values = trac.getValuesAsArray()
    moduleTrac = values[0][1] / values[0][0]
    if abs(moduleE - moduleTrac) / moduleE > 0.01:
        UTMESS("A", "MATERIAL1_5", valr=moduleTrac)


# internal methods
def check_keywords(kwargs):
    """Check for DEFI_MATERIAU keywords

    Arguments:
        kwargs (dict): User's keywords, changed in place.
    """

    if "DIS_ECRO_TRAC" in kwargs:
        kwargs["DIS_ECRO_TRAC"] = check_dis_ecro_trac(kwargs["DIS_ECRO_TRAC"])
    if "DIS_CHOC_ENDO" in kwargs:
        kwargs["DIS_CHOC_ENDO"] = check_dis_choc_endo(kwargs["DIS_CHOC_ENDO"])
    if "JONC_ENDO_PLAS" in kwargs:
        kwargs["JONC_ENDO_PLAS"] = check_dis_jvp(kwargs["JONC_ENDO_PLAS"])
    if "THER_NL" in kwargs:
        kwargs["THER_NL"] = add_enthalpy(kwargs["THER_NL"])


def check_dis_ecro_trac(keywords):
    """Check for function for DIS_ECRO_TRAC

    Arguments:
        keywords (dict): DIS_ECRO_TRAC keyword, changed in place.

    Returns:
        dict: DIS_ECRO_TRAC keyword changed in place.

        Raises '<F>' in case of error.
    """
    #
    # jean-luc.flejou@edf.fr
    #
    def _message(num, mess=""):
        UTMESS(
            "F", "DISCRETS_62", valk=("DIS_ECRO_TRAC", "FX=f(DX) | FTAN=f(DTAN)", mess), vali=num
        )

    precis = 1.0e-08
    #
    Clefs = keywords
    #
    if "FX" in Clefs:
        iffx = True
        fct = Clefs["FX"]
    elif "FTAN" in Clefs:
        iffx = False
        fct = Clefs["FTAN"]
    else:
        _message(1)
    # Les vérifications sur la fonction
    #       interpolation LIN LIN
    #       prolongée à gauche et à droite exclue
    #       paramètre 'DX' ou 'DTAN'
    OkFct = type(fct) is Function
    param = fct.Parametres()
    OkFct = OkFct and param["INTERPOL"][0] == "LIN"
    OkFct = OkFct and param["INTERPOL"][1] == "LIN"
    OkFct = OkFct and param["PROL_DROITE"] == "EXCLU"
    OkFct = OkFct and param["PROL_GAUCHE"] == "EXCLU"
    if iffx:
        OkFct = OkFct and param["NOM_PARA"] == "DX"
    else:
        OkFct = OkFct and param["NOM_PARA"] == "DTAN"
    if not OkFct:
        _message(2, "%s" % param)
    # avoir 3 points minimum ou exactement
    absc, ordo = fct.Valeurs()
    if iffx:
        OkFct = OkFct and len(absc) >= 3
    else:
        if Clefs["ECROUISSAGE"] == "ISOTROPE":
            OkFct = OkFct and len(absc) >= 3
        elif Clefs["ECROUISSAGE"] == "CINEMATIQUE":
            OkFct = OkFct and len(absc) == 3
        else:
            raise RuntimeError("Unknown value")
    #
    if not OkFct:
        _message(3, "%s" % len(absc))
    # Point n°1: (DX=0, FX=0)
    dx = absc[0]
    fx = ordo[0]
    OkFct = OkFct and dx >= 0.0 and abs(dx) <= precis
    OkFct = OkFct and fx >= 0.0 and abs(fx) <= precis
    if not OkFct:
        _message(4, "[%s %s]" % (dx, fx))
    # FX et DX sont strictement positifs, dFx >0 , dDx >0
    #   Au lieu de la boucle, on peut faire :
    #       xx=np.where(np.diff(absc) <= 0.0 or np.diff(ordo) <= 0.0)[0]
    #       if ( len(xx) != 0):
    #           message 6 : absc[xx[0]] ordo[xx[0]] absc[xx[0]+1] ordo[xx[0]+1]
    #       dfx= np.diff(ordo)/np.diff(absc)
    #       xx=np.where(np.diff(dfx) <= 0.0 )[0]
    #       if ( len(xx) != 0):
    #           message 7 : xx[0] dfx[xx[0]] dfx[xx[0]+1]
    for ii in range(1, len(absc)):
        if absc[ii] <= dx or ordo[ii] <= fx:
            _message(
                5, "Ddx, Dfx > 0 : p(i)[%s %s]  " "p(i+1)[%s %s]" % (dx, fx, absc[ii], ordo[ii])
            )
        if ii == 1:
            dfx = (ordo[ii] - fx) / (absc[ii] - dx)
            raidex = dfx
        else:
            dfx = (ordo[ii] - fx) / (absc[ii] - dx)
            if dfx > raidex:
                _message(6, "(%d) : %s > %s" % (ii, dfx, raidex))

        dx = absc[ii]
        fx = ordo[ii]
        raidex = dfx

    return Clefs


def check_dis_choc_endo(keywords):
    """Check for functions for DIS_CHOC_ENDO

    Arguments:
        keywords (dict): DIS_CHOC_ENDO keyword, changed in place.

    Returns:
        dict: DIS_CHOC_ENDO keyword changed in place. Raises '<F>' in case of
        error.
    """
    #
    # jean-luc.flejou@edf.fr
    #
    def _message(num, mess2="", mess3=""):
        UTMESS("F", "DISCRETS_63", valk=("DIS_CHOC_ENDO", mess2, mess3), vali=(num,))

    precis = 1.0e-08
    #
    Clefs = keywords

    # Conditions communes aux 3 fonctions
    #   paramètre 'DX'
    #   interpolation LIN LIN
    #   prolongée à gauche et à droite : constant ou exclue (donc pas linéaire)
    LesFcts = [Clefs["FX"], Clefs["RIGI_NOR"], Clefs["AMOR_NOR"]]
    LesFctsName = []
    for fct in LesFcts:
        OkFct = type(fct) is Function
        param = fct.Parametres()
        OkFct = OkFct and param["NOM_PARA"] == "DX"
        OkFct = OkFct and param["INTERPOL"][0] == "LIN"
        OkFct = OkFct and param["INTERPOL"][1] == "LIN"
        OkFct = OkFct and param["PROL_DROITE"] != "LINEAIRE"
        OkFct = OkFct and param["PROL_GAUCHE"] != "LINEAIRE"
        LesFctsName.append(fct.getName())
        if not OkFct:
            _message(2, fct.getName(), "%s" % param)
    # Même nombre de point, 5 points minimum
    Fxx, Fxy = Clefs["FX"].Valeurs()
    Rix, Riy = Clefs["RIGI_NOR"].Valeurs()
    Amx, Amy = Clefs["AMOR_NOR"].Valeurs()
    OkFct = len(Fxx) == len(Rix) == len(Amx)
    OkFct = OkFct and (len(Fxx) >= 5)
    if not OkFct:
        _message(3, "%s" % LesFctsName, "%d, %d, %d" % (len(Fxx), len(Rix), len(Amx)))
    # La 1ère abscisse c'est ZERO
    x1 = Fxx[0]
    if not (0.0 <= x1 <= precis):
        _message(4, "%s" % LesFctsName, "%s" % x1)
    # Même abscisses pour les 3 fonctions, à la précision près
    # Abscisses strictement croissantes, à la précision près
    #   Au lieu de la boucle, on peut faire :
    #       xx=np.where(np.abs(Fxx-Rix) + np.abs(Fxx-Amx) > precis )[0]
    #       if ( len(xx) != 0):
    #           message 5 : xx[0] Fxx[xx[0]] Rix[xx[0]] Amx[xx[0]]
    #       xx=np.where(np.diff(Fxx)<0 or np.diff(Rix)<0 or np.diff(Amx)<0)[0]
    #       if ( len(xx) != 0):
    #           message 6 : xx[0] Fxx[xx[0]] Fxx[xx[0]+1]
    xp1 = -1.0
    for ii in range(len(Fxx)):
        x1 = Fxx[ii]
        x2 = Rix[ii]
        x3 = Amx[ii]
        ddx = abs(x1 - x2) + abs(x1 - x3)
        if ddx > precis:
            _message(5, "%s" % LesFctsName, "(%d) : %s, %s, %s" % (ii + 1, x1, x2, x3))
        if xp1 >= x1:
            _message(6, "%s" % LesFctsName, "(%d) : %s, %s" % (ii + 1, xp1, x1))
        xp1 = x1
    # FX : les 2 premiers points ont la même valeur à la précision relative près
    if abs(Fxy[0] - Fxy[1]) > Fxy[0] * precis:
        _message(7, Clefs["FX"].getName())
    # RIGI_NOR : les 2 premiers points ont la même valeur à la précision relative près
    if abs(Riy[0] - Riy[1]) > Riy[0] * precis:
        _message(8, Clefs["RIGI_NOR"].getName())
    # RIGI_NOR : pente décroissante, à partir du 3ème point.
    #   Au lieu de la boucle, on peut faire :
    #       xx=np.where(np.diff(Riy[2:]) > 0.0 )[0]+2
    #       if ( len(xx) != 0):
    #           message 9 : xx[0] Riy[xx[0]] Riy[xx[0]+1]
    pente = Riy[2]
    for ii in range(3, len(Riy)):
        if Riy[ii] - pente > Riy[0] * precis:
            _message(9, "%s" % Clefs["RIGI_NOR"].getName(), "(%d) : %s %s" % (ii, pente, Riy[ii]))
        pente = Riy[ii]
    # --------------------------------------------------------------- Fin des vérifications
    # Création des fonctions
    newFx = Function()
    newFx.setResultName("Force   (plast cumulée)")
    newRi = Function()
    newRi.setResultName("Raideur (plast cumulée)")
    newAm = Function()
    newAm.setResultName("Amort.  (plast cumulée)")
    newFx.setExtrapolation("CC")
    newRi.setExtrapolation("CC")
    newAm.setExtrapolation("CC")
    newFx.setInterpolation("LIN LIN")
    newRi.setInterpolation("LIN LIN")
    newAm.setInterpolation("LIN LIN")
    newFx.setParameterName("PCUM")
    newRi.setParameterName("PCUM")
    newAm.setParameterName("PCUM")
    # Modification des abscisses pour avoir la plasticité cumulée et plus le déplacement
    pp9 = numpy.round(numpy.array(Fxx) - numpy.array(Fxy) / numpy.array(Riy), 10)
    # Le 1er point a une abscisse négative et la valeur du 2ème point
    pp9[0] = -pp9[2]
    # Vérification que p est strictement croissant.
    #   Au lieu de la boucle, on peut faire :
    #       xx=np.where(np.diff(pp9) < 0)[0]
    #       if ( len(xx) != 0):
    #           message 10 : xx[0] pp9[xx[0]] pp9[xx[0]+1]
    pp = pp9[0]
    for ii in range(1, len(pp9)):
        if pp > pp9[ii]:
            _message(10, mess3="(%d) : %s %s" % (ii, pp, pp9[ii]))
        pp = pp9[ii]
    #
    # Affectations des valeurs
    newFx.setValues(pp9, Fxy)
    newRi.setValues(pp9, Riy)
    newAm.setValues(pp9, Amy)
    # Affectation des nouvelles fonctions
    Clefs["FXP"] = newFx
    Clefs["RIGIP"] = newRi
    Clefs["AMORP"] = newAm
    #
    return Clefs


def check_dis_jvp(keywords):
    """Check for parameters in JONC_ENDO_PLAS

    Arguments:
        keywords (dict): JONC_ENDO_PLAS keyword, changed in place.

    Returns:
        dict: JONC_ENDO_PLAS keyword changed in place.

        Raises '<F>' in case of error.
    """

    def _message(num):
        UTMESS("F", "DISCRETS_65", valk=("JONC_ENDO_PLAS"), vali=num)

    #
    Clefs = keywords
    ke = Clefs["KE"]
    kp = Clefs["KP"]
    kdp = Clefs["KDP"]
    kdm = Clefs["KDM"]
    myp = Clefs["MYP"]
    mym = Clefs["MYM"]
    rdp = Clefs["RDP"]
    rdm = Clefs["RDM"]
    #
    if ke < kp:
        _message(1)
    elif kdp < kp or kdp > ke:
        _message(2)
    elif kdm < kp or kdm > ke:
        _message(3)
    elif myp < ke * rdp:
        _message(4)
    elif mym > ke * rdm:
        _message(5)
    return Clefs


def add_enthalpy(keywords):
    """Check for functions for THER_NL. Create "Beta" from "Rho_CP" if not given

    Arguments:
        keywords (dict): THER_NL keyword, changed in place.

    Returns:
        dict: THER_NL keyword changed in place. Raises '<F>' in case of
        error.
    """

    # Create "Beta" from "Rho_CP" if not given
    if keywords["BETA"] is None:
        beta = Function()
        createEnthalpy(keywords["RHO_CP"], beta)

        keywords["BETA"] = beta

    return keywords
