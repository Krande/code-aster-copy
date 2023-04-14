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

import numpy

import code_aster
from code_aster.Commands import *

test = code_aster.TestCase()

# ***********************************************************************
#
# TITRE: ESSAI DE TRACTION AVEC LA LOI DE RANKINE
# MODIFICATION DU TEST SSNV515C
#
# ***********************************************************************

# ======================================================================
DEBUT(CODE=_F(NIV_PUB_WEB="INTERNET"))

# ***********************************************************************
#
#    MAILLAGE ET MODELE
#
# ***********************************************************************

MAILLAGE = LIRE_MAILLAGE(FORMAT="MED", PARTITIONNEUR="PTSCOTCH")

MODELE = AFFE_MODELE(
    MAILLAGE=MAILLAGE, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="AXIS")
)

MAILLAGE = MODI_MAILLAGE(
    reuse=MAILLAGE,
    MAILLAGE=MAILLAGE,
    ORIE_PEAU=_F(GROUP_MA_PEAU=("DROIT", "GAUCHE", "BAS", "HAUT")),
    INFO=1,
)

# ***********************************************************************
#
#    LISTE D'INSTANTS
#
# ***********************************************************************

tarret = 10.0
tfin = 20.0
npas = 30
temps_max = 30.0

dtemps = temps_max / npas

ltemps = [dtemps * i for i in range(npas + 1)]

TEMPS = DEFI_LIST_REEL(DEBUT=0.0, INTERVALLE=(_F(JUSQU_A=20, NOMBRE=10),))
TEMPS2 = DEFI_LIST_REEL(DEBUT=20, INTERVALLE=(_F(JUSQU_A=temps_max, NOMBRE=5),))

# ***********************************************************************
#
#    DONNEES MATERIAU
#
# ***********************************************************************

# modules mecaniques [kPa]
# ------------------------
# K=516.2E6
# G=238.2E6
# # =>
# YOUNG = 9.*K*G /(3.*K+G)
# POISSON = (3.*K-2.*G) /(6.*K+2.*G)

YOUNG = 1e6
POISSON = 0.25
SIGMA_T = 1.0e3
K = YOUNG / 3.0 / (1.0 - 2.0 * POISSON)
G = YOUNG / 2.0 / (1.0 + POISSON)

SOL = DEFI_MATERIAU(ELAS=_F(E=YOUNG, NU=POISSON, ALPHA=0.0), RANKINE=_F(SIGMA_T=SIGMA_T), INFO=1)

CHMAT = AFFE_MATERIAU(MAILLAGE=MAILLAGE, AFFE=_F(TOUT="OUI", MATER=SOL))

# ***********************************************************************
#
#    CHARGEMENTS
#
# ***********************************************************************

# pression de preconsolidation [en kPa]
P0 = 5.0e3
EPZZ = 0.03

npas = 300
dtemps = temps_max / npas
linst = [dtemps * i for i in range(npas)]

SIGLAT = AFFE_CHAR_MECA(MODELE=MODELE, PRES_REP=_F(GROUP_MA=("DROIT",), PRES=P0))

DEPHAUT = AFFE_CHAR_CINE(MODELE=MODELE, MECA_IMPO=(_F(GROUP_MA=("HAUT",), DY=1.0),))

DEPL_1 = AFFE_CHAR_CINE(
    MODELE=MODELE, MECA_IMPO=(_F(GROUP_MA="BAS", DY=0.0), _F(GROUP_NO="B", DX=0.0))
)


COEF2 = DEFI_FONCTION(
    NOM_PARA="INST", PROL_DROITE="CONSTANT", VALE=(0.0, 0.0, 20, EPZZ, temps_max, 0)
)

COEF3 = DEFI_FONCTION(NOM_PARA="INST", PROL_DROITE="CONSTANT", VALE=(0.0, 1.0))


# ***********************************************************************
#
#    PRECONSOLIDATION INITIALE A 5KPA
#
# ***********************************************************************

SIG0 = CREA_CHAMP(
    INFO=2,
    TYPE_CHAM="ELGA_SIEF_R",
    OPERATION="AFFE",
    MODELE=MODELE,
    PROL_ZERO="OUI",
    AFFE=_F(GROUP_MA="BLOC", NOM_CMP=("SIXX", "SIYY", "SIZZ"), VALE=(-P0, -P0, -P0)),
)

time_init = 20.0
time_last = temps_max


def testRestart(command, restart_from, info=1):
    """unittest for restarting from a Result or from Fields.

    *Needs previously created objects from parent context.*

    Arguments:
        command (Command): Command to be used, STAT_NON_LINE or MECA_NON_LINE.
        restart_from (str): One of "result", "crea_champ", "getField".

    Returns:
        Result: Result of the second step.
    """
    first = command(
        MODELE=MODELE,
        CHAM_MATER=CHMAT,
        EXCIT=(
            _F(CHARGE=SIGLAT, FONC_MULT=COEF3),
            _F(CHARGE=DEPHAUT, FONC_MULT=COEF2),
            _F(CHARGE=DEPL_1, FONC_MULT=COEF3),
        ),
        ETAT_INIT=_F(SIGM=SIG0),
        COMPORTEMENT=_F(RELATION="RANKINE"),
        NEWTON=_F(MATRICE="TANGENTE", PREDICTION="ELASTIQUE", REAC_ITER=1),
        CONVERGENCE=_F(RESI_GLOB_RELA=1.0e-6, ITER_GLOB_MAXI=30, ARRET="OUI"),
        SOLVEUR=_F(METHODE="MUMPS", NPREC=8),
        INCREMENT=_F(LIST_INST=TEMPS),
    )

    last = first.getAccessParameters()["INST"][-1]

    if restart_from == "result":
        pass
    elif restart_from == "getField":
        depl = first.getField("DEPL", para="INST", value=last)
        sigm = first.getField("SIEF_ELGA", para="INST", value=last)
        vari = first.getField("VARI_ELGA", para="INST", value=last)
    else:
        depl = CREA_CHAMP(
            TYPE_CHAM="NOEU_DEPL_R", OPERATION="EXTR", RESULTAT=first, NOM_CHAM="DEPL", INST=last
        )
        sigm = CREA_CHAMP(
            TYPE_CHAM="ELGA_SIEF_R",
            OPERATION="EXTR",
            RESULTAT=first,
            NOM_CHAM="SIEF_ELGA",
            INST=last,
        )
        vari = CREA_CHAMP(
            TYPE_CHAM="ELGA_VARI_R",
            OPERATION="EXTR",
            RESULTAT=first,
            NOM_CHAM="VARI_ELGA",
            INST=last,
        )

    args = _F(
        MODELE=MODELE,
        CHAM_MATER=CHMAT,
        EXCIT=(
            _F(CHARGE=SIGLAT, FONC_MULT=COEF3),
            _F(CHARGE=DEPHAUT, FONC_MULT=COEF2),
            _F(CHARGE=DEPL_1, FONC_MULT=COEF3),
        ),
        COMPORTEMENT=_F(RELATION="RANKINE"),
        NEWTON=_F(MATRICE="TANGENTE", PREDICTION="ELASTIQUE", REAC_ITER=1),
        CONVERGENCE=_F(RESI_GLOB_RELA=1.0e-6, ITER_GLOB_MAXI=10, ARRET="OUI"),
        SOLVEUR=_F(METHODE="MUMPS", NPREC=8),
        INCREMENT=_F(LIST_INST=TEMPS2),
        INFO=info,
    )
    if restart_from == "result":
        args["ETAT_INIT"] = _F(EVOL_NOLI=first, INST_ETAT_INIT=time_init)  # INST=last,
    else:
        args["ETAT_INIT"] = _F(INST_ETAT_INIT=time_init, DEPL=depl, SIGM=sigm, VARI=vari)
    cont = command(**args)
    cont.userName = ("SNL" if command is STAT_NON_LINE else "MNL") + "_" + restart_from
    return cont


def compareResults(res1, res2, time, fields=("DEPL", "SIEF_ELGA", "VARI_ELGA"), label=""):
    """Compare the fields from two results at the given index.

    Arguments:
        res1 (Result): First result.
        res2 (Result): Second result.
        time (float): Time to be extracted.
        fields (list[str]): Fields to be compared.
        label (str): Label of the test.
    """
    if not label:
        label = f"{res1.userName} vs {res2.userName}: "
    for field_name in fields:
        try:
            field1 = res1.getField(field_name, para="INST", value=time)
        except ValueError as exc:
            field1 = None
            test.assertIsNotNone(field1, msg=f"{res1.userName}: {exc}")
        try:
            field2 = res2.getField(field_name, para="INST", value=time)
        except ValueError as exc:
            field2 = None
            test.assertIsNotNone(field2, msg=f"{res2.userName}: {exc}")
        if not field1 or not field2:
            continue

        field_diff = field1 - field2
        values = numpy.abs(numpy.array(field_diff.getValues()))
        diff = numpy.max(values)
        msg = f"{label}{field_name} at {time}"
        test.assertLessEqual(diff, 1.0e-10, msg=msg)


snl1 = testRestart(STAT_NON_LINE, restart_from="result")
mnl1 = testRestart(MECA_NON_LINE, restart_from="result")
compareResults(snl1, mnl1, time_init)
compareResults(snl1, mnl1, time_last)

# snl2 = testRestart(STAT_NON_LINE, restart_from="crea_champ")
# mnl2 = testRestart(MECA_NON_LINE, restart_from="crea_champ")
# compareResults(snl2, mnl2, time_init)
# compareResults(snl2, mnl2, time_last)

# compareResults(snl1, snl2, time_init)
# compareResults(snl1, snl2, time_last)
# compareResults(mnl1, mnl2, time_init)
# compareResults(mnl1, mnl2, time_last)

# # snl3 = testRestart(STAT_NON_LINE, restart_from="getField")
# mnl3 = testRestart(MECA_NON_LINE, restart_from="getField")


# compareResults(mnl1, mnl3, time_init)
# compareResults(mnl1, mnl3, time_last)
# compareResults(snl3, mnl3, time_last)

test.assertTrue(True)
test.printSummary()

FIN()
