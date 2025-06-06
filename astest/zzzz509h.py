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

import os
import numpy
from functools import lru_cache


from code_aster.Commands import *
from code_aster import CA

test = CA.TestCase()

DEBUT(CODE="OUI")

# mesh and properties extracted from forma03
mesh = LIRE_MAILLAGE(FORMAT="MED")

mesh = MODI_MAILLAGE(reuse=mesh, MAILLAGE=mesh, ORIE_PEAU=_F(GROUP_MA_PEAU="haut"))

model = AFFE_MODELE(
    MAILLAGE=mesh, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="C_PLAN")
)

curve = LIRE_FONCTION(UNITE=21, NOM_PARA="EPSI", PROL_DROITE="CONSTANT")

steel = DEFI_MATERIAU(
    ELAS=_F(E=200000.0, NU=0.3), TRACTION=_F(SIGM=curve), ECRO_LINE=_F(D_SIGM_EPSI=1930.0, SY=181.0)
)

material = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=steel))

symmetry = AFFE_CHAR_CINE(
    MODELE=model, MECA_IMPO=(_F(GROUP_MA="bas", DY=0.0), _F(GROUP_MA="gauche", DX=0.0))
)

load = AFFE_CHAR_MECA(MODELE=model, FORCE_CONTOUR=_F(GROUP_MA="haut", FY=1.0))

mult_func = DEFI_FONCTION(NOM_PARA="INST", VALE=(0.0, 0.0, 200.0, 1000.0))

time_start = 0.0
time_inter = 20.0
time_end = 30.0

time_values = DEFI_LIST_REEL(
    DEBUT=time_start, INTERVALLE=(_F(JUSQU_A=time_inter, NOMBRE=4), _F(JUSQU_A=time_end, NOMBRE=8))
)
times = DEFI_LIST_INST(DEFI_LIST=_F(LIST_INST=time_values))


@lru_cache(maxsize=16)
def testRestart(command, restart_from, time_restart, reuse=False, info=1):
    """unittest for restarting from a Result or from Fields.

    *Needs previously created objects from parent context.*

    Arguments:
        command (Command): Command to be used, STAT_NON_LINE or MECA_NON_LINE.
        restart_from (str): One of "result", "crea_champ", "getField".
        time_restart (float): Time of the assignment of the initial state
            (INST_INIT). The initial state of the second stage is always
            extracted from the last step of the first stage.
        reuse (bool): Reuse the same result object or not.
        info (int): Verbosity level.

    Returns:
        tuple: Results of the both steps.
    """
    args = _F(
        MODELE=model,
        CHAM_MATER=material,
        EXCIT=(_F(CHARGE=symmetry), _F(CHARGE=load, FONC_MULT=mult_func)),
        COMPORTEMENT=_F(RELATION="VMIS_ISOT_LINE"),
        # CONVERGENCE=_F(RESI_GLOB_RELA=1.0e-6),  # test default
        INFO=info,
    )

    first = command(INCREMENT=_F(LIST_INST=times, INST_FIN=time_inter), **args)
    first.userName = ("SNL1" if command is STAT_NON_LINE else "MNL1") + "_" + restart_from

    last = first.getLastTime()
    test.assertAlmostEqual(last, time_inter)

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

    if restart_from == "result":
        args["ETAT_INIT"] = _F(EVOL_NOLI=first)  # INST=last,
        if reuse:
            args["RESULTAT"] = first
    else:
        args["ETAT_INIT"] = _F(DEPL=depl, SIGM=sigm, VARI=vari)

    cont = command(INCREMENT=_F(LIST_INST=times, INST_INIT=time_restart), **args)
    cont.userName = ("SNL2" if command is STAT_NON_LINE else "MNL2") + "_" + restart_from
    return first, cont


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


def compare_exec(case1, case2):
    """Compare two executions.

    Arguments:
        case1 (dict): keywords defining the first execution.
        case2 (dict): keywords defining the second execution.
    """

    def _key(case):
        key = case.copy()
        key["command"] = "SNL" if case["command"] is STAT_NON_LINE else "MNL"
        key = dict([(k, str(v)) for k, v in key.items()])
        key = tuple(sorted(key.values()))
        return key

    res10, res11 = testRestart(**case1)
    res20, res21 = testRestart(**case2)

    key1 = _key(case1)
    key2 = _key(case2)
    print(f"\n\n+++ comparing\n{key1}\n{key2}\n", flush=True)
    for time_i in res10.getAccessParameters()["INST"]:
        compareResults(res10, res20, time_i)
    for time_i in res11.getAccessParameters()["INST"]:
        compareResults(res11, res21, time_i)


def check_fields(result, names):
    """Print a checksum of some fields of a result.

    Arguments:
        result (Result): Result object
        names (str): Field to be extracted.
    """
    params = result.getAccessParameters()
    for idx, time_i in zip(params["NUME_ORDRE"], params["INST"]):
        for field_name in names:
            field = result.getField(field_name, idx)
            ojb = ".VALE" if isinstance(field, CA.FieldOnNodesReal) else ".CELV"
            print("UTIMSD:", idx, time_i, field_name, flush=True)
            IMPR_CO(CHAINE=field.getName() + ojb, NIVEAU=-1, UNITE=6)


# for debugging a case, use: CASE=4 .../run_aster astest/zzzz506z.export
case = int(os.environ.get("CASE", "0"))
run_failed = case != 0

# --- from last previous step
# classic restart using 2 different results
if case in (0, 1):
    compare_exec(
        _F(command=STAT_NON_LINE, restart_from="result", reuse=False, time_restart=20.0),
        _F(command=MECA_NON_LINE, restart_from="result", reuse=False, time_restart=20.0),
    )

# classic restart using the same result
if case in (0, 2):
    compare_exec(
        _F(command=STAT_NON_LINE, restart_from="result", reuse=True, time_restart=20.0),
        _F(command=MECA_NON_LINE, restart_from="result", reuse=True, time_restart=20.0),
    )

# restart from the result or fields from CREA_CHAMP for SNL
if case in (0, 3):
    compare_exec(
        _F(command=STAT_NON_LINE, restart_from="result", reuse=False, time_restart=20.0),
        _F(command=STAT_NON_LINE, restart_from="crea_champ", reuse=False, time_restart=20.0),
    )

# restart from the result or fields from CREA_CHAMP for MNL
if case in (0, 4):
    compare_exec(
        _F(command=MECA_NON_LINE, restart_from="result", reuse=False, time_restart=20.0),
        _F(command=MECA_NON_LINE, restart_from="crea_champ", reuse=False, time_restart=20.0),
    )

# restart from fields from CREA_CHAMP or getField for SNL
if case in (0, 5) and run_failed:
    compare_exec(
        _F(command=STAT_NON_LINE, restart_from="getField", reuse=False, time_restart=20.0),
        _F(command=STAT_NON_LINE, restart_from="crea_champ", reuse=False, time_restart=20.0),
    )

# restart from fields from CREA_CHAMP or getField for MNL
if case in (0, 6):
    compare_exec(
        _F(command=MECA_NON_LINE, restart_from="getField", reuse=False, time_restart=20.0),
        _F(command=MECA_NON_LINE, restart_from="crea_champ", reuse=False, time_restart=20.0),
    )

# restart using 2 different results
if case in (0, 7) and run_failed:
    compare_exec(
        _F(command=STAT_NON_LINE, restart_from="result", reuse=False, time_restart=15.0),
        _F(command=MECA_NON_LINE, restart_from="result", reuse=False, time_restart=15.0),
    )

# restart using the same result
if case in (0, 8) and run_failed:
    compare_exec(
        _F(command=STAT_NON_LINE, restart_from="result", reuse=True, time_restart=15.0),
        _F(command=MECA_NON_LINE, restart_from="result", reuse=True, time_restart=15.0),
    )

# restart from the result or fields from CREA_CHAMP for SNL
if case in (0, 9):
    compare_exec(
        _F(command=STAT_NON_LINE, restart_from="result", reuse=False, time_restart=15.0),
        _F(command=STAT_NON_LINE, restart_from="crea_champ", reuse=False, time_restart=15.0),
    )

# restart from the result or fields from CREA_CHAMP for MNL
if case in (0, 10):
    compare_exec(
        _F(command=MECA_NON_LINE, restart_from="result", reuse=False, time_restart=15.0),
        _F(command=MECA_NON_LINE, restart_from="crea_champ", reuse=False, time_restart=15.0),
    )

# restart from fields from CREA_CHAMP or getField for MNL
if case in (0, 11):
    compare_exec(
        _F(command=MECA_NON_LINE, restart_from="getField", reuse=False, time_restart=15.0),
        _F(command=MECA_NON_LINE, restart_from="crea_champ", reuse=False, time_restart=15.0),
    )

# debugging each step of case 7:
if case == 70:
    res = testRestart(command=STAT_NON_LINE, restart_from="result", reuse=False, time_restart=15.0)
    check_fields(res[0], ["DEPL", "SIEF_ELGA", "VARI_ELGA"])
    check_fields(res[1], ["DEPL", "SIEF_ELGA", "VARI_ELGA"])

if case == 71:
    res = testRestart(command=MECA_NON_LINE, restart_from="result", reuse=False, time_restart=15.0)
    check_fields(res[0], ["DEPL", "SIEF_ELGA", "VARI_ELGA"])
    check_fields(res[1], ["DEPL", "SIEF_ELGA", "VARI_ELGA"])

FIN()
