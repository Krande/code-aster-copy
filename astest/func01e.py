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

import os

from code_aster.Commands import *
from code_aster import CA

"""On teste le bon fonctionnement de la libération de la mémoire pour l'évaluation des fonctions

1/ via l'évaluation directe : call Python à l’intérieur du C (evaluate_formula)
2/ via CALC_FONC_INTERP pour le Fortran (ffintf)
3/ via une formule dans un matériau (fointa) qu'on évaluera par l'appel à CALC_MATR_ELEM

On boucle nb_loop fois et on evalue par paquet de nb_para paramètres,
    ce qui fait un total de 3*nb_loop*nb_para evaluations
"""

long = False

if long:
    # run ~1 min
    nb_loop = 350
else:
    # run ~10 s
    nb_loop = 50
nb_para = 1000
# tole set to 2MB ~ 524288 floats
tole = 2048


CA.init("--test", ERREUR=_F(ALARME="EXCEPTION"))
test = CA.TestCase()


def get_mem():
    """return memory consumption in kB"""
    pid = os.getpid()
    proc_status_file = "/proc/%d/status" % pid
    with open(proc_status_file, "r") as f:
        status = f.read()
    return int(status[status.index("VmSize") :].split()[1])


fonc = DEFI_FONCTION(NOM_PARA="INST", VALE=(0.0, 20.0, 1.0, 30.0, 2.0, 40.0), VERIF="CROISSANT")


form = FORMULE(NOM_PARA=("INST"), VALE="FONC(INST)", FONC=fonc)

# ===
# check there is no memory leak with __call__

mem_ini = get_mem()
for i in range(nb_loop):
    for j in range(nb_para):
        form(2 * j / nb_para)
mem_end = get_mem()
print(mem_ini, mem_end)
test.assertTrue(mem_ini + tole >= mem_end)

# ===
# check there is no memory leak with ffintf

vale_para = [2 * j / nb_para for j in range(nb_para)]

mem_ini = get_mem()
for i in range(nb_loop):
    CALC_FONC_INTERP(FONCTION=form, NOM_PARA="INST", VALE_PARA=vale_para)

mem_end = get_mem()
print(mem_ini, mem_end)
test.assertTrue(mem_ini + tole >= mem_end)

# ===
# check there is no memory leak with fointa

nu = DEFI_CONSTANTE(VALE=0.3)

mater = DEFI_MATERIAU(ELAS_FO=_F(E=form, NU=nu))

mesh = LIRE_MAILLAGE()

model = AFFE_MODELE(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", PHENOMENE="MECANIQUE", MODELISATION="3D"))

chmat = AFFE_MATERIAU(MAILLAGE=mesh, AFFE=_F(TOUT="OUI", MATER=mater))

mem_ini = get_mem()
for j in range(nb_loop):
    CALC_MATR_ELEM(OPTION="RIGI_MECA", MODELE=model, INST=2 * j / nb_para, CHAM_MATER=chmat)

mem_end = get_mem()
print(mem_ini, mem_end)
test.assertTrue(mem_ini + tole >= mem_end)

CA.close()
