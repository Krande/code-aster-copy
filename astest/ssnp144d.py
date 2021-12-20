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

POURSUITE(CODE='OUI')

UTOT1.listFields()

depl = UTOT1.getFieldOnNodesReal("DEPL", 1)

# depl is real
try:
   depl_c = UTOT1.getFieldOnNodesComplex("DEPL", 1)
except IndexError:
   pass

sief_elga = UTOT1.getFieldOnCellsReal("SIEF_ELGA", 1)

# sief is real
try:
   sief_c = UTOT1.getFieldOnCellsComplex("SIEF_ELGA", 1)
except IndexError:
   pass

# sief is real
try:
   sief_i = UTOT1.getFieldOnCellsLong("SIEF_ELGA", 1)
except IndexError:
   pass


# not 150 fields
try:
   test = UTOT1.getFieldOnCellsLong("INDC_ELEM", 150)
except IndexError:
   pass

cohe_elem = UTOT1.getFieldOnCellsReal("COHE_ELEM", 1)

cont_noeu = UTOT1.getFieldOnNodesReal("CONT_NOEU", 1)

try:
   test = UTOT1.getConstantFieldOnCellsChar16("COMPORTEMENT", 150)
except IndexError:
   pass

try:
   test = UTOT1.getConstantFieldOnCellsReal("COMPORTEMENT", 1)
except IndexError:
   pass

FIN()
