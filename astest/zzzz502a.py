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

from code_aster.Commands import *
from code_aster import CA

CA.init("--test", ERREUR=_F(ERREUR_F="EXCEPTION"))

test = CA.TestCase()

# Creation du maillage
monMaillage = CA.Mesh()

# Relecture du fichier MED
monMaillage.readMedFile("zzzz502a.mmed")

# Definition du modele Aster
monModel = CA.Model(monMaillage)
monModel.addModelingOnMesh(CA.Physics.Mechanics, CA.Modelings.Tridimensional)
monModel.build()

# Definition d'un chargement de type FORCE_NODALE à partir d'une ForceReal
force = CA.ForceReal()
force.setValue(CA.PhysicalQuantityComponent.Fz, 100.0)

print(" >>>> Construction d'un chargement NodalForceReal")
CharMeca1 = CA.NodalForceReal(monModel)
# On ne peut pas imposer une force nodale sur un groupe de mailles
with test.assertRaises(RuntimeError):
    CharMeca1.setValue(force, "UP")


nameOfGroup = "A"
CharMeca1.setValue(force, nameOfGroup)
print("      sur le groupe : ", nameOfGroup)
ret = CharMeca1.build()
# CharMeca1.debugPrint()
test.assertTrue(ret)

# Definition d'un chargement de type FORCE_NODALE à partir d'un StructuralForceReal

force_pour_structure = CA.StructuralForceReal()
force_pour_structure.setValue(CA.PhysicalQuantityComponent.Mx, 10.0)
force_pour_structure.setValue(CA.PhysicalQuantityComponent.My, 20.0)
force_pour_structure.setValue(CA.PhysicalQuantityComponent.Mz, 30.0)

print(" >>>> Construction d'un chargement NodalStructuralForceReal")
print("      Ce chargement est correct pour le catalogue mais conduit à une erreur Fortran ")
CharMeca2 = CA.NodalStructuralForceReal(monModel)
nameOfGroup = "B"
CharMeca2.setValue(force_pour_structure, nameOfGroup)
print("      sur le groupe : ", nameOfGroup)

# Le ddl MX n'est pas autorisé
# fortran error
with test.assertRaises(CA.AsterError):
    CharMeca2.build()

# CharMeca2.debugPrint()

# Definition d'un chargement de type FORCE_FACE à partir d'un ForceReal
print(" >>>> Construction d'un chargement ForceOnFaceReal")

CharMeca3 = CA.ForceOnFaceReal(monModel)
nameOfGroup = "UP"
CharMeca3.setValue(force, nameOfGroup)
print("      sur le groupe : ", nameOfGroup)
ret = CharMeca3.build()
test.assertTrue(ret)


# Definition d'un chargement de type FORCE_ARETE à partir d'un ForceReal
print(" >>>> Construction d'un chargement ForceOnEdgeReal")
CharMeca4 = CA.ForceOnEdgeReal(monModel)
nameOfGroup = "UP"
CharMeca4.setValue(force, nameOfGroup)
print("      sur le groupe : ", nameOfGroup)
ret = CharMeca4.build()
test.assertTrue(ret)

# Definition d'un chargement de type FORCE_ARETE à partir d'un StructuralForceReal
print(" >>>> Construction d'un chargement StructuralForceOnEdgeReal")
# C'est bizarre, on entre un groupe qui est une face et le fortran ne détecte rien !
CharMeca5 = CA.StructuralForceOnEdgeReal(monModel)
nameOfGroup = "UP"
CharMeca5.setValue(force_pour_structure, nameOfGroup)
print("      sur le groupe : ", nameOfGroup)
ret = CharMeca5.build()
test.assertTrue(ret)

# Definition d'un chargement de type FORCE_CONTOUR à partir d'un ForceReal
print(" >>>> Construction d'un chargement LineicForceReal")

CharMeca6 = CA.LineicForceReal(monModel)
nameOfGroup = "BOTTOM"
CharMeca6.setValue(force, nameOfGroup)
print("      sur le groupe : ", nameOfGroup)
ret = CharMeca6.build()
test.assertTrue(ret)

# Definition d'un chargement de type FORCE_INTERNE à partir d'un ForceReal
print(" >>>> Construction d'un chargement InternalForceReal")

CharMeca7 = CA.InternalForceReal(monModel)
nameOfGroup = "BOTTOM"
CharMeca7.setValue(force, nameOfGroup)
print("      sur le groupe : ", nameOfGroup)
ret = CharMeca7.build()
test.assertTrue(ret)

# Definition d'un chargement de type FORCE_POUTRE à partir d'un StructuralForceReal
print(" >>>> Construction d'un chargement StructuralForceOnBeamReal")

CharMeca8 = CA.StructuralForceOnBeamReal(monModel)
nameOfGroup = "OA"
CharMeca8.setValue(force_pour_structure, nameOfGroup)
print("      sur le groupe : ", nameOfGroup)
ret = CharMeca8.build()
test.assertTrue(ret)

# Definition d'un chargement de type FORCE_POUTRE à partir d'un LocalBeamForceReal
print(" >>>> Construction d'un chargement LocalForceOnBeamReal")

fpoutre = CA.LocalBeamForceReal()
fpoutre.setValue(CA.PhysicalQuantityComponent.N, 5.0)

CharMeca9 = CA.LocalForceOnBeamReal(monModel)
nameOfGroup = "BOTTOM"
CharMeca9.setValue(fpoutre, nameOfGroup)
print("      sur le groupe : ", nameOfGroup)
ret = CharMeca9.build()
test.assertTrue(ret)

# Definition d'un chargement de type FORCE_COQUE à partir d'un StructuralForceReal
print(" >>>> Construction d'un chargement StructuralForceOnShellReal")

CharMeca10 = CA.StructuralForceOnShellReal(monModel)
nameOfGroup = "UP"
CharMeca10.setValue(force_pour_structure, nameOfGroup)
print("      sur le groupe : ", nameOfGroup)
ret = CharMeca10.build()
test.assertTrue(ret)

# Definition d'un chargement de type FORCE_COQUE à partir d'un LocalShellForceReal
print(" >>>> Construction d'un chargement LocalForceOnShellReal")

fshell = CA.LocalShellForceReal()
fshell.setValue(CA.PhysicalQuantityComponent.F1, 11.0)
fshell.setValue(CA.PhysicalQuantityComponent.F2, 12.0)
fshell.setValue(CA.PhysicalQuantityComponent.F3, 13.0)

CharMeca11 = CA.LocalForceOnShellReal(monModel)
nameOfGroup = "UP"
CharMeca11.setValue(fshell, nameOfGroup)
print("      sur le groupe : ", nameOfGroup)
ret = CharMeca11.build()
test.assertTrue(ret)

# Definition d'un chargement de type FORCE_COQUE à partir d'une PressureReal
print(" >>>> Construction d'un chargement PressureOnShellReal")

pression = CA.PressureReal()
pression.setValue(CA.PhysicalQuantityComponent.Pres, 14.0)

CharMeca12 = CA.PressureOnShellReal(monModel)
nameOfGroup = "UP"
CharMeca12.setValue(pression, nameOfGroup)
print("      sur le groupe : ", nameOfGroup)
ret = CharMeca12.build()
test.assertTrue(ret)

# Imposer une PressureReal sur un groupe de noeuds
print(" >>>> Construction d'un chargement ImposedPressureReal")
CharMeca13 = CA.ImposedPressureReal(monModel)
nameOfGroup = "O"
CharMeca13.setValue(pression, nameOfGroup)
print("      sur le groupe : ", nameOfGroup)
# fortran error
with test.assertRaises(CA.AsterError):
    CharMeca13.build()

test.printSummary()

CA.close()
