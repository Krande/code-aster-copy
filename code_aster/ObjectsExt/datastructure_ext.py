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

# person_in_charge: mathieu.courtois@edf.fr
"""
:py:class:`DataStructure` --- Base of all objects
*************************************************
"""

import libaster
from libaster import DataStructure

from ..Utilities import deprecated, import_object, injector
from ..Objects.Serialization import InternalStateBuilder


@injector(DataStructure)
class ExtendedDataStructure:
    """This class defines the base class of the DataStructures."""

    # Tell Boost that __get_state__/__set_state__ must manage __dict__
    # Search for python reference guide at https://www.boost.org/doc/libs/
    __getstate_manages_dict__ = 1
    cata_sdj = None
    ptr_class_sdj = None
    ptr_sdj = None
    internalStateBuilder = InternalStateBuilder

    orig_getName = DataStructure.getName

    def getName(self):
        """
        Overload standart getName() function to eliminate whitespace at both
        ends of the string.

        .. note:: The C++ constructor automaticaly adds a whitespace when it
            creates a new object from the result name of this overloaded
            function `getName()`.
        """
        return self.orig_getName().strip()

    def is_typco(self):
        """Tell if it is an auxiliary result."""
        return False

    def __getinitargs__(self):
        """Returns the argument required to reinitialize a derivated
        DataStructure object during unpickling.

        .. note:: This implementation does not satisfy any constructor of the
            base DataStructure. But most of the derivated class should have
            a constructor accepting the Jeveux name.
        """
        return (self.getName(),)

    def __getstate__(self):
        """Return internal state.

        Derivated *DataStructure* types should defined a dedicated *InternalStateBuilder*
        class to serialize its specific content.
        """
        return self.internalStateBuilder().save(self)

    def __setstate__(self, state):
        """Restore internal state.

        Arguments:
            state (*InternalStateBuilder*): Internal state.
        """
        assert isinstance(state, InternalStateBuilder), f"unexpected type: {state}"
        state.restore(self)

    @property
    def sdj(self):
        """Return the DataStructure catalog."""
        if self.ptr_sdj is None:
            cata_sdj = getattr(self, "cata_sdj", None)
            if not cata_sdj:
                cata_sdj = DICT_SDJ.get(self.__class__.__name__)
            assert (
                cata_sdj
            ), "The attribute 'cata_sdj' must be defined in " "the class {}".format(
                self.__class__.__name__
            )
            if self.ptr_class_sdj is None:
                self.ptr_class_sdj = import_object("code_aster." + cata_sdj)
            self.ptr_sdj = self.ptr_class_sdj(nomj=self.getName())
        return self.ptr_sdj

    def use_count(self):
        """Return the number of reference to the DataStructure.

        Warning: Use only for debugging! Supported datastructures in
        ``PythonBindings/DebugInterface.cxx``.
        """
        return libaster.use_count(self)

    # transitional functions - to remove later
    @property
    @deprecated(case=1, help="Use 'getName()' instead.")
    def nom(self):
        return self.getName()



# This dictionnary avoids to add the DataStructure "_ext.py" file just
# to define the SD definition.
DICT_SDJ = {
    "Contact": "SD.sd_contact.sd_contact",
    "CrackTip": "SD.sd_fond_fiss.sd_fond_fiss",
    "Crack": "SD.sd_fond_fissure.sd_fond_fissure",
    "HarmoGeneralizedResult": "SD.sd_dyna_gene.sd_dyna_gene",
    "MeshesMapping": "SD.sd_corresp_2_mailla.sd_corresp_2_mailla",

    "AcousticModeResult": "SD.sd_dyna_phys.sd_dyna_phys",
    "BehaviourDefinition": "SD.sd_compor.sd_compor",
    "BucklingModeResult": "SD.sd_dyna_phys.sd_dyna_phys",
    "CyclicSymmetryMode": "SD.sd_mode_cycl.sd_mode_cycl",
    "DataField": "SD.sd_champ.sd_champ",
    "ElementaryVector": "SD.sd_vect_elem.sd_vect_elem",
    "ElementaryVectorDisplacementReal": "SD.sd_vect_elem.sd_vect_elem",
    "ElementaryVectorPressureComplex": "SD.sd_vect_elem.sd_vect_elem",
    "ElementaryVectorTemperatureReal": "SD.sd_vect_elem.sd_vect_elem",
    "FiberGeometry": "SD.sd_gfibre.sd_gfibre",
    "FieldOnNodesComplex": "SD.sd_champ.sd_cham_no_class",
    "FieldOnNodesDescription" : "SD.sd_prof_chno.sd_prof_chno",
    "FiniteElementDescriptor": "SD.sd_ligrel.sd_ligrel",
    "FluidStructureInteraction": "SD.sd_type_flui_stru.sd_type_flui_stru",
    "FluidStructureModalBasis": "SD.sd_melasflu.sd_melasflu",
    "FullHarmonicAcousticResult": "SD.sd_dyna_phys.sd_dyna_phys",
    "FullHarmonicResult": "SD.sd_dyna_gene.sd_dyna_gene",
    "FullTransientResult": "SD.sd_dyna_phys.sd_dyna_phys",
    "GcpcSolver": "SD.sd_solveur.sd_solveur",
    "GeneralizedModeResult": "SD.sd_dyna_gene.sd_dyna_gene",
    "GeneralizedDOFNumbering": "SD.sd_nume_ddl_gene.sd_nume_ddl_gene",
    "GeneralizedResultComplex": "SD.sd_dyna_gene.sd_dyna_gene",
    "GeneralizedResultReal": "SD.sd_dyna_gene.sd_dyna_gene",
    "GenericModalBasis": "SD.sd_dyna_phys.sd_dyna_phys",
    "Grid": "SD.sd_grille.sd_grille",
    "InterspectralMatrix": "SD.sd_interspectre.sd_interspectre",
    "LdltSolver": "SD.sd_solveur.sd_solveur",
    "ModeResult": "SD.sd_dyna_phys.sd_dyna_phys",
    "ModeResultComplex": "SD.sd_dyna_phys.sd_dyna_phys",
    "MultFrontSolver": "SD.sd_solveur.sd_solveur",
    "MumpsSolver": "SD.sd_solveur.sd_solveur",
    "PetscSolver": "SD.sd_solveur.sd_solveur",
    "RitzBasis": "SD.sd_dyna_phys.sd_dyna_phys",
    "Skeleton": "SD.sd_squelette.sd_squelette",
    "StaticMacroElement": "SD.sd_macr_elem_stat.sd_macr_elem_stat",
    "StructureInterface": "SD.sd_interf_dyna_clas.sd_interf_dyna_clas",
    "TimeStepper": "SD.sd_list_inst.sd_list_inst",
    "TurbulentSpectrum": "SD.sd_spectre.sd_spectre",
}
