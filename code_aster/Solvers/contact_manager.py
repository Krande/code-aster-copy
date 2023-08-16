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

from .solver_features import SolverFeature
from .solver_features import SolverOptions as SOP
from ..Objects import ContactComputation, ContactPairing
from ..Utilities import no_new_attributes, profile


class ContactManager(SolverFeature):
    """Solve contact problem."""

    provide = SOP.Contact

    defi = pair = comp = None
    coef_cont = coef_frot = None
    phys_pb = None
    nb_pairing = 0
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, definition, phys_pb):
        """Initialize contact solver.

        Arguments:
            definition (ContactNew): contact definition
            phys_pb (PhysicalProblem): physical problem
        """
        super().__init__()
        self.nb_pairing = 0
        self.defi = definition
        self.phys_pb = phys_pb
        if self.defi is not None:
            self.pair = ContactPairing(self.defi)
            self.comp = ContactComputation(self.defi)
            self.coef_cont, self.coef_frot = self.comp.contactCoefficient()

    @profile
    def pairing(self, phys_pb):
        """Pairing between contact zones.

        Arguments:
            phys_pb (PhysicalProblem): physical problem
        """

        if self.enable:
            self.pair.compute()

            # Numbering does not change but connectivity yes.
            fed_defi = self.defi.getFiniteElementDescriptor()
            fed_pair = self.pair.getFiniteElementDescriptor()
            phys_pb.getListOfLoads().addContactLoadDescriptor(fed_defi, fed_pair)
            model = phys_pb.getModel()
            loads = phys_pb.getListOfLoads()
            phys_pb.getDOFNumbering().computeRenumbering(model, loads)

            self.nb_pairing += 1

    @profile
    def getPairingCoordinates(self):
        """Get the coordinates field used for pairing.

        Returns:
            MeshCoordinatesField: coordinates of nodes used for pairing:
        """

        if self.enable:
            return self.pair.getCoordinates()

        return None

    @profile
    def setPairingCoordinates(self, coor):
        """Set the coordinates field used for pairing.

        Returns:
            coor (MeshCoordinatesField): coordinates of nodes used for pairing:
        """

        if self.enable:
            self.pair.setCoordinates(coor)

    @profile
    def gap(self):
        """Compute geometric gap.

        Returns:
            (FieldOnNodesReal): geometric gap
            (FieldOnNodesReal): gap indicator
        """

        if self.enable:
            return self.comp.geometricGap(self.pair.getCoordinates())

        return None, None

    @profile
    def data(self):
        """Compute data for DiscreteComputation.

        Returns:
            (FieldOnCellsReal): data
        """

        if self.enable:
            return self.comp.contactData(
                self.pair, self.phys_pb.getMaterialField(), self.nb_pairing <= 1
            )

        return None

    @property
    def enable(self):
        """Contact is enable or not

        Returns:
         (bool): True if enable else False
        """

        if self.defi is not None:
            return True

        return False

    def update(self, phys_state):
        """Update contact solver.

        Arguments:
            phys_state (PhysicalSate): physical state
        """

        if self.enable:
            self.pair.updateCoordinates(phys_state.primal_curr)
