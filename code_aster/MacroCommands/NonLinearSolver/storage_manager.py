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

from ...Utilities import no_new_attributes, profile


class StorageManager:
    """Object that manages the storing of fields in the Result object.

    Arguments:
        result (~code_aster.Objects.Result): Result container.
    """

    class Slot:
        """Container that holds objects to be saved"""

        __slots__ = ("rank", "time", "model", "material_field", "elem_char", "load", "fields")

    result = None
    buffer = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, result):
        self.result = result
        self.buffer = []

    def getResult(self):
        """Returns the Result container.

        Returns:
            ~code_aster.Objects.Result: Result container.
        """
        return self.result

    @profile
    def storeState(self, rank, time, phys_pb, phys_state):
        """Store a new state.

        Arguments:
            rank (int): rank where to save fields
            time (float): current (pseudo)-time.
            phys_pb (PhysicalProblem): Physical problem
            phys_state (PhysicalState): Physical state
        """
        slot = StorageManager.Slot()
        slot.rank = rank
        slot.time = time
        slot.model = phys_pb.getModel()
        slot.material_field = phys_pb.getMaterialField()
        slot.elem_char = phys_pb.getElementaryCharacteristics()
        slot.load = phys_pb.getListOfLoads()
        slot.fields = phys_state.as_dict()
        slot.fields["COMPORTEMENT"] = phys_pb.getBehaviourProperty().getBehaviourField()
        self.buffer.append(slot)

    @profile
    def store(self):
        """Build result with all ranks in buffer."""

        new_size = self.result.getNumberOfRanks() + len(self.buffer)
        self.result.resize(new_size)
        for slot in self.buffer:
            rank_curr = slot.rank
            if slot.time is not None:
                self.result.setTimeValue(slot.time, rank_curr)
            if slot.model:
                self.result.setModel(slot.model, rank_curr)
            if slot.material_field:
                self.result.setMaterialField(slot.material_field, rank_curr)
            if slot.elem_char:
                self.result.setElementaryCharacteristics(slot.elem_char, rank_curr)
            if slot.load:
                self.result.setListOfLoads(slot.load, rank_curr)
            if slot.fields:
                for field_type, field in slot.fields.items():
                    if field is not None:
                        self.result.setField(field, field_type, rank_curr)

        self.buffer = []
