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

from ...Messages import UTMESS
from ...Utilities import no_new_attributes, profile, force_list


def get_index(inst, linst, prec, crit):

    assert crit in ("RELATIF", "ABSOLU")

    if crit == "RELATIF":
        min_v = inst*(1-prec)
        max_v = inst*(1+prec)
    else:
        min_v = inst-prec
        max_v = inst+prec

    min_v, max_v = sorted((min_v, max_v))
    test = tuple((i >= min_v and i <= max_v) for i in linst)
    return tuple(i for i, v in enumerate(test) if v is True)


def exists_unique(*args, **kwargs):
    return len(get_index(*args, **kwargs)) == 1


class StorageManager:
    """Object that manages the storing of fields in the Result object.

    Arguments:
        result (~code_aster.Objects.Result): Result container.
    """

    class Slot:
        """Container that holds objects to be saved"""

        __slots__ = ("index", "time", "model", "material_field",
                     "elem_char", "load", "fields", "param")

    result = None
    buffer = None
    excl_fields = set()
    crit = prec = None
    list_time = pas_arch = None
    curr_index = init_index = 0
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, result, mcf=None, **kwargs):
        """Create the storage manager object from the ARCHIVE factor keyword.

        Arguments:
            result (Result): result where store fields
            mcf (list|tuple|dict): Convenient option to pass `(kwargs,)`, the value
            returned for a factor keyword.
            kwargs (dict): Valid ARCHIVE keywords (syntax checked, with defaults).

        """
        self.result = result
        self.buffer = []

        if mcf:
            if isinstance(mcf, (list, tuple)):
                mcf = mcf[0]
            if isinstance(mcf, dict):
                kwargs = mcf

        excl_fields = kwargs.get("CHAM_EXCLU")
        if excl_fields is not None:
            for field in excl_fields:
                self.excl_fields.add(field)

        if "CRITERE" in kwargs:
            self.crit = kwargs["CRITERE"]
        if "PRECISION" in kwargs:
            self.prec = kwargs["PRECISION"]

        if "INST" in kwargs:
            self.list_time = force_list(kwargs["INST"])
        elif "LIST_INST" in kwargs:
            self.list_time = kwargs["LIST_INST"].getValues()
        elif "PAS_ARCH" in kwargs:
            self.pas_arch = kwargs["PAS_ARCH"]
        else:
            self.pas_arch = 1

        if self.list_time is not None:
            for time in self.list_time:
                assert(exists_unique(time, self.list_time, self.prec, self.crit))

    def setInitialIndex(self, index):
        """Set initial index.

        Arguments:
            index (int): initial index.
        """
        self.curr_index = index
        self.init_index = index

        if self.result.getNumberOfRanks() > 0:
            self.result.clear(self.init_index)

    def hasToBeStored(self, time):
        """To known if this time step has to be store

        Arguments:
            time (float): time step.

        Returns:
            bool: True if the time step has to be store else False
        """

        if self.pas_arch is not None:
            return (self.curr_index - self.init_index) % self.pas_arch == 0

        if self.list_time is not None:
            indexes = get_index(time, self.list_time, self.prec, self.crit)
            assert len(indexes) <= 1
            return len(indexes) == 1

        return True

    def completed(self):
        """Register the current step as completed successfully."""
        self.curr_index += 1

    def getResult(self):
        """Returns the Result container.

        Returns:
            ~code_aster.Objects.Result: Result container.
        """
        return self.result

    @profile
    def storeState(self, time, phys_pb, phys_state, param=None):
        """Store a new state.

        Arguments:
            time (float): current (pseudo)-time.
            phys_pb (PhysicalProblem): Physical problem
            phys_state (PhysicalState): Physical state
        """
        slot = StorageManager.Slot()
        slot.index = self.curr_index
        slot.time = time
        slot.param = param
        slot.model = phys_pb.getModel()
        slot.material_field = phys_pb.getMaterialField()
        slot.elem_char = phys_pb.getElementaryCharacteristics()
        slot.load = phys_pb.getListOfLoads()
        slot.fields = phys_state.as_dict()
        behav = phys_pb.getBehaviourProperty()
        if behav is not None:
            slot.fields["COMPORTEMENT"] = behav.getBehaviourField()
        self.buffer.append(slot)

        self.store()

    @profile
    def storeField(self, field, field_type, time=0.0):
        """Store a new field.

        Arguments:
            field (FieldOn***): field to store
            field_type (str) : type of the field as DEPL, SIEF_ELGA...
        """
        if field is not None and field_type not in self.excl_fields:
            self.result.setField(field, field_type, self.curr_index)
            UTMESS("I", "ARCHIVAGE_6", valk=field_type,
                   valr=time, vali=self.curr_index)

    @profile
    def store(self):
        """Build result with all ranks in buffer."""

        new_size = self.result.getNumberOfRanks() + len(self.buffer)
        self.result.resize(new_size)
        for slot in self.buffer:
            curr_index = slot.index
            if slot.time is not None:
                self.result.setTimeValue(slot.time, curr_index)
            if slot.param is not None:
                for param, value in slot.param.items():
                    self.result.setParameterValue(param, value, curr_index)
            if slot.model:
                self.result.setModel(slot.model, curr_index)
            if slot.material_field:
                self.result.setMaterialField(slot.material_field, curr_index)
            if slot.elem_char:
                self.result.setElementaryCharacteristics(
                    slot.elem_char, curr_index)
            if slot.load:
                self.result.setListOfLoads(slot.load, curr_index)
            if slot.fields:
                for field_type, field in slot.fields.items():
                    self.storeField(field, field_type, slot.time)

        self.buffer = []
