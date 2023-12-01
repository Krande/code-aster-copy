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

from ...Messages import UTMESS
from ...Utilities import SearchList, force_list, logger, no_new_attributes, profile
from ..Basics import SolverFeature
from ..Basics import SolverOptions as SOP


class StorageManager(SolverFeature):
    """Object that manages the storing of fields in the Result object.

    Arguments:
        result (~code_aster.Objects.Result): Result container.
    """

    provide = SOP.Storage

    class Slot:
        """Container that holds objects to be saved"""

        __slots__ = (
            "index",
            "store_index",
            "time",
            "model",
            "material_field",
            "elem_char",
            "load",
            "fields",
            "param",
        )

    _result = _buffer = _excl_fields = None
    _eps = _relative = None
    _timelist = _step = None
    _stor_idx = _last_idx = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, result, mcf=None, **kwargs):
        """Create the storage manager object from the ARCHIVAGE factor keyword.

        Arguments:
            result (Result): result where store fields
            mcf (list|tuple|dict): Convenient option to pass `(kwargs,)`, the value
                returned for a factor keyword.
            kwargs (dict): Valid ARCHIVAGE keywords (syntax checked, with defaults).

        """
        super().__init__()
        self._result = result
        if self._result.getNumberOfIndexes() == 0:
            self._result.resize(10)
        self._buffer = []

        if isinstance(mcf, (list, tuple)):
            mcf = mcf[0]
        if isinstance(mcf, dict):
            kwargs = mcf

        self._excl_fields = set()
        excl_fields = kwargs.get("CHAM_EXCLU")
        if excl_fields is not None:
            for field in excl_fields:
                self._excl_fields.add(field)

        times = None
        if "INST" in kwargs:
            times = force_list(kwargs["INST"])
        elif "LIST_INST" in kwargs:
            times = kwargs["LIST_INST"].getValues()
        elif "PAS_ARCH" in kwargs:
            self._step = kwargs["PAS_ARCH"]
        else:
            self._step = 1

        self._eps = 1.0e-6
        self._relative = True
        if times is not None:
            self._eps = kwargs["PRECISION"]
            self._relative = kwargs["CRITERE"] == "RELATIF"
            self._timelist = SearchList(times, self._eps, kwargs["CRITERE"])
            assert all(self._timelist.unique(t) for t in times)

        self._last_idx = -999
        self._stor_idx = 0

    def setFirstStorageIndex(self, index, first_exists=False):
        """Set initial index (index where the first step should be stored).

        Arguments:
            index (int): index to be used for the next storage.
            first_exist (bool): *True* if the first step already exists,
                *False* otherwise.
        """
        self._stor_idx = index
        if first_exists:
            self._last_idx = 0

        if self._result.getNumberOfIndexes() > 0:
            self._result.clear(self._stor_idx)

    def _to_be_stored(self, idx, time):
        """To known if this time step has to be store.

        Arguments:
            idx (int): index of the time (restarts at 0 at each new operator).
            time (float): time step.

        Returns:
            bool: *True* if the time step has to be store else *False*.
        """
        # always store first index
        if idx == 0:
            return True
        if self._step is not None:
            return idx % self._step == 0
        if self._timelist is not None:
            # FIXME use _eps/_relative, see Utilities.cmp
            return time in self._timelist
        return True

    def _was_stored(self, idx):
        """Tell if a state was already stored for this index.

        Note: It actually checks that `idx` is greater than the greatest index
        already stored.

        Arguments:
            idx (int): index of the time (restarts at 0 at each new operator).

        Returns:
            bool: *True* if the time was already stored, *False* otherwise.
        """
        return idx <= self._last_idx

    def getResult(self):
        """Returns the Result container.

        Returns:
            ~code_aster.Objects.Result: Result container.
        """
        return self._result

    def storeParam(self, idx, **kwargs):
        """Store parameters like time, model...

        Arguments:
            idx (int): index of the current (pseudo-)time
                (restarts at 0 at each new operator).
            kwargs: named parameters
        """
        # FIXME use _to_be_stored
        self._result.resize(self._result.getNumberOfIndexes() + 10)
        if "model" in kwargs:
            self._result.setModel(kwargs["model"], idx)
        if "materialField" in kwargs:
            self._result.setMaterialField(kwargs["materialField"], idx)
        if "listOfLoads" in kwargs:
            self._result.setListOfLoads(kwargs["listOfLoads"], idx)
        if "time" in kwargs:
            self._result.setTime(kwargs["time"], idx)

    @profile
    def storeState(self, idx, time, phys_pb, phys_state, param=None, force=True):
        """Store a new state.

        Arguments:
            idx (int): index of the current (pseudo-)time
                (restarts at 0 at each new operator).
            time (float): current (pseudo-)time.
            phys_pb (PhysicalProblem): Physical problem.
            phys_state (PhysicalState): Physical state.
            param (dict, optional): Dict of parameters to be stored.
            force (bool): force the storage even if storing-policy is not verified.
        """
        logger.info(
            "STORE: calc=%d, time=%f, last=%d, stor=%d", idx, time, self._last_idx, self._stor_idx
        )
        if not force and not self._to_be_stored(idx, time):
            return
        if self._was_stored(idx):
            return
        slot = StorageManager.Slot()
        slot.index = idx
        slot.store_index = self._stor_idx
        slot.time = time
        slot.param = param
        slot.model = phys_pb.getModel()
        slot.material_field = phys_pb.getMaterialField()
        slot.elem_char = phys_pb.getElementaryCharacteristics()
        slot.load = phys_pb.getListOfLoads()
        slot.fields = phys_state.as_dict()
        behav = phys_pb.getBehaviourProperty()
        if behav is not None:
            if phys_pb.isThermal():
                compor = "COMPORTHER"
            else:
                compor = "COMPORTEMENT"
            slot.fields[compor] = behav.getBehaviourField()
        self._buffer.append(slot)
        self._store()
        self.completed(idx)

    @profile
    def storeField(self, idx, field, field_type, time=0.0, force=False):
        """Store a new field.

        Arguments:
            idx (int): index of the current (pseudo-)time
                (restarts at 0 at each new operator).
            field (FieldOn***): field to store
            field_type (str): type of the field as DEPL, SIEF_ELGA...
            time (float, optional): current (pseudo-)time.
            force (bool): force the storage even if storing-policy is not verified.
        """
        logger.info(
            "STORE: calc=%d, time=%f, field=%s, last=%d, stor=%d",
            idx,
            time,
            field_type,
            self._last_idx,
            self._stor_idx,
        )
        # if not force and not self._to_be_stored(idx, time):
        #    return
        # if self._was_stored(idx):
        #    return
        self._store_field(time, field, field_type)

    def completed(self, idx):
        """Register the current step as completed successfully.

        Arguments:
            idx (int): index of the current (pseudo-)time
                (restarts at 0 at each new operator).
        """
        if idx > self._last_idx:
            self._stor_idx += 1
        self._last_idx = idx

    @profile
    def _store(self):
        """Build result with all indexes in buffer."""
        new_size = self._result.getNumberOfIndexes() + len(self._buffer)
        self._result.resize(new_size)
        for slot in self._buffer:
            idx = slot.store_index
            if slot.time is not None:
                self._result.setTime(slot.time, idx)
            if slot.param is not None:
                for param, value in slot.param.items():
                    self._result.setParameterValue(param, value, idx)
            if slot.model:
                self._result.setModel(slot.model, idx)
            if slot.material_field:
                self._result.setMaterialField(slot.material_field, idx)
            if slot.elem_char:
                self._result.setElementaryCharacteristics(slot.elem_char, idx)
            if slot.load:
                self._result.setListOfLoads(slot.load, idx)
            if slot.fields:
                for field_type, field in slot.fields.items():
                    self._store_field(slot.time, field, field_type)
        self._buffer = []

    @profile
    def _store_field(self, time, field, field_type):
        """Store a field - internal function."""
        if field is None or field_type in self._excl_fields:
            return
        self._result.setField(field, field_type, self._stor_idx)
        UTMESS("I", "ARCHIVAGE_6", valk=field_type, valr=time, vali=self._stor_idx)
