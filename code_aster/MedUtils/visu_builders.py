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

from math import isclose
from pathlib import Path
from typing import List, Dict, Tuple
from ..Utilities import medcoupling as mc, logger

import numpy as np


class VisuCutBuilder:
    def __init__(self, mesh_med: mc.MEDFileUMesh):
        self.mesh_med = mesh_med

        self.arrays: Dict[Tuple[str, int, str], np.ndarray] = {}
        self.instants: Dict[Tuple[str, int, str], float] = {}
        self.components: Dict[Tuple[str, str], List[str]] = {}
        self.nb_nodes: int = mesh_med.getNumberOfNodes()

    def __str__(self):
        return f"""
        instants: {self.instants}
        components: {self.components}
        arrays: {self.arrays}
        mesh_med: {self.mesh_med}
        nb_nodes: {self.nb_nodes}
        """

    @classmethod
    def from_aster_mesh(cls, mesh) -> "VisuCutBuilder":
        mesh_med = mesh.createMedCouplingMesh()
        return cls(mesh_med=mesh_med)

    def add_field_from_aster_result_all_timesteps(
        self, aster_result, field_name: str, nodes: list[int]
    ):
        nume_ordres = aster_result.getAccessParameters()["NUME_ORDRE"]
        instants = aster_result.getAccessParameters().get("INST", [0])

        # Create field with linearise fields
        for instant, nume_ordre in zip(instants, nume_ordres, strict=True):
            self.add_field_from_aster_result_1_timestep(
                aster_result=aster_result,
                field_name=field_name,
                nodes=nodes,
                nume_ordre=nume_ordre,
                instant=instant,
            )

    def add_field_from_aster_result_1_timestep(
        self, aster_result, field_name: str, nodes: list[int], nume_ordre: int, instant: float
    ):
        field = aster_result.getField(field_name, nume_ordre).toSimpleFieldOnNodes()
        # extract values for cut group only
        values = field.getValues()[0][nodes]
        components = field.getComponents()
        self.add_field(
            field_name=field_name,
            nodes=nodes,
            values=values,
            components=components,
            nume_ordre=nume_ordre,
            instant=instant,
        )

    def add_field(
        self,
        field_name: str,
        nodes: list[int],
        values: np.ndarray,
        components: List[str],
        nume_ordre: int,
        instant: float,
    ):
        assert len(nodes) == values.shape[0]
        stock_address: Tuple[str, int, str] = (field_name, nume_ordre, mc.ON_NODES)
        already_instanciated_field = stock_address in self.arrays
        if already_instanciated_field:
            assert isclose(self.instants[stock_address], instant)
            assert self.arrays[stock_address].shape == (self.nb_nodes, values.shape[1] + 1)
            assert self.components[(field_name, mc.ON_NODES)] == components
        else:
            self.instants[stock_address] = instant
            self.arrays[stock_address] = np.empty((self.nb_nodes, values.shape[1] + 1))
            self.components[(field_name, mc.ON_NODES)] = components

        current_array = self.arrays[stock_address]
        assert max(nodes) <= current_array.shape[0]
        current_array[nodes, 0] = nodes
        current_array[nodes, 1:] = values

    def write(self, filepath: Path):
        if not self.arrays:
            logger.warn("Cannot write visu output because no data where set")
            return
        field_file = mc.MEDFileFieldMultiTS()

        all_levels = self.mesh_med.getNonEmptyLevels()
        max_level = max(all_levels)
        for (field_name, nume_ordre, location), current_array in self.arrays.items():
            if location != mc.ON_NODES:
                raise NotImplementedError(f"location {location} is not implemented yet")
            instant = self.instants[(field_name, nume_ordre, location)]

            # Medcoupling field
            field_values = mc.DataArrayDouble(current_array[:, 1:].tolist())
            field_values.setInfoOnComponents(self.components[(field_name, location)])
            field_values.setName(field_name)

            medc_node_field = mc.MEDCouplingFieldDouble(location, mc.ONE_TIME)
            medc_node_field.setName(field_name)
            medc_node_field.setArray(field_values)
            medc_node_field.setNature(mc.IntensiveMaximum)
            if location == mc.ON_NODES:
                relative_level = max_level - 0

            medc_node_field.setMesh(self.mesh_med.getMeshAtLevel(relative_level))

            medc_node_field.checkConsistencyLight()
            medfield = mc.MEDFileField1TS()
            medfield.setFieldNoProfileSBT(medc_node_field)
            medfield.setTime(nume_ordre, 0, instant)

            print(f"medfield = {medfield}")
            field_file.pushBackTimeStep(medfield)

        print(f"field_file = {field_file}")

        self.mesh_med.write(str(filepath), 2)
        field_file.write(str(filepath), 0)
