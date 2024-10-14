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

"""
Definition of a convenient object to synchronize MEDCoupling fields.
"""

import medcoupling as MEDC
import ParaMEDMEM as PMM


class CoupledField:
    """Define the properties of an exchanged field.

    Attributes:
        name (str): Field name.
        components (list[str]): Names of the components.
        support (int): Support: ON_NODES or ON_CELLS.
        size (int): Number of values to be exchanged on the interface.
        para_field (*ParaFIELD*): *ParaMEDMEM.ParaFIELD* object.
    """

    def __init__(self, name, components, support, size, para_field):
        self.name = name
        self.components = components
        assert support in (MEDC.ON_NODES, MEDC.ON_CELLS), support
        self.support = support
        self.size = size
        self.para_field = para_field


class ExtendedDEC:
    """Object that represents a DEC and the necessary properties.

    Arguments:
        mcdec (*PMM.InterpKernelDEC*): ParaMEDMEM InterpKernelDEC object.
    """

    def __init__(self, mcdec):
        self.dec = mcdec
        self.mesh = None
        self.interf_size = 0
        self._synced = False

    def isInSourceSide(self):
        """Wrapper on DEC function."""
        return self.dec.isInSourceSide()

    def getSourceGrp(self):
        """Wrapper on DEC function."""
        return self.dec.getSourceGrp()

    def getTargetGrp(self):
        """Wrapper on DEC function."""
        return self.dec.getTargetGrp()

    def attachLocalField(self, field):
        """Wrapper on DEC function."""
        return self.dec.attachLocalField(field)

    def setMethod(self, method):
        """Wrapper on DEC function."""
        return self.dec.setMethod(method)

    @property
    def synced(self):
        """bool: Tell if the DEC has already been synced."""
        return self._synced

    def synchronize(self):
        """Wrapper on DEC function."""
        if self._synced:
            return
        self._synced = True
        return self.dec.synchronize()

    def sendData(self):
        """Wrapper on DEC function."""
        return self.dec.sendData()

    def recvData(self):
        """Wrapper on DEC function."""
        return self.dec.recvData()


class MEDCoupler:
    """Class handling the MEDCoupling related calls."""

    def __init__(self, logfunc=None):
        self.dec = {MEDC.ON_NODES: {}, MEDC.ON_CELLS: {}}
        self.xdec = {MEDC.ON_NODES: {}, MEDC.ON_CELLS: {}}
        self.interf_mesh = None
        self.exch_fields = {}
        # self.comptopo = None

        self.log = logfunc if logfunc else print

    def init_paramedmem_coupling(self, field_name, ranks1, ranks2):
        """Start ParaMEDMEM coupling and DEC.

        Arguments:
            field_name (str): Field name.
            ranks1 (list[int]): List of ranks allocated to the first application.
            ranks2 (list[int]): List of ranks allocated to the second application.
        """
        self.log("initializing ParaMEDMEM InterpKernelDEC")
        # Creating the send and recv DEC's
        for sup in (MEDC.ON_CELLS, MEDC.ON_NODES):
            sup_name = "nodes" if sup == MEDC.ON_NODES else "cells"
            key = f"{ranks1},{ranks2}"
            if key not in self.dec[sup]:
                self.log(
                    f"creating InterpKernelDEC on {sup_name} with "
                    f"ranks1={ranks1} and ranks2={ranks2}",
                    verbosity=2,
                )
                self.dec[sup][key] = PMM.InterpKernelDEC(ranks1, ranks2)
            self.xdec[sup][field_name] = ExtendedDEC(self.dec[sup][key])
            self.xdec[sup][field_name].setMethod("P0" if sup == MEDC.ON_CELLS else "P1")

    def create_coupling_mesh(self, mesh):
        """Create the ParaMEDMEM mesh, support of coupling.

        Arguments:
            mesh (*medcoupling.MEDCouplingUMesh*): The coupling mesh.
        """
        self.log("creating coupling mesh in memory", verbosity=2)
        self.interf_mesh = mesh
        for sup in (MEDC.ON_CELLS, MEDC.ON_NODES):
            for name in self.xdec[sup]:
                dec = self.xdec[sup][name]
                if dec.isInSourceSide():
                    group = dec.getSourceGrp()
                else:
                    group = dec.getTargetGrp()
                dec.mesh = PMM.ParaMESH(mesh, group, "couplingMesh")

        # self.comptopo = PMM.ComponentTopology()

    def get_coupled_field(self, name, silent=False):
        """Return a coupled field by name.

        Arguments:
            name (str): Field name.

        Returns:
            *Field*: Exchanged field or *None* if not found.
        """
        found = self.exch_fields.get(name)
        if found:
            return found
        if not silent and not found:
            msg = f"Field {name} was not defined beforehand!"
            self.log(msg)
            raise KeyError(msg)
        return found

    def define_coupled_field(self, field_name, components, field_type):
        """Add a coupled field.

        Arguments:
            field_name (str): Field name.
            components (list[str]): Components of the field.
            field_type (str): On "NODES" or "CELLS".
        """
        if not self.get_coupled_field(field_name, silent=True):
            assert field_type in ("NODES", "CELLS")
            sup = MEDC.ON_CELLS if field_type == "CELLS" else MEDC.ON_NODES
            dec = self.xdec[sup][field_name]
            if field_type == "CELLS":
                nature = MEDC.IntensiveConservation
                nb_tuples = self.interf_mesh.getNumberOfCells()
            else:
                nature = MEDC.IntensiveMaximum
                nb_tuples = self.interf_mesh.getNumberOfNodes()

            topo = PMM.ComponentTopology(len(components))
            pfield = PMM.ParaFIELD(sup, MEDC.ONE_TIME, dec.mesh, topo)
            pfield.getField().setName(field_name)
            pfield.getField().setNature(nature)
            pfield.getField().getArray().setInfoOnComponents(components)

            if not dec.interf_size:
                dec.interf_size = nb_tuples * len(components)
            assert (
                nb_tuples * len(components) == dec.interf_size
            ), "inconsistent number of components: %i vs. %i" % (
                nb_tuples * len(components),
                dec.interf_size,
            )

            field = CoupledField(field_name, components, sup, dec.interf_size, pfield)
            self.exch_fields[field_name] = field

    def sync_dec(self):
        """Sync the parallel DEC objects."""
        for sup in (MEDC.ON_CELLS, MEDC.ON_NODES):
            for name in self.xdec[sup]:
                dec = self.xdec[sup][name]
                if dec.synced:
                    continue
                cfield = [field for field in self.exch_fields.values() if field.support == sup]
                if len(cfield) > 0:
                    cfield = cfield[0]
                    pfield = cfield.para_field.getField()
                    dec.attachLocalField(pfield)
                    sup_name = "nodes" if sup == MEDC.ON_NODES else "cells"
                    self.log(
                        f"synchronize DEC #{sup}-{name} for field {name} on {sup_name}, "
                        f"interface size: {dec.interf_size}",
                        verbosity=2,
                    )
                    dec.synchronize()

    def pmm_send(self, fields):
        """Send a field to the partner code with ParaMEDMEM.

        Arguments:
            fields (dict[*ParaFIELD*]): Fields to send.
        """
        if not fields:
            return
        support = None
        for field_name in fields:
            field = fields[field_name]
            exchanged = self.get_coupled_field(field_name)
            if support is None:
                support = exchanged.support
            assert support == support, "all fields must be on the same support type."
            sup_name = "nodes" if support == MEDC.ON_NODES else "cells"
            self.log(f"sending field {field_name!r} on {sup_name}...")
            pfield = exchanged.para_field.getField()
            array = field.getArray()
            pfield.setArray(array)
            pfield.setName(field_name)
            self.log(repr(pfield), verbosity=2)
            # self.log(array, verbosity=2)
            dec = self.xdec[support][field_name]
            dec.attachLocalField(pfield)
            dec.synchronize()
            self.log(f"sendData...")
            dec.sendData()
        self.log("pmm_send: done", verbosity=2)

    def pmm_recv(self, fields_names):
        """Receive a field from the partner code with ParaMEDMEM.

        Arguments:
            fields_names (list[str]): Fields names.

        Returns:
            dict[*ParaFIELD*]: Received fields.
        """
        fields = {}
        if not fields_names:
            return fields
        support = None
        for field_name in fields_names:
            exchanged = self.get_coupled_field(field_name)
            if support is None:
                support = exchanged.support
            assert support == support, "all fields must be on the same support type."
            sup_name = "nodes" if support == MEDC.ON_NODES else "cells"
            self.log(f"waiting for field {field_name!r} on {sup_name}...")
            pfield = exchanged.para_field.getField()
            pfield.setValues([0.0] * exchanged.size)
            array = pfield.getArray()
            array.rearrange(len(exchanged.components))
            array.setInfoOnComponents(exchanged.components)
            dec = self.xdec[support][field_name]
            dec.attachLocalField(pfield)
            fields[field_name] = pfield
            dec.synchronize()
            self.log(f"recvData...", verbosity=2)
            dec.recvData()
            self.log(repr(pfield), verbosity=2)
            # self.log(array, verbosity=2)
            self.log(f"pmm_recv: done", verbosity=2)
        return fields
