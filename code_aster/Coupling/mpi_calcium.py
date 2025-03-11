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

# person_in_charge: nicolas.pignet@edf.fr

"""
This module gives common utilities for MPI communications.

Need only mpi4py package.
"""

import numpy as np

from ..Utilities.logger import logger
from ..Utilities.mpi_utils import MPI


class MPICalcium:
    """
    This class MPI is an encapsulation of CALCIUM communication.

    The same API than mpi4py is used.

    Arguments:
        comm (*MPI.Comm*): communicator between applications.
        sub_comm (*MPI.Comm*): sub-communicator for the application.
        other_root (int): root of the other application.
        log (logger): logger (default: None).
    """

    DOUBLE = MPI.DOUBLE
    INT = MPI.INT
    CHAR = MPI.CHAR
    BOOL = MPI.BOOL

    LAND = MPI.LAND

    class CommCalcium:
        def __init__(self, comm, sub_comm, other_root, log):
            self.comm = comm
            self.sub_comm = sub_comm
            self.tag = 0
            self._other_root = other_root
            self.log = log

        def recv(self, iteration, name, typ):
            """Receive a parameter (equivalent to `cs_calcium_write_xxx`).

            Arguments:
                iteration (int): Iteration number.
                name (str): Expected parameter name.
                typ (:py:class:`~MPICalcium.INT`|:py:class:`~MPICalcium.DOUBLE`): Type of MPI data.

            Returns:
                int|double: Received value of the parameter.
            """

            args = dict(source=self._other_root, tag=self.tag)
            data = None
            if self.sub_comm.rank == 0:
                self.log(
                    f"waiting for parameter {name!r} from proc #{self._other_root}...", verbosity=2
                )
                data = bytearray(128)
                self.comm.Recv((data, 128, MPICalcium.CHAR), **args)
                varname = data.decode("utf-8").strip("\x00")
                assert varname == name, f"expecting {name!r}, get {varname!r}"

                meta = np.zeros((3,), dtype=np.int32)
                self.comm.Recv((meta, 3, MPICalcium.INT), **args)
                assert meta[0] == iteration, meta
                assert meta[1] == 1, meta
                assert meta[2] == typ.size

                ctype = np.double if typ == MPICalcium.DOUBLE else np.int32
                value = np.zeros((1,), dtype=ctype)
                self.comm.Recv((value, 1, typ), **args)
                data = [varname, value[0]]
                self.log(f"Returns value is {data}", verbosity=2)

            # share the Returns, used as inputs by others
            data = self.sub_comm.bcast(data, root=0)
            self.log(f"receive parameter {name!r} (iteration {iteration}): {data[1]}")
            return data[1]

        def send(self, iteration, name, value, typ):
            """Send a parameter (equivalent to `cs_calcium_read_xxx`).

            Arguments:
                iteration (int): Iteration number.
                name (str): Parameter name.
                value (int|double): Value of the parameter.
                typ (:py:class:`~MPICalcium.INT`|:py:class:`~MPICalcium.DOUBLE`): Type of MPI data.
            """

            self.log(f"send parameter {name!r} (iteration {iteration}): {value}")
            args = dict(dest=self._other_root, tag=self.tag)
            if self.sub_comm.rank == 0:
                bname = (name + "\x00" * (128 - len(name))).encode("utf-8")
                self.comm.Send((bname, 128, MPICalcium.CHAR), **args)

                meta = np.array([iteration, 1, typ.size], dtype=np.int32)
                self.comm.Send((meta, 3, MPICalcium.INT), **args)

                ctype = np.double if typ == MPICalcium.DOUBLE else np.int32
                value = np.array(value, dtype=ctype)
                self.comm.Send((value, 1, typ), **args)

            self.sub_comm.Barrier()

        def bcast(self, root, iteration, name, value, typ):
            """Broadcast a parameter between root and receiver.

            Arguments:
                root (bool): root or not ?
                iteration (int): Iteration number.
                name (str): Parameter name.
                value (int|double): Value of the parameter.
                typ (:py:class:`~MPICalcium.INT`|:py:class:`~MPICalcium.DOUBLE`): Type of MPI data.

            Returns:
                (int|double): broadcasted value.
            """

            if root:
                self.send(iteration, name, value, typ)
            else:
                value = self.recv(iteration, name, typ)

            return value

        def allreduce(self, iteration, name, value, typ, op):
            """Allreduce a parameter between root and receiver.

            Arguments:
                root (bool): root or not ?
                iteration (int): Iteration number.
                name (str): Parameter name.
                value (int|double): Value of the parameter.
                typ (:py:class:`~MPICalcium.BOOL`): Type of MPI data.

            Returns:
                (bool): broadcasted value.
            """

            if typ == MPICalcium.BOOL:
                value = int(value)
                typ = MPICalcium.INT

                b0 = bool(self.bcast(True, iteration, name, value, typ))
                b1 = bool(self.bcast(False, iteration, name, value, typ))

                if op == MPICalcium.LAND:
                    if b0 and b1:
                        return True
                    return False
            raise NotImplementedError()

    def __init__(self, comm, sub_comm, other_root, logfunc=None):
        self.log = logfunc if logfunc else logger
        self.cpl_comm = self.CommCalcium(comm, sub_comm, other_root, self.log)

    @property
    def ASTER_COMM_WORLD(self):
        """*MPI.Comm*: (sub-)communicator for the application."""
        return self.cpl_comm.sub_comm

    @property
    def COMM_WORLD(self):
        """*MPI.Comm*: global communicator for all applications."""
        return self.cpl_comm.comm

    @property
    def COUPLING_COMM_WORLD(self):
        """*MPI.Comm*: communicator between applications."""
        return self.cpl_comm
