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

# person_in_charge: nicolas.pignet@edf.fr

"""
This module gives common utilities for MPI copuling.

Need only mpi4py package.
"""

import numpy as np

from ..Utilities.mpi_utils import MPI
from ..Utilities.logger import logger


class MPICoupling:
    """
    This class MPICoupling is an encapsulation of PLE communication.

    The same API than mpi4py is used.

    Arguments:
        comm (Comm): communicator between applications.
        sub_comm (Comm): sub-communicator for the application.
        other_root (int): root of the other code.
        log (logger): logger (default: None).
    """

    DOUBLE = MPI.DOUBLE
    INT = MPI.INT
    CHAR = MPI.CHAR
    BOOL = MPI.BOOL

    LAND = MPI.LAND

    class COMMPLE:
        def __init__(self, sub_comm, comm, other_root, log):
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
                typ (INT|DOUBLE): Type of MPI data.

            Returns:
                int|double: Received value of the parameter.
            """

            args = dict(source=self._other_root, tag=self.tag)
            data = None
            if self.comm.rank == 0:
                self.log(
                    f"waiting for parameter {name!r} from proc #{self._other_root}...", verbosity=2
                )
                data = bytearray(128)
                self.sub_comm.Recv((data, 128, MPICoupling.CHAR), **args)
                varname = data.decode("utf-8").strip("\x00")
                assert varname == name, f"expecting {name!r}, get {varname!r}"

                meta = np.zeros((3,), dtype=np.int32)
                self.sub_comm.Recv((meta, 3, MPICoupling.INT), **args)
                assert meta[0] == iteration, meta
                assert meta[1] == 1, meta
                assert meta[2] == typ.size

                ctype = np.double if typ == MPICoupling.DOUBLE else np.int32
                value = np.zeros((1,), dtype=ctype)
                self.sub_comm.Recv((value, 1, typ), **args)
                data = [varname, value[0]]
                self.log(f"Returns value is {data}", verbosity=2)

            # share the Returns, used as inputs by others
            data = self.comm.bcast(data, root=0)
            self.log(f"receive parameter {name!r} (iteration {iteration}): {data[1]}")
            return data[1]

        def send(self, iteration, name, value, typ):
            """Send a parameter (equivalent to `cs_calcium_read_xxx`).

            Arguments:
                iteration (int): Iteration number.
                name (str): Parameter name.
                value (int|double): Value of the parameter.
                typ (INT|DOUBLE): Type of MPI data.
            """

            self.log(f"send parameter {name!r} (iteration {iteration}): {value}")
            args = dict(dest=self._other_root, tag=self.tag)
            if self.comm.rank == 0:
                bname = (name + "\x00" * (128 - len(name))).encode("utf-8")
                self.sub_comm.Send((bname, 128, MPICoupling.CHAR), **args)

                meta = np.array([iteration, 1, typ.size], dtype=np.int32)
                self.sub_comm.Send((meta, 3, MPICoupling.INT), **args)

                ctype = np.double if typ == MPICoupling.DOUBLE else np.int32
                value = np.array(value, dtype=ctype)
                self.sub_comm.Send((value, 1, typ), **args)

            self.comm.Barrier()

        def bcast(self, root, iteration, name, value, typ):
            """Broadcast a parameter between root and receiver.

            Arguments:
                root (bool): root or not ?
                iteration (int): Iteration number.
                name (str): Parameter name.
                value (int|double): Value of the parameter.
                typ (MPI.INT|MPI.DOUBLE): Type of MPI data.

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
                typ (MPI.BOOL): Type of MPI data.

            Returns:
                (bool): broadcasted value.
            """

            if isinstance(type, self.BOOL):
                value = int(value)
                typ = self.INT

                b0 = bool(self.bcast(True, iteration, name, value, typ))
                b1 = bool(self.bcast(False, iteration, name, value, typ))

                if op == self.LAND:
                    if b0 and b1:
                        return True
                    else:
                        return False
                else:
                    raise NotImplemented()
            else:
                raise NotImplemented()

    def __init__(self, sub_comm, comm, other_root, logfunc=None):
        self.log = logfunc if logfunc else logger
        self.SUB_COMM = self.COMMPLE(sub_comm, comm, other_root, self.log)
        self.COMM = comm
