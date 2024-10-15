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
Definition of objects for coupled simulations with code_aster.
"""

import os

import numpy
from ple.pyple_coupler import pyple_coupler

from ..Utilities.MedUtils.med_coupler import MEDCoupler
from ..Utilities.mpi_utils import MPI

from .parameters import SchemeParams


LOGDIR = os.environ.get("CA_COUPLING_LOGDIR", "/tmp")


class ExternalCoupling:
    """Coupled simulation with code_aster and other software with PLE.

    Arguments:
        app (str): Application name (default: "code_aster").
        starter (bool): *True* if this instance starts (it sends first),
            *False* if it receives first (default: "False").
        debug (bool): Enable debugging mode (default: "False")".

    Attributes:
        ple (*ple_coupler*): PLE Coupler object.
        ctxt (dict): Execution context where code_aster objects are kept between
            the different stages of the simulation.
        tag (int): tag used for direct MPI communication (`CS_CALCIUM_MPI_TAG`
           ).
        mesh (*medcoupling.MEDCouplingUMesh*): Mesh of the interface.
        fields_in (list): List of exchanged fields (field name, number of
            components, components names).
        fields_out (list): List of exchanged fields (field name, number of
            components, components names).
        params (SchemeParams): Parameters of the coupling scheme.
    """

    def __init__(self, app="code_aster", starter=False, debug=False):
        self.whoami = app
        self.starter = starter
        self.debug = debug
        self.ple = None
        self.ctxt = {}
        self.tag = 0
        self.medcpl = None
        self.mesh = None
        self.fields_in = []
        self.fields_out = []
        self.params = SchemeParams()
        self._other_root = -1

    @property
    def comm(self):
        """Attribute that holds the sub-communicator for the application."""
        return self.ple.my_comm

    def __getstate__(self):
        """Disable pickling"""
        return None

    def sync(self, end_coupling=False):
        """Synchronization with the other application.

        Arguments:
            end_coupling (bool): *True* to stop the coupling loop.

        Returns:
            bool: *True* if the study should terminate.
        """
        return self.ple.sync_coupling_status(end_coupling=end_coupling)

    def log(self, *args, **kwargs):
        """Log messages using the logging function of the coupler object."""
        func = self.ple.log if self.ple else print
        if self.debug:
            prargs = kwargs.copy()
            prargs.pop("verbosity", None)
            prargs["flush"] = True
            print("[CPL]", *args, **prargs)
        return func(*args, **kwargs)

    def _init_metadata(self):
        if self.starter:
            self.send_parameter(0, "NBPDTM", self.params.nb_steps, MPI.INT)
            self.send_parameter(0, "NBSSIT", self.params.nb_iter, MPI.INT)
            self.send_parameter(0, "EPSILO", self.params.epsilon, MPI.DOUBLE)
            self.send_parameter(0, "TTINIT", self.params.init_time, MPI.DOUBLE)
            self.send_parameter(0, "PDTREF", self.params.delta_t, MPI.DOUBLE)
        else:
            self.params.nb_steps = self.recv_parameter(0, "NBPDTM", MPI.INT)
            self.params.nb_iter = self.recv_parameter(0, "NBSSIT", MPI.INT)
            self.params.epsilon = self.recv_parameter(0, "EPSILO", MPI.DOUBLE)
            self.params.init_time = self.recv_parameter(0, "TTINIT", MPI.DOUBLE)
            self.params.delta_t = self.recv_parameter(0, "PDTREF", MPI.DOUBLE)

    def _init_paramedmem(self, with_app):
        self.medcpl = MEDCoupler(logfunc=self.log)

        my_ranks = self.ple.get_app_ranks(self.whoami)
        other_ranks = self.ple.get_app_ranks(with_app)

        # Creating the parallel DEC
        for name, _, _ in self.fields_in + self.fields_out:
            if self.starter:
                self.medcpl.init_paramedmem_coupling(name, ranks1=my_ranks, ranks2=other_ranks)
            else:
                self.medcpl.init_paramedmem_coupling(name, ranks1=other_ranks, ranks2=my_ranks)

        # Define coupling mesh
        self.medcpl.create_coupling_mesh(self.mesh)

        # Define coupled fields
        for name, components, discr in self.fields_in + self.fields_out:
            self.medcpl.define_coupled_field(name, components, discr)

        # Sync the parallel DEC
        # self.medcpl.sync_dec()

    def recv_parameter(self, iteration, name, typ):
        """Receive a parameter (equivalent to `cs_calcium_write_xxx`).

        Arguments:
            iteration (int): Iteration number.
            name (str): Expected parameter name.
            typ (MPI.INT|MPI.DOUBLE): Type of MPI data.

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
            self.ple.base_comm.Recv((data, 128, MPI.CHAR), **args)
            varname = data.decode("utf-8").strip("\x00")
            assert varname == name, f"expecting {name!r}, get {varname!r}"

            meta = numpy.zeros((3,), dtype=numpy.int32)
            self.ple.base_comm.Recv((meta, 3, MPI.INT), **args)
            assert meta[0] == iteration, meta
            assert meta[1] == 1, meta
            assert meta[2] == typ.size

            ctype = numpy.double if typ == MPI.DOUBLE else numpy.int32
            value = numpy.zeros((1,), dtype=ctype)
            self.ple.base_comm.Recv((value, 1, typ), **args)
            data = [varname, value[0]]
            self.log(f"Returns value is {data}", verbosity=2)

        # share the Returns, used as inputs by others
        data = self.comm.bcast(data, root=0)
        self.log(f"receive parameter {name!r} (iteration {iteration}): {data[1]}")
        return data[1]

    def send_parameter(self, iteration, name, value, typ):
        """Send a parameter (equivalent to `cs_calcium_read_xxx`).

        Arguments:
            iteration (int): Iteration number.
            name (str): Parameter name.
            value (int|double): Value of the parameter.
            typ (MPI.INT|MPI.DOUBLE): Type of MPI data.
        """

        self.log(f"send parameter {name!r} (iteration {iteration}): {value}")
        args = dict(dest=self._other_root, tag=self.tag)
        if self.comm.rank == 0:
            bname = (name + "\x00" * (128 - len(name))).encode("utf-8")
            self.ple.base_comm.Send((bname, 128, MPI.CHAR), **args)

            meta = numpy.array([iteration, 1, typ.size], dtype=numpy.int32)
            self.ple.base_comm.Send((meta, 3, MPI.INT), **args)

            ctype = numpy.double if typ == MPI.DOUBLE else numpy.int32
            value = numpy.array(value, dtype=ctype)
            self.ple.base_comm.Send((value, 1, typ), **args)

        self.comm.Barrier()

    def bcast_parameter(self, root, iteration, name, value, typ):
        """Broadcast a parameter between root and receiver.

        Arguments:
            root (bool): root or not ?
            iteration (int): Iteration number.
            name (str): Parameter name.
            value (int|double): Value of the parameter.
            typ (MPI.INT|MPI.DOUBLE): Type of MPI data.

        Returns:
            int|double: Broadcasted value of the parameter.
        """

        val = value
        if root:
            self.send_parameter(iteration, name, val, typ)
        else:
            val = self.recv_parameter(iteration, name, typ)

        return val

    def recv_input_data(self):
        """Receive the inputs from the other code.

        Returns:
            data (list[*ParaFIELD*]): Data used to define the inputs at the next step.
        """
        names = [name for name, _, _ in self.fields_in]
        data = self.medcpl.pmm_recv(names)  # cs_coupler adds `.deepCopy()`
        return data

    def send_output_data(self, outputs):
        """Send the outputs to the other code.

        Arguments:
            outputs (list[*ParaFIELD*]): Result used to define the inputs of the other code.
        """
        self.medcpl.pmm_send(outputs)

    def setup(self, with_app, mesh_interf, input_fields, output_fields, **params):
        """Initialize the coupling.

        Arguments:
            with_app (str): Name of the other application to be coupled with.
            mesh_interf (MEDFileUMesh): Medcoupling mesh of the interface.
            input_fields (list): List of exchanged fields as input.
            output_fields (list): List of exchanged fields as output.
            params (dict): Parameters of the coupling scheme.
        """
        verbosity = 2 if self.debug else 1
        output = "all" if self.debug else "master"
        self.ple = pyple_coupler(verbosity=verbosity, logdir=LOGDIR, output_mode=output)
        self.ple.init_coupling(app_name=self.whoami, app_type="code_aster")

        myranks = self.ple.get_app_ranks(self.whoami)
        self.log(f"allocated ranks for {self.whoami!r}: {myranks}", verbosity=verbosity)
        other_ranks = self.ple.get_app_ranks(with_app)
        self.log(f"allocated ranks for {with_app!r}: {other_ranks}", verbosity=verbosity)
        assert other_ranks, f"Application {with_app!r} not found!"
        self._other_root = other_ranks[0]

        self.log(
            f"{self.whoami!r} coupler created from #{myranks[0]}, "
            f"{with_app!r} root proc is #{1}"
        )
        self.mesh = mesh_interf
        self.fields_in = input_fields
        self.fields_out = output_fields
        self.params.update(params)
        self._init_metadata()
        self.params.check()
        self._init_paramedmem(with_app)

    def finalize(self):
        """Finalize the coupling."""
        self.ple.finalize()

    def testing(self, exec_iteration, data=None):
        """Execute one iteration without coupling.

        Arguments:
            exec_iteration (func): Function that execute one iteration.
            data (dict[*ParaFIELD*]): fields to use in exec_iteration (default: None)

        Returns:
            dict[*ParaFIELD*]: output fields of exec_iteration function.
        """

        if data is None:
            data = {}
            for name, _, _ in self.fields_in:
                data[name] = None

        Returns = exec_iteration(0, 0.0, 1.0, data)

        return Returns

    def run(self, exec_iteration):
        """Execute the coupling loop.

        Arguments:
            exec_iteration (func): Function that execute one iteration.
        """
        # initial sync before the loop
        exit_coupling = self.sync()

        current_time = self.params.init_time
        delta_t = self.params.delta_t
        completed = False
        output_data = None
        istep = 0
        while not completed and not exit_coupling:
            istep += 1

            # ask for timestep
            delta_t = self.bcast_parameter(self.starter, istep, "DTCALC", delta_t, MPI.DOUBLE)

            current_time += delta_t
            self.log("coupling iteration #{0:d}, time: {1:f}".format(istep, current_time))

            # recv data from the other code
            if self.starter and istep == 1:
                input_data = {}
                for name, _, _ in self.fields_in:
                    input_data[name] = None
            else:
                input_data = self.recv_input_data()

            output_data = exec_iteration(istep, current_time, delta_t, input_data)

            # get convergence indicator
            icvast = self.bcast_parameter(self.starter, istep, "ICVAST", 1, MPI.INT)
            assert icvast == 1, f"icvast = {icvast}"

            completed = current_time >= self.params.final_time

            # send data to other code
            if self.starter:
                self.send_output_data(output_data)
                exit_coupling = self.sync(end_coupling=completed)
            else:
                exit_coupling = self.sync(end_coupling=completed)
                self.send_output_data(output_data)

            self.log("end of iteration status: {}".format(exit_coupling))

        self.log(
            "coupling {0} with exit status: {1}".format(
                "completed" if completed else "interrupted", exit_coupling
            )
        )
