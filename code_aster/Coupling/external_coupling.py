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

import medcoupling as MEDC
import numpy
from ple.pyple_coupler import pyple_coupler

from ..Utilities.MedUtils.med_coupler import MEDCoupler
from ..Utilities.mpi_utils import MPI
from ..CodeCommands.fin import close


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
            in code_saturne).
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

    def init_coupling(self, with_app, meshfile, input_fields, output_fields, **params):
        """Initialize the coupling.

        Arguments:
            with_app (str): Name of the other application to be coupled with.
            meshfile (str): Path to the MED file of the interface.
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
        self.mesh = MEDC.MEDFileUMesh(meshfile)[0]
        self.fields_in = input_fields
        self.fields_out = output_fields
        self.params.update(params)
        self._init_metadata()
        self.params.check()
        self._init_paramedmem(with_app)

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
        if self.starter:
            self.medcpl.init_paramedmem_coupling(ranks1=my_ranks, ranks2=other_ranks)
        else:
            self.medcpl.init_paramedmem_coupling(ranks1=other_ranks, ranks2=my_ranks)

        # Define coupling mesh
        self.medcpl.create_coupling_mesh(self.mesh)

        # Define coupled fields
        for name, components, discr in self.fields_in + self.fields_out:
            self.medcpl.define_coupled_field(name, components, discr)

        # Sync the parallel DEC
        # self.medcpl.sync_dec()

    def recv_parameter(self, iteration, cs_name, typ):
        """Receive a parameter (equivalent to `cs_calcium_write_xxx`).

        Arguments:
            iteration (int): Iteration number.
            cs_name (str): Expected parameter name in code_saturne.
            typ (MPI.INT|MPI.DOUBLE): Type of MPI data.

        Results:
            int|double: Received value of the parameter.
        """
        args = dict(source=self._other_root, tag=self.tag)
        data = None
        if self.comm.rank == 0:
            self.log(
                f"waiting for parameter {cs_name!r} from proc #{self._other_root}...", verbosity=2
            )
            data = bytearray(128)
            self.ple.base_comm.Recv((data, 128, MPI.CHAR), **args)
            varname = data.decode("utf-8").strip("\x00")
            assert varname == cs_name, f"expecting {cs_name!r}, get {varname!r}"

            meta = numpy.zeros((3,), dtype=numpy.int32)
            self.ple.base_comm.Recv((meta, 3, MPI.INT), **args)
            assert meta[0] == iteration, meta
            assert meta[1] == 1, meta
            assert meta[2] == typ.size

            ctype = numpy.double if typ == MPI.DOUBLE else numpy.int32
            value = numpy.zeros((1,), dtype=ctype)
            self.ple.base_comm.Recv((value, 1, typ), **args)
            data = [varname, value[0]]
            self.log(f"results value is {data}", verbosity=2)

        # share the results, used as inputs by others
        data = self.comm.bcast(data, root=0)
        self.log(f"receive parameter {cs_name!r} (iteration {iteration}): {data[1]}")
        return data[1]

    def send_parameter(self, iteration, cs_name, value, typ):
        """Send a parameter (equivalent to `cs_calcium_read_xxx`).

        Arguments:
            iteration (int): Iteration number.
            cs_name (str): Parameter name in code_saturne.
            value (int|double): Value of the parameter.
        """
        self.log(f"send parameter {cs_name!r} (iteration {iteration}): {value}")
        args = dict(dest=self._other_root, tag=self.tag)
        if self.comm.rank == 0:
            bname = (cs_name + "\x00" * (128 - len(cs_name))).encode("utf-8")
            self.ple.base_comm.Send((bname, 128, MPI.CHAR), **args)

            meta = numpy.array([iteration, 1, typ.size], dtype=numpy.int32)
            self.ple.base_comm.Send((meta, 3, MPI.INT), **args)

            ctype = numpy.double if typ == MPI.DOUBLE else numpy.int32
            value = numpy.array(value, dtype=ctype)
            self.ple.base_comm.Send((value, 1, typ), **args)

        self.comm.Barrier()

    def recv_input_data(self):
        """Receive the inputs from the other code.

        Returns:
            data (list[*ParaFIELD*]): Data used to define the inputs at the next step.
        """
        names = [name for name, _, _ in self.fields_in]
        data = self.medcpl.pmm_recv(names)  # cs_coupler adds `.deepCopy()`
        return data

    def send_output_result(self, outputs):
        """Send the results to the other code.

        Arguments:
            outputs (list[*ParaFIELD*]): Result used to define the inputs of the other code.
        """
        self.medcpl.pmm_send(outputs)

    def finalize(self):
        """Finalize the coupling."""
        close()
        self.ple.finalize()

    def testing(self, *args, close=True):
        """Execute the study without coupling, run the setup and one step,
        for debugging.

        Arguments:
            args (misc): Value required to execute one step, as passed to :py:meth:`exec_step`.
        """

        class _Fake:
            my_comm = 0
            log = print

        self.ple = _Fake()
        self.setup()
        self.exec_step(*args)
        if close:
            close()

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
        istep = 0
        while not completed and not exit_coupling:
            istep += 1

            # ask for timestep (propose the last used)
            self.send_parameter(istep, "DTAST", delta_t, MPI.DOUBLE)
            delta_t = self.recv_parameter(istep, "DTCALC", MPI.DOUBLE)

            current_time += delta_t
            self.log("coupling iteration #{0:d}, time: {1:f}".format(istep, current_time))

            # recv results from the other code
            data = self.recv_input_data()
            assert len(data) == 1, data

            results = exec_iteration(istep, current_time, delta_t, data)

            # get convergence indicator
            icvast = self.recv_parameter(istep, "ICVAST", MPI.INT)
            assert icvast == 1

            # send results to other code
            self.send_output_result(results)

            completed = current_time >= self.params.final_time
            exit_coupling = self.sync(end_coupling=completed)
            self.log("end of iteration status: {}".format(exit_coupling))

        self.log(
            "coupling {0} with exit status: {1}".format(
                "completed" if completed else "interrupted", exit_coupling
            )
        )
