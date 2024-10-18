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

from ..Utilities.logger import logger
from .med_coupler import MEDCoupler
from .mpi_calcium import MPI_CALCIUM
from .parameters import SchemeParams
from .ple_utils import PLE


class ExternalCoupling:
    """Coupled simulation with code_aster and other software with PLE.

    Arguments:
        app (str): Application name (default: "code_aster").
        starter (bool): *True* if this instance starts (it sends first),
            *False* if it receives first (default: "False").
        debug (bool): Enable debugging mode (default: "False")".

    Attributes:
        ple (*ple_coupler*): PLE Coupler object.
        mesh (*medcoupling.MEDCouplingUMesh*): Mesh of the interface.
        fields_in (list): List of exchanged fields (field name, number of
            components, components names).
        fields_out (list): List of exchanged fields (field name, number of
            components, components names).
        params (SchemeParams): Parameters of the coupling scheme.
        medcpl (MEDCoupler): Coupler to exchange and interpolate fields.
        MPI (MPI_CALCIUM): Wrapping of mpi4py for communication inside
            application and between applications.
    """

    def __init__(self, app="code_aster", starter=False, debug=False):
        self.whoami = app
        self.other_app = None
        self.starter = starter
        self.debug = debug
        self.ple = None
        self.MPI = None
        self.medcpl = None
        self.fields_in = []
        self.fields_out = []
        self.params = SchemeParams()

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
        func = self.ple.log if self.ple else logger
        if self.debug:
            prargs = kwargs.copy()
            prargs.pop("verbosity", None)
            prargs["flush"] = True
            print("[CPL]", *args, **prargs)
        return func(*args, **kwargs)

    def _init_paramedmem(self, with_app, interface):
        """Initialize ParaMEDMEM coupling.

        Arguments:
            with_app (str): name of coupling application
            interface (tuple(Mesh, list[str])): whole mesh and groups of the interface.
        """

        my_ranks = self.ple.get_app_ranks(self.whoami)
        other_ranks = self.ple.get_app_ranks(with_app)

        # Creating the parallel DEC
        fields_name = [name for name, _, _ in self.fields_in + self.fields_out]
        if self.starter:
            self.medcpl.init_coupling(fields_name, ranks1=my_ranks, ranks2=other_ranks)
        else:
            self.medcpl.init_coupling(fields_name, ranks1=other_ranks, ranks2=my_ranks)

        # Define coupling mesh
        self.medcpl.create_mesh_interface(interface[0], interface[1])

        # Define coupled fields
        for name, components, discr in self.fields_in + self.fields_out:
            self.medcpl.add_field(name, components, discr)

    def recv_input_fields(self):
        """Receive the input fields from the other code.

        Returns:
            (dict[*ParaFIELD*]): fields used to define the inputs at the next step.
        """
        names = [name for name, _, _ in self.fields_in]
        return self.medcpl.recv(names)

    def send_output_fields(self, outputs):
        """Send the output fields to the other code.

        Arguments:
            outputs (dict[*ParaFIELD*]): fields used to define the inputs of the other code.
        """
        assert len(outputs) == len(self.fields_out)
        self.medcpl.send(outputs)

    def set_parameters(self, params):
        """Set parameters.

        Arguments:
            params (dict): Parameters of the coupling scheme.
        """

        if self.starter:
            self.params.set_values(params)

        self.params.nb_iter = self.MPI.COUPLING_COMM_WORLD.bcast(
            self.starter, 0, "NBSSIT", self.params.nb_iter, self.MPI.INT
        )
        self.params.epsilon = self.MPI.COUPLING_COMM_WORLD.bcast(
            self.starter, 0, "EPSILO", self.params.epsilon, self.MPI.DOUBLE
        )

        if self.starter:
            if self.params.stepper:
                times = self.params.stepper._times
                nb_step = len(times)
            else:
                nb_step = 0
            self.MPI.COUPLING_COMM_WORLD.send(0, "NBPDTM", nb_step, self.MPI.INT)
            if nb_step > 0:
                self.MPI.COUPLING_COMM_WORLD.send(0, "STEP", self.params.init_time, self.MPI.DOUBLE)
                for i in range(nb_step):
                    self.MPI.COUPLING_COMM_WORLD.send(0, "STEP", times[i], self.MPI.DOUBLE)
        else:
            nb_step = self.MPI.COUPLING_COMM_WORLD.recv(0, "NBPDTM", self.MPI.INT)
            if nb_step > 0:
                times = [self.MPI.COUPLING_COMM_WORLD.recv(0, "STEP", self.MPI.DOUBLE)]
                for _ in range(nb_step):
                    times.append(self.MPI.COUPLING_COMM_WORLD.recv(0, "STEP", self.MPI.DOUBLE))

                self.params.set_values({"time_list": times})

    def init_coupling(self, with_app):
        """Initialize the coupling with an other application.

        Arguments:
            with_app (str): Name of the other application to be coupled with.
        """
        verbosity = 2 if self.debug else 1
        output = "all" if self.debug else "master"
        self.ple = PLE(verbosity=verbosity, logdir="/tmp", output_mode=output)
        self.ple.init_coupling(app_name=self.whoami, app_type="code_aster")

        self.other_app = with_app

        myranks = self.ple.get_app_ranks(self.whoami)
        self.log(f"allocated ranks for {self.whoami!r}: {myranks}", verbosity=verbosity)
        other_ranks = self.ple.get_app_ranks(with_app)
        self.log(f"allocated ranks for {with_app!r}: {other_ranks}", verbosity=verbosity)
        assert other_ranks, f"Application {with_app!r} not found!"
        self.MPI = MPI_CALCIUM(self.ple.base_comm, self.ple.my_comm, other_ranks[0], self.log)

        self.log(
            f"{self.whoami!r} coupler created from #{myranks[0]}, "
            f"{with_app!r} root proc is #{1}"
        )

        self.medcpl = MEDCoupler(logfunc=self.log)

    def setup(self, interface, input_fields, output_fields, **params):
        """Initialize the coupling.

        Arguments:
            interface (tuple(Mesh, list[str])): whole mesh and groups of the interface.
            input_fields (list): List of exchanged fields as input.
            output_fields (list): List of exchanged fields as output.
            params (dict): Parameters of the coupling scheme.
        """

        self.fields_in = input_fields
        self.fields_out = output_fields
        self.set_parameters(params)
        self._init_paramedmem(self.other_app, interface)

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

        output_data = exec_iteration(0, 0.0, 1.0, data, self.medcpl)

        return output_data

    def run(self, solver, **params):
        """Execute the coupling loop.

        Arguments:
            solver (object): Solver contains at least a method run_iteration.
            params (dict): Parameters of the coupling scheme.
        """

        # set parameters
        self.set_parameters(params)

        # initial sync before the loop
        exit_coupling = self.sync()

        stepper = self.params.stepper
        completed = False
        istep = 0

        while not stepper.isFinished():
            istep += 1
            current_time = stepper.getCurrent()
            delta_time = stepper.getIncrement()
            self.log("coupling step #{0:d}, time: {1:f}".format(istep, current_time))

            for i_iter in range(self.params.nb_iter):

                self.log("coupling iteration #{0:d}, time: {1:f}".format(i_iter, current_time))

                # recv data from the other code
                if self.starter and istep == 1:
                    input_data = {}
                    for name, _, _ in self.fields_in:
                        input_data[name] = None
                else:
                    input_data = self.recv_input_fields()

                has_cvg, output_data = solver.run_iteration(
                    i_iter, current_time, delta_time, input_data, self.medcpl
                )

                # send data to other code
                if self.starter:
                    self.send_output_fields(output_data)
                    converged = self.MPI.COUPLING_COMM_WORLD.allreduce(
                        i_iter, "ICVAST", has_cvg, self.MPI.BOOL, self.MPI.LAND
                    )
                else:
                    converged = self.MPI.COUPLING_COMM_WORLD.allreduce(
                        i_iter, "ICVAST", has_cvg, self.MPI.BOOL, self.MPI.LAND
                    )
                    self.send_output_fields(output_data)

                if converged:
                    break

                self.log("end of iteration status: {}".format(exit_coupling))

            stepper.completed()

            self.log("end of time step with status: {}".format(exit_coupling))

        # only to avoid deadlock
        if self.starter:
            input_data = self.recv_input_fields()

        self.log(
            "coupling {0} with exit status: {1}".format(
                "completed" if completed else "interrupted", exit_coupling
            )
        )

    @property
    def mesh_interface(self):
        """Mesh|ParallelMesh: mesh of the interface."""
        return self.medcpl.mesh_interf

    @property
    def mesh(self):
        """Mesh|ParallelMesh: coupled mesh."""
        return self.medcpl.mesh


class SaturneCoupling(ExternalCoupling):
    """Coupled simulation with code_aster and code_saturne with PLE.

    Arguments:
        app (str): Application name (default: "code_aster").
        starter (bool): *True* if this instance starts (it sends first),
            *False* if it receives first (default: "False").
        debug (bool): Enable debugging mode (default: "False")".

    Attributes:
        ple (*ple_coupler*): PLE Coupler object.
        tag (int): tag used for direct MPI communication (`CS_CALCIUM_MPI_TAG`
           ).
        mesh (*medcoupling.MEDCouplingUMesh*): Mesh of the interface.
        fields_in (list): List of exchanged fields (field name, number of
            components, components names).
        fields_out (list): List of exchanged fields (field name, number of
            components, components names).
        params (SchemeParams): Parameters of the coupling scheme.
    """

    def set_parameter(self, params):
        """Set parameters.

        Arguments:
            params (dict): Parameters of the coupling scheme.
        """

        raise NotImplemented()

    def run(self, solver, **params):
        """Execute the coupling loop.

        Arguments:
            solver (object): Solver contains at least a method run_iteration.
            params (dict): Parameters of the coupling scheme.
        """

        # set parameters
        self.set_parameter(params)

        # initial sync before the loop
        exit_coupling = self.sync()

        current_time = self.params.init_time
        delta_t = self.params.delta_t
        completed = False
        istep = 0

        while not completed and not exit_coupling:
            istep += 1

            for i_iter in range(self.params.nb_iter):

                self.MPI.COUPLING_COMM_WORLD.send(istep, "DTAST", delta_t, self.MPI.DOUBLE)
                delta_t = self.MPI.COUPLING_COMM_WORLD.recv(istep, "DTCALC", self.MPI.DOUBLE)

                self.log("coupling iteration #{0:d}, time: {1:f}".format(i_iter, current_time))

                # recv data from code_saturne
                current_time += delta_t
                input_data = self.recv_input_fields()

                has_cvg, output_data = solver.run_iteration(
                    i_iter, current_time, delta_time, input_data, self.medcpl
                )

                # received cvg
                converged = bool(self.recv(istep, "ICVAST", self.MPI.INT))

                # send results to code_saturne
                self.send_output_fields(output_data)

                if converged:
                    break

                exit_coupling = self.sync()
                self.log("end of iteration status: {}".format(exit_coupling))

            completed = current_time >= self.params.final_time
            exit_coupling = self.sync(end_coupling=completed)

            self.log("end of time step with status: {}".format(exit_coupling))

        self.log(
            "coupling {0} with exit status: {1}".format(
                "completed" if completed else "interrupted", exit_coupling
            )
        )
