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

import medcoupling as MEDC


from code_aster.Commands import *
from code_aster import CA

from code_aster.Coupling import ExternalCoupling


class FakeSaturne(ExternalCoupling):
    """Fake saturne Solver.

    It allows to test SaturneCoupling without code_saturne installation.

    Arguments:
        debug (bool): Enable debugging mode (default: "False")".
    """

    def __init__(self, debug=False):
        super().__init__("code_saturne", True, debug)

    def set_parameters(self, params):
        """Set parameters.

        Arguments:
            params (dict): Parameters of the coupling scheme.
        """

        self._params.set_values(params)

        nb_step = (self._params.final_time - self._params.init_time) / self._params.delta_t
        self.MPI.COUPLING_COMM_WORLD.send(0, "NBPDTM", nb_step, self.MPI.INT)

        self.MPI.COUPLING_COMM_WORLD.send(0, "NBSSIT", self._params.nb_iter, self.MPI.INT)
        self.MPI.COUPLING_COMM_WORLD.send(0, "EPSILO", self._params.epsilon, self.MPI.DOUBLE)

        self.MPI.COUPLING_COMM_WORLD.send(0, "TTINIT", self._params.init_time, self.MPI.DOUBLE)
        self.MPI.COUPLING_COMM_WORLD.send(0, "PDTREF", self._params.delta_t, self.MPI.DOUBLE)

    def run(self, solver):
        """Execute the coupling loop.

        Arguments:
            solver (object): Solver contains at least a method run_iteration.
        """

        # initial sync before the loop
        exit_coupling = self.sync()

        stepper = self._params.stepper
        first_start = self._starter
        completed = False
        input_data = None
        istep = 0

        while not stepper.isFinished():
            istep += 1
            current_time = stepper.getCurrent()
            delta_time = stepper.getIncrement()

            _ = self.MPI.COUPLING_COMM_WORLD.recv(istep, "DTAST", self.MPI.DOUBLE)
            self.MPI.COUPLING_COMM_WORLD.send(istep, "DTCALC", delta_time, self.MPI.DOUBLE)

            print("coupling step #{0:d}, time: {1:f}".format(istep, current_time))

            for i_iter in range(self._params.nb_iter):

                print("coupling iteration #{0:d}, time: {1:f}".format(i_iter, current_time))

                if first_start:
                    first_start = False
                    input_data = {}
                    for name, _, _ in self._fields_in:
                        input_data[name] = None
                elif i_iter > 0:
                    input_data = self.recv_input_fields()

                has_cvg, output_data = solver.run_iteration(
                    i_iter, current_time, delta_time, input_data
                )

                # send cvg
                self.MPI.COUPLING_COMM_WORLD.send(istep, "ICVAST", int(has_cvg), self.MPI.INT)

                # send results to code_aster
                self.send_output_fields(output_data)

                if has_cvg:
                    break

            stepper.completed()
            input_data = self.recv_input_fields()
            print(
                f"end of time step {current_time} with status: {stepper.isFinished()}", flush=True
            )
            exit_coupling = self.sync(end_coupling=stepper.isFinished())

            print("end of time step with status: {}".format(exit_coupling))

        print(
            "coupling {0} with exit status: {1}".format(
                "completed" if completed else "interrupted", exit_coupling
            )
        )


def coupled_fluid(cpl, UNITE_MA):
    """Run Fluid coupling.

    Arguments:
        cpl (ExternalCoupling): Fluid coupler
    """
    ################################################################################
    # setup the simulation
    ################################################################################
    # send signal 6 (abort) to produce a traceback

    test = CA.TestCase()

    CA.init("--test", comm=cpl.MPI.ASTER_COMM_WORLD, debug=False, ERREUR=_F(ERREUR_F="ABORT"))

    # Read the  mesh - 2 cases : 1 or several procs
    if CA.MPI.ASTER_COMM_WORLD.size > 1:
        MAFLUIDE = LIRE_MAILLAGE(FORMAT="MED", UNITE=UNITE_MA, PARTITIONNEUR="PTSCOTCH")
    else:
        MAFLUIDE = LIRE_MAILLAGE(FORMAT="MED", UNITE=UNITE_MA)

    MAFLUIDE = MODI_MAILLAGE(
        reuse=MAFLUIDE,
        MAILLAGE=MAFLUIDE,
        ORIE_PEAU=_F(GROUP_MA_PEAU=("Face1", "Face2", "Face3", "Face4", "Face5", "Face6")),
    )

    # Assign model
    MOFLUIDE = AFFE_MODELE(
        MAILLAGE=MAFLUIDE,
        AFFE=_F(
            GROUP_MA=("Face1", "Face2", "Face3", "Face4", "Face5", "Face6"),
            PHENOMENE="MECANIQUE",
            MODELISATION="3D",
        ),
    )

    FORM = FORMULE(VALE="1.E-4*INST*(X+Y+Z)", NOM_PARA=["X", "Y", "Z", "INST"])
    FORM0 = FORMULE(VALE="0.", NOM_PARA=["X", "Y", "Z", "INST"])

    ################################################################################
    # define one iteration
    ################################################################################

    class FluidSolver:
        """Class that allows to compute one iteration of the coupling.
           This class has to contains at least run() function.

        Attributes:
            result (LoadResult): result of the thermal computation.
        """

        def __init__(self, cpl):
            """cpl (ExternalCoupling): coupler."""

            self._medcpl = cpl.medcpl
            self.epsilon = 1e-7
            self.depl_prev = None
            self.result = []

        def extent(self, depl):
            fns = depl.toSimpleFieldOnNodes()

            val, desc = fns.getValuesWithDescription()

            fns.setValues([0] * 3 * MAFLUIDE.getNumberOfNodes())
            fns.setValues(desc[0], desc[1], val)

            return fns.toFieldOnNodes()

        def run_iteration(self, i_iter, current_time, delta_t, data):
            """Execute one iteration.

            Arguments:
                i_iter (int): Iteration number if the current time_step.
                current_time (float): Current time.
                delta_t (float): Time step.
                data (dict[*MEDCouplingField*]): dict of input fields.

            Returns:
                bool: True if solver has converged at the current time step, else False.
                dict[*MEDCouplingField*]: Output fields, on nodes.
            """

            assert len(data) == 1, "expecting one field"

            mc_depl = data["Displ"]
            depl = None
            if mc_depl:
                # MEDC field => .med => code_aster field
                depl = self._medcpl.import_displacement(mc_depl)

            CHINST = CA.FieldOnNodesReal(MAFLUIDE, "INST_R", {"INST": current_time})

            CHXN = CREA_CHAMP(
                OPERATION="EXTR", TYPE_CHAM="NOEU_GEOM_R", NOM_CHAM="GEOMETRIE", MAILLAGE=MAFLUIDE
            )

            if depl:
                DEPL = self.extent(depl)
                CHXN += DEPL

            PRES_F = CREA_CHAMP(
                TYPE_CHAM="NOEU_NEUT_F",
                OPERATION="AFFE",
                MAILLAGE=MAFLUIDE,
                AFFE=_F(TOUT="OUI", NOM_CMP=("X1", "X2", "X3"), VALE_F=(FORM, FORM0, FORM0)),
                INFO=1,
            )

            CHNEUT = CREA_CHAMP(
                OPERATION="EVAL", TYPE_CHAM="NOEU_NEUT_R", CHAM_F=PRES_F, CHAM_PARA=(CHXN, CHINST)
            )

            force = CHNEUT.asPhysicalQuantity("FORC_R", {"X1": "FX", "X2": "FY", "X3": "FZ"})

            force_elem = force.toFieldOnCells(MOFLUIDE.getFiniteElementDescriptor(), "ELEM")

            if i_iter == 0:
                self.result.append(force_elem)
            else:
                self.result[-1] = force_elem

            # export

            mc_pres = self._medcpl.export_field(force_elem)

            # test convergence:
            has_cvg = False
            if i_iter > 0:
                if self.depl_prev:
                    depl_incr = depl - self.depl_prev
                    resi_rela = depl_incr.norm("NORM_INFINITY") / depl.norm("NORM_INFINITY")
                    has_cvg = resi_rela < self.epsilon
                    print(f"RESI_CPL: #iter {i_iter}, #resi: {resi_rela}")

            self.depl_prev = depl

            return has_cvg, {"Forces": mc_pres}

    ################################################################################
    # loop on time steps
    ################################################################################

    fluid_solv = FluidSolver(cpl)

    cpl.setup(
        interface=(MAFLUIDE, ["Face2", "Face3", "Face4", "Face5", "Face6"]),
        input_fields=[("Displ", ["DX", "DY", "DZ"], "NODES")],
        output_fields=[("Forces", ["FX", "FY", "FZ"], "CELLS")],
        init_time=0.0,
        final_time=1.0,
        nb_step=5,
        nb_iter=10,
        epsilon=fluid_solv.epsilon,
    )

    cpl.run(fluid_solv)

    # Extract the field from the result

    test.assertEqual(len(fluid_solv.result), 5)

    ################################################################################
    # Finalize the coupled study
    ################################################################################
    cpl.finalize()

    test.printSummary()

    CA.close()
