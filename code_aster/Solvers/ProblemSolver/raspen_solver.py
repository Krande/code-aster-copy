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

from libaster import asmpi_free, asmpi_get, asmpi_set, asmpi_split

from ...Supervis import ConvergenceError
from ...Utilities import PETSc, no_new_attributes, profile
from ...Utilities.mpi_utils import MPI
from .snes_solver import SNESSolver
from .iterations_solver import IterationsSolver


def Print(*args):
    print(*args, flush=True)


class RASPENSolver(IterationsSolver):
    """Solves a step using PETSc SNES, loops on iterations."""

    _primal_incr = _resi_comp = None
    _scaling = _options = None
    local_solver = None
    __setattr__ = no_new_attributes(object.__setattr__)

    @classmethod
    def builder(cls, context):
        """Builder of RaspenSolver object.

        Args:
            context (Context): Context of the problem.

        Returns:
            instance: New object.
        """
        instance = cls()
        instance.context = context
        instance.local_solver = SNESSolver(local=False)
        instance.local_solver.context = context
        return instance

    def __init__(self):
        super().__init__()
        self._primal_incr = self._resi_comp = None
        self._scaling = self._options = None
        self.local_solver = None

    def initialize(self):
        """Initialize the object for the next step."""
        super().initialize()
        self.local_solver.initialize()

    # @profile
    def solve(self, current_matrix):
        """Solve a step with SNES.

        Raises:
            *ConvergenceError* exception in case of error.
        """

        self.oper.initialize()

        self._scaling = self.oper.getLagrangeScaling(self.matrix_type)
        self.current_matrix = self.oper.first_jacobian

        p_jac = self.current_matrix.toPetsc(local=True)

        DDPart = DomainDecomposition(self.problem.getMesh(), self.problem.getDOFNumbering())

        local_solver = self.local_solver
        local_solver.initSNES()

        self._primal_incr = self.state.primal_step.copy()

        raspen_solver = _RASPENSolver(
            local_solver,
            p_jac,
            self.state.primal_step,
            DDPart,
            None,
            0.0,
            self.current_incr,
            SolverType="RASPEN",
            BCType="Dirichlet",
            JacType="Exact",
            withFinalJac=True,
            withCoarsePb=False,
            comm=DDPart.comm,
        )
        rtol = self.get_keyword("CONVERGENCE", "RESI_GLOB_RELA")
        atol = self.get_keyword("CONVERGENCE", "RESI_GLOB_MAXI", 1.0e-24)
        maxiter = self.get_keyword("CONVERGENCE", "ITER_GLOB_MAXI")
        raspen_solver.setTolerances(rtol=rtol, atol=atol, maxiter=maxiter)

        # register the function in charge of
        # computing the nonlinear residual
        p_resi = p_jac.getVecRight()
        p_resi.set(0)

        def _monitor(snes, its, fgnorm):
            self.logManager.printConvTableRow(
                [its, " - ", fgnorm, " - ", self.local_solver.matrix_type]
            )

        raspen_solver.glbSnes.setMonitor(_monitor)

        OptDB = PETSc.Options()
        if not self._options:
            self.linear_solver.build()
            self._options = self.linear_solver.getPetscOptions()
        OptDB.insertString(self._options)
        raspen_solver.glbSnes.setFromOptions()

        # solve the nonlinear problem
        raspen_solver.solve()

        if raspen_solver.glbSnes.getConvergedReason() < 0:
            raise ConvergenceError("MECANONLINE9_7")

        self.oper.finalize()
        # delete local snes
        self.local_solver.snes = None

        return self.current_matrix


class _RASPENSolver:
    def __init__(
        self,
        LocNSL,
        locJac,
        primalStep,
        DDPart,
        Jr,
        alpha,
        current_incr,
        SolverType="RASPEN",
        BCType="Dirichlet",
        JacType="Exact",
        withFinalJac=False,
        withCoarsePb=False,
        comm=MPI.ASTER_COMM_WORLD,
    ):
        # Comm
        self.comm = comm
        # Rank
        self.rank = comm.Get_rank()
        # Size
        self.size = comm.Get_size()
        # Setting default tolerances
        self.rtol = 1e-8
        self.atol = 1e-18
        self.stol = 1e-25
        self.maxiter = 50
        # Current increment
        self.current_incr = current_incr
        # Local non-linear Solver
        self.asterLocSnl = LocNSL
        # Z0 and Z vectors
        self.ZVec = locJac.getVecRight()
        self.ZVec.set(0)
        self.Z0Vec = self.ZVec.copy()
        # aster solution
        self.asterGlbSolution = primalStep
        # Domain Decomposition partitioner
        self.DDPart = DDPart
        # Vecscatter
        self.vecscatter = self.DDPart.getVecScatter()
        # Solver
        self.SolverType = SolverType
        # Boundary conditions type
        self.BCType = BCType
        # Jacobian compute mode
        self.JacType = JacType
        # Local dofs
        self.local_dofs = self.DDPart.getLocalDofs()
        # Boundary and interior dofs
        self.locBoundaryDofs = DDPart.getLocalBoundaryDofs()
        self.locInteriorDofs = DDPart.getLocalInteriorDofs()
        # With final jacobian
        self.withFinalJac = withFinalJac
        # With coarse problem
        self.withCoarsePb = withCoarsePb
        # Global size
        self.glbSize = self.DDPart.getGlobalSize()
        # Local sizes
        self.locSize = DDPart.getLocalSize()
        self.interiorSize = len(self.locInteriorDofs)
        # The corresponding index set
        self.interiorIs = PETSc.IS().createGeneral(self.locInteriorDofs, comm=MPI.ASTER_COMM_SELF)
        # Local Jacobian
        [self.Jloc, self.Jploc, [self.JacFunc, _, _]] = self.asterLocSnl.snes.getJacobian()
        # Local Jacobian restricted on the interior
        self.Jiloc = self.Jloc.createSubMatrix(self.interiorIs, self.interiorIs)
        # Robin BC Form
        self.alpha = alpha
        # Robin BC rhs Jacobian
        self.Jr = Jr
        # A vector to hold local actions
        self.Yloc = self.ZVec.duplicate()
        # Local Rhs
        self.rhs = self.Yloc.getSubVector(self.interiorIs)
        # Global solution
        self.glbSol = PETSc.Vec().createMPI(self.glbSize, comm=self.comm)
        # Global Res
        self.Res = self.glbSol.duplicate()
        # Utility vector
        self._wrk1 = self.ZVec.duplicate()
        self._wrk2 = self.ZVec.duplicate()
        # Global Jacobian
        self.J = PETSc.Mat().createPython(self.glbSize, comm=MPI.ASTER_COMM_WORLD)
        # Global Snes
        self.glbSnes = self.setUpGlobalSnes()
        # Local snes
        self.locSnes = self.getSnes(self.asterLocSnl)
        self.locSnes.prefix = "lsnes_"
        # Local ksp
        self.locKsp = self.locSnes.getKSP()
        # Customized convergence test
        self.setGlbSnesConvergenceTest()

    def setTolerances(self, rtol=1e-6, atol=1e-24, maxiter=50):
        """
        Set Tolerances
        """
        self.rtol = rtol
        self.atol = atol
        self.maxiter = maxiter

    def setUpGlobalSnes(self):
        """
        Sets up the global snes
        """
        glbSnes = PETSc.SNES().create(comm=MPI.ASTER_COMM_WORLD)
        glbSnes.prefix = "gsnes_"
        # Snes Ksp
        SnesKsp = glbSnes.getKSP()
        # Snes type
        glbSnes.setType("newtonls")
        # Linear solver type
        SnesKsp.setType("fgmres")
        SnesKsp.setTolerances(rtol=self.rtol)
        SnesKsp.setGMRESRestart(500)
        # Linear precond type
        SnesKsp.getPC().setType("none")

        def snesksp_monitor(ksp, its, rnorm):
            if self.rank == 0:
                print(f" SNES_KSP : {its}" + f"   |    Residual Norm : {rnorm}", flush=True)

        # SnesKsp.setMonitor(snesksp_monitor)
        # Fix local boundary ddls in Dirichlet case
        self.ApplyLocalDirichlet()

        # Monitoring function
        def monitor(snes, its, fnorm):
            # This function is called at each SNES iteration
            # Computing function norm at the current iteration
            self.Fnorm = self.originalFunctionNorm()
            # self.fnorm0 = self.fnorm0 if its else self.Fnorm
            Rnorm = self.Fnorm / self.fnorm0
            snes.norm = self.Fnorm

            # if current_iteration > 0:
            #     snes.getKSP().view()
            # Outer iterations logs
            if self.rank == 0:
                print(
                    f" Outer Iteration : {its}"
                    + f"   |    Global Residual Norm : {self.Fnorm}"
                    + f"   |    Relative Residual Norm : {Rnorm}",
                    flush=True,
                )

        # Sets outer snes customized monitoring function
        glbSnes.setMonitor(monitor)

        # Defining the preconditioned function and jacobian
        def preFunc(snes, X, Y):
            # Get local solution approximation
            self.vecscatter.scatter(X, self.Z0Vec, mode="forward")
            # Set the approximation as the initial guess
            self.ZVec.getArray()[:] = self.Z0Vec.getArray()[:]
            # Synchronize code_aster's unknowns
            dx = self.asterLocSnl.snes.getSolutionUpdate()
            # Hack to circumvent the PETSc pointers' nullity if asterLocSnl hasn't start yet
            if dx.handle:
                dx.set(0)
                previous = self.asterLocSnl.snes.getSolution()
                # Print("<NIKO> asterGlbSolution.fromPetsc")
                self.asterGlbSolution.fromPetsc(previous, local=True)
                previous.copy(self._wrk1)
            # Local solve
            # -------------------------------------------------
            # - retrieve the current communicator
            mycomm = asmpi_get()
            # - split in single processes
            comm_self = asmpi_split(mycomm, MPI.ASTER_COMM_WORLD.Get_rank(), "self")
            # - change the current communicator
            MPI.ASTER_COMM_WORLD.Barrier()
            asmpi_set(comm_self)
            # - sequential solve
            self.asterLocSnl.snes.solve(None, self.ZVec)
            # Print(f"<NIKO> ZVec={self.ZVec.getArray(readonly=1)}")
            # Print(f"<NIKO> Z0Vec={self.Z0Vec.getArray(readonly=1)}")
            # -------------------------------------------------
            # Local correction
            C = self.Z0Vec - self.ZVec
            if self.SolverType in ["RAS", "RASPEN"]:
                self.DDPart.restrict(C)
            # Add corrections to get residual
            # Print("<NIKO> self.vecscatter.scatter")
            Y.set(0.0)
            self.vecscatter.scatter(-C, Y, mode="reverse", addv=True)
            # Compute final jacobian
            if self.withFinalJac or self.SolverType == "ASPIN":
                w = self.asterLocSnl.snes.getSolution()
                if self.SolverType == "ASPIN":
                    self.vecscatter.scatter(X + Y, w, mode="forward")
                self.locSnes.computeJacobian(w, self.Jloc, self.Jploc)
                self.locKsp.setOperators(self.Jloc, self.Jploc)
            # -------------------------------------------------
            # - switch back to initial communicator
            MPI.ASTER_COMM_WORLD.Barrier()
            asmpi_free(comm_self)
            asmpi_set(mycomm)
            # -------------------------------------------------
            if dx.handle:
                self.asterGlbSolution.fromPetsc(self._wrk1, local=True)
            else:
                self.asterGlbSolution *= 0

        def preJac(snes, X, J, Jp):
            # tj = time()
            self.current_incr += 1
            Jctx = JacCtx(self, snes, X)
            # Set up jacobian of type shell
            J.setPythonContext(Jctx)
            J.setUp()
            # print(" Jacobian constrcution time :",time()-tj,flush=True)

        # glbSnes.setUseMF(True)
        glbSnes.setFunction(preFunc, self.Res)
        glbSnes.setJacobian(preJac, self.J)
        glbSnes.setTolerances(rtol=self.rtol, atol=self.atol, stol=self.stol, max_it=self.maxiter)
        glbSnes.setSolution(self.glbSol)
        OptDB = PETSc.Options()
        OptDB.insertString("-gsnes_snes_linesearch_type basic  ")
        glbSnes.setFromOptions()
        return glbSnes

    def setGlbSnesConvergenceTest(self):
        """
        Sets the convergence test of the global snes
        """

        def test_on_original_function(snes, its, args):
            # residualNorm = self.originalFunctionNorm()
            if its == 0:
                return False
            cvg = self.Fnorm / self.fnorm0 < self.rtol or self.Fnorm < self.atol
            Print(f"cvg={cvg}")
            Print(f"self.rtol={self.rtol}")
            Print(f"self.atol={self.atol}")
            Print(f"self.Fnorm={self.Fnorm}")
            return cvg

        self.glbSnes.setConvergenceTest(test_on_original_function)

    def ApplyLocalDirichlet(self):
        """
        Applies Dirichlet local boundary conditions
        on each subdomain
        """
        BCDofs = self.locBoundaryDofs
        r, locFunc = self.asterLocSnl.snes.getFunction()
        A, P, locJac = self.asterLocSnl.snes.getJacobian()

        def localFunction(snes, Xl, Yl):
            locFunc[0](snes, Xl, Yl)
            Yl.getArray()[BCDofs] = 0.0
            if self.BCType == "Robin":
                Yl.getArray()[BCDofs] = self.Jr(Xl - self.Z0Vec)[BCDofs]

        def localJacobian(snes, Xl, J, Jp):
            locJac[0](snes, Xl, J, Jp)
            J.zeroRows(BCDofs)
            Jp.zeroRows(BCDofs)
            if self.BCType == "Robin":
                J[BCDofs, :] = self.Jr[BCDofs, :]
                Jp[BCDofs, :] = self.Jr[BCDofs, :]
                J.assemble()
                Jp.assemble()

        # Sets new operators
        self.asterLocSnl.snes.setFunction(localFunction, r)
        self.asterLocSnl.snes.setJacobian(localJacobian, A)

    def originalFunction(self, X, Y):
        """
        A function that do the action
        of the original problem
        """
        # Print(f"<NIKO> X={X.getArray(readonly=1)}")
        self.vecscatter.scatter(X, self._wrk1)
        # Print(f"<NIKO> work={work.getArray(readonly=1)}")
        self.asterLocSnl.snes.computeFunction(self._wrk1, self._wrk2)
        # Print(f"<NIKO> y={y.getArray(readonly=1)}")
        self.DDPart.restrict(self._wrk2)
        Y.set(0.0)
        self.vecscatter.scatter(self._wrk2, Y, mode="reverse", addv=True)
        # Print(f"<NIKO> Y={Y.getArray(readonly=1)}")

    def originalFunctionNorm(self):
        """
        Computes the norm of the original
        function
        """
        Y = self.glbSol.duplicate()
        # Print(f"<NIKO> glbSol={self.glbSol.getArray(readonly=1)}")
        self.originalFunction(self.glbSol, Y)
        return Y.norm()

    def setUpLocalKsp(self):
        """
        Sets up the local KSP to solves the linear
        problem on the interior of the subdomain
        """
        ksp = PETSc.KSP().create(MPI.ASTER_COMM_SELF)
        if self.BCType == "Robin":
            ksp.setOperators(self.Jloc, self.Jloc)
        else:
            ksp.setOperators(self.Jiloc, self.Jiloc)
        ksp.setFromOptions()
        ksp.setType("preonly")
        ksp.getPC().setType("lu")
        ksp.setUp()
        return ksp

    def getSnes(self, Locnsl):
        """
        Gets the local snes from the non-linear solver class
        """
        return self.asterLocSnl.snes

    def getGlbSnes(self):
        """
        Gets global snes
        """
        return self.glbSnes

    def getLocalKsp(self):
        """
        Gets the local ksp
        """
        return self.locKsp

    def getSolution(self):
        """
        Gets the solution petsc vector
        """
        return self.glbSol

    def solve(self, rhs=None):
        """
        RASPEN Solve
        """
        # Scattering data
        self.glbSol.set(0.0)
        self.DDPart.restrict(self.Z0Vec)
        self.vecscatter.scatter(self.Z0Vec, self.glbSol, mode="reverse", addv=True)
        self.fnorm0 = self.originalFunctionNorm()
        # Solving
        self.glbSnes.solve(rhs, self.glbSol)


class JacCtx:
    """
    Jacobian context class
    """

    def __init__(self, Sl: RASPENSolver, snes, Xs):
        self.F = snes.getFunction()[0]
        self.Sl = Sl
        self.snes = snes
        self.Xs = Xs
        self.Jloc = Sl.Jloc
        self.locVec = self.Sl.Yloc.duplicate()

    def mult(self, mat, X, Y):
        """
        Mat-Vec multiplication
        """
        # tm = time()
        if self.Sl.JacType == "Fd":
            # Finite difference
            X1, X2 = X.duplicate(), X.duplicate()
            eps = 1e-4
            self.Xs.copy(X1)
            X1.axpy(eps, X)
            self.snes.computeFunction(X1, X2)
            ((1.0 / eps) * (X2 - self.F)).copy(Y)
            # print("Jacobian multiplication time", time()-tm,flush=True)
        else:
            if self.Sl.SolverType in ["RASPEN", "ASPIN"]:
                # Scattering global direction to local vecs
                self.Sl.vecscatter.scatter(X, self.locVec, mode="forward")
                self.Sl.Jloc.mult(self.locVec, self.Sl.Yloc)
                self.Sl.Yloc.getArray()[self.Sl.locBoundaryDofs] = 0.0
                # tksp = time()
                # Print("************ SOLVE *************")
                self.Sl.locKsp.solve(self.Sl.Yloc, self.Sl.Yloc)
                #         "Exec Time :",time()-tksp,flush=True)
                if self.Sl.SolverType in ["RAS", "RASPEN"]:
                    self.Sl.DDPart.restrict(self.Sl.Yloc)
                Y.set(0.0)
                self.Sl.vecscatter.scatter(self.Sl.Yloc, Y, mode="reverse", addv=True)
                # print("Jacobian multiplication time", time()-tm,flush=True)
            else:
                X.copy(Y)


class DomainDecomposition:
    def __init__(self, mesh, nume_ddl, overlap=1):
        """
        DDM constructor
        """

        # Setting inputs
        self.mesh = mesh
        self.dofNbg = nume_ddl
        numeq = nume_ddl.getEquationNumbering()
        self.overlap = overlap
        self.comm = MPI.ASTER_COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()
        # in global numbering
        LGMap = nume_ddl.localToGlobalDOF
        neq = nume_ddl.getNumberOfDOFs(local=True)
        self.local_dofs = [LGMap(i) for i in range(neq)]
        # in global numbering
        self.non_overlap_dofs = numeq.getNoGhostDOFs(local=False)
        # all in ovp with ghosts in global numbering
        self.on_overlap_dofs = numeq.getGhostDOFs(local=False)
        # all in ovlp in local numbering
        self.local_ovlp_dofs = numeq.getGhostDOFs(local=True)
        # ghost in local numbering
        self.local_boundary_dofs = numeq.getGhostDOFs(local=True)
        # all except ghost in local numbering
        self.local_interior_dofs = numeq.getNoGhostDOFs(local=True)
        # Local size
        self.localSize = nume_ddl.getNumberOfDOFs(local=True)
        self.nonOvlpLocalSize = self.localSize - len(self.on_overlap_dofs)
        # Initializing DD variables
        self.createVecscatters()

    def createVecscatters(self):
        """
        Creates the vector scatters
        """
        # Create global and local vectors
        # Global Vec
        GlbVec = PETSc.Vec().create(comm=self.comm)
        GlbVec.setSizes(self.getGlobalSize())
        GlbVec.setFromOptions()
        GlbVec.setUp()
        # Local Vec
        LocVec = PETSc.Vec().create(comm=MPI.ASTER_COMM_SELF)
        LocVec.setSizes(self.localSize)
        LocVec.setFromOptions()
        LocVec.setUp()

        # Non-overlapping local Vec
        NonOvlpLocVec = PETSc.Vec().create(comm=MPI.ASTER_COMM_SELF)
        NonOvlpLocVec.setSizes(self.nonOvlpLocalSize)
        NonOvlpLocVec.setFromOptions()
        NonOvlpLocVec.setUp()

        # Creating Index Sets
        localIs = PETSc.IS().createGeneral(self.local_dofs, comm=self.comm)
        nonOvlpLocIs = PETSc.IS().createGeneral(self.non_overlap_dofs, comm=self.comm)
        locRangeIs = PETSc.IS().createStride(self.localSize, step=1)
        nonOvrlpRangeIs = PETSc.IS().createStride(self.nonOvlpLocalSize, step=1)

        # Create the scatter context for local to global mapping
        self.Vecscatter = PETSc.Scatter().create(GlbVec, localIs, LocVec, locRangeIs)
        self.nonOvlpVecScatter = PETSc.Scatter().create(
            GlbVec, nonOvlpLocIs, NonOvlpLocVec, nonOvrlpRangeIs
        )

    def restrict(self, petscVec):
        """
        Restricts a local vector on its local non-overlap
        """
        # Local size
        size = petscVec.getSize()
        # Assertions
        assert size == self.localSize
        # Restriction on the non-overlapping part
        petscVec.getArray()[self.local_ovlp_dofs] = 0.0

    def createVecOfWeights(self, type="restrict"):
        """
        Creates weights to apply on corrections
        """
        assert type in ["restrict", "average"]
        # Initializing vector of weights
        WeightsVec = PETSc.Vec().createSeq(self.localSize, comm=MPI.ASTER_COMM_SELF)
        WeightsVec.set(0.0)
        if type == "restrict":
            IsOnlocal = [list(self.local_dofs).index(dof) for dof in self.non_overlap_dofs]
            WeightsVec.getArray()[IsOnlocal] = 1.0
        #  TODO ("average" option)
        return WeightsVec

    def getVecScatter(self):
        """
        Gets vecscatter from global to local
        """

        if not self.Vecscatter:
            raise (RuntimeError("createVecscatters method must be " + "called before this routine"))

        return self.Vecscatter

    def getGlobalSize(self):
        """
        Gets the global size of the mesh
        """
        return self.dofNbg.getNumberOfDOFs(local=False)

    def getLocalSize(self):
        """
        Gets the local size of the mesh
        """
        return self.localSize

    def getLocalDofs(self):
        """
        Gets overlapping set of local DOFs in
        Global numbering representation
        """
        return self.local_dofs

    def getNonOverlappingLocalDofs(self):
        """
        Gets non-overlapping set of local DOFs in
        Global numbering representation
        """
        return self.non_overlap_dofs

    def getOnOverlapDofs(self):
        """
        Gets the set of dofs on the overlap
        """
        return self.on_overlap_dofs

    def getOnOverlapLocalDofs(self):
        """
        Gets the set of local dofs on the overlap
        """
        return self.local_ovlp_dofs

    def getLocalBoundaryDofs(self):
        """
        Gets local boundary dofs on local numbering
        """
        return self.local_boundary_dofs

    def getLocalInteriorDofs(self):
        """
        Gets local interior dofs on local numbering
        """
        return self.local_interior_dofs
