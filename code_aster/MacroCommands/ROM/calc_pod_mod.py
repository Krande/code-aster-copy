# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
import abc
import scipy.sparse
import numpy as np
from ...Utilities import PETSc, SLEPc
from mpi4py import MPI

comm = MPI.COMM_WORLD
global_size = comm.Get_size()

if global_size > 1:
    assert 1 == 0


# FUNCTIONNALITIES TO EXTRACT SNAPSHOT MATRIX FROM RESULT
def findIndexCHAM(lst, s):
    try:
        return lst.index(s)
    except ValueError:
        return None


def extractSnapshotsFromResult(result, chamName, indexSteps=None):
    """
    Extraction of snapshots from a SD RESULTAT.

    Parameters
    ----------
    result : SD RESULTAT
        code_aster result in which we seek snapshots
    chamName : str
        Name of the field in the result (example: DEPL or SIEF_ELGA)
    indexSteps : list or None
        List of the indices of the snapshots we seek to keep

    Returns
    -------
    snapshots : numpy.ndarray
        Snapshot array. Each column contains a given snapshot (size = number of dofs * number of snapshots).
    """
    # - Checks for the extraction procedure
    fieldsNames = result.getFieldsNames()
    if not fieldsNames:
        raise ValueError("Error in extraction procedure: No fields available in RESULTAT")
    chamIndex = findIndexCHAM(fieldsNames, chamName)
    if chamIndex is None:
        raise ValueError("Error in extraction procedure: Couldn't find the asked field in RESULTAT")
    # - Get proper data-structure for indexSteps : either None, int or list
    # if indexSteps is an integer, should be modified to be a list
    if indexSteps is None:
        # if None, all the timesteps are taken into account
        indStepsList = result.getIndexes()
        if isinstance(indStepsList, int):
            indStepsList = [indStepsList]
    else:
        if isinstance(indexSteps, int):
            indStepsList = [indexSteps]
        else:
            indStepsList = indexSteps
    # - Extraction of the snapshots
    snapshots = []

    for idx in indStepsList:
        if idx not in result.getIndexes():
            raise ValueError(
                f"Error in extraction procedure: Timestep index {idx} is not available in RESULTAT"
            )
        cham = result.getField(chamName, idx)
        values = cham.getValues()
        snapshots.append(np.array(values))

    return np.column_stack(snapshots)


def transferSnapshotsToPETSC(snapshots):
    if global_size == 1:
        snapshots_petsc = PETSc.Mat().createDense(
            snapshots.shape, array=snapshots, comm=PETSc.COMM_SELF
        )
        snapshots_petsc.assemble()
    else:
        raise ValueError("MPI version not implemented yet")
        # mat_petsc = PETSc.Mat().create(comm)
        # mat_petsc.setSizes(((None, numpy_matrix.shape[0]), (None, numpy_matrix.shape[1])))
        # mat_petsc.setType('dense')
        # mat_petsc.setUp()

        # # each process should properly copy its part
        # rstart, rend = mat_petsc.getOwnershipRange()
        # local_view = mat_petsc.getDenseArray()
        # if rend > rstart:
        #     local_view[:, :] = numpy_matrix[rstart:rend, :]
        # mat_petsc.restoreDenseArray(local_view)
        # mat_petsc.assemble()
    return snapshots_petsc


# POST-TREATMENT FUNCTIONNALITIES


def computeProjectionErrors(Phi, snapshots):
    """Compute the projection errors on snaphots knowing a reduced order basis

    Parameters
    ----------
    Phi : numpy.ndarray
        Reduced order basis.
    snapshots : numpy.ndarray
        Snapshot array. Each column contains a given snapshot (size = number of dofs * number of snapshots).

    Returns
    -------
    abs_errors_arr : numpy.ndarray
        Array of absolute projection error
    rel_errors_arr : numpy.ndarray
        Array of relative projection error
    """
    abs_errors_arr = []
    rel_errors_arr = []

    for i in range(snapshots.shape[1]):
        u = snapshots[:, i]
        ## - Projection and reconstruction
        u_proj = Phi @ (Phi.T @ u)
        ## - Compute errors
        norm_u = np.linalg.norm(u)
        abs_error = np.subtract(u, u_proj)

        if norm_u > 1e-12:
            rel_error = np.linalg.norm(abs_error) / norm_u
        else:
            rel_error = 0.0

        abs_errors_arr.append(abs_error)
        rel_errors_arr.append(rel_error)

    return np.array(abs_errors_arr), np.array(rel_errors_arr)


def computeProjectionErrors_petsc(Phi_petsc, snapshots_petsc):
    # - Get the current communicator
    comm = snapshots_petsc.getComm()
    n_dim, n_snapshots = snapshots_petsc.getSize()
    # - Create temporay vectors
    u_proj_vec, projected_coords_vec = Phi_petsc.createVecs()
    # # - PETSC matrix to store the differences for error computations
    # abs_errors_mat = PETSc.Mat().createDense([n_dim, n_snapshots], comm=comm)
    # abs_errors_mat.setUp()
    # - List to store
    abs_errors_list = []
    rel_errors_list = []

    for i in range(n_snapshots):
        u_vec = snapshots_petsc.getColumnVector(i)
        Phi_petsc.multTranspose(u_vec, projected_coords_vec)
        Phi_petsc.mult(projected_coords_vec, u_proj_vec)

        # - compute u - u_proj
        abs_error_vec = u_vec.copy()
        abs_error_vec.axpy(-1, u_proj_vec)

        # - compute errors
        norm_u = u_vec.norm(PETSc.NormType.NORM_2)
        abs_error = abs_error_vec.norm(PETSc.NormType.NORM_2)
        if norm_u > 1e-12:
            rel_error = abs_error / norm_u
        else:
            rel_error = 0.0

        # abs_errors_mat.setColumnVector(i, abs_error_vec)
        abs_errors_list.append(abs_error)
        rel_errors_list.append(rel_error)
    # Clean the temporary vectors
    u_proj_vec.destroy()
    projected_coords_vec.destroy()
    return np.array(abs_errors_list), np.array(rel_errors_list)


# CLASS DEFINITION FOR A POD ANALYSIS
POD_VALID_METHOD = ["SVD", "snapshot", "GS-classical", "GS-modified"]
POD_METHOD_WITHOUT_CRIT = ["GS-classical", "GS-modified"]
assert all(item in POD_VALID_METHOD for item in POD_METHOD_WITHOUT_CRIT)
POD_CRITERION_METHOD = ["energy", "nbModes"]
INCR_POD_VALID_METHOD = ["HPOD", "HAPOD"]
GS_METHOD = ["classical", "modified"]


class PODAnalysisBase(abc.ABC):
    """
    Abstract class for building a base incrementally using POD.

    This class is designed to build a base from a set of snapshots.
    """

    def __init__(
        self, snapshots, method="SVD", criterion="energy", tolerance=None, nbModes=None, CorrOp=None
    ):
        self._snapshots = snapshots  # Each column is a snapshot
        ## - Initialisation parameters
        self._methodCompress = None
        self._criterionModes = None
        self._crit_tolerance = tolerance
        self._crit_nbModes = nbModes
        self._tol_num = 1e-12
        ## - Prepare the operators for POD Analysis
        self.setInfosSnapshots()
        self.setCompressionMethod(method)
        self.setCorrelationOperator(CorrOperator=CorrOp)
        self.setCriterionModes(criterion)
        ## - Tests
        self.correctionSnapshots()
        self.runCompatibilityTests()
        assert self._methodCompress is not None
        assert self._criterionModes is not None

    def runCompatibilityTests(self):
        """Testing compatibility between options (for arguments)"""
        if self._methodCompress in ["GS-classical", "GS-modified"] and self._crit_tolerance is None:
            raise ValueError("If method is of GS type, user should provide a tolerance")
        if self._criterionModes == "energy" and self._crit_tolerance is None:
            raise ValueError("If criterion = energy, user should provide a tolerance")
        if self._criterionModes == "nbModes" and self._crit_nbModes is None:
            raise ValueError("If criterion = nbModes, user should provide a number of modes")

    def setCompressionMethod(self, method):
        """Set method for the compression method"""
        if method in POD_VALID_METHOD:
            self._methodCompress = method
        else:
            raise ValueError(
                f"PODAnalysis: Method '{method}' is not valid. Choose method in {POD_VALID_METHOD}."
            )

    def setCriterionModes(self, criterion):
        """Set method for the criterion for selecting modes"""
        if criterion in POD_CRITERION_METHOD:
            self._criterionModes = criterion
        else:
            raise ValueError(
                f"PODAnalysis: Method '{criterion}' is not valid. Choose method in {POD_CRITERION_METHOD}."
            )

    def getCompressionMethod(self, method):
        """Get method for the criterion for selecting modes"""
        return self._methodCompress

    def getCriterionModes(self, method):
        """Get method for the compression method"""
        return self._criterionModes

    def computePODBasis(self, option=1):
        """Method to construct a reduced order basis by POD
        using the stored snapshots

        Arguments
        ----------
        option : int
            Changes the outputs of the function. If option=1, only reduced order basis.
            If option=2, returns reduced order basis and singular values.

        """
        return self.POD(self._snapshots, option=option)

    @abc.abstractmethod
    def setInfosSnapshots(self):
        pass

    @abc.abstractmethod
    def setCorrelationOperator(self, CorrOperator=None):
        pass

    @abc.abstractmethod
    def correctionSnapshots(self):
        pass

    def selectModes(self, Phi, singval, tolerance=None, nbModes=None):
        ## - Test
        if tolerance is None and nbModes is None:
            raise ValueError(
                "selectModes: Either tolerance or nbModes should be provided to the method"
            )
        ## - Choice of the number of modes to keep
        if self._criterionModes == "energy":
            s_squared = singval**2
            sum_i = 0
            i = 0
            while i < len(singval) and sum_i / np.sum(s_squared) < (1 - tolerance):
                sum_i += s_squared[i]
                i += 1
            nbModes_v = i
        elif self._criterionModes == "nbModes":
            if nbModes is None or nbModes > len(singval):
                raise ValueError("selectModes: nbModes should be given or is too big!")
            nbModes_v = nbModes
        else:
            raise ValueError(
                f"PODAnalysis: Method '{self._criterionModes}' is not valid. Choose method in {POD_CRITERION_METHOD}."
            )
        ## - Truncation of the POD Basis
        Phi_v = self.correctionModesAfterSelection(Phi, nbModes_v)
        singval_v = singval[:nbModes_v]
        return Phi_v, singval_v

    @abc.abstractmethod
    def correctionModesAfterSelection(self, Phi, nbModes_v):
        pass

    def POD(self, matS, option):
        """Method to construct a reduced order basis by POD

        Arguments
        ----------
        matS : format depends on class
            Matrix of snapshots on which the POD operator should be applied
        option : int
            Changes the outputs of the function. If option=1, only reduced order basis.
            If option=2, returns reduced order basis and singular values.

        Returns
        -------
        Phi : format depends on class
            Reduced order basis
        singval : numpy.ndarray
            Singular values (only if option=2)
        """
        ## - Compression step
        if self._methodCompress == "snapshot":
            Phi, singval = self.snapshotMethod(matS)
        elif self._methodCompress == "SVD":
            Phi, singval = self.SVDMethod(matS)
        elif self._methodCompress == "GS-classical":
            Phi, singval = self.GSMethod(matS, self._crit_tolerance, "classical")
        elif self._methodCompress == "GS-modified":
            Phi, singval = self.GSMethod(matS, self._crit_tolerance, "modified")
        else:
            raise ValueError(
                f"PODAnalysis: Method '{self._methodCompress}' is not implemented yet."
            )
        if self._methodCompress in POD_METHOD_WITHOUT_CRIT:
            Phi_t, singval_t = Phi, singval
        else:
            ## - Apply truncation
            Phi_t, singval_t = self.selectModes(
                Phi, singval, self._crit_tolerance, self._crit_nbModes
            )
        ## - Return outputs
        if option == 1:
            return Phi_t
        elif option == 2:
            return Phi_t, singval_t
        else:
            raise ValueError("PODAnalysis: computePODBasis should be 1 or 2.")

    @abc.abstractmethod
    def SVDMethod(self, matS, verbose=True):
        pass

    @abc.abstractmethod
    def snapshotMethod(self, matS):
        pass

    @abc.abstractmethod
    def updateGStype(self, Phi, singval, snapshot_new, tole, methodGS):
        pass

    @abc.abstractmethod
    def computeGSprojection(self, Phi, s_new, methodGS):
        pass

    @abc.abstractmethod
    def GSMethod(self, matS, tole, methodGS):
        pass

    @abc.abstractmethod
    def computePODBasisIncremental(self, Phi, singval=None, method="HPOD"):
        pass

    def computeDecayRate(self, singval):
        """Compute a decay rate of a list of singular values

        Parameters
        ----------
        singval : numpy.ndarray
            Singular values (order in a decreasing manner)

        Returns
        -------
        Phi : float
            Decay rate
        """
        nbSing = len(singval)

        N = np.log(np.arange(1, nbSing + 1))
        Y = np.log(singval)

        sum_N = np.sum(N)
        sum_Y = np.sum(Y)
        sum2_N = np.sum(N**2)

        decayRate = (nbSing * np.dot(Y.T, N) - sum_N * sum_Y) / (nbSing * sum2_N - sum_N**2)
        return decayRate


class PODAnalysisNumpy(PODAnalysisBase):
    def setInfosSnapshots(self):
        """Store information about the size of the snapshot matrix"""
        self._numberOfDOFs = self._snapshots.shape[0]
        self._numberOfSnapshots = self._snapshots.shape[1]

    def setCorrelationOperator(self, CorrOperator=None):
        """Set method for the correlation operator
        
       Arguments
        ----------
        CorrOperator : scipy.sparse.csr_matrix or None
            Correlation operator when using the method of snapshots
        """
        if CorrOperator is None:
            self._CorrOperator = scipy.sparse.identity(self._numberOfDOFs, format="csr")
        else:
            self._CorrOperator = CorrOperator

    def correctionSnapshots(self):
        """Correct the snapshots by removing null values"""
        norms = np.linalg.norm(self._snapshots, axis=0)
        self._snapshots = self._snapshots[:, norms > self._tol_num]

    def correctionModesAfterSelection(self, Phi, nbModes_v):
        """Construction of a reduced basis by selecting
        the most important modes (the first nbModes_v)

        Arguments
        ----------
        Phi : numpy.ndarray
            Basis
        nbModes_v : int
            number of modes to keep (the first ones)
        """
        return Phi[:, :nbModes_v]

    def computePODBasisIncremental(self, Phi, singval=None, method="HPOD"):
        """Method to enrich a reduced order basis with the stored snapshots

        Arguments
        ----------
        Phi : numpy.ndarray
            Reduced order basis which has been previously computed
        method : str
            Name of the incremental approach to use
        """
        assert method in INCR_POD_VALID_METHOD
        if method == "HPOD":
            matS = self._snapshots
            projS = np.zeros(np.shape(matS))
            for i in range(matS.shape[1]):
                projS[:, i] = self.computeGSprojection(Phi, matS[:, i], "modified")
            Phi_new = self.POD(projS, option=1)
            return np.column_stack((Phi, Phi_new))
        elif method == "HAPOD":
            assert singval is not None
            mPhi = singval * Phi
            assert np.shape(Phi) == np.shape(mPhi)
            mS = np.column_stack((mPhi, self._snapshots))
            return self.POD(mS, option=1)
        else:
            raise ValueError(
                f"PODAnalysis: Method '{method}' is not valid. Choose method in {INCR_POD_VALID_METHOD}."
            )

    def updateGStype(self, Phi, singval, snapshot_new, tole, methodGS):
        """Update a basis with a new snapshot method using a Gram-Schmidt process

        Arguments
        ----------
        Phi : numpy.ndarray
            Basis to enrich
        singval : numpy.ndarray
            Singular values
        snapshot_new : numpy.ndarray
            Snapshot to add
        tole : float
            Tolerance used for the test when adding new snapshot
        methodGS : str
            Should be "classical" or "modified" = GS method applied

        Returns
        -------
        Phi : numpy.ndarray
            Reduced order basis
        singval : numpy.ndarray
            Singular values
        """
        assert methodGS in GS_METHOD
        ## - Check that the added snapshot is 1D
        s_new = snapshot_new.flatten()
        s_new_norm = np.linalg.norm(snapshot_new)
        ## - GS orthogonalisation
        s_new = self.computeGSprojection(Phi, s_new, methodGS)
        s_new_perp_norm = np.linalg.norm(s_new)
        if s_new_perp_norm > tole * s_new_norm:
            Phi = np.column_stack((Phi, s_new / s_new_perp_norm))
            singval = np.hstack([singval, s_new_perp_norm])
        return Phi, singval

    def GSMethod(self, matS, tole, methodGS):
        """Compression method using a Gram-Schmidt process

        Arguments
        ----------
        matS : numpy.ndarray
            Matrix of snapshots on which the POD operator should be applied
        tole : float
            Tolerance used for the test when adding new snapshot
        methodGS : str
            Should be "classical" or "modified" = GS method applied

        Returns
        -------
        Phi : numpy.ndarray
            Reduced order basis
        singval : numpy.ndarray
            Singular values
        """
        s_0_norm = np.linalg.norm(matS[:, 0])
        Phi = matS[:, 0:1] / s_0_norm
        singval = np.array([s_0_norm])
        n_snap = matS.shape[1]
        for i in range(1, n_snap):
            Phi, singval = self.updateGStype(Phi, singval, matS[:, i : i + 1], tole, methodGS)
        return Phi, singval

    def computeGSprojection(self, Phi, s_new, methodGS):
        """Compute a Gram-Schmidt projection

        Arguments
        ----------
        Phi : numpy.ndarray
            Basis
        s_new : numpy.ndarray
            New snapshot considered
        methodGS : str
            Should be "classical" or "modified" = GS method applied

        Returns
        -------
        s_proj : numpy.ndarray
            Projected snapshot
        """
        assert methodGS in GS_METHOD
        ## - GS orthogonalisation
        n_modes = Phi.shape[1]
        for kp in range(2):  # Kahan-Parlett process
            if methodGS == "classical":
                s_new_loc = s_new
                for k in range(n_modes):
                    s_new = s_new - np.dot(s_new_loc, Phi[:, k]) * Phi[:, k]
            if methodGS == "modified":
                for k in range(n_modes):
                    s_new = s_new - np.dot(s_new, Phi[:, k]) * Phi[:, k]
        return s_new

    def SVDMethod(self, matS, verbose=True):
        """Compression method using SVD on the snapshot matrix

        Arguments
        ----------
        matS : numpy.ndarray
            Matrix of snapshots on which the POD operator should be applied
        Returns
        -------
        Phi : numpy.ndarray
            Reduced order basis
        singval : numpy.ndarray
            Singular values
        """
        ## - Apply SVD directly on the snapshot matrix
        U, sigma, _ = np.linalg.svd(matS, full_matrices=False)
        ## - Order eigenvalues and compute basis
        n = np.where(sigma == 0)[0]
        if n.size == 0:
            n = len(sigma)
        else:
            n = n[0]
        return U[:, :n], sigma[:n]

    def snapshotMethod(self, matS):
        """Compression method using the snapshot method on a correlation matrix

        Arguments
        ----------
        matS : numpy.ndarray
            Matrix of snapshots on which the POD operator should be applied

        Returns
        -------
        Phi : numpy.ndarray
            Reduced order basis
        singval : numpy.ndarray
            Singular values
        """
        ## - Compute correlation matrix
        corrMatrix = matS.T @ self._CorrOperator @ matS
        ## - Solve eigenproblem
        eigenvalues, eigenvectors = np.linalg.eigh(corrMatrix)
        ## - Order eigenvalues and compute basis
        idx = np.argsort(eigenvalues)[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvalues = np.where(eigenvalues < 0, 0, eigenvalues)
        eigenvectors = eigenvectors[:, idx]

        n = np.where(eigenvalues == 0)[0]
        if n.size == 0:
            n = len(eigenvalues)
        else:
            n = n[0]
        # - Return reduced order basis and singular values
        return np.dot(matS, eigenvectors[:, :n]) / np.sqrt(eigenvalues[:n]), np.sqrt(
            eigenvalues[:n]
        )

class PODAnalysisPetsc(PODAnalysisBase):

    def __init__(self, *args, **kwargs):
        "Initialization of MatrixScaler"
        super().__init__(*args, **kwargs)
        # - Add communicator to handle the parallel version
        self._comm = self._snapshots.getComm()

    def setInfosSnapshots(self):
        """Store information about the size of the snapshot matrix"""
        self._numberOfDOFs, self._numberOfSnapshots = self._snapshots.getSize()

    def setCorrelationOperator(self, CorrOperator=None):
        """Set method for the correlation operator

        Arguments
        ----------
        CorrOperator : PETSc.Mat or None
            Correlation operator when using the method of snapshots
        """
        if CorrOperator is None:
            I = PETSc.Mat().create(self._comm)
            I.setSizes(((None, self._numberOfDOFs), (None, self._numberOfDOFs)))
            I.setType("aij")
            I.setPreallocationNNZ(1)
            # diag_vec = PETSc.Vec().create(self._comm)
            # diag_vec.setSizes((None, self._numberOfDOFs))
            # diag_vec.setUp()
            diag_vec = I.createVecLeft()
            diag_vec.set(1.0)

            I.setDiagonal(diag_vec)
            diag_vec.destroy()
            I.assemble()
            self._CorrOperator = I
        else:
            self._CorrOperator = CorrOperator

    def correctionSnapshots(self):
        """Correct the snapshots by removing null values"""
        col_vec = self._snapshots.createVecLeft()
        indices_to_keep = []
        for j in range(self._numberOfSnapshots):
            self._snapshots.getColumnVector(j, col_vec)
            if col_vec.norm() > self._tol_num:
                indices_to_keep.append(j)

        if len(indices_to_keep) < self._numberOfSnapshots:
            is_rows = PETSc.IS().createStride(M, 0, 1)
            is_cols = PETSc.IS().createGeneral(indices_to_keep)
            self._snapshots = self._snapshots.createSubMatrix(is_rows, is_cols)

    def correctionModesAfterSelection(self, Phi, nbModes_v):
        """Construction of a reduced basis by selecting
        the most important modes (the first nbModes_v)

        Arguments
        ----------
        Phi : PETSc.Mat
            Basis
        nbModes_v : int
            number of modes to keep (the first ones)
        """
        num_rows = Phi.getSize()[0]
        is_rows = PETSc.IS().createStride(num_rows, 0, 1)
        is_cols = PETSc.IS().createStride(nbModes_v, 0, 1)
        Phi_v = Phi.createSubMatrix(is_rows, is_cols)
        return Phi_v

    def _stack_matrices_horizontally(self, mat_a, mat_b):
        """Compute a new matrix [mat_a, mat_b]

        Arguments
        ----------
        mat_a : PETSc.Mat
        mat_b : PETSc.Mat
        """
        M, N1 = mat_a.getSize()
        _, N2 = mat_b.getSize()

        new_mat = PETSc.Mat().createDense([M, N1 + N2], comm=mat_a.getComm())

        col_vec = mat_a.createVecLeft()
        for j in range(N1):
            mat_a.getColumnVector(j, col_vec)
            new_mat.setColumnVector(j, col_vec)
        for j in range(N2):
            mat_b.getColumnVector(j, col_vec)
            new_mat.setColumnVector(N1 + j, col_vec)

        return new_mat

    def computePODBasisIncremental(self, Phi, singval=None, method="HPOD"):
        """
        Method to enrich a reduced order basis with the stored snapshots

        Arguments
        ----------
        Phi : PETSc.Mat
            Reduced order basis which has been previously computed
        method : str
            Name of the incremental approach to use
        """
        assert method in INCR_POD_VALID_METHOD

        if method == "HPOD":
            matS = self._snapshots
            M, N = matS.getSize()

            projS = matS.duplicate(copy=False)
            projS.zeroEntries()

            # Temporary vectors
            col_s = matS.createVecLeft()
            col_proj = matS.createVecLeft()

            for i in range(N):
                matS.getColumnVector(i, col_s)
                col_proj = self.computeGSprojection(Phi, col_s, "modified")
                projS.setColumnVector(i, col_proj)

            Phi_new = self.POD(projS, option=1)
            return self._stack_matrices_horizontally(Phi, Phi_new)

        elif method == "HAPOD":
            assert singval is not None
            mPhi = Phi.duplicate(copy=True)

            col_vec = mPhi.createVecLeft()
            for j in range(len(singval)):
                mPhi.getColumnVector(j, col_vec)
                col_vec.scale(singval[j])
                mPhi.setColumnVector(j, col_vec)

            mS = self._stack_matrices_horizontally(mPhi, self._snapshots)
            return self.POD(mS, option=1)

        else:
            raise ValueError(
                f"PODAnalysis: Method '{method}' is not valid. Choose method in {INCR_POD_VALID_METHOD}."
            )

    def updateGStype(self, Phi, singval, snapshot_new, tole, methodGS):
        """Update of a reduced order basis using a Gram-Schmidt process

        Arguments
        ----------
        Phi : PETSc.Mat
            Basis
        singval : numpy.ndarray
            Singular values
        snapshot_new : PETSc.Vec
            New snapshot considered
        tole : float
            Tolerance used for the test when adding new snapshot
        methodGS : str
            Should be "classical" or "modified" = GS method applied

        Returns
        -------
        Phi : PETSc.Mat
            Reduced order basis
        singval : numpy.ndarray
            Singular values
        """
        assert methodGS in GS_METHOD

        s_new_norm = snapshot_new.norm()

        s_new_proj = self.computeGSprojection(Phi, snapshot_new, methodGS)
        s_new_perp_norm = s_new_proj.norm()

        if s_new_perp_norm > tole * s_new_norm:
            M, N = Phi.getSize()
            Phi_new = PETSc.Mat().createDense([M, N + 1], comm=Phi.getComm())

            # Copy old columns
            col_vec = Phi.createVecLeft()
            for j in range(N):
                Phi.getColumnVector(j, col_vec)
                Phi_new.setColumnVector(j, col_vec)

            # Add new column
            s_new_proj.scale(1.0 / s_new_perp_norm)
            Phi_new.setColumnVector(N, s_new_proj)

            # Update the singular values
            singval_new = np.hstack([singval, s_new_perp_norm])

            return Phi_new, singval_new
        else:
            # If nothing is needed, return basis and singular values unchanged
            return Phi, singval

    def GSMethod(self, matS, tole, methodGS):
        """Compression method using a Gram-Schmidt process

        Arguments
        ----------
        matS : PETSc.Mat
            Matrix of snapshots on which the POD operator should be applied
        tole : float
            Tolerance used for the test when adding new snapshot
        methodGS : str
            Should be "classical" or "modified" = GS method applied

        Returns
        -------
        Phi : PETSc.Mat
            Reduced order basis
        singval : numpy.ndarray
            Singular values
        """
        col_vec = matS.createVecLeft()
        matS.getColumnVector(0, col_vec)
        s_0_norm = col_vec.norm()

        singval = np.array([s_0_norm])
        Phi = PETSc.Mat().createDense([M, 1], comm=self._comm)
        if s_0_norm > self._tol_num:
            col_vec.scale(1.0 / s_0_norm)
        Phi.setColumnVector(0, col_vec)

        for i in range(1, n_snap):
            matS.getColumnVector(i, col_vec)
            Phi, singval = self.updateGStype(Phi, singval, col_vec, tole, methodGS)
        return Phi, singval

    def computeGSprojection(self, Phi, s_new, methodGS):
        """Compute a Gram-Schmidt projection

        Arguments
        ----------
        Phi : PETSc.Mat
            Basis
        s_new : PETSc.Vec
            New snapshot considered
        methodGS : str
            Should be "classical" or "modified" = GS method applied

        Returns
        -------
        s_proj : PETSc.Vec
            Projected snapshot
        """
        assert methodGS in GS_METHOD

        s_proj = s_new.copy()

        M, n_modes = Phi.getSize()
        phi_k = Phi.createVecLeft()

        for kp in range(2):
            if methodGS == "classical":
                s_loc = s_proj.copy()
                for k in range(n_modes):
                    Phi.getColumnVector(k, phi_k)
                    alpha = s_loc.dot(phi_k)
                    # s_proj = s_proj - alpha * phi_k
                    s_proj.axpy(-alpha, phi_k)
            elif methodGS == "modified":
                for k in range(n_modes):
                    Phi.getColumnVector(k, phi_k)
                    alpha = s_proj.dot(phi_k)
                    # s_proj = s_proj - alpha * phi_k
                    s_proj.axpy(-alpha, phi_k)
        return s_proj

    def SVDMethod(self, matS, verbose=True):
        """Compression method using SVD on the snapshot matrix

        Arguments
        ----------
        matS : PETSc.Mat
            Matrix of snapshots on which the POD operator should be applied
        Returns
        -------
        Phi : PETSc.Mat
            Reduced order basis
        singval : numpy.ndarray
            Singular values
        """
        tol = 1e-12
        comm = matS.getComm()
        Print = PETSC.Sys.Print

        # - Creation and configuration of the SVD solver
        svd = SLEPc.SVD().create(comm=comm)
        svd.setOperator(matS)
        S.setType(S.Type.TRLANCZOS)

        # - Ask to compute all the needed singular values
        _, n_cols = matS.getSize()
        svd.setDimensions(nsv=n_cols)
        S.setFromOptions()
        # - Resolution
        svd.solve()

        n_conv = svd.getConverged()
        if verbose and comm.rank == 0:
            Print("******************************")
            Print("*** SLEPc SVD Solution Results ***")
            Print("******************************\n")
            svd_type = svd.getType()
            Print(f"Solution method: {svd_type}")
            its = svd.getIterationNumber()
            Print(f"Number of iterations: {its}")
            nsv, _, _ = svd.getDimensions()
            Print(f"Number of requested singular values: {nsv}")
            tol, maxit = svd.getTolerances()
            Print(f"Stopping condition: tol={tol:.4g}, maxit={maxit}")
            Print(f"Number of converged singular triplets: {n_conv}\n")

        # - Construction of the outputs
        singular_values = []
        singular_vectors = []

        if n_conv > 0:
            u_tmp, v_tmp = matS.createVecs()

            if verbose and comm.rank == 0:
                Print("    sigma       residual norm ")
                Print("-------------  ---------------")

            for i in range(n_conv):
                sigma_i = svd.getSingularTriplet(i, u_tmp, v_tmp)
                singular_values.append(sigma_i)
                singular_vectors.append(u_tmp.copy())

                if verbose and comm.rank == 0:
                    error = svd.computeError(i)
                    Print(f"   {sigma_i:6f}     {error:12g}")

            if verbose and comm.rank == 0:
                Print()

            u_tmp.destroy()
            v_tmp.destroy()

        sigma = np.array(singular_values)
        zero_indices = np.where(sigma < tol)[0]
        n = zero_indices[0] if zero_indices.size > 0 else len(sigma)

        n_rows, _ = matS_petsc.getSize()
        Phi_petsc = PETSc.Mat().createDense([n_rows, n], comm=comm)
        Phi_petsc.setUp()

        for i in range(n):
            Phi_petsc.setColumnVector(i, singular_vectors[i])
        Phi_petsc.assemble()

        # - Clean vectors
        for vec in singular_vectors:
            vec.destroy()
        svd.destroy()

        return Phi_petsc, sigma[:n]

    def snapshotMethod(self, matS):
        """Compression method using the snapshot method on a correlation matrix

        Arguments
        ----------
        matS : PETSc.Mat
            Matrix of snapshots on which the POD operator should be applied

        Returns
        -------
        Phi : PETSc.Mat
            Reduced order basis
        singval : numpy.ndarray
            Singular values
        """
        # - Compute correlation matrix
        corrMatrix = self._CorrOperator.ptap(matS)
        # - Solve eigenproblem using SLEPc
        eps = SLEPc.EPS().create(comm=matS.getComm())
        eps.setOperators(corrMatrix)
        eps.setProblemType(SLEPc.EPS.ProblemType.HEP)
        # - Extract eigenvalues and compute basis
        nconv = eps.getConverged()
        eigenvalues = []
        basis_vectors = []
        
        _, v = corrMatrix.createVecs()
        phi_i = matS.createVecRight()
        for i in range(nconv):
            k = eps.getEigenvalue(i)

            if k.real <= self._tol_num:
                break

            eigenvalues.append(k.real)

            # - Get eigenvector and compute the mode
            eps.getEigenvector(i, v)
            sing_val = np.sqrt(k.real)

            matS.mult(v, phi_i)
            phi_i.scale(1.0/sing_val)

            basis_vectors.append(phi_i.copy())
            
        # - Assembly of the final bais  matrix from all the vectors computed before
        n_modes = len(basis_vectors)

        if n_modes == 0:
            Phi = PETSc.Mat().createDense([matS.getSize()[0], 0], comm=matS.getComm())
            Phi.assemble()
            return Phi, np.array([])
        else:
            Phi = PETSc.Mat().createDense([matS.getSize()[0], n_modes], comm=matS.getComm())
            Phi.setUp()
            for i, vec in enumerate(basis_vectors):
                Phi.setValues(range(matS.getSize()[0]), [i], vec.getArray(readonly=True), PETSc.InsertMode.INSERT_VALUES)
            Phi.assemble()
            singular_values = np.sqrt(np.array(eigenvalues))
            return Phi, singular_values

class PODAnalysis:
    def __init__(self, snapshots, **kwargs):
        """
        The constructor acts as an internal factory.
        It selects and instantiates the correct implementation."""
        if isinstance(snapshots, np.ndarray):
            self._impl = PODAnalysisNumpy(snapshots, **kwargs)
        elif isinstance(snapshots, PETSc.Mat):
            self._impl = PODAnalysisPetsc(snapshots, **kwargs)
        else:
            raise TypeError(f"The snapshot type is not supported : {type(snapshots).__name__}")

    def __getattr__(self, name):
        """Delegates all method calls to the implementation object."""
        return getattr(self._impl, name)


# class PODAnalysis:
#     """
#     Class for building a base incrementally using POD.

#     This class is designed to build a base from a set of snapshots.
#     """

#     def __init__(
#         self, snapshots, method="SVD", criterion="energy", tolerance=None, nbModes=None, CorrOp=None
#     ):
#         """
#         Initializes a PODAnalysis.

#         Arguments
#         ----------
#         snapshots : numpy.ndarray
#             Snapshot array. Each column contains a given snapshot (size = number of dofs * number of snapshots).
#         method : str
#             Data compression method used.
#         criterion : str
#             Criteria for selecting the number of modes used (energy or nbModes)
#         tolerance : str or NoneType
#             POD compression tolerance for an energy criterion (criterion=energy).
#         nbModes : str or NoneType
#             Number of modes used for a criterion where the number of modes is provided (criterion=nbModes).
#         CorrOp : numpy.ndarray or NoneType
#             Correlation operator, provided in matrix form for a snapshot approach (method=snapshot).
#         """
#         self._snapshots = snapshots  # Each column is a snapshot
#         ## - Initialisation parameterss
#         self._methodCompress = None
#         self._criterionModes = None
#         self._crit_tolerance = tolerance
#         self._crit_nbModes = nbModes
#         ## - Prepare the operators for POD Analysis
#         self.setInfosSnapshots()
#         self.setCompressionMethod(method)
#         self.setCorrelationOperator(CorrOperator=CorrOp)
#         self.setCriterionModes(criterion)
#         ## - Tests
#         self._correctionSnapshots()
#         self._runCompatibilityTests()
#         assert self._methodCompress is not None
#         assert self._criterionModes is not None

#     def _runCompatibilityTests(self):
#         """Testing compatibility between options (for arguments)"""
#         if self._methodCompress in ["GS-classical", "GS-modified"] and self._crit_tolerance is None:
#             raise ValueError("If method is of GS type, user should provide a tolerance")
#         if self._criterionModes == "energy" and self._crit_tolerance is None:
#             raise ValueError("If criterion = energy, user should provide a tolerance")
#         if self._criterionModes == "nbModes" and self._crit_nbModes is None:
#             raise ValueError("If criterion = nbModes, user should provide a number of modes")

#     def setInfosSnapshots(self):
#         """Store information about the size of the snapshot matrix"""
#         self._numberOfDOFs = self._snapshots.shape[0]
#         self._numberOfSnapshots = self._snapshots.shape[1]

#     def setCorrelationOperator(self, CorrOperator=None):
#         """Set method for the correlation operator"""
#         if CorrOperator is None:
#             self._CorrOperator = scipy.sparse.identity(self._numberOfDOFs, format="csr")
#         else:
#             self._CorrOperator = CorrOperator

#     def setCompressionMethod(self, method):
#         """Set method for the compression method"""
#         if method in POD_VALID_METHOD:
#             self._methodCompress = method
#         else:
#             raise ValueError(
#                 f"PODAnalysis: Method '{method}' is not valid. Choose method in {POD_VALID_METHOD}."
#             )

#     def setCriterionModes(self, criterion):
#         """Set method for the criterion for selecting modes"""
#         if criterion in POD_CRITERION_METHOD:
#             self._criterionModes = criterion
#         else:
#             raise ValueError(
#                 f"PODAnalysis: Method '{criterion}' is not valid. Choose method in {POD_CRITERION_METHOD}."
#             )

#     def _correctionSnapshots(self):
#         """Correct the snapshots by removing null values"""
#         tol = 1e-12
#         norms = np.linalg.norm(self._snapshots, axis=0)
#         self._snapshots = self._snapshots[:, norms > 0]

#     def getCompressionMethod(self, method):
#         """Get method for the criterion for selecting modes"""
#         return self._methodCompress

#     def getCriterionModes(self, method):
#         """Get method for the compression method"""
#         return self._criterionModes

#     def selectModes(self, Phi, singval, tolerance=None, nbModes=None):
#         """Method for selecting the number of modes given a basis and previously calculated singular values

#         Arguments
#         ----------
#         Phi : numpy.ndarray
#             Previously calculated reduced order basis (size = number of DOFs * number of modes).
#             Here number of modes should be close to the number of snapshots. Only redundant information can be removed.
#         singval : numpy.ndarray
#             Previously calculated singular values.
#         tolerance : str or NoneType
#             POD compression tolerance for an energy criterion (criterion=energy).
#         nbModes : str or NoneType
#             Number of modes used for a criterion where the number of modes is provided (criterion=nbModes).

#         Returns
#         -------
#         Phi : numpy.ndarray
#             Reduced order basis after truncation
#         singval : numpy.ndarray
#             Singular values after truncation
#         """
#         ## - Test
#         if tolerance is None and nbModes is None:
#             raise ValueError(
#                 "selectModes: Either tolerance or nbModes should be provided to the method"
#             )
#         ## - Choice of the number of modes to keep
#         if self._criterionModes == "energy":
#             s_squared = singval**2
#             sum_i = 0
#             i = 0
#             while i < len(singval) and sum_i / np.sum(s_squared) < (1 - tolerance):
#                 sum_i += s_squared[i]
#                 i += 1
#             nbModes_v = i
#         elif self._criterionModes == "nbModes":
#             if nbModes is None or nbModes > len(singval):
#                 raise ValueError("selectModes: nbModes should be given or is too big!")
#             nbModes_v = nbModes
#         else:
#             raise ValueError(
#                 f"PODAnalysis: Method '{self._criterionModes}' is not valid. Choose method in {POD_CRITERION_METHOD}."
#             )
#         ## - Truncation of the POD Basis
#         Phi_v = Phi[:, :nbModes_v]
#         singval_v = singval[:nbModes_v]
#         return Phi_v, singval_v

#     def POD(self, matS, option):
#         """Method to construct a reduced order basis by POD

#         Arguments
#         ----------
#         matS : numpy.ndarray
#             Matrix of snapshots on which the POD operator should be applied
#         option : int
#             Changes the outputs of the function. If option=1, only reduced order basis.
#             If option=2, returns reduced order basis and singular values.

#         Returns
#         -------
#         Phi : numpy.ndarray
#             Reduced order basis
#         singval : numpy.ndarray
#             Singular values (only if option=2)
#         """
#         ## - Compression step
#         if self._methodCompress == "snapshot":
#             Phi, singval = self.snapshotMethod(matS)
#         elif self._methodCompress == "SVD":
#             Phi, singval = self.SVDMethod(matS)
#         elif self._methodCompress == "GS-classical":
#             Phi, singval = self.GSMethod(matS, self._crit_tolerance, "classical")
#         elif self._methodCompress == "GS-modified":
#             Phi, singval = self.GSMethod(matS, self._crit_tolerance, "modified")
#         else:
#             raise ValueError(
#                 f"PODAnalysis: Method '{self._methodCompress}' is not implemented yet."
#             )
#         if self._methodCompress in POD_METHOD_WITHOUT_CRIT:
#             Phi_t, singval_t = Phi, singval
#         else:
#             ## - Apply truncation
#             Phi_t, singval_t = self.selectModes(
#                 Phi, singval, self._crit_tolerance, self._crit_nbModes
#             )
#         ## - Return outputs
#         if option == 1:
#             return Phi_t
#         elif option == 2:
#             return Phi_t, singval_t
#         else:
#             raise ValueError("PODAnalysis: computePODBasis should be 1 or 2.")

#     def computePODBasis(self, option=1):
#         """Method to construct a reduced order basis by POD
#         using the stored snapshots

#         Arguments
#         ----------
#         option : int
#             Changes the outputs of the function. If option=1, only reduced order basis.
#             If option=2, returns reduced order basis and singular values.

#         """
#         return self.POD(self._snapshots, option=option)

#     def computePODBasisIncremental(self, Phi, singval=None, method="HPOD"):
#         """Method to enrich a reduced order basis with the stored snapshots

#         Arguments
#         ----------
#         Phi : numpy.ndarray
#             Reduced order basis which has been previously computed
#         method : str
#             Name of the incremental approach to use
#         """
#         assert method in INCR_POD_VALID_METHOD
#         if method == "HPOD":
#             matS = self._snapshots
#             projS = np.zeros(np.shape(matS))
#             for i in range(matS.shape[1]):
#                 projS[:, i] = self.computeGSprojection(Phi, matS[:, i], "modified")
#             Phi_new = self.POD(projS, option=1)
#             return np.column_stack((Phi, Phi_new))
#         elif method == "HAPOD":
#             assert singval is not None
#             mPhi = singval * Phi
#             assert np.shape(Phi) == np.shape(mPhi)
#             mS = np.column_stack((mPhi, self._snapshots))
#             return self.POD(mS, option=1)
#         else:
#             raise ValueError(
#                 f"PODAnalysis: Method '{method}' is not valid. Choose method in {INCR_POD_VALID_METHOD}."
#             )

#     def SVDMethod(self, matS):
#         """Compression method using SVD on the snapshot matrix

#         Arguments
#         ----------
#         matS : numpy.ndarray
#             Matrix of snapshots on which the POD operator should be applied
#         Returns
#         -------
#         Phi : numpy.ndarray
#             Reduced order basis
#         singval : numpy.ndarray
#             Singular values
#         """
#         ## - Apply SVD directly on the snapshot matrix
#         U, sigma, _ = np.linalg.svd(matS, full_matrices=False)
#         ## - Order eigenvalues and compute basis
#         n = np.where(sigma == 0)[0]
#         if n.size == 0:
#             n = len(sigma)
#         else:
#             n = n[0]
#         return U[:, :n], sigma[:n]

#     def snapshotMethod(self, matS):
#         """Compression method using the snapshot method on a correlation matrix

#         Arguments
#         ----------
#         matS : numpy.ndarray
#             Matrix of snapshots on which the POD operator should be applied

#         Returns
#         -------
#         Phi : numpy.ndarray
#             Reduced order basis
#         singval : numpy.ndarray
#             Singular values
#         """
#         ## - Compute correlation matrix
#         corrMatrix = matS.T @ self._CorrOperator @ matS
#         ## - Solve eigenproblem
#         eigenvalues, eigenvectors = np.linalg.eigh(corrMatrix)
#         ## - Order eigenvalues and compute basis
#         idx = np.argsort(eigenvalues)[::-1]
#         eigenvalues = eigenvalues[idx]
#         eigenvalues = np.where(eigenvalues < 0, 0, eigenvalues)
#         eigenvectors = eigenvectors[:, idx]

#         n = np.where(eigenvalues == 0)[0]
#         if n.size == 0:
#             n = len(eigenvalues)
#         else:
#             n = n[0]
#         # - Return reduced order basis and singular values
#         return np.dot(matS, eigenvectors[:, :n]) / np.sqrt(eigenvalues[:n]), np.sqrt(
#             eigenvalues[:n]
#         )

#     def GSMethod(self, matS, tole, methodGS):
#         """Compression method using a Gram-Schmidt process

#         Arguments
#         ----------
#         matS : numpy.ndarray
#             Matrix of snapshots on which the POD operator should be applied
#         tole : float
#             Tolerance used for the test when adding new snapshot
#         methodGS : str
#             Should be "classical" or "modified" = GS method applied

#         Returns
#         -------
#         Phi : numpy.ndarray
#             Reduced order basis
#         singval : numpy.ndarray
#             Singular values
#         """
#         s_0_norm = np.linalg.norm(matS[:, 0])
#         Phi = matS[:, 0:1] / s_0_norm
#         singval = np.array([s_0_norm])
#         n_snap = matS.shape[1]
#         for i in range(1, n_snap):
#             Phi, singval = self.updateGStype(Phi, singval, matS[:, i : i + 1], tole, methodGS)
#         return Phi, singval

#     def computeGSprojection(self, Phi, s_new, methodGS):
#         assert methodGS in GS_METHOD
#         ## - GS orthogonalisation
#         n_modes = Phi.shape[1]
#         for kp in range(2):  # Kahan-Parlett process
#             if methodGS == "classical":
#                 s_new_loc = s_new
#                 for k in range(n_modes):
#                     s_new = s_new - np.dot(s_new_loc, Phi[:, k]) * Phi[:, k]
#             if methodGS == "modified":
#                 for k in range(n_modes):
#                     s_new = s_new - np.dot(s_new, Phi[:, k]) * Phi[:, k]
#         return s_new

#     def updateGStype(self, Phi, singval, snapshot_new, tole, methodGS):
#         """Update a basis with a new snapshot method using a Gram-Schmidt process

#         Arguments
#         ----------
#         Phi : numpy.ndarray
#             Basis to enrich
#         singval : numpy.ndarray
#             Singular values
#         snapshot_new : numpy.ndarray
#             Snapshot to add
#         tole : float
#             Tolerance used for the test when adding new snapshot
#         methodGS : str
#             Should be "classical" or "modified" = GS method applied

#         Returns
#         -------
#         Phi : numpy.ndarray
#             Reduced order basis
#         singval : numpy.ndarray
#             Singular values
#         """
#         assert methodGS in GS_METHOD
#         ## - Check that the added snapshot is 1D
#         s_new = snapshot_new.flatten()
#         s_new_norm = np.linalg.norm(snapshot_new)
#         ## - GS orthogonalisation
#         s_new = self.computeGSprojection(Phi, s_new, methodGS)
#         s_new_perp_norm = np.linalg.norm(s_new)
#         if s_new_perp_norm > tole * s_new_norm:
#             Phi = np.column_stack((Phi, s_new / s_new_perp_norm))
#             singval = np.hstack([singval, s_new_perp_norm])
#         return Phi, singval

#     def computeDecayRate(self, singval):
#         """Compute a decay rate of a list of singular values

#         Parameters
#         ----------
#         singval : numpy.ndarray
#             Singular values (order in a decreasing manner)

#         Returns
#         -------
#         Phi : float
#             Decay rate
#         """
#         nbSing = len(singval)

#         N = np.log(np.arange(1, nbSing + 1))
#         Y = np.log(singval)

#         sum_N = np.sum(N)
#         sum_Y = np.sum(Y)
#         sum2_N = np.sum(N**2)

#         decayRate = (nbSing * np.dot(Y.T, N) - sum_N * sum_Y) / (nbSing * sum2_N - sum_N**2)
#         return decayRate


# """
# On utilise notamment https://slepc.upv.es/release/slepc4py/demo/ex10.html
# """


# class PODAnalysisPETSC:
#     """
#     Class for building a base incrementally using POD.

#     This class is designed to build a base from a set of snapshots.
#     """

#     def __init__(
#         self, snapshots, method="SVD", criterion="energy", tolerance=None, nbModes=None, CorrOp=None
#     ):
#         """
#         Initializes a PODAnalysis.

#         Arguments
#         ----------
#         snapshots : numpy.ndarray
#             Snapshot array. Each column contains a given snapshot (size = number of dofs * number of snapshots).
#         method : str
#             Data compression method used.
#         criterion : str
#             Criteria for selecting the number of modes used (energy or nbModes)
#         tolerance : str or NoneType
#             POD compression tolerance for an energy criterion (criterion=energy).
#         nbModes : str or NoneType
#             Number of modes used for a criterion where the number of modes is provided (criterion=nbModes).
#         CorrOp : numpy.ndarray or NoneType
#             Correlation operator, provided in matrix form for a snapshot approach (method=snapshot).
#         """
#         self._snapshots = snapshots  # Each column is a snapshot
#         self._comm = self._snapshots.getComm()
#         ## - Initialisation parameterss
#         self._methodCompress = None
#         self._criterionModes = None
#         self._crit_tolerance = tolerance
#         self._crit_nbModes = nbModes
#         ## - Prepare the operators for POD Analysis
#         self.setInfosSnapshots()
#         self.setCompressionMethod(method)
#         self.setCorrelationOperator(CorrOperator=CorrOp)
#         self.setCriterionModes(criterion)
#         ## - Tests
#         self._correctionSnapshots()
#         self._runCompatibilityTests()
#         assert self._methodCompress is not None
#         assert self._criterionModes is not None

#     def _runCompatibilityTests(self):
#         """Testing compatibility between options (for arguments)"""
#         if self._methodCompress in ["GS-classical", "GS-modified"] and self._crit_tolerance is None:
#             raise ValueError("If method is of GS type, user should provide a tolerance")
#         if self._criterionModes == "energy" and self._crit_tolerance is None:
#             raise ValueError("If criterion = energy, user should provide a tolerance")
#         if self._criterionModes == "nbModes" and self._crit_nbModes is None:
#             raise ValueError("If criterion = nbModes, user should provide a number of modes")

#     def setInfosSnapshots(self):
#         """Store information about the size of the snapshot matrix"""
#         self._numberOfDOFs = self._snapshots.shape[0]
#         self._numberOfSnapshots = self._snapshots.shape[1]

#     def setCorrelationOperator(self, CorrOperator=None):
#         """Set method for the correlation operator"""
#         if CorrOperator is None:
#             I = PETSc.Mat().create(self._comm)
#             I.setSizes(((None, self._numberOfDOFs), (None, self._numberOfDOFs)))
#             I.setType("aij")
#             I.setPreallocationNNZ(1)
#             # diag_vec = PETSc.Vec().create(self._comm)
#             # diag_vec.setSizes((None, self._numberOfDOFs))
#             # diag_vec.setUp()
#             diag_vec = I.createVecLeft()
#             diag_vec.set(1.0)

#             I.setDiagonal(diag_vec)
#             diag_vec.destroy()
#             I.assemble()
#             self._CorrOperator = I
#         else:
#             self._CorrOperator = CorrOperator

#     def setCompressionMethod(self, method):
#         """Set method for the compression method"""
#         if method in POD_VALID_METHOD:
#             self._methodCompress = method
#         else:
#             raise ValueError(
#                 f"PODAnalysis: Method '{method}' is not valid. Choose method in {POD_VALID_METHOD}."
#             )

#     def setCriterionModes(self, criterion):
#         """Set method for the criterion for selecting modes"""
#         if criterion in POD_CRITERION_METHOD:
#             self._criterionModes = criterion
#         else:
#             raise ValueError(
#                 f"PODAnalysis: Method '{criterion}' is not valid. Choose method in {POD_CRITERION_METHOD}."
#             )

#     def _correctionSnapshots(self):
#         """Correct the snapshots by removing null values"""
#         tol = 1e-12
#         norms = np.linalg.norm(self._snapshots, axis=0)
#         self._snapshots = self._snapshots[:, norms > 0]

#     def getCompressionMethod(self, method):
#         """Get method for the criterion for selecting modes"""
#         return self._methodCompress

#     def getCriterionModes(self, method):
#         """Get method for the compression method"""
#         return self._criterionModes

#     def selectModes(self, Phi, singval, tolerance=None, nbModes=None):
#         """Method for selecting the number of modes given a basis and previously calculated singular values

#         Arguments
#         ----------
#         Phi : numpy.ndarray
#             Previously calculated reduced order basis (size = number of DOFs * number of modes).
#             Here number of modes should be close to the number of snapshots. Only redundant information can be removed.
#         singval : numpy.ndarray
#             Previously calculated singular values.
#         tolerance : str or NoneType
#             POD compression tolerance for an energy criterion (criterion=energy).
#         nbModes : str or NoneType
#             Number of modes used for a criterion where the number of modes is provided (criterion=nbModes).

#         Returns
#         -------
#         Phi : numpy.ndarray
#             Reduced order basis after truncation
#         singval : numpy.ndarray
#             Singular values after truncation
#         """
#         ## - Test
#         if tolerance is None and nbModes is None:
#             raise ValueError(
#                 "selectModes: Either tolerance or nbModes should be provided to the method"
#             )
#         ## - Choice of the number of modes to keep
#         if self._criterionModes == "energy":
#             s_squared = singval**2
#             sum_i = 0
#             i = 0
#             while i < len(singval) and sum_i / np.sum(s_squared) < (1 - tolerance):
#                 sum_i += s_squared[i]
#                 i += 1
#             nbModes_v = i
#         elif self._criterionModes == "nbModes":
#             if nbModes is None or nbModes > len(singval):
#                 raise ValueError("selectModes: nbModes should be given or is too big!")
#             nbModes_v = nbModes
#         else:
#             raise ValueError(
#                 f"PODAnalysis: Method '{self._criterionModes}' is not valid. Choose method in {POD_CRITERION_METHOD}."
#             )
#         ## - Truncation of the POD Basis
#         Phi_v = Phi[:, :nbModes_v]
#         singval_v = singval[:nbModes_v]
#         return Phi_v, singval_v

#     def POD(self, matS, option):
#         """Method to construct a reduced order basis by POD

#         Arguments
#         ----------
#         matS : numpy.ndarray
#             Matrix of snapshots on which the POD operator should be applied
#         option : int
#             Changes the outputs of the function. If option=1, only reduced order basis.
#             If option=2, returns reduced order basis and singular values.

#         Returns
#         -------
#         Phi : numpy.ndarray
#             Reduced order basis
#         singval : numpy.ndarray
#             Singular values (only if option=2)
#         """
#         ## - Compression step
#         if self._methodCompress == "snapshot":
#             Phi, singval = self.snapshotMethod(matS)
#         elif self._methodCompress == "SVD":
#             Phi, singval = self.SVDMethod(matS)
#         elif self._methodCompress == "GS-classical":
#             Phi, singval = self.GSMethod(matS, self._crit_tolerance, "classical")
#         elif self._methodCompress == "GS-modified":
#             Phi, singval = self.GSMethod(matS, self._crit_tolerance, "modified")
#         else:
#             raise ValueError(
#                 f"PODAnalysis: Method '{self._methodCompress}' is not implemented yet."
#             )
#         if self._methodCompress in POD_METHOD_WITHOUT_CRIT:
#             Phi_t, singval_t = Phi, singval
#         else:
#             ## - Apply truncation
#             Phi_t, singval_t = self.selectModes(
#                 Phi, singval, self._crit_tolerance, self._crit_nbModes
#             )
#         ## - Return outputs
#         if option == 1:
#             return Phi_t
#         elif option == 2:
#             return Phi_t, singval_t
#         else:
#             raise ValueError("PODAnalysis: computePODBasis should be 1 or 2.")

#     def computePODBasis(self, option=1):
#         """Method to construct a reduced order basis by POD
#         using the stored snapshots

#         Arguments
#         ----------
#         option : int
#             Changes the outputs of the function. If option=1, only reduced order basis.
#             If option=2, returns reduced order basis and singular values.

#         """
#         return self.POD(self._snapshots, option=option)

#     def computePODBasisIncremental(self, Phi, singval=None, method="HPOD"):
#         """Method to enrich a reduced order basis with the stored snapshots

#         Arguments
#         ----------
#         Phi : numpy.ndarray
#             Reduced order basis which has been previously computed
#         method : str
#             Name of the incremental approach to use
#         """
#         assert method in INCR_POD_VALID_METHOD
#         if method == "HPOD":
#             matS = self._snapshots
#             projS = np.zeros(np.shape(matS))
#             for i in range(matS.shape[1]):
#                 projS[:, i] = self.computeGSprojection(Phi, matS[:, i], "modified")
#             Phi_new = self.POD(projS, option=1)
#             return np.column_stack((Phi, Phi_new))
#         elif method == "HAPOD":
#             assert singval is not None
#             mPhi = singval * Phi
#             assert np.shape(Phi) == np.shape(mPhi)
#             mS = np.column_stack((mPhi, self._snapshots))
#             return self.POD(mS, option=1)
#         else:
#             raise ValueError(
#                 f"PODAnalysis: Method '{method}' is not valid. Choose method in {INCR_POD_VALID_METHOD}."
#             )

#     # def SVDMethod(self, matS):
#     #     """Compression method using SVD on the snapshot matrix

#     #     Arguments
#     #     ----------
#     #     matS : numpy.ndarray
#     #         Matrix of snapshots on which the POD operator should be applied
#     #     Returns
#     #     -------
#     #     Phi : numpy.ndarray
#     #         Reduced order basis
#     #     singval : numpy.ndarray
#     #         Singular values
#     #     """
#     #     ## - Apply SVD directly on the snapshot matrix
#     #     U, sigma, _ = np.linalg.svd(matS, full_matrices=False)
#     #     ## - Order eigenvalues and compute basis
#     #     n = np.where(sigma == 0)[0]
#     #     if n.size == 0:
#     #         n = len(sigma)
#     #     else:
#     #         n = n[0]
#     #     return U[:, :n], sigma[:n]

#     def SVDMethod(self, matS, verbose=True):
#         tol = 1e-12
#         comm = matS.getComm()
#         Print = PETSC.Sys.Print

#         # - Creation and configuration of the SVD solver
#         svd = SLEPc.SVD().create(comm=comm)
#         svd.setOperator(matS)
#         S.setType(S.Type.TRLANCZOS)

#         # - Ask to compute all the needed singular values
#         _, n_cols = matS.getSize()
#         svd.setDimensions(nsv=n_cols)
#         S.setFromOptions()
#         # - Resolution
#         svd.solve()

#         n_conv = svd.getConverged()
#         if verbose and comm.rank == 0:
#             Print("******************************")
#             Print("*** SLEPc SVD Solution Results ***")
#             Print("******************************\n")
#             svd_type = svd.getType()
#             Print(f"Solution method: {svd_type}")
#             its = svd.getIterationNumber()
#             Print(f"Number of iterations: {its}")
#             nsv, _, _ = svd.getDimensions()
#             Print(f"Number of requested singular values: {nsv}")
#             tol, maxit = svd.getTolerances()
#             Print(f"Stopping condition: tol={tol:.4g}, maxit={maxit}")
#             Print(f"Number of converged singular triplets: {n_conv}\n")

#         # - Construction of the outputs
#         singular_values = []
#         singular_vectors = []

#         if n_conv > 0:
#             u_tmp, v_tmp = matS.createVecs()

#             if verbose and comm.rank == 0:
#                 Print("    sigma       residual norm ")
#                 Print("-------------  ---------------")

#             for i in range(n_conv):
#                 sigma_i = svd.getSingularTriplet(i, u_tmp, v_tmp)
#                 singular_values.append(sigma_i)
#                 singular_vectors.append(u_tmp.copy())

#                 if verbose and comm.rank == 0:
#                     error = svd.computeError(i)
#                     Print(f"   {sigma_i:6f}     {error:12g}")

#             if verbose and comm.rank == 0:
#                 Print()

#             u_tmp.destroy()
#             v_tmp.destroy()

#         sigma = np.array(singular_values)
#         zero_indices = np.where(sigma < tol)[0]
#         n = zero_indices[0] if zero_indices.size > 0 else len(sigma)

#         n_rows, _ = matS_petsc.getSize()
#         Phi_petsc = PETSc.Mat().createDense([n_rows, n], comm=comm)
#         Phi_petsc.setUp()

#         for i in range(n):
#             Phi_petsc.setColumnVector(i, singular_vectors[i])
#         Phi_petsc.assemble()

#         # - Clean vectors
#         for vec in singular_vectors:
#             vec.destroy()
#         svd.destroy()

#         return Phi_petsc, sigma[:n]

#     def snapshotMethod(self, matS):
#         """Compression method using the snapshot method on a correlation matrix

#         Arguments
#         ----------
#         matS : numpy.ndarray
#             Matrix of snapshots on which the POD operator should be applied

#         Returns
#         -------
#         Phi : numpy.ndarray
#             Reduced order basis
#         singval : numpy.ndarray
#             Singular values
#         """
#         ## - Compute correlation matrix
#         corrMatrix = matS.T @ self._CorrOperator @ matS
#         ## - Solve eigenproblem
#         eigenvalues, eigenvectors = np.linalg.eigh(corrMatrix)
#         ## - Order eigenvalues and compute basis
#         idx = np.argsort(eigenvalues)[::-1]
#         eigenvalues = eigenvalues[idx]
#         eigenvalues = np.where(eigenvalues < 0, 0, eigenvalues)
#         eigenvectors = eigenvectors[:, idx]

#         n = np.where(eigenvalues == 0)[0]
#         if n.size == 0:
#             n = len(eigenvalues)
#         else:
#             n = n[0]
#         # - Return reduced order basis and singular values
#         return np.dot(matS, eigenvectors[:, :n]) / np.sqrt(eigenvalues[:n]), np.sqrt(
#             eigenvalues[:n]
#         )

#     def snapshotMethod(self, matS):
#         ## - Compute correlation matrix
#         corrMatrix = self._CorrOperator.ptap(matS)
#         ## - Solve eigenproblem using SLEPc
#         eps = SLEPc.EPS().create(comm=matS.getComm())
#         eps.setOperators(corrMatrix)
#         eps.setProblemType(SLEPc.EPS.ProblemType.HEP)

#     def GSMethod(self, matS, tole, methodGS):
#         """Compression method using a Gram-Schmidt process

#         Arguments
#         ----------
#         matS : numpy.ndarray
#             Matrix of snapshots on which the POD operator should be applied
#         tole : float
#             Tolerance used for the test when adding new snapshot
#         methodGS : str
#             Should be "classical" or "modified" = GS method applied

#         Returns
#         -------
#         Phi : numpy.ndarray
#             Reduced order basis
#         singval : numpy.ndarray
#             Singular values
#         """
#         s_0_norm = np.linalg.norm(matS[:, 0])
#         Phi = matS[:, 0:1] / s_0_norm
#         singval = np.array([s_0_norm])
#         n_snap = matS.shape[1]
#         for i in range(1, n_snap):
#             Phi, singval = self.updateGStype(Phi, singval, matS[:, i : i + 1], tole, methodGS)
#         return Phi, singval

#     def computeGSprojection(self, Phi, s_new, methodGS):
#         assert methodGS in GS_METHOD
#         ## - GS orthogonalisation
#         n_modes = Phi.shape[1]
#         for kp in range(2):  # Kahan-Parlett process
#             if methodGS == "classical":
#                 s_new_loc = s_new
#                 for k in range(n_modes):
#                     s_new = s_new - np.dot(s_new_loc, Phi[:, k]) * Phi[:, k]
#             if methodGS == "modified":
#                 for k in range(n_modes):
#                     s_new = s_new - np.dot(s_new, Phi[:, k]) * Phi[:, k]
#         return s_new

#     def updateGStype(self, Phi, singval, snapshot_new, tole, methodGS):
#         """Update a basis with a new snapshot method using a Gram-Schmidt process

#         Arguments
#         ----------
#         Phi : numpy.ndarray
#             Basis to enrich
#         singval : numpy.ndarray
#             Singular values
#         snapshot_new : numpy.ndarray
#             Snapshot to add
#         tole : float
#             Tolerance used for the test when adding new snapshot
#         methodGS : str
#             Should be "classical" or "modified" = GS method applied

#         Returns
#         -------
#         Phi : numpy.ndarray
#             Reduced order basis
#         singval : numpy.ndarray
#             Singular values
#         """
#         assert methodGS in GS_METHOD
#         ## - Check that the added snapshot is 1D
#         s_new = snapshot_new.flatten()
#         s_new_norm = np.linalg.norm(snapshot_new)
#         ## - GS orthogonalisation
#         s_new = self.computeGSprojection(Phi, s_new, methodGS)
#         s_new_perp_norm = np.linalg.norm(s_new)
#         if s_new_perp_norm > tole * s_new_norm:
#             Phi = np.column_stack((Phi, s_new / s_new_perp_norm))
#             singval = np.hstack([singval, s_new_perp_norm])
#         return Phi, singval

#     def computeDecayRate(self, singval):
#         """Compute a decay rate of a list of singular values

#         Parameters
#         ----------
#         singval : numpy.ndarray
#             Singular values (order in a decreasing manner)

#         Returns
#         -------
#         Phi : float
#             Decay rate
#         """
#         nbSing = len(singval)

#         N = np.log(np.arange(1, nbSing + 1))
#         Y = np.log(singval)

#         sum_N = np.sum(N)
#         sum_Y = np.sum(Y)
#         sum2_N = np.sum(N**2)

#         decayRate = (nbSing * np.dot(Y.T, N) - sum_N * sum_Y) / (nbSing * sum2_N - sum_N**2)
#         return decayRate
