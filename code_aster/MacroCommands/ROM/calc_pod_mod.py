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
import scipy.sparse
import numpy as np


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


# CLASS DEFINITION FOR A POD ANALYSIS
POD_VALID_METHOD = ["SVD", "snapshot", "GS-classical", "GS-modified"]
POD_METHOD_WITHOUT_CRIT = ["GS-classical", "GS-modified"]
assert all(item in POD_VALID_METHOD for item in POD_METHOD_WITHOUT_CRIT)
POD_CRITERION_METHOD = ["energy", "nbModes"]


class PODAnalysis:
    """
    Class for building a base incrementally using POD.

    This class is designed to build a base from a set of snapshots.
    """

    def __init__(
        self, snapshots, method="SVD", criterion="energy", tolerance=None, nbModes=None, CorrOp=None
    ):
        """
        Initializes a PODAnalysis.

        Arguments
        ----------
        snapshots : numpy.ndarray
            Snapshot array. Each column contains a given snapshot (size = number of dofs * number of snapshots).
        method : str
            Data compression method used.
        criterion : str
            Criteria for selecting the number of modes used (energy or nbModes)
        tolerance : str or NoneType
            POD compression tolerance for an energy criterion (criterion=energy).
        nbModes : str or NoneType
            Number of modes used for a criterion where the number of modes is provided (criterion=nbModes).
        CorrOp : numpy.ndarray or NoneType
            Correlation operator, provided in matrix form for a snapshot approach (method=snapshot).
        """
        self._snapshots = snapshots  # Each column is a snapshot
        ## - Initialisation parameterss
        self._methodCompress = None
        self._criterionModes = None
        self._crit_tolerance = tolerance
        self._crit_nbModes = nbModes
        ## - Prepare the operators for POD Analysis
        self.setInfosSnapshots()
        self.setCompressionMethod(method)
        self.setCorrelationOperator(CorrOperator=CorrOp)
        self.setCriterionModes(criterion)
        ## - Tests
        self._correctionSnapshots()
        self._runCompatibilityTests()
        assert self._methodCompress is not None
        assert self._criterionModes is not None

    def _runCompatibilityTests(self):
        """Testing compatibility between options (for arguments)"""
        if self._methodCompress in ["GS-classical", "GS-modified"] and self._crit_tolerance is None:
            raise ValueError("If method is of GS type, user should provide a tolerance")
        if self._criterionModes == "energy" and self._crit_tolerance is None:
            raise ValueError("If criterion = energy, user should provide a tolerance")
        if self._criterionModes == "nbModes" and self._crit_nbModes is None:
            raise ValueError("If criterion = nbModes, user should provide a number of modes")

    def setInfosSnapshots(self):
        """Store information about the size of the snapshot matrix"""
        self._numberOfDOFs = self._snapshots.shape[0]
        self._numberOfSnapshots = self._snapshots.shape[1]

    def setCorrelationOperator(self, CorrOperator=None):
        """Set method for the correlation operator"""
        if CorrOperator is None:
            self._CorrOperator = scipy.sparse.identity(self._numberOfDOFs, format="csr")
        else:
            self._CorrOperator = CorrOperator

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

    def _correctionSnapshots(self):
        tol = 1e-12
        norms = np.linalg.norm(self._snapshots, axis=0)
        self._snapshots = self._snapshots[:, norms > 0]

    def getCompressionMethod(self, method):
        """Get method for the criterion for selecting modes"""
        return self._methodCompress

    def getCriterionModes(self, method):
        """Get method for the compression method"""
        return self._criterionModes

    def selectModes(self, Phi, singval, tolerance=None, nbModes=None):
        """Method for selecting the number of modes given a basis and previously calculated singular values

        Arguments
        ----------
        Phi : numpy.ndarray
            Previously calculated reduced order basis (size = number of DOFs * number of modes).
            Here number of modes should be close to the number of snapshots. Only redundant information can be removed.
        singval : numpy.ndarray
            Previously calculated singular values.
        tolerance : str or NoneType
            POD compression tolerance for an energy criterion (criterion=energy).
        nbModes : str or NoneType
            Number of modes used for a criterion where the number of modes is provided (criterion=nbModes).

        Returns
        -------
        Phi : numpy.ndarray
            Reduced order basis after truncation
        singval : numpy.ndarray
            Singular values after truncation
        """
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
        Phi_v = Phi[:, :nbModes_v]
        singval_v = singval[:nbModes_v]
        return Phi_v, singval_v

    def computePODBasis(self, option=1):
        """Method to construct a reduced order basis

        Arguments
        ----------
        option : int
            Changes the outputs of the function. If option=1, only reduced order basis.
            If option=2, returns reduced order basis and singular values.

        Returns
        -------
        Phi : numpy.ndarray
            Reduced order basis
        singval : numpy.ndarray
            Singular values (only if option=2)
        """
        ## - Compression step
        if self._methodCompress == "snapshot":
            Phi, singval = self.snapshotMethod()
        elif self._methodCompress == "SVD":
            Phi, singval = self.SVDMethod()
        elif self._methodCompress == "GS-classical":
            Phi, singval = self.GSmethod(self._crit_tolerance, "classical")
        elif self._methodCompress == "GS-modified":
            Phi, singval = self.GSmethod(self._crit_tolerance, "modified")
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

    def SVDMethod(self):
        """Compression method using SVD on the snapshot matrix

        Returns
        -------
        Phi : numpy.ndarray
            Reduced order basis
        singval : numpy.ndarray
            Singular values (only if option=2)
        """
        S = self._snapshots
        ## - Apply SVD directly on the snapshot matrix
        U, sigma, _ = np.linalg.svd(S, full_matrices=False)
        ## - Order eigenvalues and compute basis
        n = np.where(sigma == 0)[0]
        if n.size == 0:
            n = len(sigma)
        else:
            n = n[0]
        return U[:, :n], sigma[:n]

    def snapshotMethod(self):
        """Compression method using the snapshot method on a correlation matrix

        Returns
        -------
        Phi : numpy.ndarray
            Reduced order basis
        singval : numpy.ndarray
            Singular values (only if option=2)
        """
        ## - Compute correlation matrix
        S = self._snapshots
        corrMatrix = S.T @ self._CorrOperator @ S
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
        return np.dot(S, eigenvectors[:, :n]) / np.sqrt(eigenvalues[:n]), np.sqrt(eigenvalues[:n])

    def GSmethod(self, tole, methodGS):
        s_0_norm = np.linalg.norm(self._snapshots[:, 0])
        Phi = self._snapshots[:, 0:1] / s_0_norm
        singval = np.array([s_0_norm])
        n_snap = self._snapshots.shape[1]
        for i in range(1, n_snap):
            Phi, singval = self.update_GStype(
                Phi, singval, self._snapshots[:, i : i + 1], tole, methodGS
            )
        return Phi, singval

    def update_GStype(self, Phi, singval, snapshot_new, tole, methodGS):
        ## - Check that the added snapshot is 1D
        s_new = snapshot_new.flatten()
        s_new_norm = np.linalg.norm(snapshot_new)
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
        ## - Check the relevance of the new information
        s_new_perp_norm = np.linalg.norm(s_new)
        if s_new_perp_norm > tole * s_new_norm:
            Phi = np.column_stack((Phi, s_new / s_new_perp_norm))
            singval = np.hstack([singval, s_new_perp_norm])
        return Phi, singval

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
