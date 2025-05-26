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

import copy
import re
from math import sqrt

from ...Objects import DiscreteComputation
from ...Utilities import MPI, logger, no_new_attributes, profile
from ..Basics import ContextMixin


class ConvergenceManager(ContextMixin):
    """Object that decides about the convergence status."""

    __needs__ = ("problem", "state", "keywords")

    _param = _residual_reference = None
    __setattr__ = no_new_attributes(object.__setattr__)

    class UnDefined:
        """Used to identify undefined values."""

        def __repr__(self):
            return "UNDEF"

    undef = UnDefined()

    class Parameter:
        """A parameter to be checked for convergence.

        Arguments:
            reference (float|int): Reference value.
        """

        match = None
        _refe = _value = _minValue = None
        __setattr__ = no_new_attributes(object.__setattr__)

        @classmethod
        def factory(cls, name, reference):
            """Build a parameter object.

            Arguments:
                name (str): Parameter name.
                reference (float|int): Reference value.
            """
            for subclass in cls.__subclasses__():
                if subclass.match(name):
                    return subclass(reference)
            raise TypeError(f"unknown parameter: {name!r}")

        def __init__(self, reference):
            self._refe = reference
            self._value = ConvergenceManager.undef
            self._minValue = ConvergenceManager.undef

        def __repr__(self) -> str:
            minV = ""
            if self.minSet():
                minV = f", min: {self._minValue}"
            return f"{self._value} (ref.: {self._refe}{minV})"

        def isSet(self):
            """Tell if the parameter value has been assigned.

            Returns:
                bool: *True* if the parameter has a value, *False* if not.
            """
            return self._value is not ConvergenceManager.undef

        def minSet(self):
            """Tell if a minimal value has been assigned.

            Returns:
                bool: *True* if the parameter has a minimum, *False* if not.
            """
            return self._minValue is not ConvergenceManager.undef

        def hasRef(self):
            """Tell if a reference value is defined for the parameter.

            Returns:
                bool: *True* if the parameter has a reference value, *False* if not.
            """
            return self._refe is not ConvergenceManager.undef

        def reset(self):
            """Reset the parameter value."""
            self._value = ConvergenceManager.undef

        @property
        def reference(self):
            """float|int: Reference value of the parameter."""
            return self._refe

        @property
        def value(self):
            """float|int: Current value of the parameter."""
            return self._value

        @value.setter
        def value(self, value):
            """Set the parameter value.

            Arguments:
                value (float|int): Parameter value.
            """
            self._value = value

        @property
        def minValue(self):
            """float|int: Current minimal value of the parameter."""
            return self._minValue

        @minValue.setter
        def minValue(self, value):
            """Set the parameter value.

            Arguments:
                value (float|int): Parameter value.
            """
            self._minValue = value

        def isConverged(self):
            """Tell if the current value is converged.

            Returns:
                bool: *True* if the value is converged, *False* otherwise.
            """
            raise NotImplementedError("must be subclassed!")

        def isFinished(self):
            """Tell if the current parameter should stop the calculation.

            Returns:
                bool: *True* if the calculation should be stopped, *False* otherwise.
            """
            raise NotImplementedError("must be subclassed!")

    class ResidualParameter(Parameter):
        """Type of *Parameter* for a residual.

        Arguments:
            name (str): Parameter name.
            reference (float|int): Reference value.
        """

        match = re.compile("^RESI_").search

        def __init__(self, reference):
            super().__init__(reference)
            self.minValue = 0.0

        def isConverged(self):
            """Tell if the current value is converged.

            Returns:
                bool: *True* if the value is converged, *False* otherwise.
            """
            if not self.hasRef():
                return True
            if not self.isSet():
                return not self.hasRef()
            checkMin = not self.minSet() or self._minValue <= self._value
            return checkMin and self._value <= self._refe

        def isFinished(self):
            """Tell if the current parameter should stop the calculation.

            Returns:
                bool: *True* if the calculation should be stopped, *False* otherwise.
            """
            return self.isConverged()

    class IterationParameter(Parameter):
        """Type of *Parameter* for a number of iteration.

        Arguments:
            name (str): Parameter name.
            reference (float|int): Reference value.
        """

        match = re.compile("^ITER_").search

        def isConverged(self):
            """The number of iteration is not a convergence criteria.
            but it can nullify the convergence if it is less than a minimal value.

            Returns:
                bool: *True*.
            """

            if not self.hasRef() or not self.isSet() or not self.minSet():
                return True
            return self._minValue <= self._value

        def isPrediction(self):
            """Return True if self._value<1

            Returns:
                bool: *True*.
            """
            if not self.hasRef() or not self.isSet() or not self.minSet():
                return False
            return self._minValue > self._value

        def isFinished(self):
            """Tell if the current parameter should stop the calculation.

            Returns:
                bool: *True* if the calculation should be stopped, *False* otherwise.
            """
            return self.isSet() and self._value >= self._refe

    class ReferenceParameter(Parameter):
        """reerences for RESI_REFE_RELA"""

        match = re.compile("_REFE$").search

        def isConverged(self):
            return True

        def isFinished(self):
            return True

    @classmethod
    def builder(cls, context):
        """Default builder for :py:class:`ContextMixin` object.
        Should be subclassed for non trivial constructor.

        Args:
            context (Context): Context of the problem.

        Returns:
            instance: New object.
        """
        instance = super().builder(context)
        for crit in ("RESI_GLOB_RELA", "RESI_GLOB_MAXI", "ITER_GLOB_MAXI", "RESI_REFE_RELA"):
            value = instance.get_keyword("CONVERGENCE", crit)
            if value is not None:
                instance.setdefault(crit, value)
        if instance.get_keyword("CONVERGENCE", "RESI_REFE_RELA"):
            for crit in [
                x for x in instance.get_keyword("CONVERGENCE").keys() if x.endswith("_REFE")
            ]:
                value = instance.get_keyword("CONVERGENCE", crit)
                instance.setdefault(crit, value)
        if instance.get_keyword("CONTACT"):
            instance.setdefault("RESI_GEOM", instance.get_keyword("CONTACT", "RESI_GEOM"))
        return instance

    def __init__(self):
        super().__init__()
        self._param = {}

    def initialize(self, *mandatory):
        """Initialize the object for a new iteration.

        Mandatory parameters must be defined before checking their convergency
        (here they are assigned to a negative value, not converged).

        Arguments:
            mandatory (tuple): Name of parameters that are mandatory.
        """
        for para in self._param.values():
            para.reset()
        for name in mandatory:
            para = self._param.get(name)
            if para:
                para.value = -1

    def hasResidual(self):
        """Tell if there is at least one residual convergence parameter.

        Returns:
            bool: *True* if at least one parameter is defined, *False* otherwise.
        """
        return bool(self._residuals)

    def setdefault(self, name, reference=undef):
        """Add a convergence parameter if it does not yet exist.

        Arguments:
            name (str): Name of the parameter.
            reference (float): Reference value of the parameter.

        Returns:
            *Parameter*: The existing parameter if it already exists or a newly
            created one.
        """
        if name not in self._param:
            self._param[name] = ConvergenceManager.Parameter.factory(name, reference)
        return self._param[name]

    def get(self, name):
        """Return the value of a convergence parameter.

        Arguments:
            name (str): Name of the parameter.

        Returns:
            float|int: Parameter value or *undef* it the parameter is not defined.
        """
        if name not in self._param:
            return ConvergenceManager.undef
        return self._param.get(name).value

    @property
    def _residuals(self):
        return [
            (name, para)
            for name, para in self._param.items()
            if isinstance(para, ConvergenceManager.ResidualParameter)
        ]

    @property
    def _iterations(self):
        return [
            (name, para)
            for name, para in self._param.items()
            if isinstance(para, ConvergenceManager.IterationParameter)
        ]

    @property
    def _references(self):
        return [
            (name, para)
            for name, para in self._param.items()
            if isinstance(para, ConvergenceManager.ReferenceParameter)
        ]

    def getParameters(self):
        """Return a copy of the parameters with their current values.

        Returns:
            dict: Dict of parameters
        """
        return copy.deepcopy(self._param)

    @profile
    def getDirichletResidual(self, residual):
        """Return the residual with Dirichlet imposed values.

        Arguments:
            residual (FieldOnNodesReal): Residual.

        Returns:
            FieldOnNodesReal: Residual changed in place.
        """
        loads = self.problem.getListOfLoads()

        # maybe not really efficient
        if loads.hasDirichletBC():
            disc_comp = DiscreteComputation(self.problem)
            diriBCs = disc_comp.getIncrementalDirichletBC(
                self.state.time_curr, self.state.primal_curr
            )
            eliminatedDofs = self.problem.getDirichletBCDOFs()
            nbElimination = len(eliminatedDofs)
            assert residual.size() == nbElimination

            residual.updateValuePointers()
            diriBCs.updateValuePointers()
            for ieq in range(nbElimination):
                if eliminatedDofs[ieq] == 1:
                    residual[ieq] = diriBCs[ieq]

        return residual

    @profile
    def getResidualReference(self):
        if self._residual_reference is not None:
            return self._residual_reference
        disc_comp = DiscreteComputation(self.problem)
        vale_by_name = {}
        for name, reference in self._references:
            vale_by_name[name] = reference.reference
        self._residual_reference = disc_comp.getResidualReference(vale_by_name)
        return self._residual_reference

    @profile
    def getRelativeScaling(self, residuals):
        """Returns the scaling factor to compute the relative error

        Arguments:
            residuals (Residuals): Collections of residuals.

        Returns:
            float: scaling factor.
        """
        scaling = 0.0
        eliminatedDofs = self.problem.getDirichletBCDOFs()

        residuals.update()

        nume_equa = self.problem.getDOFNumbering().getEquationNumbering()
        cmp2dof = nume_equa.getDOFFromNodeAndComponent()

        disc_comp = DiscreteComputation(self.problem)
        varc = None

        computeVarc = (
            self.state.externVar
            and self.problem.getMaterialField().hasExternalStateVariableForLoad()
            and self.problem.isMechanical()
        )
        if computeVarc:
            varc = disc_comp.getExternalStateVariablesForces(
                self.state.time_curr,
                self.state.externVar,
                self.state.getState(-1).externVar,
                self.state.internVar,
                self.state.getState(-1).stress,
            ).getValues()

        for [iNode, cmp], ieq in cmp2dof.items():
            f_int = 0.0
            f_ext = 0.0
            f_cont = 0.0
            f_mass = 0.0
            if varc:
                f_varc = varc[ieq]
            else:
                f_varc = 0.0

            if cmp in ("LAGR_C", "LAGR_F1", "LAGR_F2"):
                continue

            if eliminatedDofs[ieq] == 1:
                f_int = -residuals.resi_int[ieq]
            else:
                f_int = residuals.resi_dual[ieq]
                f_cont = residuals.resi_cont[ieq]
                f_ext = residuals.resi_ext[ieq]
                if residuals.resi_mass:
                    f_mass = residuals.resi_mass[ieq]

            value = abs(f_int + f_cont + f_mass - f_ext) + abs(f_varc)

            if scaling < value:
                scaling = value

        return MPI.ASTER_COMM_WORLD.allreduce(scaling, MPI.MAX)

    @profile
    def evalNormResidual(self, residuals):
        """Evaluate global residual.

        Arguments:
            residuals (Residuals): Collections of residuals.
        """
        residual = self.getDirichletResidual(residuals.resi)
        resi_maxi = self.setdefault("RESI_GLOB_MAXI")
        resi_rela = self.setdefault("RESI_GLOB_RELA")
        resi_refe = self.setdefault("RESI_REFE_RELA")

        resi_maxi.value = residual.norm("NORM_INFINITY")

        # idx = np.abs(np.asarray(residual.getValues())).argmax()
        # info = residual.getValuesWithDescription()[1]
        # print(f"MaxAbs for node {info[0][idx]+1} dof {info[1][idx]}")
        scaling = self.getRelativeScaling(residuals)
        residual_rela = residual.copy()
        if scaling == 0.0:
            resi_rela.value = -1.0
            residual_rela.setValues([-1] * residual_rela.size())
        else:
            resi_rela.value = resi_maxi.value / scaling
            residual_rela = residual / scaling

        if resi_refe.hasRef():
            residual_refe = residual.copy()
            residual_refe /= self.getResidualReference()
            resi_refe.value = residual_refe.norm("NORM_INFINITY")

        resi_fields = {"RESI_NOEU": residual, "RESI_RELA_NOEU": residual_rela}

        return resi_fields

    @profile
    def evalGeometricResidual(self, displ_delta):
        """Evaluate geometric residual.

        Arguments:
            displ_dela (FieldOnNodesReal): variation of displacement.
        """
        # scaling with diagonal of bounding box
        TABG = self.problem.getMesh().getTable("CARA_GEOM")
        x_diag = TABG["X_MAX", 1] - TABG["X_MIN", 1]
        y_diag = TABG["Y_MAX", 1] - TABG["Y_MIN", 1]
        z_diag = TABG["Z_MAX", 1] - TABG["Z_MIN", 1]
        diag = sqrt(pow(x_diag, 2) + pow(y_diag, 2) + pow(z_diag, 2))

        resi_geom = self.setdefault("RESI_GEOM")
        if diag == 0.0:
            resi_geom.value = -1.0
        else:
            resi_geom.value = displ_delta.norm("NORM_INFINITY", ["DX", "DY", "DZ"]) / diag

    # @with_loglevel()
    def isConverged(self):
        """Tell if the convergence parameters are verified.

        Returns:
            bool: *True* if converged, *False* otherwise.
        """
        logger.debug("isConverged ? %r", self._param)
        defined = [para for para in self._param.values() if para.isSet() and para.hasRef()]
        if not defined:
            logger.debug("no parameter set: not converged")
            return False
        for name in self._param:
            para = self._param[name]
            if not para.isConverged():
                logger.debug("parameter %s is not converged", name)
                return False
        return True

    def isPrediction(self):
        """Tell if the current Nuewton iteration is the prediction
        iteration.

        Returns:
            bool: *True* if prediction, *False* otherwise.
        """
        name = "ITER_GLOB_MAXI"
        para = self._param[name]
        return para.isPrediction()

    # @with_loglevel()
    def isFinished(self):
        """Tell if a parameter is stopping the calculation.

        Returns:
                bool: *True* if the calculation should be stopped, *False* otherwise.
        """
        # one iteration parameter can stop
        for name, para in self._iterations:
            if para.isFinished():
                logger.debug("parameter %s would finish", name)
                return True
            if not para.isConverged():
                logger.debug("parameter %s is not converged", name)
                return False
        return self.isConverged()
