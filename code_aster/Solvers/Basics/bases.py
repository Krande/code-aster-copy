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

"""
Base objects used to solve generic non linear problems.
"""

from ...Cata.Language.SyntaxObjects import _F
from ...Utilities import force_list, logger, no_new_attributes


class KeywordsStore:
    """Container that stores and gives access to some user keywords."""

    _keywords = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, keywords):
        self._keywords = keywords

    def get(self, keyword, parameter=None, default=None):
        """ "Return a keyword value.

        Args:
            keyword (str): Simple or factor keyword.
            parameter (str|None): Simple keyword under the factor keyword, or *None*.
            default (*misc*): Default value if the keyword is undefined.

        Returns:
            *misc*: Keyword value.
        """
        args = self._keywords
        if parameter is not None:
            if args.get(keyword) is None:
                return default
            return _F(args[keyword])[0].get(parameter, default)

        return args.get(keyword, default)


class ContextMixin:
    """Mixin object that store the objects required by the command or
    the informations to build the right object when needed.

    Attributes:
        problem: PhysicalProblem object
        state: PhysicalState object
        result: Result object (NonLinearResult, ThermalResult)
        problem_type: ProblemType enum value
        keywords: Part of the user keywords
        oper: Operators object
        contact: ContactManager object
        linear_solver: LinearSolver object
    """

    # FIXME: @dataclass with python>=3.7
    class Data:
        _problem = _type = _state = _keywords = _result = None
        _oper = _contact = _linsolv = None
        __setattr__ = no_new_attributes(object.__setattr__)

        def __init__(self):
            self._problem = None
            # FIXME: check if it's needed
            self._type = None
            self._state = None
            self._keywords = KeywordsStore()
            self._result = None
            # FIXME: to be removed and created deeper/later from '.keywords'?
            # FIXME: with TimeScheme.Multiple? several '.oper'? one '.oper' per StepSolver?
            self._oper = None
            self._contact = None
            self._linsolv = None

    _ctxt = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self):
        self._ctxt = ContextMixin.Data()

    @property
    def context(self):
        """Data: Context attached to the object."""
        return self._ctxt

    def useProblemContext(self, parent):
        """Set context content from the parent one."""
        self._ctxt = parent._ctxt

    # convenient shortcuts properties
    @property
    def problem(self):
        """PhysicalProblem: current problem description."""
        return self._ctxt._problem

    @problem.setter
    def problem(self, problem):
        assert not self._ctxt._problem, "must be set only once!"
        self._ctxt._problem = problem

    @property
    def state(self):
        """PhysicalState: current state."""
        return self._ctxt._state

    @state.setter
    def state(self, state):
        self._ctxt._state = state

    @property
    def result(self):
        """Result: Attribute that holds the result object."""
        return self._ctxt._result

    @result.setter
    def result(self, value):
        self._ctxt._result = value

    # FIXME: à voir : on garde certains mots-clés ou uniquement quelques infos
    @property
    def keywords(self):
        """Dict: Attribute that holds the keywords object."""
        return self._ctxt._keywords

    @keywords.setter
    def keywords(self, value_dict):
        assert isinstance(value_dict, dict), f"unsupported type: {type(value_dict)}"
        self._ctxt._keywords = KeywordsStore(value_dict)

    def get_keyword(self, keyword, parameter=None, default=None):
        """ "Return a keyword value.

        Args:
            keyword (str): Simple or factor keyword.
            parameter (str|None): Simple keyword under the factor keyword, or *None*.
            default (*misc*): Default value if the keyword is undefined.

        Returns:
            *misc*: Keyword value.
        """
        return self._ctxt._keywords.get(keyword, parameter, default)

    @property
    def problem_type(self):
        """ProblemType: Attribute that holds the type of problem."""
        return self._ctxt._type

    @problem_type.setter
    def problem_type(self, value):
        self._ctxt._type = value

    @property
    def oper(self):
        """Operators: Objects that adapts operators for each type of problem."""
        return self._ctxt._oper

    @oper.setter
    def oper(self, value):
        self._ctxt._oper = value

    @property
    def contact(self):
        """ContactManager: Objects to solve contact conditions"""
        logger.debug("CTXT: ctxt.contact: %s", self._ctxt._contact)
        return self._ctxt._contact

    @contact.setter
    def contact(self, value):
        self._ctxt._contact = value

    @property
    def linear_solver(self):
        """LinearSolver: Attribute that holds the linear solver."""
        return self._linsolv

    @linear_solver.setter
    def linear_solver(self, value):
        self.linsolv = value


def convert_F(input, max=1):
    """Return a normalized value of a factor keyword.
    If max=1, a `_F` object is returned.
    Otherwise, a list of `_F` is returned.

    Args:
        input (list[dict|_F]): A dict or a `_F`, or a list of.
        max (int|'**'): Value of the 'max' attribute in the keyword description.

    Returns:
        `_F`|list[`_F`]: one `_F` object if max=1 or a list of `_F` objects.
    """
    if max == 1:
        return _F(force_list(input)[0])
    return [_F(i) for i in force_list(input)]
