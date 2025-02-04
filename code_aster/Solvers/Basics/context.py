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

from ...Utilities import no_new_attributes


class ContextMixin:
    """Mixin object that store the objects required by the command or
    the informations to build the right object when needed."""

    class Data:
        _problem = _type = _state = _keywords = _result = _contact = None
        __setattr__ = no_new_attributes(object.__setattr__)

        def __init__(self):
            self._problem = None
            self._type = None
            self._state = None
            self._keywords = {}
            self._result = None
            self._contact = None  # FIXME: remove? create deeper?

    _ctxt = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self):
        self._ctxt = ContextMixin.Data()

    def set_context(self, context):
        self._ctxt = context

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

    @property
    def contact(self):
        """ContactManager: Objects to solve contact conditions"""
        return self._ctxt._contact

    @contact.setter
    def contact(self, value):
        self._ctxt._contact = value

    # FIXME: à voir : on garde certains mots-clés ou uniquement quelques infos
    @property
    def keywords(self):
        """Dict: Attribute that holds the keywords object."""
        return self._ctxt._keywords

    @keywords.setter
    def keywords(self, value):
        self._ctxt._keywords = value

    @property
    def problem_type(self):
        """ProblemType: Attribute that holds the type of problem."""
        return self._ctxt._type

    @problem_type.setter
    def problem_type(self, value):
        self._ctxt._type = value
