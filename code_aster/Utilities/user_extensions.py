# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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


class WithEmbeddedObjects:
    """Base object for user classes that embed code_aster objects in their
    attributes.

    code_aster datastructures can directly be stored as attribute, in
    attributes of types list, tuple or dict or in such nested objects.

    Attributes:
        aster_embedded (list[str]): The list of the names of the attributes
            that contain code_aster objects.

    .. note:: User classes must inherit from this class to be serialized with
        other code_aster datastructures.
    """
    aster_embedded = []
