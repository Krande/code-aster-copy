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

import numpy as np


def coef_pts_B2(ix1, ix2, iy1, iy2):

    coef = np.array(
        [
            [
                [0.558, 0.558, 0.558, 0.558, 0.558],
                [0.558, 0.544, 0.531, 0.518, 0.506],
                [0.559, 0.533, 0.511, 0.490, 0.471],
                [0.560, 0.517, 0.482, 0.452, 0.426],
                [0.562, 0.498, 0.449, 0.411, 0.381],
                [0.565, 0.487, 0.432, 0.391, 0.359],
                [0.568, 0.480, 0.421, 0.378, 0.346],
                [0.573, 0.472, 0.408, 0.364, 0.332],
                [0.578, 0.468, 0.401, 0.357, 0.324],
            ],
            [
                [0.657, 0.657, 0.657, 0.657, 0.657],
                [0.658, 0.640, 0.622, 0.606, 0.590],
                [0.660, 0.626, 0.596, 0.569, 0.545],
                [0.662, 0.606, 0.560, 0.521, 0.487],
                [0.669, 0.583, 0.519, 0.470, 0.431],
                [0.675, 0.570, 0.497, 0.444, 0.403],
                [0.681, 0.562, 0.484, 0.428, 0.388],
                [0.692, 0.553, 0.468, 0.411, 0.370],
                [0.701, 0.549, 0.460, 0.402, 0.361],
            ],
            [
                [0.695, 0.695, 0.695, 0.695, 0.695],
                [0.698, 0.678, 0.659, 0.642, 0.625],
                [0.701, 0.665, 0.632, 0.603, 0.576],
                [0.706, 0.645, 0.593, 0.550, 0.514],
                [0.716, 0.621, 0.550, 0.496, 0.453],
                [0.726, 0.608, 0.527, 0.469, 0.424],
                [0.734, 0.601, 0.513, 0.452, 0.408],
                [0.749, 0.592, 0.497, 0.434, 0.389],
                [0.761, 0.589, 0.489, 0.425, 0.380],
            ],
            [
                [0.709, 0.709, 0.709, 0.709, 0.709],
                [0.712, 0.691, 0.672, 0.654, 0.636],
                [0.716, 0.679, 0.645, 0.615, 0.587],
                [0.723, 0.659, 0.606, 0.562, 0.524],
                [0.734, 0.636, 0.562, 0.506, 0.462],
                [0.745, 0.623, 0.539, 0.478, 0.433],
                [0.755, 0.616, 0.525, 0.462, 0.416],
                [0.772, 0.608, 0.510, 0.444, 0.397],
                [0.786, 0.605, 0.502, 0.435, 0.388],
            ],
            [
                [0.713, 0.713, 0.713, 0.713, 0.713],
                [0.717, 0.696, 0.676, 0.658, 0.640],
                [0.721, 0.684, 0.650, 0.619, 0.591],
                [0.728, 0.664, 0.611, 0.566, 0.528],
                [0.741, 0.641, 0.567, 0.510, 0.466],
                [0.752, 0.628, 0.543, 0.482, 0.436],
                [0.762, 0.621, 0.529, 0.466, 0.419],
                [0.781, 0.615, 0.514, 0.448, 0.400],
                [0.797, 0.613, 0.507, 0.439, 0.391],
            ],
            [
                [0.702, 0.702, 0.702, 0.702, 0.702],
                [0.703, 0.682, 0.662, 0.643, 0.625],
                [0.705, 0.666, 0.632, 0.601, 0.573],
                [0.709, 0.645, 0.591, 0.546, 0.508],
                [0.722, 0.622, 0.548, 0.491, 0.446],
                [0.733, 0.609, 0.524, 0.463, 0.416],
                [0.744, 0.603, 0.510, 0.446, 0.399],
                [0.766, 0.597, 0.496, 0.429, 0.381],
                [0.785, 0.597, 0.490, 0.421, 0.373],
            ],
        ]
    )

    return coef[ix1, iy1, :], coef[ix2, iy1, :], coef[ix1, iy2, :], coef[ix2, iy2, :]
