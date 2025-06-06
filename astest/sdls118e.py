# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
from code_aster.Commands import DEFI_FONCTION


def F_FONC():

    ACCE1 = DEFI_FONCTION(
        NOM_PARA="INST",
        VALE=(
            0.0,
            0.0,
            0.013672,
            -0.01765,
            0.027344,
            -0.020471,
            0.041016,
            -0.017756,
            0.054688,
            -0.013225,
            0.06836,
            -0.01394,
            0.082032,
            -0.020436,
            0.095704,
            -0.021801,
            0.109376,
            -0.018725,
            0.123048,
            -0.0081087,
            0.13672,
            0.0095333,
            0.150392,
            0.01939,
            0.164064,
            0.015589,
            0.177736,
            0.0081624,
            0.191408,
            0.0052446,
            0.20508,
            0.012322,
            0.218752,
            0.0061115,
            0.232424,
            -0.014552,
            0.246096,
            -0.017576,
            0.259768,
            -0.0030371,
            0.27344,
            0.021728,
            0.287112,
            0.032462,
            0.300784,
            0.030741,
            0.314456,
            0.030607,
            0.328128,
            0.025483,
            0.3418,
            0.015453,
            0.355472,
            -0.015674,
            0.369144,
            -0.044803,
            0.382816,
            -0.032032,
            0.396488,
            -0.0012607,
            0.41016,
            0.0079026,
            0.423832,
            -0.0063606,
            0.437504,
            -0.015367,
            0.451176,
            0.01177,
            0.464848,
            0.053128,
            0.47852,
            0.056559,
            0.492192,
            0.027427,
            0.505864,
            0.0099848,
            0.519536,
            0.018467,
            0.533208,
            0.034255,
            0.54688,
            0.023207,
            0.560552,
            -0.013605,
            0.574224,
            -0.038973,
            0.587896,
            -0.042627,
            0.601568,
            -0.029954,
            0.61524,
            -0.041828,
            0.628912,
            -0.073322,
            0.642584,
            -0.022112,
            0.656256,
            0.069766,
            0.669928,
            0.087114,
            0.6836,
            0.049358,
            0.697272,
            0.0045632,
            0.710944,
            -0.0033541,
            0.724616,
            -0.00082904,
            0.738288,
            -0.019078,
            0.75196,
            -0.041723,
            0.765632,
            -0.054787,
            0.779304,
            -0.0056487,
            0.792976,
            0.06785,
            0.806648,
            0.097993,
            0.82032,
            0.077278,
            0.833992,
            0.023088,
            0.847664,
            -0.01176,
            0.861336,
            -0.0050909,
        ),
        PROL_DROITE="CONSTANT",
        PROL_GAUCHE="CONSTANT",
    )

    ACCE2 = DEFI_FONCTION(
        NOM_PARA="INST",
        VALE=(
            0.0,
            0.0,
            0.013672,
            0.01765,
            0.027344,
            0.020471,
            0.041016,
            0.017756,
            0.054688,
            0.013225,
            0.06836,
            0.01394,
            0.082032,
            0.020436,
            0.095704,
            0.021801,
            0.109376,
            0.018725,
            0.123048,
            0.0081087,
            0.13672,
            -0.0095333,
            0.150392,
            -0.01939,
            0.164064,
            -0.015589,
            0.177736,
            -0.0081624,
            0.191408,
            -0.0052446,
            0.20508,
            -0.012322,
            0.218752,
            -0.0061115,
            0.232424,
            0.014552,
            0.246096,
            0.017576,
            0.259768,
            0.0030371,
            0.27344,
            -0.021728,
            0.287112,
            -0.032462,
            0.300784,
            -0.030741,
            0.314456,
            -0.030607,
            0.328128,
            -0.025483,
            0.3418,
            -0.015453,
            0.355472,
            0.015674,
            0.369144,
            0.044803,
            0.382816,
            0.032032,
            0.396488,
            0.0012607,
            0.41016,
            -0.0079026,
            0.423832,
            0.0063606,
            0.437504,
            0.015367,
            0.451176,
            -0.01177,
            0.464848,
            -0.053128,
            0.47852,
            -0.056559,
            0.492192,
            -0.027427,
            0.505864,
            -0.0099848,
            0.519536,
            -0.018467,
            0.533208,
            -0.034255,
            0.54688,
            -0.023207,
            0.560552,
            0.013605,
            0.574224,
            0.038973,
            0.587896,
            0.042627,
            0.601568,
            0.029954,
            0.61524,
            0.041828,
            0.628912,
            0.073322,
            0.642584,
            0.022112,
            0.656256,
            -0.069766,
            0.669928,
            -0.087114,
            0.6836,
            -0.049358,
            0.697272,
            -0.0045632,
            0.710944,
            0.0033541,
            0.724616,
            0.00082904,
            0.738288,
            0.019078,
            0.75196,
            0.041723,
            0.765632,
            0.054787,
            0.779304,
            0.0056487,
            0.792976,
            -0.06785,
            0.806648,
            -0.097993,
            0.82032,
            -0.077278,
            0.833992,
            -0.023088,
            0.847664,
            0.01176,
            0.861336,
            0.0050909,
        ),
        PROL_DROITE="CONSTANT",
        PROL_GAUCHE="CONSTANT",
    )

    ACCE3 = DEFI_FONCTION(
        NOM_PARA="INST",
        VALE=(
            0.0,
            0.0,
            0.013672,
            -0.01765,
            0.027344,
            -0.020471,
            0.041016,
            -0.017756,
            0.054688,
            -0.013225,
            0.06836,
            -0.01394,
            0.082032,
            -0.020436,
            0.095704,
            -0.021801,
            0.109376,
            -0.018725,
            0.123048,
            -0.0081087,
            0.13672,
            0.0095333,
            0.150392,
            0.01939,
            0.164064,
            0.015589,
            0.177736,
            0.0081624,
            0.191408,
            0.0052446,
            0.20508,
            0.012322,
            0.218752,
            0.0061115,
            0.232424,
            -0.014552,
            0.246096,
            -0.017576,
            0.259768,
            -0.0030371,
            0.27344,
            0.021728,
            0.287112,
            0.032462,
            0.300784,
            0.030741,
            0.314456,
            0.030607,
            0.328128,
            0.025483,
            0.3418,
            0.015453,
            0.355472,
            -0.015674,
            0.369144,
            -0.044803,
            0.382816,
            -0.032032,
            0.396488,
            -0.0012607,
            0.41016,
            0.0079026,
            0.423832,
            -0.0063606,
            0.437504,
            0.0,
            0.451176,
            -0.01765,
            0.464848,
            -0.020471,
            0.47852,
            -0.017756,
            0.492192,
            -0.013225,
            0.505864,
            -0.01394,
            0.519536,
            -0.020436,
            0.533208,
            -0.021801,
            0.54688,
            -0.018725,
            0.560552,
            -0.0081087,
            0.574224,
            0.0095333,
            0.587896,
            0.01939,
            0.601568,
            0.015589,
            0.61524,
            0.0081624,
            0.628912,
            0.0052446,
            0.642584,
            0.012322,
            0.656256,
            0.0061115,
            0.669928,
            -0.014552,
            0.6836,
            -0.017576,
            0.697272,
            -0.0030371,
            0.710944,
            0.021728,
            0.724616,
            0.032462,
            0.738288,
            0.030741,
            0.75196,
            0.030607,
            0.765632,
            0.025483,
            0.779304,
            0.015453,
            0.792976,
            -0.015674,
            0.806648,
            -0.044803,
            0.82032,
            -0.032032,
            0.833992,
            -0.0012607,
            0.847664,
            0.0079026,
            0.861336,
            -0.0063606,
        ),
        PROL_DROITE="CONSTANT",
        PROL_GAUCHE="CONSTANT",
    )

    return ACCE1, ACCE2, ACCE3
