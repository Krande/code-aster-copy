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

"""
    extract_propnames.py [options] file.mfront

This script extracts the properties names from a '.mfront' file using
the Python API of MFront.

The path to the MFront file is passed as first argument.
The properties names are printed on stdout.

NB : Using an external script avoids to interrupt the build process in case of
fatal error (not raised by exceptions).
"""

import argparse
import os
import os.path as osp


def getMaterialProperties(filename, suffix):
    """Return the list of material properties

    Args:
        filename (str): Path of the '.mfront' file.
        suffix (str): Installation suffix (TFELSUFFIX).

    Returns:
        list[str]: List of the names of the properties.
    """
    suffix = suffix.replace("-", "_").replace(".", "_")
    tfel = __import__("tfel" + suffix + ".material", fromlist=[])
    mfront = __import__("mfront" + suffix)
    ModellingHypothesis = tfel.material.ModellingHypothesis

    os.environ["MFRONT_INCLUDE_PATH"] = osp.dirname(filename)
    props = []

    dsl = mfront.getDSL(filename)
    dsl.analyseFile(filename, [])
    # behaviour description
    desc = dsl.getBehaviourDescription()

    # Should also check "RequiresThermalExpansionCoefficientTensor" and "Symmetry"
    # (see old bibcxx/mfront/MFrontBehaviour.cxx up to 13.3.13)
    req = desc.getBooleanAttribute("requiresStiffnessTensor", False)
    if req:
        props.extend(["YoungModulus{0}".format(i) for i in (1, 2, 3)])
        props.extend(["PoissonRatio{0}".format(i) for i in (12, 23, 13)])
        props.extend(["ShearModulus{0}".format(i) for i in (12, 23, 13)])

    data = desc.getBehaviourData(ModellingHypothesis.TRIDIMENSIONAL)
    for var in data.getMaterialProperties():
        name = data.getExternalName(var.name)
        if var.arraySize <= 1:
            props.append(name)
        else:
            for i in range(var.arraySize):
                props.append("{0}_{1}".format(name, i))

    return props


if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("--suffix", action="store", default="", help="TFEL libraries suffix")
    parser.add_argument("filename", metavar="file.mfront", help="MFront file")
    args = parser.parse_args()

    props = getMaterialProperties(args.filename, args.suffix)
    print(" ".join(props))
