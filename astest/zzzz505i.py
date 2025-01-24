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

import os
from pathlib import Path

import numpy as np

from code_aster import CA
from code_aster.CA import MPI
from code_aster.Helpers.debugging import DebugChrono
from code_aster.Utilities import useHPCMode

CA.init("--test")
comm = MPI.ASTER_COMM_WORLD
test = CA.TestCase()

refinement = int(os.environ.get("ZZZZ505_REFINEMENT", 1))
create_mesh = False

meshdir = CA.basedir / f"cube_refined_x{refinement}-{comm.size}p"
local_file = meshdir / f"{comm.rank}.med"

if comm.size > 1:
    if create_mesh:
        filename = CA.basedir / f"cube_refined_x{refinement}.med"
        if comm.rank == 0:
            mesh = CA.Mesh.buildCube(refine=refinement)
            mesh.printMedFile(filename)
        comm.Barrier()

        mesh = CA.ParallelMesh()
        mesh.readMedFile(filename, partitioned=False, deterministic=True, verbose=1)
        meshdir.mkdir(exist_ok=True)
        mesh.printMedFile(local_file, local=True)

        test.assertTrue(local_file.exists())
        CA.close(exit=True)

    test.assertTrue(local_file.exists())
    mesh = CA.ParallelMesh()
    mesh.readMedFile(local_file, partitioned=True, verbose=1)
else:
    mesh = CA.Mesh.buildCube(refine=refinement)

model = CA.Model(mesh)
model.addModelingOnMesh(CA.Physics.Mechanics, CA.Modelings.Tridimensional)
model.build()

# extract coordinates as a SimpleFieldOnNodes
with DebugChrono.measure("getCoordinates"):
    chs = mesh.getCoordinatesAsSimpleFieldOnNodes()

names = chs.getComponents()
print("components:", names)

with DebugChrono.measure("global"):
    cumsize = 0

    # use False to force 'getComponentOnNodes' to reset the cache
    before = chs.getValues(copy=False)[0]
    copied = before.copy()

    with DebugChrono.measure("chs.X"):
        valx = chs.getComponentOnNodes("X")
    print(repr(valx))

    with DebugChrono.measure("chs.Y"):
        valy = chs.Y

    # total number of values in the field (globally)
    size = valx.size

    with DebugChrono.measure("1 add"):
        xp1 = 1.0 + valx
    with DebugChrono.measure("1 iadd"):
        xp1 += 1.0
    with DebugChrono.measure("2 sums"):
        s1, s0 = xp1.sum(), valx.sum()
    test.assertAlmostEqual(s1, s0 + 2 * size, 8, msg="add, radd, iadd")

    with DebugChrono.measure("3 subs"):
        opp = -valx
        xm1 = valx - 1.0
        xrm1 = 1.0 - valx

    with DebugChrono.measure("4 sums"):
        so, s0, sm1, srm1 = opp.sum(), valx.sum(), xm1.sum(), xrm1.sum()
    test.assertAlmostEqual(so, -s0, 8, msg="neg")
    test.assertAlmostEqual(sm1, s0 - size, 8, msg="sub")
    test.assertAlmostEqual(sm1, -srm1, 8, msg="rsub")
    cumsize += opp.sizeof() + xm1.sizeof() + xrm1.sizeof()
    del opp, xm1, xrm1

    with DebugChrono.measure("1 isub"):
        xp1 -= 1.0
    with DebugChrono.measure("1 mul + 1 add"):
        x2 = valx * 0.0 + 2.0
    test.assertAlmostEqual(x2.sum(), 2 * size, 8, msg="mul, add")

    with DebugChrono.measure("1 div"):
        div2 = valx / x2
    test.assertAlmostEqual(div2.sum(), valx.sum() * 0.5, 8, msg="div")
    cumsize += div2.sizeof()
    del div2

    with DebugChrono.measure("4 muls + 2 adds + 2 negs"):
        comb = 2.0 * valx + valy * 3
        zero = comb - valx * 2 - 3.0 * valy
        null = abs(zero)

    with DebugChrono.measure("1 min + 1 max + 1 mean + 1 sum"):
        vmin, vmax, vmean, vsum = null.min(), null.max(), null.mean(), null.sum()
    test.assertAlmostEqual(vmin, 0.0, 8, msg="vmin")
    test.assertAlmostEqual(vmax, 0.0, 8, msg="vmax")
    test.assertAlmostEqual(vmean, 0.0, 8, msg="mean")
    test.assertAlmostEqual(vsum, 0.0, 8, msg="sum")
    cumsize += comb.sizeof() + zero.sizeof()
    del comb, zero

    with DebugChrono.measure("4 muls + 1 iadd + 2 isubs + 1 abs"):
        comb = 2.0 * valx
        comb += valy * 3
        comb -= valx * 2
        comb -= 3.0 * valy
        null = abs(comb)
    test.assertAlmostEqual(null.max(), 0.0, 8, msg="null comb + abs with i-ops")
    cumsize += comb.sizeof() + null.sizeof()
    del comb, null

    with DebugChrono.measure("1 div"):
        one = x2 / 2.0
    with DebugChrono.measure("1 rdiv"):
        one = 2.0 / x2
    test.assertAlmostEqual(one.sum(), size, 8, msg="div cst")
    test.assertAlmostEqual(one.sum(), size, 8, msg="rdiv cst")
    cumsize += one.sizeof()
    del one

    with DebugChrono.measure("1 pow (int)"):
        x2 = valx**2
    with DebugChrono.measure("1 pow (float)"):
        x0 = x2**0.5
    test.assertAlmostEqual(x2.max(), valx.max() ** 2, 8, msg="pow int")
    test.assertAlmostEqual(abs(x0 - valx).max(), 0.0, 8, msg="pow float")

    with DebugChrono.measure("setComponentValues"):
        chs.setComponentValues("X", xp1)

    after = chs.getValues()[0]

    added = len(x0.values)  # 'getValues()' returns all values (inner and outer nodes)
    test.assertAlmostEqual(
        np.sum(after), np.sum(copied) + added, 8, msg="setComponentValues values"
    )
    # 'before' points on the same data, so it now contains 'xp1' as first component
    test.assertAlmostEqual(
        np.sum(after), np.sum(before), 8, msg="setComponentValues: unchanged origin"
    )
    cumsize += x0.sizeof() + x2.sizeof() + xp1.sizeof()
    del x2, xp1
    del after, before, copied

    # restore initial values
    chs.setComponentValues("X", x0)

    test.assertIsNone(valx.restr, msg="on all nodes")
    facez0 = mesh.getNodes(["N1", "N2", "N3", "N4"])
    if useHPCMode():
        facez0 = list(set(facez0).intersection(mesh.getInnerNodes()))
        nbnodes = comm.allreduce(len(facez0), op=MPI.SUM)
    else:
        nbnodes = len(facez0)
    test.assertEqual(nbnodes, 4, msg="getNodes")

    with DebugChrono.measure("restrict"):
        valx.restrict(facez0)

    test.assertTrue(np.all(facez0 == valx.restr), msg="restrict: nodes")
    test.assertEqual(len(valx.values), len(facez0), msg="restrict: values")
    test.assertEqual(valx._descr._nbval, len(valy.values), msg="restrict: nbval")
    test.assertEqual(valx._descr._nbnodes, valy._descr._nbnodes, msg="restrict: nbnodes")

    with DebugChrono.measure("expand"):
        xe = valx.expand()
    test.assertAlmostEqual(sum(xe.values), sum(valx.values), 8, msg="expand")
    cumsize += xe.sizeof()
    del xe

    with DebugChrono.measure("setComponentValues"):
        chs.setComponentValues("X", valx)
    cumsize += valx.sizeof() + valy.sizeof()
    del valx

    # create a field with 1.0 on corners, 0.0 elsewhere
    one_f = CA.FieldOnNodesReal(mesh, "GEOM_R", ("X", "Y", "Z"))
    one_f.setValues(1.0)
    one_s = one_f.toSimpleFieldOnNodes()
    one = one_s.X
    # change mask to test different supports
    mask = np.ones(one.size, dtype=bool)
    mask[facez0] = False
    one.values.mask = mask
    test.assertAlmostEqual(one.sum(), 4.0, 6, msg="sum of non-masked values")

    with test.assertRaises(IndexError):
        one.onSupportOf(valy, strict=True)

    with DebugChrono.measure("onSupportOf, prol needed"):
        vy = valy.onSupportOf(one, strict=False)

    test.assertEqual(vy.sum(), 2.0, msg="onSupportOf: only 4 non-masked values")
    cumsize += vy.sizeof()
    del valy, vy

    with DebugChrono.measure("filtering"):
        vz = chs.Z
        maxi = vz.max() * 0.99999
        nodes = vz.filterByValues(maxi, 1.0, strict_maxi=False)
        expect = (2**refinement + 1) ** 2
        filtered = comm.allreduce(len(nodes), op=MPI.SUM)
        test.assertEqual(filtered, expect, msg="number filtered nodes")
        cumsize += vz.sizeof()
        del vz

# check for description plot
dest = Path("one_descr.png")
ret = one.plot_descr(filename=dest)
test.assertTrue(not ret or dest.exists(), msg="figure png")

test.printSummary()
print(f"cumulative size of created ComponentOnNodes objects: {cumsize:,.0f} bytes")

# for perf statistics
fpick = os.environ.get("ZZZZ505_PICKLE")
if fpick:
    DebugChrono.save(fpick)

CA.close()
