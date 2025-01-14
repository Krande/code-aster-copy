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
import time
from contextlib import contextmanager
from pathlib import Path

import numpy as np

from code_aster.Commands import CALC_CHAM_ELEM
from code_aster import CA


@contextmanager
def chrono(title):
    """Measure elapsed time."""
    t0 = time.time()
    yield
    tf = time.time()
    print(f"elapsed time: {title}: {tf-t0:.6f}")


CA.init("--test")

test = CA.TestCase()

refinement = int(os.environ.get("ZZZZ505_REFINEMENT", 0))

mesh = CA.Mesh.buildCube(refine=refinement)
model = CA.Model(mesh)
model.addModelingOnMesh(CA.Physics.Mechanics, CA.Modelings.Tridimensional)
model.build()

ch0 = CALC_CHAM_ELEM(MODELE=model, OPTION="COOR_ELGA")
names = ch0.getComponents()
print("components:", names)

with chrono("toSimpleFieldOnCells"):
    chs = ch0.toSimpleFieldOnCells()

cumsize = 0

# use False to force 'getComponentValues' to reset the cache
before = chs.getValues(copy=False)[0]
copied = before.copy()

with chrono("chs.X"):
    valx = chs.getComponentValues("X")
with chrono("chs.Y"):
    valy = chs.Y

size = valx.size

with chrono("1 add"):
    xp1 = 1.0 + valx
with chrono("1 iadd"):
    xp1 += 1.0
with chrono("2 sums"):
    s1, s0 = xp1.sum(), valx.sum()
test.assertAlmostEqual(s1, s0 + 2 * size, 8, msg="add, radd, iadd")

with chrono("3 subs"):
    opp = -valx
    xm1 = valx - 1.0
    xrm1 = 1.0 - valx

with chrono("4 sums"):
    so, s0, sm1, srm1 = opp.sum(), valx.sum(), xm1.sum(), xrm1.sum()
test.assertAlmostEqual(so, -s0, 8, msg="neg")
test.assertAlmostEqual(sm1, s0 - size, 8, msg="sub")
test.assertAlmostEqual(sm1, -srm1, 8, msg="rsub")
cumsize += opp.sizeof() + xm1.sizeof() + xrm1.sizeof()
del opp, xm1, xrm1

with chrono("1 isub"):
    xp1 -= 1.0
with chrono("1 mul + 1 add"):
    x2 = valx * 0.0 + 2.0
test.assertAlmostEqual(x2.sum(), 2 * size, 8, msg="mul, add")

with chrono("1 div"):
    div2 = valx / x2
test.assertAlmostEqual(div2.sum(), valx.sum() * 0.5, 8, msg="div")
cumsize += div2.sizeof()
del div2

with chrono("4 muls + 2 adds + 2 negs"):
    comb = 2.0 * valx + valy * 3
    zero = comb - valx * 2 - 3.0 * valy
    null = abs(zero)

with chrono("1 min + 1 max + 1 mean + 1 sum"):
    vmin, vmax, vmean, vsum = null.min(), null.max(), null.mean(), null.sum()
test.assertAlmostEqual(vmin, 0.0, 8, msg="vmin")
test.assertAlmostEqual(vmax, 0.0, 8, msg="vmax")
test.assertAlmostEqual(vmean, 0.0, 8, msg="mean")
test.assertAlmostEqual(vsum, 0.0, 8, msg="sum")
cumsize += comb.sizeof() + zero.sizeof()
del comb, zero

with chrono("4 muls + 1 iadd + 2 isubs + 1 abs"):
    comb = 2.0 * valx
    comb += valy * 3
    comb -= valx * 2
    comb -= 3.0 * valy
    null = abs(comb)
test.assertAlmostEqual(null.max(), 0.0, 8, msg="null comb + abs with i-ops")
cumsize += comb.sizeof() + null.sizeof()
del comb, null

with chrono("1 div"):
    one = x2 / 2.0
with chrono("1 rdiv"):
    one = 2.0 / x2
test.assertAlmostEqual(one.sum(), size, 8, msg="div cst")
test.assertAlmostEqual(one.sum(), size, 8, msg="rdiv cst")
cumsize += one.sizeof()
del one

with chrono("1 pow (int)"):
    x2 = valx**2
with chrono("1 pow (float)"):
    x0 = x2**0.5
test.assertAlmostEqual(x2.max(), valx.max() ** 2, 8, msg="pow int")
test.assertAlmostEqual(abs(x0 - valx).max(), 0.0, 8, msg="pow float")

with chrono("setComponentValues"):
    chs.setComponentValues("X", xp1)

after = chs.getValues()[0]
test.assertAlmostEqual(np.sum(after), np.sum(copied) + x0.size, 8, msg="setComponentValues values")
# 'before' points on the same data, so it now contains 'xp1' as first component
test.assertAlmostEqual(np.sum(after), np.sum(before), 8, msg="setComponentValues: unchanged origin")
cumsize += x0.sizeof() + x2.sizeof() + xp1.sizeof()
del x0, x2, xp1
del after, before, copied

test.assertIsNone(valx._cells, msg="on all cells")
faces = mesh.getCells(["TOP", "FRONT", "RIGHT", "S68", "S26"])
test.assertEqual(len(faces), 3 * 2 ** (refinement * 2) + 2 * 2**refinement, msg="getCells")

with chrono("restrict"):
    valx.restrict(faces)

test.assertTrue(np.all(faces == valx._cells), msg="restrict: cells")
test.assertEqual(valx.size, valx._descr._discr.prod(axis=1).sum(), msg="restrict: values")
test.assertEqual(valx._descr._nbval, valy.size, msg="restrict: nbval")
test.assertEqual(valx._descr._nbcells, valy._descr._nbcells, msg="restrict: nbcells")

with chrono("expand"):
    xe = valx.expand()
test.assertAlmostEqual(sum(xe.values), sum(valx.values), 8, msg="expand")
cumsize += xe.sizeof()
del xe

with chrono("setComponentValues"):
    chs.setComponentValues("X", valx)
cumsize += valx.sizeof() + valy.sizeof()
del valx, valy

one_f = CA.FieldOnCellsReal(model, "ELGA", "SIEF_R")
one_f.setValues(1.0)
one_s = one_f.toSimpleFieldOnCells()
one = one_s.SIXX  # only assigned on VOLUME
weight = chs.W
zone = mesh.getCells(["VOLUME", "FRONT"])
vol = np.array(mesh.getCells(["VOLUME"]))

with test.assertRaises(IndexError):
    one.onSupportOf(weight, strict=True)

with chrono("onSupportOf, prol needed"):
    zerone = one.onSupportOf(weight, strict=False)  # 1.0 on VOLUME, 0.0 elsewhere
test.assertEqual(zerone.size, weight.size, msg="onSupportOf: extended by zero")
test.assertAlmostEqual(sum(zerone.values), sum(one.values), 8, msg="onSupportOf: values")
cumsize += zerone.sizeof()
del zerone

with chrono("onSupportOf"):
    hexa_weight = weight.onSupportOf(one)
test.assertEqual(hexa_weight.size, one.size, msg="onSupportOf")
test.assertAlmostEqual(hexa_weight.sum(), 1.0, 8, msg="onSupportOf: values")
cumsize += one.sizeof()
del one_s, one

with chrono("getValuesByCells"):
    cells_sizes = weight.getValuesByCells().sum(axis=1)

quad_size = 1.0 / 2 ** (refinement * 2)
test.assertAlmostEqual(
    (cells_sizes[mesh.getCells("FRONT")]).mean(), quad_size, 8, msg="size of quad"
)
hexa_size = 1.0 / 2 ** (refinement * 3)
test.assertAlmostEqual(
    (cells_sizes[mesh.getCells("VOLUME")]).mean(), hexa_size, 8, msg="size of hexa"
)

# \int_{0}^{1} y.dy = 0.5
vy = chs.Y
vy.restrict(vol)
wvol = weight.onSupportOf(vy)
integr = (vy * wvol).sum()
test.assertAlmostEqual(integr, 0.5, msg="integr")

# check for getValuesByCells on a restricted component
restr = wvol.getValuesByCells()
v_all = hexa_weight.getValuesByCells()
test.assertAlmostEqual(restr.sum(), v_all.sum(), 8, msg="getValuesByCells/restrict")
test.assertEqual(wvol.getNumberOfCells(), hexa_weight.getNumberOfCells(), msg="getNumberOfCells")
del hexa_weight, restr, v_all

# only relevant with refinement > 0
vz = chs.Z
vz.restrict(vol)
vz.maximum(vy)
maxi = vz.max() * 0.99999
cells = vz.filterByValues(maxi, 1.0, strict_mini=False)
expect = (2**refinement + (2**refinement - 1)) * 2**refinement  # last term for x direction
test.assertEqual(len(cells), expect, msg="number filtered cells")
cumsize += vz.sizeof()
del vz

# check for description plot
dest = Path("weight_descr.png")
weight.plot_descr(filename=dest)
test.assertTrue(dest.exists(), msg="figure png")

test.printSummary()
print(f"cumulative size of created ComponentOnCells objects: {cumsize:,.0f} bytes")

CA.close()
