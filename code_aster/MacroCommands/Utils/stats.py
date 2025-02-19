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

# extract from pystat module
#
from math import exp, log, pi, pow, sqrt

import numpy

# CDF
#
# -------------------------------------------------------------------


def normcdf(X):
    # Cumulative normal distribution

    (a1, a2, a3, a4, a5) = (0.31938153, -0.356563782, 1.781477937, -1.821255978, 1.330274429)
    L = numpy.absolute(X)
    K = 1.0 / (1.0 + 0.2316419 * L)
    w = 1.0 - 1.0 / sqrt(2 * pi) * exp(-L * L / 2.0) * (
        a1 * K + a2 * K * K + a3 * pow(K, 3) + a4 * pow(K, 4) + a5 * pow(K, 5)
    )
    if X < 0:
        w = 1.0 - w
    return w


# Inverse CDF


def normicdf(v):
    if v > 0.5:
        r = -1.0
    else:
        r = 1.0
    xp = 0.0
    lim = 1.0e-20
    p = [-0.322232431088, -1.0, -0.342242088547, -0.0204231210245, -0.453642210148e-4]
    q = [0.0993484626060, 0.588581570495, 0.531103462366, 0.103537752850, 0.38560700634e-2]

    if v < lim or v == 1:
        return -1.0 / lim
    elif v == 0.5:
        return 0
    elif v > 0.5:
        v = 1.0 - v
    y = sqrt(log(1.0 / v ** 2.0))
    xp = y + ((((y * p[4] + p[3]) * y + p[2]) * y + p[1]) * y + p[0]) / (
        (((y * q[4] + q[3]) * y + q[2]) * y + q[1]) * y + q[0]
    )
    if v < 0.5:
        xp *= -1.0
    return xp * r


# --


def linregress(x, y=None):
    """
    Calculate a regression line

    This computes a least-squares regression for two sets of measurements.

    Parameters
    ----------
    x, y : array_like
    two sets of measurements. Both arrays should have the same length.
    If only x is given (and y=None), then it must be a two-dimensional
    array where one dimension has length 2. The two sets of measurements
    are then found by splitting the array along the length-2 dimension.

    Returns
    -------
    slope : float
    slope of the regression line
    intercept : float
    intercept of the regression line
    r-value : float
    correlation coefficient
    p-value : float
    two-sided p-value for a hypothesis test whose null hypothesis is
    that the slope is zero.
    stderr : float
    Standard error of the estimate


    Examples
    --------
    >>> from scipy import stats
    >>> import numpy as np
    >>> x = np.random.random(10)
    >>> y = np.random.random(10)
    >>> slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

    # To get coefficient of determination (r_squared)

    >>> print "r-squared:", r_value**2
    r-squared: 0.15286643777
    """
    TINY = 1.0e-20
    if y is None:  # x is a (2, N) or (N, 2) shaped array_like
        x = numpy.asarray(x)
        if x.shape[0] == 2:
            x, y = x
        elif x.shape[1] == 2:
            x, y = x.T
        else:
            msg = (
                "If only `x` is given as input, it has to be of shape (2, N) \
or (N, 2), provided shape was %s"
                % str(x.shape)
            )
            raise ValueError(msg)
    else:
        x = numpy.asarray(x)
        y = numpy.asarray(y)
    n = len(x)
    xmean = numpy.mean(x, None)
    ymean = numpy.mean(y, None)

    # average sum of squares:
    ssxm, ssxym, ssyxm, ssym = numpy.cov(x, y, bias=1).flat
    r_num = ssxym
    r_den = numpy.sqrt(ssxm * ssym)
    if r_den == 0.0:
        r = 0.0
    else:
        r = r_num / r_den
        # test for numerical error propagation
        if r > 1.0:
            r = 1.0
        elif r < -1.0:
            r = -1.0

    df = n - 2
    t = r * numpy.sqrt(df / ((1.0 - r + TINY) * (1.0 + r + TINY)))
    #    prob = distributions.t.sf(numpy.abs(t),df)*2
    slope = r_num / ssxm
    intercept = ymean - slope * xmean
    sterrest = numpy.sqrt((1 - r * r) * ssym / ssxm / df)
    pred = intercept + slope * x
    sigma = numpy.sqrt(1.0 / (len(x) - 1) * numpy.sum((y - pred) ** 2))
    return slope, intercept, sigma
