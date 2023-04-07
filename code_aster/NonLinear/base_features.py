# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
Useful objects used to build operators.
"""

from functools import wraps
from math import log


class BaseFeaturesOptions:
    """Options that describes the features an object can provide."""

    Nothing = 0

    @classmethod
    def name(cls, idx):
        """Return the features names from an option.

        Arguments:
            idx (int): Enabled options.

        Returns:
            str: Features names, "&" separated.
        """
        attrs = dir(cls)
        lab = []
        values = {1: ""}
        for attr in attrs:
            if isinstance(getattr(cls, attr), int):
                values[getattr(cls, attr)] = attr
        last = round(log(max(values.keys())) / log(2))
        for expo in range(last + 1):
            if idx & 2**expo:
                lab.append(values[2**expo])
        return "|".join(lab)


class FeatureMeta(type):
    """Metaclass to set the list of supported features."""

    def __new__(cls, name, base, defcl):
        obj = super().__new__(cls, name, base, defcl)
        # compute '_supported' attribute value
        obj._supported = getattr(base, "_supported", 0)
        for feat in defcl.get("required_features", []) + defcl.get("optional_features", []):
            obj._supported |= feat
        return obj


class BaseFeature(metaclass=FeatureMeta):
    """A `BaseFeature` is defined by an object that does the job and the kind of
    services it provides. A feature *uses* other features.
    """

    options = BaseFeaturesOptions
    provide = BaseFeaturesOptions.Nothing
    required_features = []
    optional_features = []
    # for no_new_attributes
    _use = _checked = None

    def __init__(self):
        self._use = []
        self._checked = False

    def use(self, obj, provide=BaseFeaturesOptions.Nothing, replace=False):
        """Add a feature to be used.

        The provided services are defined by the
        'provide' attribute if `obj` is a `BaseFeature` + the `provide` argument.

        Arguments:
            obj (BaseFeature|misc): Object that to be used.
            provide (FeaturesOptions, optional): Features provided by the object.
        """
        if not obj:
            return
        try:
            provide |= obj.provide
        except AttributeError:
            pass
        if not provide:
            raise ValueError("This object provides no feature, use 'provide'.")
        if not provide & self._supported:
            raise TypeError(
                f"{self.__class__.__name__} does not support {self.options.name(provide)!r}"
            )
        self._use.append((obj, provide))

    def check_features(self):
        """Check that required features are defined."""
        for feat in self.required_features:
            if not self.has_feature(feat):
                name = self.options.name(feat)
                raise TypeError(f"{self.__class__.__name__} requires the {name!r} feature")
        self._checked = True

    def has_feature(self, feature):
        """Tell if the feature is available.

        Arguments:
            feature (FeatureOptions): Resquested feature.

        Returns:
            bool: *True* if the feature is available, *False* otherwise.
        """
        return bool(self.get_features(feature))

    def discard(self, worker):
        """Remove a registered worker if it is a member.

        Arguments:
            worker (BaseFeature): Feature object (not the services it provides).
        """
        use = []
        for obj, provide in self._use:
            if obj is not worker:
                use.append((obj, provide))
        self._use = use

    def get_features(self, feature):
        """Get the features that provide a feature.

        Arguments:
            feature (FeatureOptions): Resquested feature.

        Returns:
            list[object]: List of feature objects.
        """
        return [obj for obj, provide in self._use if provide & feature == feature]

    def get_feature(self, feature, optional=False):
        """Get the feature that provide a feature, to be used when exactly one
        feature is expected.

        Arguments:
            feature (FeatureOptions): Resquested feature.
            optional (bool): if *True* and there is no feature, returns *None*.

        Returns:
            object: Worker object.
        """
        features = self.get_features(feature)
        if optional and not features:
            return None
        if len(features) != 1:
            name = self.options.name(feature)
            raise ValueError(f"expecting one {name!r} feature, found {features!r}")
        return features[0]

    @staticmethod
    def check_once(method):
        """Decorator that checks that the required features are registered
        before calling a method.

        Arguments:
            method (func): Method of a *BaseFeature* instance.
        """

        @wraps(method)
        def wrapper(inst, *args, **kwds):
            """wrapper"""
            if not inst._checked:
                inst.check_features()
            return method(inst, *args, **kwds)

        return wrapper
