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


class FeatureMeta(type):
    """Metaclass to set the list of supported features."""

    def __new__(mcs, name, base, defcl):
        obj = super().__new__(mcs, name, base, defcl)
        # obj.required_features = defcl.get("required_features", getattr(base, "required_features"))
        # obj.optional_features = defcl.get("optional_features", getattr(base, "optional_features"))
        # compute '_supported' attribute value
        obj._supported = getattr(base, "_supported", 0)
        for feat in obj.required_features + obj.optional_features:
            obj._supported |= feat
        return obj


class BaseFeature(metaclass=FeatureMeta):
    """A `BaseFeature` is defined by an object that does the job and the kind of
    services it provides. A feature *uses* other features.
    """

    provide = 0
    _supported = 0
    required_features = []
    optional_features = []
    # for no_new_attributes
    _use = _checked = None

    def __init__(self):
        super().__init__()
        self._use = []
        self._checked = False

    def use(self, obj, provide=0):
        """Add a feature to be used.

        The provided services are defined by the
        'provide' attribute if `obj` is a `BaseFeature` + the `provide` argument.

        Arguments:
            obj (BaseFeature|misc): Object that to be used.
            provide (FeaturesOptions, optional): Features provided by the object.
        """
        if type(obj) in (list, tuple):
            assert provide, "provide is required using list of features"
            return [self.use(obj_i, provide) for obj_i in obj]
        if not obj or obj in [reg for reg, _ in self._use]:
            return
        try:
            provide |= obj.provide
        except AttributeError:
            pass
        if not provide:
            raise ValueError("This object provides no feature, use 'provide'.")
        if not provide & self._supported:
            raise TypeError(f"{self.__class__.__name__} does not support {provide!r}")
        self._use.append((obj, provide))

    def check_features(self):
        """Check that required features are defined."""
        for feat in self.required_features:
            if not self.has_feature(feat):
                raise TypeError(f"{self.__class__.__name__} requires the {feat!r} feature")
        self._checked = True

    def has_feature(self, feature):
        """Tell if the feature is available.

        Arguments:
            feature (FeatureOptions): Requested feature.

        Returns:
            bool: *True* if the feature is available, *False* otherwise.
        """
        return bool(self.get_features(feature))

    def undefined(self):
        """Returns the list of the undefined features.

        For each option, a bool indicates if the feature is required or not.

        Returns:
            list[int, bool]: List of undefined options. bool is *True* for required
            options, *False* otherwise.
        """
        undef = [(feat, True) for feat in self.required_features if not self.has_feature(feat)]
        undef.extend(
            [(feat, False) for feat in self.optional_features if not self.has_feature(feat)]
        )
        return undef

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
            feature (FeatureOptions): Requested feature.

        Returns:
            list[object]: List of feature objects.
        """
        return [obj for obj, provide in self._use if provide & feature == feature]

    def get_feature(self, feature, optional=False):
        """Get the feature that provide a feature, to be used when exactly one
        feature is expected.

        Arguments:
            feature (FeatureOptions): Requested feature.
            optional (bool): if *True* and there is no feature, returns *None*.

        Returns:
            object: Worker object.
        """
        features = self.get_features(feature)
        if optional and not features:
            return None
        if len(features) != 1:
            raise ValueError(f"expecting one {feature!r} feature, found {features!r}")
        return features[0]

    def get_childs(self, feature):
        """Get the features that provide a feature from all known features from childs.

        Arguments:
            feature (FeatureOptions): Requested feature.

        Returns:
            list[object]: List of feature objects. If *feature* is None, it returns
            list of tuple (object, provide value).
        """
        all_feat = self._use[:]
        for obj, _ in self._use:
            if not hasattr(obj, "get_childs"):
                continue
            all_feat.extend(obj.get_childs(None))
        if not feature:
            return all_feat
        return [obj for obj, provide in all_feat if provide & feature == feature]

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


class Observer:
    """The Observer interface declares the `notify` method, used by events."""

    def notify(self, event):
        """Receive notification from event.

        Arguments:
            event (EventSource): Object that sends the notification.
        """
        # calls event.get_state()
        raise NotImplementedError("must be subclassed")


class EventSource:
    """The EventSource interface declares a set of methods for managing observers."""

    # for no_new_attributes
    _observers = None

    def __init__(self) -> None:
        super().__init__()
        self._observers = []

    def add_observer(self, observer):
        """Attach an observer to the event.

        Arguments:
            observer (Observer): Observer object to be added.
        """
        self._observers.append(observer)

    def remove_observer(self, observer):
        """Detach an observer from the event.

        Arguments:
            observer (Observer): Observer object to be removed.
        """
        self._observers.remove(observer)

    def notifyObservers(self):
        """Notify all observers about an event."""
        for obs in self._observers:
            obs.notify(self)

    def get_state(self):
        """Returns the current state to be shared with observers."""
        raise NotImplementedError("must be subclassed")
