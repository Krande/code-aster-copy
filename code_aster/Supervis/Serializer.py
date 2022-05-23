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
:py:mod:`Serializer` --- Serialization of code_aster objects
************************************************************

code_aster objects are saved and reloaded using the *pickle* protocol.

:func:`saveObjects` does the saving of objects which are available in the user
context (in which :func:`saveObjects` is called).
The command :func:`~code_aster.Commands.FIN` automatically calls this function.

Objects are reloaded by the function :func:`loadObjects` called just after the
Jeveux database has been reloaded by :func:`~code_aster.Commands.debut.init`
if the ``--continue`` argument is passed.
The command :func:`~code_aster.Commands.debut.POURSUITE` does also the same.

To be properly pickled, the classes defined by the user must be inherit
from :class:`~code_aster.Objects.user_extensions.WithEmbeddedObjects`.
"""

import gc
import os.path as osp
import pickle
import re
import traceback
import types
from hashlib import sha256
from io import IOBase

import libaster
import numpy

from .. import Objects
from ..Objects import DataStructure, InternalStateBuilder, ResultNaming, WithEmbeddedObjects
from ..Utilities import (
    DEBUG,
    MPI,
    ExecutionParameter,
    Options,
    get_caller_context,
    logger,
    no_new_attributes,
)

ARGS = "_MARK_DS_ARGS_"
STATE = "_MARK_DS_STATE_"
LIST = "_MARK_LIST_"
DICT = "_MARK_DICT_"

# same values in op9999
class FinalizeOptions:
    """Options for closure."""

    SaveBase = 1
    InfoResu = 2
    Repack = 4
    OnlyProc0 = 8


class Serializer(object):

    """This class manages 'save & reload' feature.

    Arguments:
        context (dict): The context to be saved or in which the loaded objects
            have to been set.

    Attributes:
        _ctxt (dict): Working context.
        _pick_filename (str): Filename of the pickle file.
        _base (str): Filename of the Jeveux database file.
        _sha_filename (str): Filename containing the SHA of the both previous
            files.
    """

    _sha_filename = "pick.code_aster.sha"
    _pick_filename = "pick.code_aster.objects"
    _info_filename = "pick.code_aster.infos"
    _base = "glob.1"
    _ctxt = None

    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, context=None):
        """Initialization
        :param context: context to save or in which loading objects
        :type context: dict
        """
        self._ctxt = context

    @classmethod
    def canRestart(cls):
        """Tell if a restart is possible.
        This means that glob & pickled files are consistent.

        Returns:
            bool: *True* if the previous execution can be continued.
        """
        for fname in (cls._base, cls._pick_filename, cls._info_filename, cls._sha_filename):
            if not osp.exists(fname):
                logger.error("Can not restart, no such file: %s", fname)
                return False

        sign = read_signature(cls._sha_filename)
        if len(sign) != 3:
            logger.error("Invalid sha file: %r", cls._sha_filename)
            return False
        ref_pick, ref_info, ref_base = sign
        sign_pick = file_signature(cls._pick_filename)
        if sign_pick != ref_pick:
            logger.error("Current pickled file: %s", sign_pick)
            logger.error("Expected signature  : %s", ref_pick)
            logger.error("The %r file is not the expected one.", cls._pick_filename)
            return False

        sign_info = file_signature(cls._info_filename)
        if sign_info != ref_info:
            logger.error("Current info file : %s", sign_info)
            logger.error("Expected signature: %s", ref_info)
            logger.error("The %r file is not the expected one.", cls._info_filename)
            return False

        sign_base = file_signature(cls._base, 0, 8000000)
        if sign_base != ref_base:
            logger.error("Current base file : %s", sign_base)
            logger.error("Expected signature: %s", ref_base)
            logger.error("The %r file is not the expected one.", cls._base)
            return False
        return True

    def save(self):
        """Save objects of the context.

        Returns:
            list[str]: List of the names of the DataStructure objects actually
            saved.
        """
        assert self._ctxt is not None, "context is required"
        ctxt = _filteringContext(self._ctxt)
        saved = []
        with open(self._pick_filename, "wb") as pick:
            pickler = AsterPickler(pick)
            # ordered list of objects names
            logger.info("Saving objects...")
            objList = []
            for name, obj in ctxt.items():
                if name == "CO" or obj is logger:
                    continue
                try:
                    logger.info("%-24s %s", name, type(obj))
                    pickler.save_one(obj)
                    objList.append(name)
                except Exception:
                    logger.warning("object can not be pickled: %s", name)
                    logger.debug(traceback.format_exc())
                    continue
                if isinstance(obj, DataStructure):
                    saved.append(name)

        logger.debug("Objects saved: %s", objList)
        with open(self._info_filename, "wb") as pick:
            # add management objects on the stack
            pickle.dump(objList, pick)
            pickle.dump(ResultNaming.getCurrentName(), pick)
        return saved

    def sign(self):
        """Sign the saved files and store their SHA signatures."""
        with open(self._sha_filename, "wb") as pick:
            pickler = pickle.Pickler(pick)

            sign_pick = file_signature(self._pick_filename)
            logger.info("Signature of pickled file   : %s", sign_pick)
            sign_info = file_signature(self._info_filename)
            logger.info("Signature of info file      : %s", sign_info)
            sign_base = file_signature(self._base, 0, 8000000)
            logger.info("Signature of Jeveux database: %s", sign_base)

            pickler.dump(sign_pick)
            pickler.dump(sign_info)
            pickler.dump(sign_base)

    def load(self):
        """Load objects into the context."""
        assert self._ctxt is not None, "context is required"
        with open(self._info_filename, "rb") as pick:
            # add management objects on the stack
            objList = pickle.load(pick)
            lastId = int(pickle.load(pick), 16)
        # restore the objects counter
        ResultNaming.initCounter(lastId)

        should_fail = ExecutionParameter().option & Options.StrictUnpickling
        pool = objList[:]
        logger.debug("Objects pool: %s", pool)
        with open(self._pick_filename, "rb") as pick:
            unpickler = AsterUnpickler(pick)
            # load all the objects
            objects = []
            names = []
            try:
                while True:
                    name = pool.pop(0) if pool else None
                    logger.debug("loading: %s...", name)
                    try:
                        obj = unpickler.load_one()
                        logger.debug("object restored: %s %s...", name, type(obj))
                    except Exception as exc:
                        if isinstance(exc, EOFError):
                            raise
                        logger.info(traceback.format_exc())
                        logger.info("can not restore object: %s", name)
                        if should_fail:
                            raise
                        continue
                    names.append(name)
                    objects.append(obj)
            except EOFError:
                pass

        not_read = set(objList).difference(names)
        if not_read:
            logger.warning("These objects have not been reloaded: %s", tuple(not_read))
        logger.info("Restored objects:")
        for name, obj in zip(names, objects):
            logger.debug("restoring %s...", name)
            try:
                obj = _restore(name, obj)
            except Exception:
                if should_fail:
                    raise
                logger.warning(traceback.format_exc())
                logger.error("can not restore object: %s <%s>", name, obj)
                continue
            self._ctxt[name] = obj
            if not name:
                logger.warning("restoring %s with name 'None'!", type(obj))
                continue
            if isinstance(obj, DataStructure):
                obj.userName = name
            logger.info("%-24s %s", name, type(obj))
            assert not isinstance(obj, AsterUnpickler.BufferObject)


def _restore(name, obj):
    """Build instance from BufferObject.

    Arguments:
        name (str): Object name, only for debugging.
        obj (*misc*): Object in which *DataStructure* objects will be instanciated.

    Returns:
        *misc*: Same object with *DataStructure* objects.
    """
    if isinstance(obj, list):
        return [_restore(name, i) for i in obj]
    if isinstance(obj, tuple):
        return tuple(_restore(name, list(obj)))
    if isinstance(obj, dict):
        return obj.__class__([(i, _restore(name, obj[i])) for i in obj])
    if isinstance(obj, InternalStateBuilder):
        obj._st = _restore(f"{name}._st", obj._st)
        return obj
    if isinstance(obj, AsterUnpickler.BufferObject):
        return obj.instance
    if isinstance(obj, WithEmbeddedObjects):
        logger.info("restoring user object %s, attrs: %s", name, obj.aster_embedded)
        for attr in obj.aster_embedded:
            logger.debug("attr: %s, was: %s", attr, getattr(obj, attr))
            setattr(obj, attr, _restore(name, getattr(obj, attr)))
            logger.debug("new: %s", getattr(obj, attr))
        return obj
    return obj


def saveObjects(level=1, delete=True, options=0):
    """Save objects of the caller context in a file.

    Arguments:
        level (int): Number of frames to go back to find the user context.
        delete (bool): If *True* the saved objects are deleted from the context.
        options (*FinalizeOptions*): Options for finalization.
    """
    gc.collect()
    options |= FinalizeOptions.SaveBase
    rank = MPI.COMM_WORLD.Get_rank()
    if options & FinalizeOptions.OnlyProc0 and rank != 0:
        logger.info("Objects not saved on processor #{0}".format(rank))
        libaster.jeveux_finalize(FinalizeOptions.OnlyProc0)
        return

    context = get_caller_context(level)

    if options & FinalizeOptions.InfoResu:
        for name, obj in context.items():
            if hasattr(obj, "printInfo"):
                libaster.write("\n ======> " + name)
                obj.printInfo()

    # if ExecutionParameter().option & Options.Debug:
    #     libaster.debugJeveuxContent("Saved jeveux objects:")

    # orig = logger.getEffectiveLevel()
    # logger.setLevel(DEBUG)
    pickler = Serializer(context)
    saved = pickler.save()

    # close Jeveux files (should not be done before pickling)
    libaster.jeveux_finalize(options)
    pickler.sign()
    # logger.setLevel(orig)

    if delete:
        # Remove the objects from the context
        for name in saved:
            context[name] = None


def loadObjects(level=1):
    """Load objects from a file in the caller context.

    Arguments:
        level (int): Number of frames to go back to find the user context.
    """
    context = get_caller_context(level)

    # if ExecutionParameter().option & Options.Debug:
    #     libaster.debugJeveuxContent("Reloaded jeveux objects:")
    # orig = logger.getEffectiveLevel()
    # logger.setLevel(DEBUG)
    Serializer(context).load()
    # logger.setLevel(orig)


def contains_datastructure(sequence):
    """Tell if a sequence contains a DataStructure.

    Arguments:
        sequence (*iterable*): List-like object.

    Returns:
        bool: *True* if *sequence* contains a *DataStructure*,
        *False* otherwise.
    """
    for item in sequence:
        if isinstance(item, DataStructure):
            return True
    return False


class AsterPickler(pickle.Pickler):

    """Adapt pickling of DataStructure objects.

    In the Python namespace, DataStructures are wrappers on *shared_ptr* through
    *pybind11* instances. So there are several *pointers* for the same instance.
    Standard pickling creates new objects for each *pointers* and during
    unpickling this creates new *pybind11* instance for each Python wrapper.
    To avoid that, the pickling step only saves arguments (returned by
    :py:meth:`__getinitargs__`), a state (created by :py:meth:`__getstate__`
    as an instance of
    :py:class:`~code_aster.Objects.Serialization.InternalStateBuilder`)
    and an identifier of the DataStructure (its Jeveux name).

    See :py:class:`.AsterUnpickler` for unpickling phase.

    See *Pickling and unpickling external objects* from the :py:mod:`pickle`
    documentation.
    """

    def __init__(self, *args, **kwargs):
        pickle.Pickler.__init__(self, *args, **kwargs)
        self._memods = set()
        self._depth = 0

    def save_one(self, obj):
        """Save one object.

        Arguments:
            obj (*misc*): Object to save.
        """
        self._depth += 1
        logger.debug("save_one: %s / %s", self._depth, obj)
        if isinstance(obj, (list, tuple)):
            if obj and contains_datastructure(obj):
                self.save_one(LIST)
                self.save_one(len(obj))
                for item in obj:
                    self.save_one(item)
        elif isinstance(obj, dict):
            if obj and contains_datastructure(obj.values()):
                self.save_one(DICT)
                self.save_one(len(obj))
                for item in obj.values():
                    self.save_one(item)
        elif isinstance(obj, InternalStateBuilder):
            self.save_one(STATE)
            flat = obj.flat()
            self.save_one(len(flat))
            for item in flat:
                self.save_one(item)
        elif isinstance(obj, DataStructure):
            ds_id = unique_id(obj)
            if ds_id not in self._memods:
                self._memods.add(ds_id)
                name = obj.getName()
                # save initial arguments
                if hasattr(obj, "__getinitargs__"):
                    init_args = obj.__getinitargs__()
                    assert isinstance(init_args, tuple), (name, init_args)
                else:
                    init_args = ()
                self.save_one(ARGS)
                self.save_one(len(init_args))
                logger.debug("saved initargs: len %d: %s", len(init_args), init_args)
                self.save_one(name)
                for item in init_args:
                    self.save_one(item)
                # save state
                state = obj.__getstate__()
                assert isinstance(state, InternalStateBuilder), state
                logger.debug("saved state: %s", state)
                self.save_one(state)
            else:
                logger.debug("skip object %s", ds_id)

        self.dump(obj)
        self._depth -= 1

    def persistent_id(self, obj):
        """Compute a persistent id for DataStructure.

        Returns:
            str: Identifier containing " mark, class name, object name".
        """
        if isinstance(obj, DataStructure):
            pers_id = unique_id(obj)
            logger.debug("persistent id: %s", pers_id)
            return pers_id

        # other objects, pickled as usual
        return None


def unique_id(obj):
    """Compute a unique id for DataStructure, used as *persistent_id* for
    pickle.

    Returns:
        str: Identifier containing " mark, class name, object name".
    """
    class_ = type(obj).__name__
    pers_id = ("DataStructure", class_, obj.getName().strip())
    return str(pers_id)


class AsterUnpickler(pickle.Unpickler):

    """Adapt unpickling of DataStructure objects.

    See :py:class:`.AsterPickler` for pickling phase.

    During unpickling, only :py:class:`.BufferObject` are created from the
    persistent identifiers that are reloaded.
    Only after this step, all *BufferObjects* are reloaded.
    The instances are created on demand. The parent instances are
    also automatically and recursively created when needed.

    See *Pickling and unpickling external objects* from the :py:mod:`pickle`
    documentation.
    """

    class BufferObject(object):

        """This class defines a temporary object that is created instead of
        an instance of DataStructure.

        Attributes:
            _name (str): *Jeveux* name of the object.
            _args (tuple[misc]): Initial arguments to pass to the constructor.
            _state (tuple[misc]): Arguments pass to the ``__setstate__`` method
                if it exists.
            _class (str): Class name of the *DataStructure* (in Objects module).
            _inst (*DataStructure*): Cached object.
        """

        def __init__(self, name):
            self._name = name
            self._args = None
            self._state = None
            self._class = None
            self._inst = None

        def __repr__(self):
            return "Buffer:{0._name!r}<{0._class}>".format(self)

        @property
        def args(self):
            """Arguments to pass to the DataStructure's ``__init__``."""
            assert self._args is not None, self._name
            return self._args

        @args.setter
        def args(self, args_):
            """Register the initial arguments."""
            self._args = args_

        @property
        def state(self):
            """Arguments to pass to the DataStructure's ``__setstate__``."""
            assert self._state is not None, self._name
            return self._state

        @state.setter
        def state(self, state_):
            """Register the object state to restore."""
            self._state = state_

        @property
        def classname(self):
            """Class name of the DataStructure to create."""
            assert self._class is not None, self._name
            return self._class

        @classname.setter
        def classname(self, classname_):
            """Register the name of the class to create."""
            self._class = classname_

        @property
        def instance(self):
            """Return the instance, build it if it does not yet exist.

            If DataStructure objects are present in values returned by
            :py:meth:`__getinitargs__` they must bedirectly present in the
            returned tuple (not in sub-objects).

            Returns:
                misc: DataStructure object.
            """
            if self._inst is None:
                logger.debug("building %r of type %s", self._name, self._class)
                # DataStructure must be in args (not in sub-objects)
                args = [i.instance if isinstance(i, type(self)) else i for i in self.args]
                logger.debug("initargs: %s", args)
                self._inst = getattr(Objects, self.classname)(*args)
                logger.debug("new object: %s", self._inst)
                logger.debug("setting state: %s", self.state)
                _restore(f"{self._name} .state", self.state)
                getattr(self._inst, "__setstate__")(self.state)
                if hasattr(self._inst, "update"):
                    self._inst.update()
            return self._inst

    class BufferStack(object):

        """Helper object to store *BufferObject* objects."""

        def __init__(self):
            self._store = {}

        def __iter__(self):
            for name in self._store:
                yield name

        def buffer(self, name):
            name = name.strip()
            if name not in self._store:
                self._store[name] = AsterUnpickler.BufferObject(name)
            return self._store[name]

    def __init__(self, fileobj):
        pickle.Unpickler.__init__(self, fileobj)
        self._stack = AsterUnpickler.BufferStack()
        self._depth = 0

    def load_one(self):
        """Load one object.

        Returns:
            *misc*: Loaded object.
        """
        obj = self.load()
        self._depth += 1
        logger.debug("load_one: %s / %s", self._depth, obj)
        if not isinstance(obj, str):
            self._depth -= 1
            return obj
        if obj == LIST:
            size = self.load_one()
            for _ in range(size):
                self.load_one()
            logger.debug("list: %s / %s", self._depth, obj)
            obj = self.load_one()
        if obj == DICT:
            size = self.load_one()
            for _ in range(size):
                self.load_one()
            logger.debug("dict: %s / %s", self._depth, obj)
            obj = self.load_one()
        if obj == ARGS:
            nbobj = self.load_one()
            name = self.load_one()
            init_args = []
            for _ in range(nbobj):
                init_args.append(self.load_one())
            buffer = self._stack.buffer(name)
            buffer.args = init_args
            logger.debug("loaded initargs: %s", init_args)
            # expecting the STATE mark
            mark = self.load_one()
            assert mark == STATE, mark
            nbobj = self.load_one()
            for _ in range(nbobj):
                self.load_one()
            try:
                buffer.state = self.load_one()
            except:
                logger.info(traceback.format_exc())
                logger.info("internal state can not be loaded")
                buffer.state = None
                raise pickle.PicklingError("internal state can not be loaded")
            # 'load' will call 'persistent_load'
            obj = self.load_one()
            logger.debug("loaded DataStructure %r", obj)
        self._depth -= 1
        return obj

    def recover_ds(self, class_id, key_id):
        """Create a new *BufferObject*.

        Arguments:
            class_id (str): Name of the DataStructure type.
            key_id (str): Name of the *Jeveux* object.

        Returns:
            *BufferObject*
        """
        buffer = self._stack.buffer(key_id)
        buffer.classname = class_id
        return buffer

    def persistent_load(self, pers_id):
        """Action called when a persistent id is reloaded.

        Arguments:
            str: Identifier of the DataStructure.

        Returns:
            *BufferObject*: New object to store DataStructure informations.
        """
        decoded_id = eval(pers_id)
        logger.debug("persistent id loaded: %s", decoded_id)
        type_tag, class_id, key_id = decoded_id
        if type_tag == "DataStructure":
            return self.recover_ds(class_id, key_id)
        else:
            raise pickle.UnpicklingError("unsupported persistent object")


def _filteringContext(context):
    """Return a context by filtering the input objects by excluding:
    - modules,
    - code_aster objects,
    - ...

    Arguments:
        context (dict): Context to be filtered.

    Returns:
        dict: New cleaned context.
    """
    # functions to be ignored
    ignored = ("code_aster", "DETRUIRE", "FIN", "VARIABLE")
    re_system = re.compile("^__.*__$")
    ctxt = {}
    for name, obj in list(context.items()):
        if not name or name in ignored or re_system.search(name):
            continue
        if getattr(numpy, name, None) is obj:  # see issue29282
            continue
        # check attr needed for python<=3.6
        if hasattr(obj, "__class__") and isinstance(obj, IOBase):
            continue
        if type(obj) in (
            types.ModuleType,
            type,
            types.MethodType,
            types.FunctionType,
            types.BuiltinMethodType,
            types.BuiltinFunctionType,
        ):
            continue
        # aster-legacy try 'dumps' before keeping the object
        ctxt[name] = obj
    return ctxt


def read_signature(sha_file):
    """Read the signatures from the file containing the SHA strings.

    Arguments:
        sha_file (str): Filename of the pickled file to read.

    Returns:
        list[str]: List of the two signatures as SHA256 strings that identify
            the pickled file and the (first) Jeveux database file.
    """
    sign = []
    try:
        with open(sha_file, "rb") as pick:
            sign.append(pickle.Unpickler(pick).load())
            sign.append(pickle.Unpickler(pick).load())
            sign.append(pickle.Unpickler(pick).load())
    except Exception:
        traceback.print_exc()
    logger.debug("pickled signatures: %s", sign)
    return sign


def file_signature(filename, offset=0, bufsize=-1):
    """Compute a signature of the file.

    Arguments:
        filename (str): File to sign.
        offset (int): Offset before reading content (default to 0).
        bufsize (int): Size of the content to read (default to -1, content is
            read up to EOF).

    Returns:
        str: Signature as SHA256 string to identify the file.
    """
    try:
        with open(filename, "rb") as fobj:
            fobj.seek(offset, 0)
            sign = sha256(fobj.read(bufsize)).hexdigest()
    except Exception:
        traceback.print_exc()
    return sign
