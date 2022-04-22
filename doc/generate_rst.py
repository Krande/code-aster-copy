#!/usr/bin/env python3
# coding: utf-8

"""
Generate automodule blocks and fake libaster.
"""

EPILOG = """
EXAMPLES:
    Generate a fake libaster as python file::

        python3 generate_rst.py --libaster

    Generate a file like ``supervis.rst``::

        python3 generate_rst.py --manual code_aster/__init__.py code_aster/Supervis/*.py

    Generate files for *DataStructure* and derivated subclasses::

        python3 generate_rst.py --objects
"""

import argparse
import os
import os.path as osp
from collections import OrderedDict

from pydoc_pyrenderer import get_python_code

automodule_block = """.. automodule:: {0}
   :show-inheritance:
   :members:
   :special-members: __init__
""".format

autoclass_block = """.. autoclass:: code_aster.Objects.{0}
   :show-inheritance:
   :members:
""".format

auto_documentation = """.. AUTOMATICALLY CREATED BY generate_rst.py - DO NOT EDIT MANUALLY!

.. _devguide-{link}:

{intro}

{content}
""".format

title_ds = """
********************************************************************************
:py:class:`~code_aster.Objects.{0}` subclasses
********************************************************************************
""".format

title_ds_alone = """
********************************************************************************
:py:class:`~code_aster.Objects.{0}` object
********************************************************************************
""".format

subtitle = """
================================================================================
:py:class:`~code_aster.Objects.{0}` object
================================================================================
""".format


def automodule(filename):
    """Return a block to document a module."""
    name = osp.splitext(filename)[0].replace("/", ".")
    return automodule_block(name)


def all_objects(destdir):
    """Generate sphinx blocks for all libaster objects."""
    import code_aster.Objects as OBJ

    pyb_instance = OBJ.DataStructure.mro()[1]
    pyb_enum = OBJ.Physics.mro()[1]

    # sections: directly derivated from pybind11 instance
    sections = [OBJ.DataStructure, OBJ.GenericMaterialProperty]
    addsect = []
    for name, obj in list(OBJ.__dict__.items()):
        if not isinstance(obj, type):
            continue
        # if pyb_instance in obj.mro():
        # some objects are missed? or included via dependencies?
        # check AcousticDirichletBC for example
        if obj.mro()[1] is pyb_instance:
            # print("1:", name, obj.mro())
            addsect.append((name, obj))
    for _, obj in sorted(addsect):
        sections.append(obj)
    sections.append(pyb_enum)
    sections.append(Exception)
    # print(len(sections), "sections")

    # dict of subclasses
    dictobj = OrderedDict([(i, []) for i in sections])
    for name, obj in list(OBJ.__dict__.items()):
        # if obj is not OBJ.Material:
        #     continue
        if not isinstance(obj, type) or issubclass(
            obj, (OBJ.OnlyParallelObject, OBJ.InternalStateBuilder, OBJ.WithEmbeddedObjects)
        ):
            continue
        found = False
        for subtyp in sections:
            if issubclass(obj, subtyp) or obj is subtyp:
                dictobj[subtyp].append(name)
                found = True
                # print("Found:", name ,">>>", subtyp)
                break
        if not found and not issubclass(obj, OBJ.PyDataStructure):
            raise KeyError("pybind11 class not found: {0}".format(obj.mro()))

    dicttext = OrderedDict()
    for subtyp, objs in list(dictobj.items()):
        typename = subtyp.__name__
        lines = []
        # put subclass first
        try:
            objs.remove(typename)
        except ValueError:
            if subtyp not in (pyb_enum, Exception):
                print(subtyp, typename)
                print(objs)
                raise
        objs.sort()
        if subtyp not in (pyb_enum, Exception):
            objs.insert(0, typename)

        if len(objs) > 1:
            lines.append(title_ds(typename))
        else:
            lines.append(title_ds_alone(typename))
        for name in objs:
            if typename in ("DataStructure", "GenericMaterialProperty"):
                lines.append(subtitle(name))
            lines.append(autoclass_block(name))

        dicttext[typename] = os.linesep.join(lines)

    # generate a page for each of the first two classes
    with open(osp.join(destdir, "objects_datastructure.rst"), "w") as fobj:
        params = dict(link="objects_datastructure", content=dicttext["DataStructure"], intro="")
        fobj.write(auto_documentation(**params))

    with open(osp.join(destdir, "objects_materialbehaviour.rst"), "w") as fobj:
        params = dict(
            link="objects_materialbehaviour", content=dicttext["GenericMaterialProperty"], intro=""
        )
        fobj.write(auto_documentation(**params))

    # generate a page for all other classes
    with open(osp.join(destdir, "objects_others.rst"), "w") as fobj:
        params = dict(
            link="objects_others",
            content=os.linesep.join(list(dicttext.values())[2:]),
            intro="""
####################################
Index of all other available objects
####################################

Documentation of all other types.
""",
        )
        fobj.write(auto_documentation(**params))


def build_pylibaster(filename):
    """Create a fake libaster as Python file with only signatures and docstrings.

    Arguments:
        filename (str): Destination file
    """
    # do not import code_aster not to extend objects
    import libaster

    pyb_instance = libaster.DataStructure.mro()[1]
    blocks = []
    for name, obj in list(libaster.__dict__.items()):
        export = True
        if isinstance(obj, type):
            export = pyb_instance in obj.mro() or Exception in obj.mro()
        else:
            if "builtin_function_or_method" not in repr(type(obj)):
                export = False
        if export:
            blocks.append(get_python_code("libaster." + name))
        else:
            print("not exported:", name, obj)

    with open(filename, "w") as flib:
        flib.write("\n".join(blocks))


def main():
    default_dest = osp.join(osp.dirname(__file__), "devguide")
    parser = argparse.ArgumentParser(
        description=__doc__, epilog=EPILOG, formatter_class=argparse.RawTextHelpFormatter
    )
    # libaster and objects can be run in the same process
    parser.add_argument(
        "--libaster",
        action="store_const",
        dest="action",
        const="libaster",
        help="build fake libaster as a pure Python file",
    )
    parser.add_argument(
        "--objects",
        action="store_const",
        dest="action",
        const="objects",
        help="for C++ only objects (needs to import libaster)",
    )
    parser.add_argument(
        "--manual", action="store_const", dest="action", const="manual", help="for only few objects"
    )
    parser.add_argument(
        "-d",
        "--destdir",
        action="store",
        metavar="DIR",
        default=default_dest,
        help="directory where `rst` files will be written",
    )
    parser.add_argument("file", metavar="FILE", nargs="*", help="file to analyse")
    args = parser.parse_args()

    pylib = osp.join(osp.dirname(__file__), "_fake", "libaster.py")
    if args.action == "libaster":
        build_pylibaster(pylib)
    if args.action == "objects":
        all_objects(args.destdir)
    elif args.action == "manual":
        for name in args.file:
            print(automodule(name))


if __name__ == "__main__":
    main()
