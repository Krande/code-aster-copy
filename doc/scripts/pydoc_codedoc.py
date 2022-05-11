#!/usr/bin/env python3

"""This module extracts docstrings from pybind11 objects to generate static
Python files with the same docstrings.

It defines a specialized renderer based on `pydoc` ones.
"""

# pydoc_codedoc Copyright 2022 Mathieu Courtois <mathieu.courtois.at.gmail.com>
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.

# Python and pydoc are provided under the terms of the Python License, Version 2,
# Copyright (c) 2001, 2002, 2003, 2004, 2005, 2006 Python Software Foundation.

import builtins
import inspect
import re
import sys
from collections import deque
from pydoc import (
    TextDoc,
    _is_bound_method,
    _split_list,
    classify_class_attrs,
    classname,
    getdoc,
    plaintext,
    render_doc,
    sort_attributes,
    visiblename,
)


# The code of 'CodeDoc' is a copy from 'TextDoc' with the less changes as
# possible to easily be updated with newer versions.
#
# Changes are commented by '# CodeDoc:'


class CodeDoc(TextDoc):
    """Subclass `pydoc.TextDoc` renderer to generate equivalent Python code
    (only signature and docstring) from pybind11 objects.

    This is useful to publish use Sphinx autodoc capabilities on static Python
    files, without compiling huge pybind11 extension.
    """

    # CodeDoc: as for _PlainTextDoc
    def bold(self, text):
        return text

    def docclass(self, object, name=None, mod=None, *ignored):
        """Produce text documentation for a given class object."""
        realname = object.__name__
        name = name or realname
        bases = object.__bases__

        def makename(c, m=object.__module__):
            return classname(c, m)

        # CodeDoc: declare class + hide inheriting from pybind11
        title = "class " + self.bold(realname)
        if bases:
            parents = map(makename, bases)
            parents = [i.replace("builtins.", "") for i in parents if not i.startswith("pybind11")]
            if parents:
                title = title + "(%s)" % ", ".join(parents)
        title += ":"

        contents = []
        push = contents.append

        try:
            signature = inspect.signature(object)
        except (ValueError, TypeError):
            signature = None
        if signature:
            argspec = str(signature)
            if argspec and argspec != "()":
                push(name + argspec + "\n")

        # CodeDoc: as multilines string
        doc = _format_doc(getdoc(object) or "")
        if doc:
            push(doc + "\n")

        # List the mro, if non-trivial.
        mro = deque(inspect.getmro(object))
        if len(mro) > 2:
            # CodeDoc: keep as comment
            push("# Method resolution order:")
            for base in mro:
                push("#     " + makename(base))
            push("")

        # CodeDoc: keep as comment
        # Cute little class to pump out a horizontal rule between sections.
        class HorizontalRule:
            def __init__(self):
                self.needone = 0

            def maybe(self):
                if self.needone:
                    push("#" + "-" * 70)
                self.needone = 1

        hr = HorizontalRule()

        def spill(msg, attrs, predicate):
            ok, attrs = _split_list(attrs, predicate)
            # CodeDoc: ignore inherited objects
            if ok and not _inherited_from(msg):
                hr.maybe()
                # CodeDoc: keep as comment
                push("# " + msg)
                for name, kind, homecls, value in ok:
                    try:
                        value = getattr(object, name)
                    except Exception:
                        # Some descriptors may meet a failure in their __get__.
                        # (bug #1785)
                        push(self._docdescriptor(name, value, mod))
                    else:
                        push(self.document(value, name, mod, object))
            return attrs

        def spilldescriptors(msg, attrs, predicate):
            ok, attrs = _split_list(attrs, predicate)
            # CodeDoc: ignore inherited objects
            if ok and not _inherited_from(msg):
                hr.maybe()
                push("# " + msg)
                for name, kind, homecls, value in ok:
                    push(self._docdescriptor(name, value, mod))
            return attrs

        def spilldata(msg, attrs, predicate):
            ok, attrs = _split_list(attrs, predicate)
            if ok:
                hr.maybe()
                # CodeDoc: keep as comment
                push("# " + msg)
                for name, kind, homecls, value in ok:
                    if callable(value) or inspect.isdatadescriptor(value):
                        doc = getdoc(value)
                    else:
                        doc = None
                    try:
                        obj = getattr(object, name)
                    except AttributeError:
                        obj = homecls.__dict__[name]
                    # CodeDoc: do not export data, except pybind11 Enums
                    if hasattr(obj, "name") and hasattr(obj, "value"):
                        if obj.name not in builtins.__dict__:
                            push(obj.name + " = " + str(obj.value) + "\n")
            return attrs

        attrs = [
            (name, kind, cls, value)
            for name, kind, cls, value in classify_class_attrs(object)
            if visiblename(name, obj=object)
        ]

        while attrs:
            if mro:
                thisclass = mro.popleft()
            else:
                thisclass = attrs[0][2]
            attrs, inherited = _split_list(attrs, lambda t: t[2] is thisclass)

            if thisclass is builtins.object:
                attrs = inherited
                continue
            elif thisclass is object:
                tag = "defined here"
            else:
                tag = "inherited from %s" % classname(thisclass, object.__module__)

            sort_attributes(attrs, object)

            # Pump out the attrs, segregated by kind.
            attrs = spill("Methods %s:\n" % tag, attrs, lambda t: t[1] == "method")
            attrs = spill("Class methods %s:\n" % tag, attrs, lambda t: t[1] == "class method")
            attrs = spill("Static methods %s:\n" % tag, attrs, lambda t: t[1] == "static method")
            attrs = spilldescriptors(
                "Data descriptors %s:\n" % tag, attrs, lambda t: t[1] == "data descriptor"
            )
            attrs = spilldata(
                "Data and other attributes %s:\n" % tag, attrs, lambda t: t[1] == "data"
            )

            assert attrs == []
            attrs = inherited

        contents = "\n".join(contents)
        if not contents:
            return title + "\n"
        # CodeDoc: only spaces
        return title + "\n" + self.indent(contents.rstrip(), "    ") + "\n"

    def docroutine(self, object, name=None, mod=None, cl=None):
        """Produce text documentation for a function or method object."""
        realname = object.__name__
        name = name or realname
        note = ""
        skipdocs = 0
        if _is_bound_method(object):
            imclass = object.__self__.__class__
            if cl:
                if imclass is not cl:
                    note = " from " + classname(imclass, mod)
            else:
                if object.__self__ is not None:
                    note = " method of %s instance" % classname(object.__self__.__class__, mod)
                else:
                    note = " unbound %s method" % classname(imclass, mod)
        # CodeDoc: skipped
        # exception for "method of builtins.PyCapsule" as created by pybind
        if note and "method of builtins.PyCapsule" not in note:
            return ""
        note = ""

        # CodeDoc: just keep the name
        if True or name == realname:
            title = self.bold(realname)
        else:
            if cl and inspect.getattr_static(cl, realname, []) is object:
                skipdocs = 1
            title = self.bold(name) + " = " + realname
        # CodeDoc: declare routine
        title = "def " + title
        argspec = None

        if inspect.isroutine(object):
            try:
                signature = inspect.signature(object)
            except (ValueError, TypeError):
                signature = None
            if signature:
                argspec = str(signature)
                if realname == "<lambda>":
                    title = self.bold(name) + " lambda "
                    # XXX lambda's won't usually have func_annotations['return']
                    # since the syntax doesn't support but it is possible.
                    # So removing parentheses isn't truly safe.
                    argspec = argspec[1:-1]  # remove parentheses
        if not argspec:
            argspec = "(...)"

        # CodeDoc: search for argspec from pybind11 docstring
        if skipdocs:
            doc = ""
        else:
            doc = getdoc(object) or ""

        if "/," in argspec and sys.version_info < (3, 8):
            argspec = argspec.replace("/,", "")
        if argspec == "(...)":
            argspec, doc = pybind_args(realname, doc)

        doc = self.indent(_format_doc(doc)).rstrip() + "\n"

        argspec += ":"
        decl = title + argspec + note
        return decl + "\n" + doc

    def _docdescriptor(self, name, value, mod):
        results = []
        push = results.append

        # CodeDoc: define properties
        # why docproperty is not called?
        if name:
            if name in ("__weakref__",):
                return ""
            assert isinstance(value, property), "unsupported: %s: %s" % (name, value)
            push("@property" + "\n")
            push("def " + self.bold(name) + "(self):")
            # push(self.bold(name))
            push("\n")
        doc = _format_doc(getdoc(value) or "")
        if doc:
            push(self.indent(doc))
            push("\n")
        return "".join(results)


# CodeDoc: add helper functions
def _inherited_from(msg):
    """use to ignore inherited methods"""
    return "inherited from" in msg


def pybind_args(name, doc, remove_sign=True):
    """Extract argspec from pybind docstring.

    Arguments:
        name (str): routine name.
        doc (str): docstring.
        remove_sign (bool): if *True*, the pybind signature is removed from
            the docstring.

    Returns:
        str: representation of the arguments.
    """
    if "Overloaded function" in doc:
        # remove all before 'Overloaded function'
        doc = doc[doc.index("Overloaded function") :]
        remove_sign = False
        # check if it is a function or method
        argstr = ("self, " if "self:" in doc else "") + "*args, **kwargs"
    else:
        expr = re.compile(name + r"\((.*?)\) *->.*$", re.M)
        mat = expr.search(doc)
        if not mat:
            return "(self)", doc
        argstr = mat.group(1)
        # remove template arguments
        argstr = _remove_enclosed(argstr, "<")
        argstr = re.sub(": .*?(?P<sep>[=,]|$)", r"\g<sep>", argstr)
        # default value has been removed
        argstr = re.sub("= *(?P<sep>, *|$)", r"\g<sep>", argstr)
    if remove_sign:
        doc = expr.sub("", doc).lstrip()
    return "(" + argstr + ")", doc


def _remove_enclosed(string, mark):
    end = ">]"["<[".index(mark)]
    contains = re.compile("{beg}.*{end}".format(beg=mark, end=end), re.M)
    expr = re.compile("{beg}([^{beg}{end}]+){end}".format(beg=mark, end=end), re.M)
    while contains.search(string):
        string = expr.sub("", string)
    return string


def _format_doc(doc):
    if doc:
        doc = '"""' + doc + "\n" + '"""'
    else:
        doc = "pass"
    return doc


renderer = CodeDoc()


def get_python_code(name):
    """Return the Python code for a pybind11 object.

    Arguments:
        name (str): object name.

    Returns:
        str: Python code as a string.
    """
    return render_doc(name, title="# %s", renderer=renderer)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--python",
        action="store_const",
        dest="fmt",
        const="python",
        default="python",
        help="use python format",
    )
    parser.add_argument(
        "--text", action="store_const", dest="fmt", const="text", help="use text format"
    )
    parser.add_argument("object", nargs="+", help="Name of the object to extract")
    args = parser.parse_args(sys.argv[1:])

    if args.fmt == "python":
        func = get_python_code
    else:

        def func(name):
            return render_doc(name, title="# %s", renderer=plaintext)

    for obj in args.object:
        print(func(obj))
