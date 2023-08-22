"""
Generate fake libaster.
"""

import os.path as osp
import sys
from subprocess import PIPE, CalledProcessError, run
from pathlib import Path

sys.path.insert(0, osp.dirname(__file__))

from pydoc_codedoc import get_python_code


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

    lines = [i.rstrip() for i in "\n".join(blocks).splitlines()]

    path = Path(filename)
    path.parent.mkdir(exist_ok=True)
    with path.open("w") as flib:
        flib.write("\n".join(lines))
    print("output written to:", filename)
    try:
        cmd = ["black", filename]
        run(cmd, check=True, stderr=PIPE)
        print("black formatter passed with success")
    except CalledProcessError:
        print("black formatter failed!")


if __name__ == "__main__":
    assert len(sys.argv[1:]) == 1, "exactly one argument is expected!"
    build_pylibaster(sys.argv[1])
