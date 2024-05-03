# must install dlltracer using p

import dlltracer
import sys

from config import TMP_DIR

sys.path.insert(0, TMP_DIR.as_posix())
with dlltracer.Trace(out=sys.stdout):
    import code_aster
    from code_aster import CA
