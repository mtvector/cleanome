import importlib.abc
import os
import sys


class _BlockPyArrowFinder(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "pyarrow" or fullname.startswith("pyarrow."):
            raise ImportError(
                "cleanome blocked pyarrow import to avoid broken binary compatibility"
            )
        return None


def block_broken_pyarrow():
    if os.environ.get("CLEANOME_ALLOW_PYARROW"):
        return
    for finder in sys.meta_path:
        if isinstance(finder, _BlockPyArrowFinder):
            return
    sys.meta_path.insert(0, _BlockPyArrowFinder())
