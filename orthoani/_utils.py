"""Some metaprogramming helpers used in the main `orthoani` module.
"""

import contextlib
import collections.abc
import tempfile

from pathlib import Path
from typing import Iterator


class ExitStack(contextlib.ExitStack):
    """A `contextlib.ExitStack` aliasing `enter_context` to the ``<<`` operator.

    Inspired by the `contexter <https://pypi.org/project/contexter/>`_ library,
    which is not developed anymore, since most of its features are now covered
    by the built-in `contexlib` module.
    """

    def __lshift__(self, r):  # noqa: D105
        return self.enter_context(r)


class BlockIterator(collections.abc.Iterator):
    """An iterator that yields even-sized blocks from a sliceable input.

    Example:
        >>> for x in BlockIterator("abcde", 2):
        ...     print(x)
        ab
        cd

    """

    def __init__(self, data, blocksize):  # noqa: D105, D107
        self.data = data[:]
        self.blocksize = blocksize

    def __iter__(self):  # noqa: D105
        return self

    def __next__(self):  # noqa: D105
        if len(self.data) > self.blocksize:
            block = self.data[: self.blocksize]
            self.data = self.data[self.blocksize :]
            return block
        raise StopIteration


@contextlib.contextmanager
def temppath(suffix: str = "") -> Iterator[Path]:
    """Get a context manager that returns a path to a temporary file.

    The file is created when the context is entered, and deleted when the
    context is exited.
    """
    try:
        temp = tempfile.NamedTemporaryFile(mode="w", suffix=suffix)
        yield Path(temp.name)
    finally:
        temp.close()
