"""Some metaprogramming helpers used in the main `orthoani` module.
"""

import contextlib
import tempfile
from typing import Any, ContextManager, Generic, Iterator, Sequence, TypeVar

from pathlib import Path
from typing import Iterator


_T = TypeVar("_T")
_S = TypeVar("_S", bound=Sequence[Any])


class ExitStack(contextlib.ExitStack):
    """A `contextlib.ExitStack` aliasing `enter_context` to the ``<<`` operator.

    Inspired by the `contexter <https://pypi.org/project/contexter/>`_ library,
    which is not developed anymore, since most of its features are now covered
    by the built-in `contexlib` module.
    """

    def __lshift__(self, r: ContextManager[_T]) -> _T:  # noqa: D105
        return self.enter_context(r)


class BlockIterator(Generic[_S], Iterator[_S]):
    """An iterator that yields even-sized blocks from a sliceable input.

    Example:
        >>> for x in BlockIterator("abcde", 2):
        ...     print(x)
        ab
        cd

    """

    def __init__(self, data: _S, blocksize: int):  # noqa: D105, D107
        self.data = data
        self.blocksize = blocksize
        self.cursor = 0

    def __iter__(self) -> "BlockIterator[_S]":  # noqa: D105
        return self

    def __next__(self) -> _S:  # noqa: D105
        if self.cursor + self.blocksize < len(self.data):
            block = self.data[self.cursor : self.cursor + self.blocksize]
            self.cursor += self.blocksize
            return block  # type: ignore
        raise StopIteration


@contextlib.contextmanager
def temppath(prefix: str = "orthoani", suffix: str = "") -> Iterator[Path]:
    """Get a context manager that returns a path to a temporary file.

    The file is created when the context is entered, and deleted when the
    context is exited.
    """
    try:
        temp = tempfile.NamedTemporaryFile(mode="w", prefix=prefix, suffix=suffix)
        yield Path(temp.name)
    finally:
        temp.close()
