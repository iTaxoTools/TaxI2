from __future__ import annotations

from typing import Callable, Generic, Iterable, Iterator, NamedTuple, TypeVar

from itaxotools.common.types import Type, TypeMeta  # noqa

Item = TypeVar("Item", bound=NamedTuple)


class Container(Generic[Item]):
    """Container that can be iterated multiple times"""

    def __init__(
        self,
        source: Iterable[Item] | Callable[..., iter[Item]],
        *args,
        **kwargs,
    ):
        """The `source` argument is either an iterable or a callable"""
        self.iterable = None
        self.callable = None
        self.args = []
        self.kwargs = {}
        if callable(source):
            self.callable = source
            self.args = args
            self.kwargs = kwargs
        else:  # iterable
            self.iterable = source
            if args or kwargs:
                raise TypeError("Cannot pass arguments to iterable source")

    def __iter__(self) -> Iterator[Item]:
        if self.callable:
            return self.callable(*self.args, **self.kwargs)
        return iter(self.iterable)

    def __len__(self):
        return sum(1 for _ in self)


class Percentage(float):
    def __str__(self):
        return f"{100*self:.2f}%"
