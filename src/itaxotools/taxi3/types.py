from __future__ import annotations

from typing import Callable, Generic, Iterable, NamedTuple, TypeVar


class TypeMeta(type):
    _inheritors = dict()

    def __new__(cls, name, bases, attrs):
        obj = super().__new__(cls, name, bases, attrs)
        cls._inheritors[obj] = dict()
        obj._parent = None
        for base in bases:
            if issubclass(base, Type):
                cls._inheritors[base][name] = obj
                obj._parent = base
        return obj

    def __dir__(self):
        return super().__dir__() + list(self._inheritors[self].keys())

    def __getattr__(self, attr):
        if attr in self._inheritors[self]:
            return self._inheritors[self][attr]
        raise AttributeError(f'Type {repr(self.__name__)} has no subtype {repr(attr)}')

    def __iter__(self):
        return iter(self._inheritors[self].values())


class Type(metaclass=TypeMeta):
    """All subclasses are added as class attributes"""

    def __repr__(self):
        if self._parent:
            return f'<{self._parent.__name__}.{type(self).__name__}>'
        return super().__repr__()

    def __eq__(self, other):
        return type(self) is type(other)

    @property
    def type(self):
        return type(self)


Item = TypeVar('Item', bound=NamedTuple)


class Container(Generic[Item]):
    """Container that can be iterated multiple times"""

    def __init__(
        self, source: Iterable[Item] | Callable[..., iter[Item]],
        *args, **kwargs,
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
                raise TypeError('Cannot pass arguments to iterable source')

    def __iter__(self) -> iter[Item]:
        if self.callable:
            return self.callable(*self.args, **self.kwargs)
        return iter(self.iterable)

    def __len__(self):
        return sum(1 for _ in self)


class Percentage(float):
    def __str__(self):
        return f'{100*self:.2f}%'
