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

    @property
    def type(self):
        return type(self)
