from __future__ import annotations

import pytest

from itaxotools.taxi2.types import Type


def test_type_inheritance() -> None:
    class Parent(Type):
        pass

    class Child_A(Parent):
        pass

    class Child_B(Parent):
        pass

    class GrandChild_A(Child_A):
        pass

    class GrandChild_B(Child_A, Parent):
        pass

    assert Child_A in Parent
    assert Child_B in Parent

    assert GrandChild_A in Child_A
    assert GrandChild_A not in Parent

    assert GrandChild_B in Child_A
    assert GrandChild_B in Parent

    assert Child_A() not in Parent
    with pytest.raises(TypeError):
        assert Child_A() not in Parent()
    with pytest.raises(TypeError):
        assert Child_A not in Parent()
