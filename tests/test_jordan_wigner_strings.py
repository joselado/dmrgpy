"""The Jordan-Wigner strings of multioperatortk/jordanwigner.py are built
in one step; they used to be built by one MultiOperator product per F
factor, which is quadratic in the span of a term. The old recipe is kept
here as the reference, and the new terms must be the same, coefficient
(its float value and sign) and factor order included."""

import pytest

from dmrgpy import multioperator
from dmrgpy.multioperatortk import jordanwigner as jw


def _o(a):
    return multioperator.obj2MO([a])


def _old_CdagC(i, j):
    if i == j:
        return _o(["Adag", i])*_o(["A", j])
    elif i < j:
        m = _o(["Adag", i])
        for k in range(i, j-1):
            m = m*_o(["F", k+1])
        return m*_o(["A", j])
    return -1*_old_CCdag(j, i)


def _old_CCdag(i, j):
    if i == j:
        return _o(["A", i])*_o(["Adag", j])
    elif i < j:
        m = -1*_o(["A", i])
        for k in range(i, j-1):
            m = m*_o(["F", k+1])
        return m*_o(["Adag", j])
    return -_old_CdagC(j, i)


def _old_C(i):
    m = 1.0
    for k in range(i):
        m = _o(["F", k])*m
    return m*_o(["A", i])


def _old_Cdag(i):
    m = _o(["Adag", i])
    for k in range(i):
        m = _o(["F", k])*m
    return m


def _old_CC(i, j):
    if i == j:
        return _o(["A", i])*_o(["A", i])
    elif i < j:
        m = -1*_o(["A", i])
        for k in range(i, j-1):
            m = m*_o(["F", k+1])
        return m*_o(["A", j])
    return -1*_old_CC(j, i)


def _same(a, b):
    assert len(a.op) == len(b.op) == 1
    ta, tb = a.op[0], b.op[0]
    assert ta[0] == tb[0] and type(ta[0]) is type(tb[0])
    assert [list(f) for f in ta[1:]] == [list(f) for f in tb[1:]]


@pytest.mark.parametrize("name", ["CdagC", "CCdag", "CC", "CdagCdag"])
def test_two_point_strings_are_unchanged(name):
    new = getattr(jw, name)
    old = {"CdagC": _old_CdagC, "CCdag": _old_CCdag, "CC": _old_CC,
           "CdagCdag": lambda i, j: _old_CC(j, i).get_dagger()}[name]
    for i in range(9):
        for j in range(9):
            _same(new(i, j), old(i, j))


def test_one_point_strings_are_unchanged():
    for i in range(9):
        _same(jw.C(i), _old_C(i))
        _same(jw.Cdag(i), _old_Cdag(i))
