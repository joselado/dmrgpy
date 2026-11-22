---
paths:
  - "src/dmrgpy/multioperator.py"
  - "src/dmrgpy/multioperatortk/**"
  - "tests/test_multioperator*.py"
  - "tests/test_audit_2026_09_24c_hermiticity.py"
  - "tests/test_audit_2026_09_24_operators.py"
  - "examples/algebra/operator_canonical_form/**"
---

# The canonical form of a `MultiOperator` and the Hermiticity proof

This is the paragraph CLAUDE.md used to carry under "Operator
representation: `MultiOperator`"; it loads when `multioperator.py` or a
file under `multioperatortk/` is read or edited.

**Canonical form and the Hermiticity proof.**
`multioperatortk/canonical.py` rewrites every term into one spelling
(factors sorted by site, stable so same-site factors keep their order, the
coefficient carrying the sign of every fermionic exchange the sort makes,
identity factors dropped) so that terms which are the same operator have
the same signature and collect in a dict.
`MultiOperator.simplify()`/`is_zero()`/`is_hermitian()`/`is_antihermitian()`
all read off it. It replaced a sympy round trip (`sympymultioperator.py`,
deleted with it, taking sympy out of the dependencies) that collected only
terms already spelled identically: 0.834 s on an n=20 spin-1/2 chain to
report `False` for an ordinary Heisenberg Hamiltonian, now proven Hermitian
in 0.5 ms. The rewrite is exact, so the proof is **one-sided**: an empty
canonical form of `H-H^dagger` means Hermitian, full stop, while a nonempty
one means not proven, since a same-site identity (`Sx Sx = 1/4`) or a name
alias (`Sp = Sx + i Sy`, `Sz = (Nup-Ndn)/2`) needs the site type, which
lives on the chain and not on the operator. Two factors sharing a site are
reordered only when both are diagonal (`_DIAGONAL`), which is what makes the
Hubbard `Nup[i]*Ndn[i]` provable while leaving `Sx[i]*Sz[i]` alone; do NOT
widen that to a general same-site sort, those two are both even and do not
commute. Nor is a term reordered when it names both a pre-transform
`C`-type name and a bare `A`-type one (`_LOCAL_LADDER`): with
`C_j = F_0...F_{j-1} A_j`, `A_i` anticommutes with `C_j` for j>i and
commutes for j<i, so their exchange sign depends on the site order and is
not a parity (2026-09-24 audit finding 1). Do NOT fix that by giving the
`A`-type names no parity: they are also the boson ladder operators, and the
boson and Jordan-Wigner-transformed proofs need them even. So
`Many_Body_Chain.is_hermitian()` takes the proof when it lands and falls
back to its random-witness probe otherwise; do not turn that fallback into
a rejection. Only the names in `canonical.py`'s `_PARITY` table are
reordered; anything else (parafermionic `Sig`/`Tau`, which reorder with a
Z_n phase, or a caller's own name) is left spelled as written, and the
Hermiticity proof refuses outright when any name is off the table, since
`get_dagger()` leaves an unrecognized name untouched (and carries a phase
for a name whose adjoint is a multiple of itself, `_dagger_phase`, so far
only `ISy` -> `-ISy`) and such an operator would otherwise cancel against
its own "dagger". See `docs/development.md` 4.2a,
`tests/test_multioperator_canonical.py` and
`examples/algebra/operator_canonical_form`.
