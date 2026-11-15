# ITensorNHDMRG.jl, a copy

This folder is a copy of the `src/` folder of ITensorNHDMRG.jl
(<https://github.com/tipfom/ITensorNHDMRG.jl>), commit
`c7f107c55cc2a19939e396cb5ab2596bdfe2b26b` of 2025-12-29, under its MIT
license (`LICENSE`). The files are not changed.

We keep a copy because the package declares `ITensorMPS = "0.3.19"` in its
compat, so as a dependency it blocks ITensorMPS 0.4. The code itself runs
unchanged on ITensorMPS 0.4, since 0.4 only removes `TimeDependentSum` and
`to_vec`, which the package does not use.

`juliasession.py` includes `ITensorNHDMRG.jl` as a module of `Main` before
any other Julia file, and `get_gs.jl` and `nhdmrg.jl` load it with
`using .ITensorNHDMRG`. Its one dependency outside ITensors, ITensorMPS,
KrylovKit and the standard library is VectorInterface, which is listed in
`../../juliapkg.json`. To update the copy, replace the `.jl` files by the
`src/` of a newer commit and change the hash above.
