using KrylovKit: bieigsolve, eigsolve, BiArnoldi, linsolve
using LinearAlgebra, KrylovKit
using VectorInterface: VectorInterface

function nhproblemsolver!(
    ::Algorithm"twosided",
    PH,
    Θl,
    Θr;
    eigsolve_tol,
    eigsolve_krylovdim,
    eigsolve_maxiter,
    eigsolve_verbosity,
    eigsolve_which_eigenvalue,
)
    fA = x -> product(PH, x)
    fAH = x -> adjointproduct(PH, x)
    f = (fA, fAH)
    
    vals, (V, W), info = bieigsolve(
        f,
        Θr,
        Θl,
        1,
        eigsolve_which_eigenvalue,
        BiArnoldi(;
            tol=eigsolve_tol,
            krylovdim=eigsolve_krylovdim,
            maxiter=eigsolve_maxiter,
            verbosity=eigsolve_verbosity,
        ),
    )
    
    return first(vals), first(W), first(V)
end

function nhproblemsolver!(
    ::Algorithm"onesided",
    PH,
    Θl,
    Θr;
    eigsolve_tol,
    eigsolve_krylovdim,
    eigsolve_maxiter,
    eigsolve_verbosity,
    eigsolve_which_eigenvalue,
)
    fA = x -> product(PH, x)
    fAH = x -> adjointproduct(PH, x)

    valsH, vecsH = eigsolve(
        fAH,
        Θl,
        1,
        eigsolve_which_eigenvalue;
        ishermitian=false,
        tol=eigsolve_tol,
        krylovdim=eigsolve_krylovdim,
        maxiter=eigsolve_maxiter,
        verbosity=eigsolve_verbosity,
    )

    vals, vecs = eigsolve(
        fA,
        Θr,
        1,
        eigsolve_which_eigenvalue;
        ishermitian=false,
        tol=eigsolve_tol,
        krylovdim=eigsolve_krylovdim,
        maxiter=eigsolve_maxiter,
        verbosity=eigsolve_verbosity,
    )

    return first(vals), noprime(first(vecsH)), noprime(first(vecs))
end

function selecteigenvalue(f::Function, x)
    return f(x)
end

function selecteigenvalue(f::Number, x)
    return f
end

function selecteigenvalue(f, x)
    return x
end

function iteraterayleighquotient(
    PH,
    left,
    right;
    eigsolve_tol,
    eigsolve_krylovdim,
    eigsolve_maxiter,
    eigsolve_verbosity,
    eigsolve_which_eigenvalue,
)
    ρ = inner(left, product(PH, right)) / inner(left, right)

    λt = selecteigenvalue(eigsolve_which_eigenvalue, ρ)

    nr, info = linsolve(
        x -> product(PH, x),
        right,
        -λt;
        tol=eigsolve_tol,
        krylovdim=eigsolve_krylovdim,
        maxiter=eigsolve_maxiter,
        verbosity=eigsolve_verbosity,
    )
    nl, info = linsolve(
        x -> adjointproduct(PH, x),
        left,
        -conj(λt);
        tol=eigsolve_tol,
        krylovdim=eigsolve_krylovdim,
        maxiter=eigsolve_maxiter,
        verbosity=eigsolve_verbosity,
    )

    return nl, nr
end

function nhproblemsolver!(
    ::Algorithm"twosidedinverse",
    PH,
    Θl,
    Θr;
    eigsolve_tol,
    eigsolve_krylovdim,
    eigsolve_maxiter,
    eigsolve_verbosity,
    eigsolve_which_eigenvalue,
    eigsolve_refinements=2,
)
    # this function uses the two-sided inverse iteration and two-sided rayleigh quotient 
    # to anchor to the closest eigenvalue, see: https://doi.org/10.1002/nla.1945

    # by supplying a function to eigsolve_which_eigenvalue, the selected eigenvalue may be altered

    for _ in 1:eigsolve_refinements
        Θl, Θr = iteraterayleighquotient(
            PH,
            Θl,
            Θr;
            eigsolve_tol,
            eigsolve_krylovdim,
            eigsolve_maxiter,
            eigsolve_verbosity,
            eigsolve_which_eigenvalue,
        )

        normalize!(Θl)
        normalize!(Θr)
        ov = complex(inner(Θl, Θr))
        Θl /= conj(sqrt(ov))
        Θr /= sqrt(ov)
    end

    fx = inner(Θl, product(PH, Θr)) / inner(Θl, Θr)

    return fx, Θl, Θr
end

function selecttau(s::Symbol)
    if s == :SR 
        return -1
    elseif s == :LR
        return +1 
    elseif s == :SI
        return -1im 
    elseif s == :LI 
        return +1im 
    end
    return -1
end
