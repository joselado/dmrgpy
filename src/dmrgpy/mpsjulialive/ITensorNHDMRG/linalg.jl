function biorthoblocktransform(
    M::ITensor, leftinds, rightinds; checknormal=false, kwargs...
)
    # Linds, Rinds may not have the correct directions
    Lis = ITensors.indices(leftinds)
    Ris = ITensors.indices(rightinds)

    # Ensure the indices have the correct directions,
    # QNs, etc.
    # First grab the indices in A, then permute them
    # correctly.
    Lis = permute(commoninds(M, Lis), Lis)
    Ris = permute(commoninds(M, Ris), Ris)

    for (l, r) in zip(Lis, Ris)
        if space(l) != space(r)
            error("In transform, indices must come in pairs with equal spaces.")
        end
        if hasqns(M)
            if dir(l) == dir(r)
                error("In transform, indices must come in pairs with opposite directions")
            end
        end
    end

    CL = combiner(Lis...)
    CR = combiner(dag(Ris)...)

    cL = combinedind(CL)
    cR = dag(combinedind(CR))

    M = M * CL * dag(CR)
    if inds(M) != (cL, cR)
        M = permute(M, cL, cR)
    end
    Mtensor = NDTensors.Tensor(M)

    if checknormal && isapprox(Mtensor * dag(Mtensor), dag(Mtensor) * Mtensor)
        # This is Technique 1 in the paper
        @info "using normality condition"
        K = (Mtensor + dag(Mtensor)) / 2
        D, U, spec = eigen(K; kwargs...)

        return D, U, dag(U), spec
    else
        Bi, Yi, Ybari, spec = biorthoblocktransform(Mtensor; kwargs...)

        Yi = Yi * dag(CL)
        Ybari = Ybari * CR

        return Bi, Yi, Ybari, spec
    end
end

function biorthoblocktransform(
    M::ITensors.Tensor{ElT,2,<:ITensors.Dense}; kwargs...
) where {ElT<:Union{Real,Complex}}
    # transforms the matrix M according to the procedure outlined in App. C in 2401.15000
    lB, lY, lYbar, spec = biorthoblocktransform(matrix(M); kwargs...)
    B = first(lB)
    Y = first(lY)
    Ybar = first(lYbar)

    link = Index(size(B, 1))

    Yi = itensor(Y, inds(M)[1], link)
    Ybari = itensor(Ybar, inds(M)[2], link')
    Bi = itensor(B, link, link')

    return Bi, Yi, Ybari, spec
end

function biorthoblocktransform(
    M::ITensors.Tensor{ElT,2,<:ITensors.BlockSparse}; kwargs...
) where {ElT<:Union{Real,Complex}}
    # transforms the matrix M according to the procedure outlined in App. C in 2401.15000

    for b in eachnzblock(M)
        all(==(b[1]), b) ||
            error("The transformation currently only supports block diagonal matrices.")
    end

    Ms = [collect(M[b]) for b in eachnzblock(M)]

    # transform all the blocks
    lB, lY, lYbar, spec = biorthoblocktransform(Ms...; kwargs...)

    # find all the blocks for which the truncation reduced the size to zero
    dropblocks = Int64[]
    for i in 1:nnzblocks(M)
        if size(lB[i], 1) == 0
            # the block has shrunk to zero 
            push!(dropblocks, i)
        end
    end

    # Get the list of blocks of M that are not dropped
    nzblocksM = nzblocks(M)
    deleteat!(nzblocksM, dropblocks)
    deleteat!(lB, dropblocks)
    deleteat!(lY, dropblocks)
    deleteat!(lYbar, dropblocks)

    # The number of blocks of T remaining
    nnzblocksM = nnzblocks(M) - length(dropblocks)

    # reconstruct ITensor
    i1, i2 = inds(M)
    l = sim(i1)

    lkeepblocks = Int[b[1] for b in nzblocksM]
    ldropblocks = setdiff(1:nblocks(l), lkeepblocks)
    deleteat!(l, ldropblocks)

    # l may have too many blocks
    (nblocks(l) > nnzblocksM) && error("New index l in transform has too many blocks")

    # Truncation may have changed some block sizes
    for n in 1:nnzblocksM
        ITensors.setblockdim!(l, minimum(size(lB[n])), n)
    end

    indsB = (l', dag(l))
    indsY = (i1, dag(l)')
    indsYbar = (i2, l)

    nzblocksB = Vector{ITensors.Block{2}}(undef, nnzblocksM)
    nzblocksY = Vector{ITensors.Block{2}}(undef, nnzblocksM)
    nzblocksYbar = Vector{ITensors.Block{2}}(undef, nnzblocksM)
    for n in 1:nnzblocksM
        blockM = nzblocksM[n]

        blockB = (n, n)
        nzblocksB[n] = blockB

        blockY = (blockM[1], n)
        nzblocksY[n] = blockY
        nzblocksYbar[n] = blockY
    end

    ElB = eltype(first(lB))
    ElY = eltype(first(lY))

    B = ITensors.BlockSparseTensor(
        ITensors.set_eltype(ITensors.unwrap_array_type(M), ElB), undef, nzblocksB, indsB
    )
    Y = ITensors.BlockSparseTensor(
        ITensors.set_eltype(ITensors.unwrap_array_type(M), ElY), undef, nzblocksY, indsY
    )
    Ybar = ITensors.BlockSparseTensor(
        ITensors.set_eltype(ITensors.unwrap_array_type(M), ElY),
        undef,
        nzblocksYbar,
        indsYbar,
    )

    # copy the data 
    for n in 1:nnzblocksM
        blockB = nzblocksB[n]
        copyto!(ITensors.blockview(B, blockB), lB[n])

        blockY = nzblocksY[n]
        copyto!(ITensors.blockview(Y, blockY), lY[n])
        copyto!(ITensors.blockview(Ybar, blockY), lYbar[n])
    end

    return itensor(B), itensor(Y), itensor(Ybar), spec
end

function biorthoblocktransform(
    Ms::Matrix{ElT}...;
    maxdim,
    mindim,
    cutoff,
    biorthonormalize=true,
    unitarize=true,
    verbosity=1,
    use_absolute_cutoff=false,
    use_relative_cutoff=true,
) where {ElT<:Union{Real,Complex}}
    # transforms the matrix M according to the procedure outlined in App. C in 2401.15000
    cumdims = cumsum([size(M, 1) for M in Ms])

    Fs = Schur[]
    sizehint!(Fs, length(Ms))

    vals = zeros(complex(ElT), cumdims[end])
    for (i, M) in enumerate(Ms)
        F = schur(M)
        push!(Fs, F)
        # We have that M = F.vectors * F.Schur * F.vectors'

        vals[(i == firstindex(Ms) ? 1 : cumdims[i - 1] + 1):cumdims[i]] = F.values
    end

    # permute the Schur decomposition such that (0 ... keep) are in the beginning and the discarded 
    # part of M is at (keep+1 ... end)
    valsperm = sortperm(vals; by=abs2, rev=true)
    select = zeros(Bool, size(vals))
    minvalkept = min(firstindex(valsperm) + mindim, lastindex(valsperm))
    maxvalkept = min(firstindex(valsperm) + maxdim - 1, lastindex(valsperm))
   
    @inbounds for i in firstindex(valsperm):minvalkept
        select[valsperm[i]] = true
    end
    
    if use_absolute_cutoff
        relative_cutoff = use_relative_cutoff ? abs(vals[1]) * cutoff : cutoff
        @inbounds for i in minvalkept+1:maxvalkept
            abs(vals[valsperm[i]]) <= relative_cutoff && break
            select[valsperm[i]] = true
        end
    else
        relative_cutoff = use_relative_cutoff ? sum(abs, vals) * cutoff : cutoff
        itrunc = maxvalkept
        truncatedweight = sum(i -> abs(vals[i]), valsperm[itrunc:end])
        @inbounds while itrunc > minvalkept && truncatedweight < relative_cutoff
            itrunc -= 1
            truncatedweight += abs(vals[valsperm[itrunc]]) 
        end
        @inbounds for i in minvalkept+1:itrunc
            select[valsperm[i]] = true
        end
    end

    @assert sum(select) <= maxdim

    lB = Matrix[]
    lY = Matrix[]
    lYbar = Matrix[]

    eigvalskept = eltype(vals)[]
    truncerr = 0.0

    for (i, F) in enumerate(Fs)
        selecti = select[(i == firstindex(Ms) ? 1 : cumdims[i - 1] + 1):cumdims[i]]
        F = ordschur!(F, selecti)
        keepi = sum(selecti)

        if keepi == 0
            push!(lB, Matrix{eltype(F.Schur)}(undef, 0, 0))
            push!(lY, Matrix{eltype(F.vectors)}(undef, 0, 0))
            push!(lYbar, Matrix{eltype(F.vectors)}(undef, 0, 0))
            continue
        end

        if eltype(F) <: Real &&
            keepi < length(F.values) &&
            !iszero(F.Schur[keepi + 1, keepi])
            @info "Increasing keep in block $i because of 2x2 block in the Schur decomposition"
            keepi += 1
        end
        append!(eigvalskept, F.values[1:keepi])
        if keepi >= length(F.values)
            push!(lB, F.Schur)
            push!(lY, F.vectors)
            push!(lYbar, F.vectors)
            continue
        end

        truncerr += sum(abs, @views F.values[(keepi + 1):end])

        @views A = F.Schur[1:keepi, 1:keepi]
        @views B = F.Schur[(keepi + 1):end, 1:keepi]
        @views C = F.Schur[(keepi + 1):end, (keepi + 1):end]
        @views D = F.Schur[1:keepi, (keepi + 1):end]

        @assert all(iszero, B) "matrix B contains non-zero elements $B"

        Y = copy(F.vectors)
        @views Ys = Y[:, 1:keepi]
        @views Yd = Y[:, (keepi + 1):end]
        Ybars = copy(Ys) 
        Ybard = copy(Yd)

        # solve the Sylvester equation AX - XC = D -> A X + X (-C) + (-D) = 0
        # X = sylvester(A, -C, -D)
        try
            X = LAPACK.trsyl!('N', 'N', A, C, D, -1)
        catch e
            if e isa LAPACKException && e.info == 1
                if verbosity > 0
                    @warn "Lapack call to trsyl! had to alter the eigenvalues indicating almost degenerate eigenvalues in the matrix"
                end
            else
                throw(e)
            end
        end
        X = D

        mul!(Ybars, Ybard, Adjoint(X), 1.0, 1.0)
        mul!(Yd, Ys, D, -1.0, 1.0)

        push!(lB, A)
        push!(lY, Ys)
        push!(lYbar, Ybars)
    end

    for i in eachindex(lB)
        size(lB[i], 1) == 0 && continue
        # Gram-Schmidt algorithm on the columns of Y and Ybar
        # This is Technique 3 in the paper

        Y = lY[i]
        Ybar = lYbar[i]

        for k in axes(Y, 2)
            for j in firstindex(Y, 1):(k - 1)
                # this should always be one
                n = adjoint(Ybar[:, j]) * Y[:, j]
                # @assert isapprox(n, one(n); rtol=1e-4) "norm $n deviates significantly from one"

                Y[:, k] -= ((adjoint(Y[:, k]) * Ybar[:, j]) / n) * Y[:, j]
                Ybar[:, k] -= ((adjoint(Ybar[:, k]) * Y[:, j]) / n) * Ybar[:, j]
            end

            # normalize the remainder 
            n = sqrt(adjoint(Ybar[:, k]) * Y[:, k])
            Y[:, k] ./= n
            Ybar[:, k] ./= n
        end

        @assert adjoint(Ybar) * Y ≈ I
    end


    if unitarize
        for i in eachindex(lB)
            size(lB[i], 1) == 0 && continue
            # Unitarize both Y and Ybar
            # This is Technique 4 in the paper
            F = svd(lY[i])
            lY[i] = F.U * F.Vt
            lYbar[i] = conj(lY[i])
        end
    end

    return lB, lY, lYbar, Spectrum(eigvalskept, truncerr)
end

function invertitensor(A::ITensor, finv, Linds...)
    Lis = commoninds(A, ITensors.indices(Linds...))
    Ris = uniqueinds(A, Lis)

    Cr = combiner(Ris...)
    Cl = combiner(Lis...)

    A = A * Cr * Cl

    Minv = nothing

    return if hasqns(A)
        Minv = deepcopy(A)
        for b in eachnzblock(A)
            Amat = matrix(A[b])
            Amatinv = finv(Amat)

            Minv[b] .= adjoint(Amatinv)
        end
        Minv * dag(Cr) * dag(Cl)
    else
        M = matrix(A, combinedind(Cr), combinedind(Cl))
        Minv = itensor(adjoint(finv(M)), combinedind(Cr), combinedind(Cl); tol=1e-8)
        Minv * dag(Cr) * dag(Cl)
    end
end

function LinearAlgebra.inv(A::ITensor, Linds...)
    return invertitensor(A, inv, Linds...)
end

function LinearAlgebra.pinv(A::ITensor, Linds...; kwargs...)
    return invertitensor(A, x -> pinv(x; kwargs...), Linds...)
end

function LinearAlgebra.inv(A::ITensor)
    Ris = filterinds(A; plev=0)
    Lis = Ris'
    return LinearAlgebra.inv(A, Lis)
end

function LinearAlgebra.pinv(A::ITensor; kwargs...)
    Ris = filterinds(A; plev=0)
    Lis = Ris'
    return LinearAlgebra.pinv(A, Lis; kwargs...)
end