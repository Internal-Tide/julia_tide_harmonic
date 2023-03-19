"""
Robust MLR via iteratively reweighted least squares.

"""

using LinearAlgebra
using Distributions
using Einsum
# Weighting functions:
function andrews(r)
    r = abs.(r)
    r = max(sqrt(eps()), r)
    w = (r .< pi) .* sin.(r) ./ r
    return w
end


function bisquare(r)
    r = abs.(r)
    w = (r .< 1) .* (1 - r .^ 2) .^ 2
    return w
end


function cauchy(r)
    r = abs.(r)
    w = 1 ./ (1 .+ r .^ 2)
    return w
end


function fair(r)
    w = 1 ./ (1 .+ abs.(r))
    return w
end


function huber(r)
    w = 1 ./ max.(1, abs.(r))
    return w
end


function logistic(r)
    r = abs.(r)
    r = max(sqrt(1.0), r)
    w = tanh.(r) ./ r
    return w
end


function ols(r)
    w = ones(length(r))
    return w
end


function talwar(r)
    w = (abs.(r) .< 1) .* 1.0
    return w
end


function welsch(r)
    r = abs.(r)
    w = exp.(-(r .^ 2))
    return w
end

wfuncdict = Dict(
    "andrews" => andrews,
    "bisquare" => bisquare,
    "cauchy" => cauchy,
    "fair" => fair,
    "huber" => huber,
    "logistic" => logistic,
    "ols" => ols,
    "talwar" => talwar,
    "welsch" => welsch,
)

tune_defaults = Dict(
    "andrews" => 1.339,
    "bisquare" => 4.685,
    "cauchy" => 2.385,
    "fair" => 1.400,
    "huber" => 1.345,
    "logistic" => 1.205,
    "ols" => 1,
    "talwar" => 2.795,
    "welsch" => 2.985,
)


function sigma_hat(x)
    """
    Robust estimate of standard deviation based on medians.
    """
    # The center could be based on the mean or some other function.
    return median(abs.(x .- median(x))) / 0.6745
end


function leverage(x)
    """
    Calculate leverage as the diagonal of the "Hat" matrix of the
    model matrix, x.
    """

    # The Hat is x times its pseudo-inverse.
    # In einum, the diagonal is calculated for each row of x
    # and column of pinv as the dot product of column j of x.T
    # and column j of pinv; hence the 'j' in the output means
    # *don't* sum over j.

    hdiag = Einsum._einsum("ij, ij -> j", x', pinv(x))
    # This should be real and positive, but with floating point
    # arithmetic the imaginary part is not exactly zero.
    return abs.(hdiag)
end


function r_normed(R, rfac)
    """
    Normalized residuals from raw residuals and a multiplicative factor.
    """
    return rfac .* R ./ sigma_hat(R)
end

#
function lstsq(A::AbstractMatrix{T}, B::AbstractMatrix{T}; rcond::Real=eps(T)) where {T<:AbstractFloat}
    m, n = size(A)
    k = min(m, n)
    if m < n
        A = A'
        m, n = n, m
    end
    if !isequal(T, eltype(B))
        B = convert(Matrix{T}, B)
    end
    if size(B, 1) != m
        throw(DimensionMismatch("B has wrong dimensions"))
    end
    if rcond < 0
        throw(ArgumentError("rcond must be non-negative"))
    end
    if rcond == 0
        rcond = eps(T)
    end
    if m == n
        # square case
        if issymmetric(A)
            # symmetric case
            F = cholesky(A, rcond=rcond)
            X = F \ (F' \ B)
            R = B - A * X
            S = nothing
            rank = F.rank
            svdvals = nothing
        else
            # nonsymmetric case
            F = qr(A, Val(true), rcond=rcond)
            X = F \ B
            R = B - A * X
            S = nothing
            rank = F.rank
            svdvals = nothing
        end
    else
        # rectangular case
        F = qr(A, Val(true), rcond=rcond)
        X = F \ B
        R = B - A * X
        S = svd(A, Val(true), rcond=rcond).S
        rank = F.rank
        svdvals = nothing
    end
    return X, R, S, rank, svdvals
end
#


function robustfit(X, y; weight_function="bisquare", tune=nothing, rcond=1, tol=0.001, maxit=50)
    if tune === nothing
        tune = tune_defaults[weight_function]
    end
    _wfunc = wfuncdict[weight_function]
    if ndims(X) == 1
        X = reshape(X, (size(X), 1))
    end
    n, p = size(X)
    lev = leverage(X)
    out = Dict(
        "weight_function" => weight_function,
        "tune" => tune,
        "rcond" => rcond,
        "tol" => tol,
        "maxit" => maxit,
        "leverage" => lev
    )
    rfac = 1 / (tune * sqrt(1 - lev))
    oldrsumsq = nothing
    oldrmeansq = nothing
    oldlstsq = nothing
    oldw = nothing
    iterations = 0
    w = ones(size(y))
    for i = 1:maxit
        wX = w .* reshape(X, (n, p))
        wy = w .* y
        b, rsumsq, rank, sing = linalg.lstsq(wX, wy, rcond)
        rsumsq = rsumsq[1]
        if i == 1
            rms_resid = sqrt(rsumsq / n)
            out["ols_b"] = b
            out["ols_rms_resid"] = rms_resid
        end
        rmeansq = rsumsq / sum(w)
        if oldrsumsq != nothing
            improvement = (oldrmeansq - rmeansq) / oldrmeansq
            if improvement < 0
                b, rsumsq, rank, sing = oldlstsq
                w = oldw
                iterations = i
                break
            end
            if improvement < tol
                iterations = i + 1
                break
            end
        end
        oldlstsq = (b, rsumsq, rank, sing)
        oldw = w
        oldrsumsq = rsumsq
        oldrmeansq = rmeansq
        resid = y - X * b
        w = _wfunc(r_normed(resid, rfac))
        if iterations == 0
            iterations = maxit
            rms_resid = sqrt(mean(abs.(resid) .^ 2))
            out["iterations"] = iterations
            out["b"] = b
            out["s"] = sing
            out["w"] = w
            out["rank"] = rank
            out["rms_resid"] = rms_resid
        end
    end
    return out
end

using StatsBase








###########
using LinearAlgebra

function robustfit(X, y; weight_function="bisquare", tune=nothing, rcond=1, tol=0.001, maxit=50)
    n, p = size(X)
    X = [ones(n) X]
    w = zeros(n)
    iter = 0
    if tune === nothing
        if weight_function == "bisquare"
            tune = 4.685
        else
            tune = 6.0
        end
    end
    while iter < maxit
        iter += 1
        Xw = X .* sqrt.(w)
        yw = y .* sqrt.(w)
        b = ldiv!(Xw, yw)
        res = y - X * b
        s = median(abs.(res))
        if s == 0
            break
        end
        u = res / s
        if weight_function == "bisquare"
            w = ((1 - u .^ 2) .^ 2) .* ((abs.(u) .<= 1) .* (abs.(u) .> 0))
            w = w .* ((n * sum(w .^ 2)) / sum(w)^2)
            w = min.(w, tune)
        elseif weight_function == "huber"
            w = abs.(u) .<= 1
            w = w .* ((n * sum(w .^ 2)) / sum(w)^2)
            w = min.(w, tune)
        else
            w = 1 ./ (abs.(u) .+ 1)
            w = w .* ((n * sum(w .^ 2)) / sum(w)^2)
        end
        if maximum(abs.(w .- sqrt.(w))) <= tol
            break
        end
    end
    b = ldiv!(X .* sqrt.(w), y .* sqrt.(w))
    y_hat = X * b
    residuals = y .- y_hat
    df_resid = n - p
    tvalues = b ./ sqrt.(diag(inv(X' * (diagm(w) * X)) * rcond))
    pvalues = 2 .* (1 .- cdf.(tdist(df_resid), abs.(tvalues)))
    return b, residuals, w, tvalues, pvalues
end




using LinearAlgebra


n = 10000
x = collect(0:n-1)
x0 = ones(n)
x1 = exp.(1im .* x ./ 9)
x2 = exp.(1im .* x ./ 7)
y = (1 + 1im) .* x1 .+ (2 - 1im) .* x2 .+ (0.1 .* randn(n) .+ 0.1im .* randn(n))
y[1:10:end] .= randn(n)[1:10:end] .+ 1im .* randn(n)[1:10:end]

y[11] = 3
y[21] = 2im
y[31] = -2 - 3im

A = hcat(x0, x1, x2)
c = A \ y
println("OLS: ", c)

rf1 = robustfit(A, y)

println("robust: ", rf1.b)
println("another test: a very short real series")

x = collect(1:20)
x0 = ones(x)
xx = hcat(x0, x)

# Signal for the model: linear trend.
y = 2 .* x

# Some outliers.
y[1] = 1.5
y[3] = -2
y[5] = 9.6

# Use a sine as the "noise" component; not part of the model.
y = y .+ 0.1 .* sin.(x)

rf2 = robustfit(xx, y)
println(xx \ y)
println(rf2.b)

using GLM
using StatsModels


# Generate some data
n = 100
p = 2
X = randn(n, p)
y = 2 * X[:, 1] - 3 * X[:, 2] + 0.5 * randn(n)

# Create a weight vector
weights = ones(n)

# Fit a linear model using IRLS
gm1 = fit(
    @formula(y ~ X[:, 1] + X[:, 2]),
    (X=X, y=y),
    Weights(weights),
    Normal(),
    LogLink()
)

