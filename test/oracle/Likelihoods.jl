# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# The independent oracle: pure-Julia log-likelihoods for the hurdle NB and the
# ZINB, written from the definitions in docs/method-conditions/, and samplers
# for the same two models.
#
# Nothing in this file calls R, `pscl` or any of the package's own code. That is
# the whole point of it: it exists to be an *independent* evaluation of the
# likelihood at the parameters R reports, so a mistake in the R call (a sign, a
# missing offset, an off-by-one in the truncation) shows up as a disagreement
# rather than cancelling out.
#
# The oracle uses the natural parameterisation (means and probabilities), the
# package's R bridge uses the model-matrix parameterisation; converting between
# them is the caller's job, in the tests.

module Likelihoods

export nb_logpdf,
    nb_pmf,
    tnb_logpdf,
    tnb_pmf,
    hurdle_loglik,
    hurdle_zero_loglik,
    zinb_loglik,
    loggamma,
    invlogit,
    logit,
    rand_nb,
    rand_tnb,
    rand_hurdle,
    rand_zinb,
    mean_of_nb

using Random: Random, AbstractRNG, randn, rand, randexp

# --- log-gamma ----------------------------------------------------------------
# Lanczos (g = 7, n = 9), written out rather than taken from a special-functions
# package so that the oracle shares no code with the code it checks.
const LANCZOS_G = 7.0
const LANCZOS_C = (
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7,
)

"""
    loggamma(x) -> Float64

Natural logarithm of the gamma function, by the Lanczos approximation, with
reflection for `x < 1/2`. Accurate to a few ulps over the range this oracle
uses (counts and size parameters below a few hundred).
"""
function loggamma(x::Real)
    xf = Float64(x)
    if xf < 0.5
        return log(π / sin(π * xf)) - loggamma(1 - xf)
    end
    z = xf - 1
    acc = LANCZOS_C[1]
    for i = 2:9
        acc += LANCZOS_C[i] / (z + i - 1)
    end
    t = z + LANCZOS_G + 0.5
    return 0.5 * log(2π) + (z + 0.5) * log(t) - t + log(acc)
end

# --- links ---------------------------------------------------------------------
logit(p::Real) = log(p / (1 - p))
invlogit(η::Real) = 1 / (1 + exp(-η))

# --- negative binomial ---------------------------------------------------------
"""
    nb_logpdf(y, μ, θ) -> Float64

Log density of the negative binomial with mean `μ` and size (dispersion) `θ`:

    log Γ(y + θ) - log Γ(θ) - log Γ(y + 1)
        + θ log(θ/(θ + μ)) + y log(μ/(θ + μ))
"""
function nb_logpdf(y::Real, μ::Real, θ::Real)
    return loggamma(y + θ) - loggamma(θ) - loggamma(y + 1) +
           θ * log(θ / (θ + μ)) +
           y * log(μ / (θ + μ))
end

"""
    nb_pmf(y, μ, θ) -> Float64

Density of the negative binomial at `y` (exp of [`nb_logpdf`](@ref)).
"""
nb_pmf(y::Real, μ::Real, θ::Real) = exp(nb_logpdf(y, μ, θ))

"""
    mean_of_nb(μ, θ) -> Float64

The mean of the negative binomial, `μ`. Named rather than hard-wired so a test
can state that it is checking the parameterisation and not the name.
"""
mean_of_nb(μ::Real, θ::Real) = μ

# --- zero-truncated negative binomial ------------------------------------------
"""
    tnb_logpdf(y, μ, θ) -> Float64

Log density of the negative binomial conditioned on `y > 0`:
`log f(y) - log(1 - f(0))`. `y` must be at least 1.
"""
function tnb_logpdf(y::Real, μ::Real, θ::Real)
    y >= 1 || throw(ArgumentError("the zero-truncated NB has support y >= 1"))
    return nb_logpdf(y, μ, θ) - log1p(-nb_pmf(0, μ, θ))
end

"""
    tnb_pmf(y, μ, θ) -> Float64

Density of the zero-truncated negative binomial at `y >= 1`.
"""
tnb_pmf(y::Real, μ::Real, θ::Real) = exp(tnb_logpdf(y, μ, θ))

# --- the two models ------------------------------------------------------------
"""
    hurdle_zero_loglik(y, p) -> Float64

The hurdle zero part's contribution: `log(1 - p)` where nothing was seen,
`log p` where something was, with `p` the probability that the count is
positive.
"""
function hurdle_zero_loglik(y::Real, p::Real)
    return y == 0 ? log(1 - p) : log(p)
end

"""
    hurdle_loglik(y, μ, p, θ) -> Float64

Log-likelihood of the hurdle model:

* `y == 0` contributes `log(1 - p_i)`;
* `y > 0` contributes `log p_i + log f_TNB(y_i; μ_i, θ)`.

`μ` and `p` are vectors of the same length as `y`.
"""
function hurdle_loglik(y::AbstractVector, μ::AbstractVector, p::AbstractVector, θ::Real)
    length(y) == length(μ) == length(p) ||
        throw(ArgumentError("y, μ and p must have the same length"))
    total = 0.0
    for i in eachindex(y)
        if y[i] == 0
            total += log(1 - p[i])
        else
            total += log(p[i]) + tnb_logpdf(y[i], μ[i], θ)
        end
    end
    return total
end

"""
    zinb_loglik(y, μ, π, θ) -> Float64

Log-likelihood of the zero-inflated negative binomial, with `π` the structural
zero weight:

* `y == 0` contributes `log(π_i + (1 - π_i) f_NB(0; μ_i, θ))`;
* `y > 0` contributes `log(1 - π_i) + log f_NB(y_i; μ_i, θ)`.
"""
function zinb_loglik(y::AbstractVector, μ::AbstractVector, π::AbstractVector, θ::Real)
    length(y) == length(μ) == length(π) ||
        throw(ArgumentError("y, μ and π must have the same length"))
    total = 0.0
    for i in eachindex(y)
        if y[i] == 0
            total += log(π[i] + (1 - π[i]) * nb_pmf(0, μ[i], θ))
        else
            total += log(1 - π[i]) + nb_logpdf(y[i], μ[i], θ)
        end
    end
    return total
end

# --- samplers ------------------------------------------------------------------
"""
    rand_gamma(rng, shape, scale) -> Float64

Gamma draw by Marsaglia–Tsang (2000). Exact for all positive shapes.
"""
function rand_gamma(rng::AbstractRNG, shape::Float64, scale::Float64)
    if shape < 1
        return rand_gamma(rng, shape + 1, scale) * rand(rng)^(1 / shape)
    end
    d = shape - 1 / 3
    c = 1 / sqrt(9 * d)
    while true
        x = randn(rng)
        v = (1 + c * x)^3
        v <= 0 && continue
        u = rand(rng)
        if u < 1 - 0.0331 * x^4 || log(u) < 0.5 * x^2 + d * (1 - v + log(v))
            return d * v * scale
        end
    end
end

"""
    rand_poisson(rng, λ) -> Int

Poisson draw by inversion, `p(k)` built up by the recurrence
`p(k+1) = p(k) λ/(k+1)`. Exact; `O(λ)` per draw, which is what the simulations
in this suite can afford and is impossible to get subtly wrong.
"""
function rand_poisson(rng::AbstractRNG, λ::Float64)
    λ <= 0 && return 0
    target = rand(rng) * exp(-λ)
    k = 0
    p = exp(-λ)
    acc = p
    while acc < target
        k += 1
        p *= λ / k
        acc += p
    end
    return k
end

"""
    rand_nb(rng, μ, θ) -> Int

Negative binomial draw, as the Poisson–gamma mixture: `λ ~ Gamma(θ, μ/θ)`, then
`y ~ Poisson(λ)`.
"""
function rand_nb(rng::AbstractRNG, μ::Real, θ::Real)
    λ = rand_gamma(rng, Float64(θ), Float64(μ) / Float64(θ))
    return rand_poisson(rng, λ)
end

"""
    rand_tnb(rng, μ, θ) -> Int

Zero-truncated negative binomial draw: rejection sampling from the negative
binomial, discarding zeros.
"""
function rand_tnb(rng::AbstractRNG, μ::Real, θ::Real)
    while true
        y = rand_nb(rng, μ, θ)
        y >= 1 && return y
    end
end

"""
    rand_hurdle(rng, μ, p, θ) -> Int

Draw one observation from the hurdle model: a positive count with probability
`p`, drawn from the zero-truncated negative binomial, and zero otherwise.
"""
function rand_hurdle(rng::AbstractRNG, μ::Real, p::Real, θ::Real)
    return rand(rng) < p ? rand_tnb(rng, μ, θ) : 0
end

"""
    rand_zinb(rng, μ, π, θ) -> Int

Draw one observation from the ZINB: a structural zero with probability `π`,
otherwise a negative binomial draw.
"""
function rand_zinb(rng::AbstractRNG, μ::Real, π::Real, θ::Real)
    return rand(rng) < π ? 0 : rand_nb(rng, μ, θ)
end

end # module Likelihoods
