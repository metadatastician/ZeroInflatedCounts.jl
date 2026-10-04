# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Test-side helpers: simulated data with planted coefficients, and the
# conversion from the model-matrix parameterisation R uses to the natural
# parameterisation the oracle uses.
#
# These live outside test/cases/ deliberately: runtests.jl includes every .jl
# file under test/cases/ as a case, and these are not cases.

using Random: MersenneTwister

# `Likelihoods` is included into Main by test/runtests.jl before this file, so
# the oracle's names are reached through it directly.

"""
    Simulation

A simulated table: `counts` (samples × 1 taxon column, but shaped for the
package's API), `groups`, `size_factors`, and the planted parameters.
"""
struct Simulation
    counts::Matrix{Int}
    groups::Vector{String}
    size_factors::Vector{Float64}
    beta::Vector{Float64}      # count part: intercept, group, (log size factor is an offset)
    gamma::Vector{Float64}     # zero part: intercept, group, log size factor
    theta::Float64
    kind::Symbol
end

"""
    simulate_table(kind, rng; n_per_group=20, beta=(1.2, 0.9), gamma=(0.4, 0.5, 0.8),
                   theta=2.5, ref="control", contrast="treated") -> Simulation

Simulate one taxon from the same model the package fits, using the oracle's
samplers and the same formulas:

* count part: `μ_i = exp(β₁ + β₂·g_i + log s_i)`, with `g_i = 1` in the
  contrast group;
* hurdle zero part: `p_i = logistic(γ₁ + γ₂·g_i + γ₃·log s_i)`, the probability
  that the count is positive;
* ZINB zero part: `π_i = logistic(γ₁ + γ₂·g_i + γ₃·log s_i)`, the structural
  zero weight.

`counts` is returned as an `n × 1` matrix so it can be handed to the package
unchanged.
"""
function simulate_table(kind::Symbol, rng; n_per_group::Int = 20,
    beta::Tuple = (1.2, 0.9), gamma::Tuple = (0.4, 0.5, 0.8),
    theta::Real = 2.5, ref::String = "control", contrast::String = "treated",
    size_factor_spread::Real = 0.4, size_factors::Union{Nothing,AbstractVector} = nothing,
    groups_in::Union{Nothing,AbstractVector} = nothing)
    n = 2 * n_per_group
    groups = groups_in === nothing ? vcat(fill(ref, n_per_group), fill(contrast, n_per_group)) :
             string.(groups_in)
    n = length(groups)
    g = Float64.(groups .!= ref)
    s = size_factors === nothing ? exp.(size_factor_spread .* randn(rng, n)) :
        Float64.(size_factors)
    log_s = log.(s)

    μ = exp.(beta[1] .+ beta[2] .* g .+ log_s)
    η = gamma[1] .+ gamma[2] .* g .+ gamma[3] .* log_s
    rate = 1 ./ (1 .+ exp.(-η))

    y = Vector{Int}(undef, n)
    if kind === :hurdle_nb
        for i in 1:n
            y[i] = Likelihoods.rand_hurdle(rng, μ[i], rate[i], theta)
        end
    elseif kind === :zinb
        for i in 1:n
            y[i] = Likelihoods.rand_zinb(rng, μ[i], rate[i], theta)
        end
    else
        throw(ArgumentError("unknown kind $(kind)"))
    end

    return Simulation(reshape(y, n, 1), groups, s, [beta[1], beta[2]],
        [gamma[1], gamma[2], gamma[3]], Float64(theta), kind)
end

"""
    group_indicator(groups, ref) -> Vector{Float64}

`1.0` for the contrast group, `0.0` for the reference group: the contrast
coding both the package's R formula and the oracle's formulas use.
"""
group_indicator(groups::AbstractVector, ref::String) = Float64.(groups .!= ref)

"""
    hurdle_probabilities(sim, gamma) -> Vector{Float64}

`p_i = logistic(γ₁ + γ₂ g_i + γ₃ log s_i)`, the probability the count is
positive, from a coefficient vector in the order the package reports.
"""
function hurdle_probabilities(log_s::AbstractVector, g::AbstractVector, gamma::AbstractVector)
    return Likelihoods.invlogit.(gamma[1] .+ gamma[2] .* g .+ gamma[3] .* log_s)
end

"""
    count_means(log_s, g, beta) -> Vector{Float64}

`μ_i = exp(β₁ + β₂ g_i + log s_i)`, the count part's mean.
"""
count_means(log_s::AbstractVector, g::AbstractVector, beta::AbstractVector) =
    exp.(beta[1] .+ beta[2] .* g .+ log_s)

"""
    wald_interval(estimate, se; level=0.95) -> (lo, hi)

Normal-approximation Wald interval. Used by the simulation-recovery case, which
is exactly the claim the method-conditions document makes about intervals: not
a proof, a measured rate.
"""
function wald_interval(estimate::Real, se::Real; z::Real = 1.959963984540054)
    return (estimate - z * se, estimate + z * se)
end

"""
    std_error(xs) -> Float64

Monte-Carlo standard error of the mean of `xs`. Used by the simulation-recovery
case to size its tolerance from the observed spread rather than from a formula
that assumes the estimator's variance is known.
"""
function std_error(xs::AbstractVector)
    n = length(xs)
    n < 2 && return 0.0
    m = sum(xs) / n
    return sqrt(sum((x .- m)^2 for x in xs) / (n - 1) / n)
end
