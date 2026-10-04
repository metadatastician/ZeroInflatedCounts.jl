# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# The oracle checked against values that can be written down by hand.
#
# Nothing here calls R. These are the identities the method-conditions document
# names as proofs' territory — separability, the truncation normaliser, the
# ZINB zero probability — checked numerically here so that the Agda development
# in verification/proofs/agda/ is checking the same statements the code relies
# on.

using Random: MersenneTwister

@testset "oracle — known answers" begin
    # NB(μ = 3, θ = 2), by hand: f(y) = Γ(y+θ)/(Γ(θ)Γ(y+1)) · (2/5)^θ · (3/5)^y.
    #   f(0) = (2/5)^2                        = 0.16
    #   f(1) = (2/1) · 0.16 · 0.6             = 0.192
    #   f(2) = (6/2) · 0.16 · 0.36            = 0.1728
    @test Likelihoods.nb_pmf(0, 3.0, 2.0) ≈ 0.16 atol = 1.0e-12
    @test Likelihoods.nb_pmf(1, 3.0, 2.0) ≈ 0.192 atol = 1.0e-12
    @test Likelihoods.nb_pmf(2, 3.0, 2.0) ≈ 0.1728 atol = 1.0e-12

    # The zero-truncated version renormalises by exactly 1 - f(0) = 0.84.
    @test Likelihoods.tnb_pmf(1, 3.0, 2.0) ≈ 0.192 / 0.84 atol = 1.0e-12
    @test Likelihoods.tnb_pmf(2, 3.0, 2.0) ≈ 0.1728 / 0.84 atol = 1.0e-12

    # log Γ against known factorials.
    @test Likelihoods.loggamma(1.0) ≈ 0.0 atol = 1.0e-12
    @test Likelihoods.loggamma(2.0) ≈ 0.0 atol = 1.0e-12
    @test isapprox(Likelihoods.loggamma(5.0), log(24.0); rtol = 1.0e-12)
    @test isapprox(Likelihoods.loggamma(0.5), log(sqrt(π)); rtol = 1.0e-12)

    # Both densities sum to 1 over the support (the truncation normaliser is
    # 1 - f(0), so the truncated one must too).
    @test sum(Likelihoods.nb_pmf(k, 3.0, 2.0) for k in 0:2000) ≈ 1.0 atol = 1.0e-10
    @test sum(Likelihoods.tnb_pmf(k, 3.0, 2.0) for k in 1:2000) ≈ 1.0 atol = 1.0e-10
end

@testset "oracle — hurdle separability, numerically" begin
    y = [0, 1, 4, 2, 0, 7, 3]
    μ = [1.5, 2.0, 3.5, 2.5, 1.0, 4.0, 2.2]
    p = [0.3, 0.6, 0.8, 0.7, 0.2, 0.9, 0.5]
    θ = 2.5

    zero_part = sum(Likelihoods.hurdle_zero_loglik(y[i], p[i]) for i in eachindex(y))
    count_part = sum(Likelihoods.tnb_logpdf(y[i], μ[i], θ) for i in eachindex(y) if y[i] > 0)

    @test isapprox(Likelihoods.hurdle_loglik(y, μ, p, θ), zero_part + count_part; rtol = 1.0e-14)

    # Moving a count-part parameter leaves the zero part alone: the two parts
    # are separate sums, which is what makes the joint MLE the pair of
    # separately fitted MLEs.
    μ2 = [1.5, 2.0, 3.5, 2.5, 1.0, 4.0, 2.2] .* 3
    moved = Likelihoods.hurdle_loglik(y, μ2, p, θ) - Likelihoods.hurdle_loglik(y, μ, p, θ)
    expected = sum(Likelihoods.tnb_logpdf(y[i], μ2[i], θ) -
                   Likelihoods.tnb_logpdf(y[i], μ[i], θ)
                   for i in eachindex(y) if y[i] > 0)
    @test isapprox(moved, expected; rtol = 1.0e-12)
end

@testset "oracle — ZINB zero probability stays in [0, 1]" begin
    for π in 0.0:0.1:1.0, f0 in 0.0:0.1:1.0
        q = π + (1 - π) * f0
        @test 0.0 <= q <= 1.0
    end
    # And the mixture probability matches a direct sum of the two components.
    y = [0, 0, 2, 5]
    μ = [1.0, 2.0, 3.0, 4.0]
    π = [0.2, 0.7, 0.1, 0.9]
    θ = 1.5
    direct = 0.0
    for i in eachindex(y)
        if y[i] == 0
            direct += log(π[i] + (1 - π[i]) * Likelihoods.nb_pmf(0, μ[i], θ))
        else
            direct += log(1 - π[i]) + Likelihoods.nb_logpdf(y[i], μ[i], θ)
        end
    end
    @test isapprox(Likelihoods.zinb_loglik(y, μ, π, θ), direct; rtol = 1.0e-14)
end

@testset "oracle — samplers reproduce the model means" begin
    rng = MersenneTwister(20261004)
    draws = 200_000
    μ, θ = 3.0, 2.5
    p = 0.6

    nb_total = 0
    for _ in 1:draws
        nb_total += Likelihoods.rand_nb(rng, μ, θ)
    end
    nb_mean = nb_total / draws
    h = Likelihoods.nb_pmf(0, μ, θ)
    @test isapprox(nb_mean, μ; atol = 4 * sqrt(μ * (1 + μ / θ) / draws))

    f0 = Likelihoods.nb_pmf(0, μ, θ)
    hurdle_total = 0
    for _ in 1:draws
        hurdle_total += Likelihoods.rand_hurdle(rng, μ, p, θ)
    end
    hurdle_mean = hurdle_total / draws
    theoretical = p * μ / (1 - f0)
    @test isapprox(hurdle_mean, theoretical; rtol = 0.05)

    zinb_total = 0
    for _ in 1:draws
        zinb_total += Likelihoods.rand_zinb(rng, μ, 0.3, θ)
    end
    @test isapprox(zinb_total / draws, 0.7 * μ; rtol = 0.05)
end
