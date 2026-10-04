# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Simulation recovery and a null simulation.
#
# The method-conditions document is explicit that whether an estimator recovers
# the truth is a TESTED claim, never a proved one. This is that test:
#
#   * data simulated from each model with planted coefficients, the planted
#     values required to fall inside the reported 95% intervals at the stated
#     rate across repeats;
#   * a null simulation, in which the Benjamini–Hochberg false discovery rate is
#     measured and reported.
#
# The repetition counts are small on purpose — this runs inside `Pkg.test()` —
# so the coverage floor is deliberately loose (0.70 for a nominal 0.95): it is
# a regression alarm for a broken estimator, not a calibration certificate.

using Random: MersenneTwister

# NOTE: no `continue`/`break` inside a `for` inside `@testset` — Test's
# loop-form detection re-parses the loop body at top level (see runtests.jl).
# The guards below are `if` blocks for that reason.

@testset "simulation recovery — hurdle NB" begin
    rng = MersenneTwister(20261005)
    repeats = 30
    planted_lfc = 0.9
    planted_zero = 0.5
    covered_lfc = 0
    covered_zero = 0
    lfc = Float64[]
    zero = Float64[]
    for _ in 1:repeats
        sim = simulate_table(:hurdle_nb, rng; n_per_group = 20,
            beta = (1.2, planted_lfc), gamma = (0.4, planted_zero, 0.8), theta = 2.5)
        res = hurdle_nb(sim.counts, sim.groups, sim.size_factors; ref = "control")
        f = fitted(res)
        if length(f) == 1
            lo, hi = wald_interval(f[1].log_fold_change, f[1].se_log_fold_change)
            covered_lfc += (lo <= planted_lfc <= hi)
            push!(lfc, f[1].log_fold_change)
            if isfinite(f[1].se_zero_log_odds)
                zlo, zhi = wald_interval(f[1].zero_log_odds, f[1].se_zero_log_odds)
                covered_zero += (zlo <= planted_zero <= zhi)
                push!(zero, f[1].zero_log_odds)
            end
        end
    end

    @test length(lfc) >= repeats - 2
    @test covered_lfc / repeats >= 0.70
    @test isapprox(sum(lfc) / length(lfc), planted_lfc; atol = 4 * std_error(lfc))
    # The hurdle zero part reads as P(count > 0), so the planted positive-vs-zero
    # log odds γ₂ is on the same scale as the reported coefficient.
    @test length(zero) >= repeats - 2
    @test covered_zero / repeats >= 0.70
end

@testset "simulation recovery — ZINB" begin
    rng = MersenneTwister(20261006)
    repeats = 30
    planted_lfc = 0.9
    planted_zero = 0.5
    covered_lfc = 0
    covered_zero = 0
    lfc = Float64[]
    zero = Float64[]
    theta = Float64[]
    for _ in 1:repeats
        sim = simulate_table(:zinb, rng; n_per_group = 20,
            beta = (1.2, planted_lfc), gamma = (0.4, planted_zero, 0.8), theta = 2.5)
        res = zinb(sim.counts, sim.groups, sim.size_factors; ref = "control")
        f = fitted(res)
        if length(f) == 1
            lo, hi = wald_interval(f[1].log_fold_change, f[1].se_log_fold_change)
            covered_lfc += (lo <= planted_lfc <= hi)
            push!(lfc, f[1].log_fold_change)
            if isfinite(f[1].se_zero_log_odds)
                # A positive planted γ₂ means more structural zeros in the
                # contrast group, and the reported log odds is on that scale.
                zlo, zhi = wald_interval(f[1].zero_log_odds, f[1].se_zero_log_odds)
                covered_zero += (zlo <= planted_zero <= zhi)
                push!(zero, f[1].zero_log_odds)
                push!(theta, f[1].theta)
            end
        end
    end

    @test length(lfc) >= repeats - 2
    @test covered_lfc / repeats >= 0.70
    @test isapprox(sum(lfc) / length(lfc), planted_lfc; atol = 4 * std_error(lfc))
    @test length(zero) >= repeats - 2
    @test sum(zero) / length(zero) > 0
    @test covered_zero / repeats >= 0.70
    # θ is recovered within a factor, not to three digits: it is the hardest
    # parameter here and this is a smoke check on the sign of the mistake.
    @test 1.0 <= sum(theta) / length(theta) <= 8.0
end

@testset "null simulation — the reported false discovery rate" begin
    rng = MersenneTwister(20261007)
    repeats = 20
    n_taxa = 10
    discoveries = 0
    tested = 0
    for _ in 1:repeats
        seed_taxon = simulate_table(:hurdle_nb, rng; n_per_group = 20,
            beta = (1.2, 0.0), gamma = (0.4, 0.0, 0.8), theta = 2.5)
        groups = seed_taxon.groups
        sizes = seed_taxon.size_factors
        columns = Vector{Vector{Int}}(undef, n_taxa)
        columns[1] = seed_taxon.counts[:, 1]
        for j in 2:n_taxa
            columns[j] = simulate_table(:hurdle_nb, rng; n_per_group = 20,
                beta = (1.2, 0.0), gamma = (0.4, 0.0, 0.8), theta = 2.5,
                size_factors = sizes, groups_in = groups).counts[:, 1]
        end
        counts = reduce(hcat, columns)
        res = hurdle_nb(counts, groups, sizes)
        for f in fitted(res)
            tested += 1
            discoveries += (f.p_adjusted < 0.05)
        end
    end

    fdr = tested == 0 ? 0.0 : discoveries / tested
    @info "null simulation: Benjamini–Hochberg false discovery rate" fdr = fdr tested = tested
    @test tested > 0
    @test fdr <= 0.15
end
