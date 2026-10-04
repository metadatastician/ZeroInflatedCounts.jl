# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# The independent reference: the same data fitted by a second, independently
# written `pscl` call, and the log-likelihood evaluated by the pure-Julia
# oracle at R's own coefficients.
#
# The R below is NOT src/core/pscl.jl's helper. If the package's bridge
# mis-signs a coefficient, drops the offset, or compares the wrong two
# log-likelihoods, a shared implementation would agree with itself. This one is
# derived from the method-conditions document instead.

using RCall
using Random: MersenneTwister

# Returns, in order:
#   count intercept, count g, count g SE,
#   zero intercept, zero g, zero log(s), zero g SE,
#   theta, logLik(full), logLik(null).
const R_REFERENCE = raw"""
zic_reference <- function(y, g, s, kind, levels) {
    d <- data.frame(y = y, g = factor(g, levels = levels), s = s)
    if (kind == "hurdle") {
        full <- pscl::hurdle(y ~ g + offset(log(s)) | g + log(s), data = d,
                             dist = "negbin", zero.dist = "binomial")
        null <- pscl::hurdle(y ~ offset(log(s)) | log(s), data = d,
                             dist = "negbin", zero.dist = "binomial")
    } else {
        full <- pscl::zeroinfl(y ~ g + offset(log(s)) | g + log(s), data = d,
                               dist = "negbin")
        null <- pscl::zeroinfl(y ~ offset(log(s)) | log(s), data = d,
                               dist = "negbin")
    }
    cc <- coef(full, model = "count")
    zc <- coef(full, model = "zero")
    sc <- summary(full)$coefficients$count
    sz <- summary(full)$coefficients$zero
    c(cc["(Intercept)"], cc[grep("^g", names(cc))[1]],
      sc[grep("^g", rownames(sc))[1], 2],
      zc["(Intercept)"], zc[grep("^g", names(zc))[1]], zc["log(s)"],
      sz[grep("^g", rownames(sz))[1], 2],
      full$theta,
      as.numeric(logLik(full)), as.numeric(logLik(null)))
}
"""

@testset "cross-check: an independently written pscl call" begin
    R"eval(parse(text = $R_REFERENCE))"

    for kind in (:hurdle_nb, :zinb)
        rng = MersenneTwister(kind === :hurdle_nb ? 20261011 : 20261012)
        sim = simulate_table(kind, rng; n_per_group = 25, beta = (1.3, 0.8),
            gamma = (0.5, 0.6, 0.7), theta = 3.0, size_factor_spread = 0.5)
        res = kind === :hurdle_nb ?
              hurdle_nb(sim.counts, sim.groups, sim.size_factors; ref = "control") :
              zinb(sim.counts, sim.groups, sim.size_factors; ref = "control")
        f = fitted(res)
        @test length(f) == 1
        fit = f[1]

        y = Float64.(sim.counts[:, 1])
        s = sim.size_factors
        gg = sim.groups
        kind_string = kind === :hurdle_nb ? "hurdle" : "zeroinfl"
        levels = ["control", "treated"]
        v = rcopy(R"zic_reference($y, $gg, $s, $kind_string, $levels)")
        @test length(v) == 10

        # The package's numbers against the independent call's.
        @test isapprox(fit.log_fold_change, v[2]; rtol = 1.0e-8)
        @test isapprox(fit.se_log_fold_change, v[3]; rtol = 1.0e-8)
        @test isapprox(fit.zero_log_odds, v[5]; rtol = 1.0e-8)
        @test isapprox(fit.se_zero_log_odds, v[7]; rtol = 1.0e-8)
        @test isapprox(fit.theta, v[8]; rtol = 1.0e-8)
        @test isapprox(fit.loglik_full, v[9]; rtol = 1.0e-9)

        # The test recomputed from the two independent log-likelihoods.
        lr = 2 * (v[9] - v[10])
        @test isapprox(fit.lr_statistic, lr; rtol = 1.0e-8)
        @test fit.df == 2
        @test isapprox(fit.pvalue, exp(-lr / 2); rtol = 1.0e-8)

        # The Julia oracle, at R's own coefficients. This is the independent
        # evaluation of the likelihood: a wrong sign, a missing offset or a
        # mis-normalised truncation would show up here as a disagreement.
        log_s = log.(s)
        gind = group_indicator(gg, "control")
        mu = count_means(log_s, gind, [v[1], v[2]])
        eta_zero = v[4] .+ v[5] .* gind .+ v[6] .* log_s
        zero_prob = Likelihoods.invlogit.(eta_zero)
        if kind === :hurdle_nb
            @test isapprox(Likelihoods.hurdle_loglik(y, mu, zero_prob, v[8]), v[9];
                rtol = 1.0e-9)
        else
            @test isapprox(Likelihoods.zinb_loglik(y, mu, zero_prob, v[8]), v[9];
                rtol = 1.0e-9)
        end
    end
end

@testset "the zero part's reading is measured, not assumed" begin
    # Each model is simulated so that the contrast group has the named
    # direction, and the sign of the reported group coefficient is then the
    # measurement.
    rng = MersenneTwister(20261013)

    # ZINB, contrast with MORE structural zeros: the log-odds of a structural
    # zero must come back positive.
    simz = simulate_table(:zinb, rng; n_per_group = 60, beta = (1.5, 0.2),
        gamma = (-0.5, 1.5, 0.3), theta = 4.0)
    fitz = fitted(zinb(simz.counts, simz.groups, simz.size_factors; ref = "control"))[1]
    zeros_ref = sum(simz.counts[simz.groups.=="control"] .== 0) / 60
    zeros_con = sum(simz.counts[simz.groups.=="treated"] .== 0) / 60
    @test zeros_con > zeros_ref
    @test fitz.zero_log_odds > 0

    # Hurdle, contrast with MORE positives (fewer zeros): pscl's hurdle zero
    # part reads as P(count > 0), so its group coefficient must also come back
    # positive here.
    simh = simulate_table(:hurdle_nb, rng; n_per_group = 60, beta = (1.5, 0.2),
        gamma = (0.5, 1.5, 0.3), theta = 4.0)
    fith = fitted(hurdle_nb(simh.counts, simh.groups, simh.size_factors; ref = "control"))[1]
    zeros_ref_h = sum(simh.counts[simh.groups.=="control"] .== 0) / 60
    zeros_con_h = sum(simh.counts[simh.groups.=="treated"] .== 0) / 60
    @test zeros_con_h < zeros_ref_h
    @test fith.zero_log_odds > 0

    # And the reading that makes those signs coherent is the one the fitted
    # probabilities satisfy: for the hurdle, `plogis` of the zero part tracks
    # the observed positives, sample by sample.
    y = Float64.(simh.counts[:, 1])
    s = simh.size_factors
    gg = simh.groups
    kind_string = "hurdle"
    levels = ["control", "treated"]
    v = rcopy(R"zic_reference($y, $gg, $s, $kind_string, $levels)")
    p_hat = Likelihoods.invlogit.(v[4] .+ v[5] .* Float64.(gg .!= "control") .+ v[6] .* log.(s))
    @test isapprox(sum(p_hat) / length(p_hat), sum(y .> 0) / length(y); atol = 0.05)
end
