# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# The negative controls: every refusal in the method-conditions document, and
# the rule that a refused taxon is outside the multiple-testing family.

@testset "refusals happen for the stated reasons" begin
    # 8 samples, 4 per group. Columns:
    #  1 "separated"      — only zeros in the control group
    #  2 "few-positives"  — 2 positive counts in the control group
    #  3 "fine"           — fittable in both models
    counts = [
        0 3 1;
        0 0 2;
        0 0 0;
        0 2 4;
        2 1 0;
        3 0 2;
        1 4 1;
        4 3 3
    ]
    groups = [
        "control",
        "control",
        "control",
        "control",
        "treated",
        "treated",
        "treated",
        "treated",
    ]
    sizes = ones(8)
    taxa = ["separated", "few-positives", "fine"]

    res = hurdle_nb(counts, groups, sizes; taxa = taxa)
    refused_names = [r.taxon for r in refused(res)]
    @test "separated" in refused_names
    @test "few-positives" in refused_names
    @test [f.taxon for f in fitted(res)] == ["fine"]
    reason_of(r, name) = r.taxon == name ? r.reason : ""
    @test occursin("zero part is separated", reason_of(refused(res)[1], "separated"))
    @test occursin(
        "at least 3 positive counts",
        reason_of(refused(res)[2], "few-positives"),
    )

    # A refused taxon has no p-value: the family is exactly the fitted taxa.
    @test length(pvalues(res)) == 1
    @test length(adjusted_pvalues(res)) == 1
    @test all(isfinite, adjusted_pvalues(res))

    # The same separation rule applies to the ZINB; the hurdle-only floor does
    # not. Whether pscl can fit 6 parameters on 8 samples is pscl's business, so
    # this asserts the *refusal* set and that every taxon is accounted for, not
    # that a particular fit succeeds.
    resz = zinb(counts, groups, sizes; taxa = taxa)
    @test [r.taxon for r in refused(resz)] == ["separated"]
    @test length(fitted(resz)) + length(failed(resz)) == 2
    @test all(f -> isempty(f.message) == false, failed(resz))
    @test length(pvalues(resz)) == length(fitted(resz))
    @test all(isfinite, pvalues(resz))

    # ZINB is fitted jointly with 6 parameters: 6 or fewer samples is refused
    # for every taxon, and nothing is fitted.
    small = counts[1:6, :]
    res_small = zinb(small, groups[1:6], sizes[1:6]; taxa = taxa)
    @test length(fitted(res_small)) == 0
    @test length(refused(res_small)) == 3
    @test all(r -> occursin("6 or fewer samples", r.reason), refused(res_small))

    # The hurdle zero part needs at least 4 samples.
    tiny_rows = [1, 2, 5] # Keep both groups so input validation reaches the sample floor.
    res_tiny =
        hurdle_nb(counts[tiny_rows, :], groups[tiny_rows], sizes[tiny_rows]; taxa = taxa)
    @test length(fitted(res_tiny)) == 0
    @test length(refused(res_tiny)) == 3
    @test all(r -> occursin("at least 4 samples", r.reason), refused(res_tiny))

    # The caller's prevalence floor is applied, and says so. This table passes
    # every earlier check — 8 samples, no separation, 3 positives per group — so
    # the prevalence floor is the only reason left to refuse it (6/8 = 0.75).
    prev_counts = [
        0 1 1 1;
        1 0 1 1;
        1 1 0 1;
        1 1 1 0;
        0 1 1 1;
        1 0 1 1;
        1 1 0 1;
        1 1 1 0
    ]
    prev_taxa = ["p1", "p2", "p3", "p4"]
    res_prev = hurdle_nb(prev_counts, groups, sizes; taxa = prev_taxa, min_prevalence = 0.9)
    @test length(fitted(res_prev)) == 0
    @test length(refused(res_prev)) == 4
    @test all(r -> occursin("prevalence", r.reason), refused(res_prev))

    # Refusals are reported per taxon, in input order, and carry no p-value.
    @test [r.taxon for r in res_prev.fits] == prev_taxa
end

@testset "the test's own arithmetic" begin
    # exp(-x/2) is the chi-squared(2) upper tail: the 95% point is 5.991464...,
    # so 5.991464547107979 must come back as 0.05.
    @test ZeroInflatedCounts.chi2_2_sf(0.0) == 1.0
    @test isapprox(ZeroInflatedCounts.chi2_2_sf(5.991464547107979), 0.05; rtol = 1.0e-12)

    # Benjamini–Hochberg by the definition, including the step-up.
    @test ZeroInflatedCounts.bh_adjust(Float64[]) == Float64[]
    @test ZeroInflatedCounts.bh_adjust([0.4]) == [0.4]
    @test isapprox(ZeroInflatedCounts.bh_adjust([0.01, 0.02, 0.03]), [0.03, 0.03, 0.03])
    @test isapprox(ZeroInflatedCounts.bh_adjust([0.001, 0.5, 0.6]), [0.003, 0.6, 0.6])
    # An inconclusive p-value never comes back below 1.
    @test all(p -> p <= 1.0, ZeroInflatedCounts.bh_adjust([0.9, 0.95, 1.0]))

    # The θ bound notes, and the collapse rule, on their own branch.
    @test ZeroInflatedCounts.theta_note_for(2.5) === nothing
    @test ZeroInflatedCounts.theta_note_for(1.0e8) !== nothing
    @test ZeroInflatedCounts.theta_note_for(1.0e-9) !== nothing
    @test !ZeroInflatedCounts.zero_inflation_detectable([1.0e-9, 1.0e-8])
    @test ZeroInflatedCounts.zero_inflation_detectable([1.0e-9, 1.0e-4])
    # No probabilities at all is not "detectable": nothing was fitted, so
    # nothing may be reported as tested.
    @test !ZeroInflatedCounts.zero_inflation_detectable(Float64[])
end
