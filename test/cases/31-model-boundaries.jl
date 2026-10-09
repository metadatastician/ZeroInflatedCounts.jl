# SPDX-License-Identifier: MPL-2.0
# Deterministic model decisions: no optimiser or simulated data is needed.

@testset "prefit floors and precedence" begin
    check = ZeroInflatedCounts.check_inputs
    refusal = ZeroInflatedCounts.prefit_refusal
    groups = repeat(["control", "treated"], inner = 4)
    y = [0, 1, 2, 3, 0, 4, 5, 6]
    for kind in (:hurdle_nb, :zinb)
        # The observed prevalence equals the requested floor: equality is accepted.
        inputs = check(reshape(y, :, 1), groups, ones(8); min_prevalence = 0.75)
        @test refusal(kind, inputs, 1) === nothing
        above = check(reshape(y, :, 1), groups, ones(8); min_prevalence = nextfloat(0.75))
        @test occursin("prevalence", refusal(kind, above, 1))

        for grp in ("control", "treated"), value in (0, 1)
            separated = copy(y)
            separated[groups .== grp] .= value
            inputs = check(reshape(separated, :, 1), groups, ones(8))
            reason = refusal(kind, inputs, 1)
            @test occursin("separated in group $grp", reason)
            @test occursin(value == 0 ? "0 of 4" : "4 of 4", reason)
        end
    end

    # Two positives in either group fails the hurdle floor but is allowed for ZINB.
    for row in (2, 6)
        few = copy(y)
        few[row] = 0
        inputs = check(reshape(few, :, 1), groups, ones(8))
        @test occursin("at least 3 positive counts", refusal(:hurdle_nb, inputs, 1))
        @test refusal(:zinb, inputs, 1) === nothing
    end

    for n in (6, 7)
        inputs = check(reshape(y[1:n], :, 1), groups[1:n], ones(n))
        if n == 6
            @test occursin("6 or fewer samples", refusal(:zinb, inputs, 1))
        else
            @test refusal(:zinb, inputs, 1) === nothing
        end
    end
    for n in (3, 4)
        inputs = check(reshape([0, 1, 0, 1][1:n], :, 1), ["a", "a", "b", "b"][1:n], ones(n))
        expected = n == 3 ? "at least 4 samples" : "at least 3 positive counts"
        @test occursin(expected, refusal(:hurdle_nb, inputs, 1))
    end

    # Sample size wins even when separation and prevalence would also refuse.
    tiny = check(zeros(Int, 3, 1), ["a", "a", "b"], ones(3); min_prevalence = 1)
    @test occursin("at least 4 samples", refusal(:hurdle_nb, tiny, 1))
    @test occursin("6 or fewer samples", refusal(:zinb, tiny, 1))
end

@testset "statistical boundaries and unsorted multiple testing" begin
    # Hand-calculated BH values in input order, with ties and a step-up plateau.
    p = [0.04, 0.001, 0.03, 0.03, 1.0, 0.0]
    original = copy(p)
    expected = [0.048, 0.003, 0.045, 0.045, 1.0, 0.0]
    @test ZeroInflatedCounts.bh_adjust(p) ≈ expected
    @test p == original
    permutation = [6, 4, 1, 5, 2, 3]
    @test ZeroInflatedCounts.bh_adjust(view(p, permutation)) ≈ expected[permutation]
    @test ZeroInflatedCounts.bh_adjust([0, 1]) == [0.0, 1.0]
    @test ZeroInflatedCounts.bh_adjust([0.9, 0.95, 1.0]) == ones(3)
    @test ZeroInflatedCounts.chi2_2_sf(Inf) == 0.0

    @test occursin("upper bound", ZeroInflatedCounts.theta_note_for(1.0e7))
    @test ZeroInflatedCounts.theta_note_for(prevfloat(1.0e7)) === nothing
    @test occursin("lower bound", ZeroInflatedCounts.theta_note_for(1.0e-8))
    @test ZeroInflatedCounts.theta_note_for(nextfloat(1.0e-8)) === nothing
    @test !ZeroInflatedCounts.zero_inflation_detectable([0.0, prevfloat(1.0e-6)])
    @test ZeroInflatedCounts.zero_inflation_detectable([0.0, 1.0e-6])
    @test ZeroInflatedCounts.zero_inflation_detectable([nextfloat(1.0e-6), 0.0])
end
