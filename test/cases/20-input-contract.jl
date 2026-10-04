# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# The data contract, refused rather than repaired.
#
# Every one of these must throw before R is consulted: a table that violates
# the contract is a caller's error, not a per-taxon state, and the whole call is
# rejected.

@testset "input contract — refused, not rounded" begin
    counts = [0 4 1; 3 0 2; 5 2 0; 0 1 6]
    groups = ["a", "a", "b", "b"]
    sizes = [1.0, 1.0, 1.0, 1.0]

    # Non-integer counts: nothing is rounded to make them fit.
    @test_throws ArgumentError hurdle_nb(Float64.(counts), groups, sizes)
    @test_throws ArgumentError zinb(Float64.(counts), groups, sizes)

    # Negative counts.
    @test_throws ArgumentError hurdle_nb(counts .- 10, groups, sizes)

    # Row counts that do not line up.
    @test_throws ArgumentError hurdle_nb(counts, groups[1:3], sizes)
    @test_throws ArgumentError hurdle_nb(counts, groups, sizes[1:2])

    # Anything other than exactly two groups.
    @test_throws ArgumentError hurdle_nb(counts, ["a", "a", "b", "c"], sizes)
    @test_throws ArgumentError zinb(counts, ["a", "a", "a", "a"], sizes)

    # Size factors must be finite and strictly positive.
    @test_throws ArgumentError hurdle_nb(counts, groups, [-1.0, 1.0, 1.0, 1.0])
    @test_throws ArgumentError hurdle_nb(counts, groups, [0.0, 1.0, 1.0, 1.0])
    @test_throws ArgumentError hurdle_nb(counts, groups, [Inf, 1.0, 1.0, 1.0])

    # A reference group that does not exist, and an impossible prevalence floor.
    @test_throws ArgumentError hurdle_nb(counts, groups, sizes; ref = "c")
    @test_throws ArgumentError hurdle_nb(counts, groups, sizes; min_prevalence = 1.5)

    # A named taxa vector of the wrong length.
    @test_throws ArgumentError hurdle_nb(counts, groups, sizes; taxa = ["only-one"])

    # The valid call does not throw (it is refused per taxon at most).
    res = hurdle_nb(counts, groups, sizes)
    @test res isa ZIFit
    @test res.provenance.method === :hurdle_nb
    @test length(res.fits) == size(counts, 2)
end
