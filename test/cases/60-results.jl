# SPDX-License-Identifier: MPL-2.0
# Result inspection and reporting are independent of R and model fitting.

@testset "result accessors, adjustment, and human-readable summaries" begin
    original = TaxonFit(
        "tested", 0.75, 0.25, -0.5, 0.125, 2.5, -10.0, -13.0,
        6.0, 2, 0.05, NaN, "theta note", ["fit warning"],
    )
    adjusted = ZeroInflatedCounts.set_adjusted_p(original, 0.1)
    @test adjusted.p_adjusted == 0.1
    @test isnan(original.p_adjusted)
    for field in fieldnames(TaxonFit)
        if field !== :p_adjusted
            @test isequal(getfield(adjusted, field), getfield(original, field))
        end
    end
    refusal = TaxonRefused("refused", "too few positives")
    failure = TaxonFailed("failed", "singular covariance")

    for kind in (:hurdle_nb, :zinb)
        provenance = Provenance(
            kind, "count formula", "zero formula", "null count", "null zero",
            "offset method", "fixture-R", "fixture-pscl", "fixture-MASS",
            1, 1, 1, ["tested: fit warning"],
        )
        res = ZIFit(kind, ["refused", "tested", "failed"],
                    ZeroInflatedCounts.TaxonResult[refusal, adjusted, failure], provenance)
        @test fitted(res) == [adjusted]
        @test refused(res) == [refusal]
        @test failed(res) == [failure]
        @test pvalues(res) == [0.05]
        @test adjusted_pvalues(res) == [0.1]
        @test sprint(show, res) == "ZIFit($kind: 1 tested, 1 refused, 1 failed)"

        lines = split(summary_table(res), '\n')
        @test lines[1] == "method: $kind"
        interpretation = kind === :hurdle_nb ? "P(count > 0)" : "P(structural zero)"
        @test occursin(interpretation, lines[2])
        # Tokenise columns so whitespace padding is not part of the contract.
        @test split(lines[4])[1:7] == ["refused", "-", "-", "-", "-", "-", "-"]
        @test endswith(lines[4], "refused: too few positives")
        @test split(lines[5]) == ["tested", "0.75", "0.25", "-0.5", "2.5", "0.05", "0.1", "tested"]
        @test split(lines[6])[1:7] == ["failed", "-", "-", "-", "-", "-", "-"]
        @test endswith(lines[6], "failed: singular covariance")
        @test lines[7] == "R fixture-R; pscl fixture-pscl; MASS fixture-MASS"

        empty_provenance = Provenance(
            kind, "", "", "", "", "", "R", "pscl", "MASS", 0, 0, 0, String[],
        )
        empty_result = ZIFit(kind, String[], ZeroInflatedCounts.TaxonResult[], empty_provenance)
        @test fitted(empty_result) isa Vector{TaxonFit}
        @test refused(empty_result) isa Vector{TaxonRefused}
        @test failed(empty_result) isa Vector{TaxonFailed}
        @test pvalues(empty_result) == Float64[]
        @test adjusted_pvalues(empty_result) == Float64[]
        @test length(split(summary_table(empty_result), '\n')) == 4
    end
    @test occursin("LFC=0.75 (SE 0.25)", sprint(show, adjusted))
    @test occursin("p_adj=0.1", sprint(show, adjusted))
    @test sprint(show, refusal) == "TaxonRefused(refused: too few positives)"
    @test sprint(show, failure) == "TaxonFailed(failed: singular covariance)"
end
