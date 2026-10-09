# SPDX-License-Identifier: MPL-2.0
# Control the R bridge response to exercise Julia's decisions without relying
# on an optimiser to produce a particular failure. Restore every replaced R
# binding even when a test throws, so the later pscl cross-checks remain real.

using RCall

function with_pscl_fixture(f, fit_function; available = true)
    saved_fit = R"zic_fit"
    saved_available = R"zic_has_pscl"
    saved_versions = R"zic_versions"
    try
        R"zic_fit <- $fit_function"
        R"zic_has_pscl <- local({ value <- $available; function() value })"
        R"""
        zic_versions <- function() list(r_version = "fixture-R",
            pscl_version = "fixture-pscl", mass_version = "fixture-MASS")
        """
        return f()
    finally
        R"zic_fit <- $saved_fit"
        R"zic_has_pscl <- $saved_available"
        R"zic_versions <- $saved_versions"
    end
end

@testset "controlled fit responses" begin
    counts = reshape([0, 1, 2, 3, 0, 4, 5, 6], :, 1)
    groups = repeat(["control", "treated"], inner = 4)
    sizes = [0.5, 1.0, 1.5, 2.0, 0.75, 1.25, 1.75, 2.25]
    base_response = R"""
    list(ok = TRUE, count_coef = 0.7, count_se = 0.2,
         zero_coef = -0.4, zero_se = 0.3, theta = 2.5,
         loglik_full = -10, loglik_null = -13, converged = TRUE,
         zero_part_prob = rep(0.25, 8), warnings = c("first warning", "second warning"))
    """

    for (kind, fit, r_kind) in
        ((:hurdle_nb, hurdle_nb, "hurdle"), (:zinb, zinb, "zeroinfl"))
        @testset "$kind forwards inputs and reports a successful fit" begin
            fit_function = R"""
            local({
                answer <- $base_response
                expected_y <- as.numeric($counts)
                expected_g <- $groups
                expected_s <- $sizes
                expected_kind <- $r_kind
                function(y, g, s, kind, levels) {
                    stopifnot(identical(y, expected_y), identical(g, expected_g),
                              identical(s, expected_s), identical(kind, expected_kind),
                              identical(levels, c("treated", "control")))
                    answer
                }
            })
            """
            with_pscl_fixture(fit_function) do
                res = fit(counts, groups, sizes; ref = "treated", taxa = ["named"])
                @test res.method === kind
                @test res.taxa == ["named"]
                @test length(fitted(res)) == 1
                f = only(fitted(res))
                @test f.taxon == "named"
                @test (f.log_fold_change, f.se_log_fold_change) == (0.7, 0.2)
                @test (f.zero_log_odds, f.se_zero_log_odds) == (-0.4, 0.3)
                @test f.theta == 2.5
                @test (f.loglik_full, f.loglik_null, f.lr_statistic, f.df) ==
                      (-10, -13, 6, 2)
                @test f.pvalue ≈ exp(-3)
                @test f.p_adjusted == f.pvalue
                @test f.theta_note === nothing
                @test f.warnings == ["first warning", "second warning"]
                p = res.provenance
                @test p.method === kind
                @test (p.n_tested, p.n_refused, p.n_failed) == (1, 0, 0)
                @test (p.r_version, p.pscl_version, p.mass_version) ==
                      ("fixture-R", "fixture-pscl", "fixture-MASS")
                @test p.count_formula == "count ~ group + offset(log(size_factor))"
                @test p.zero_formula == "~ group + log(size_factor)"
                @test p.null_count_formula == "count ~ offset(log(size_factor))"
                @test p.null_zero_formula == "~ log(size_factor)"
                @test occursin("as an offset", p.offset_method)
                @test occursin("as a covariate", p.offset_method)
                @test p.warnings == ["named: first warning", "named: second warning"]
            end
        end

        @testset "$kind rejects unusable fit responses" begin
            cases = (
                ("loglik_full", NaN, "non-finite log-likelihood"),
                ("loglik_null", -Inf, "non-finite log-likelihood"),
                ("theta", Inf, "size parameter"),
                ("converged", false, "non-convergence"),
                ("count_se", NaN, "standard error"),
                ("zero_se", Inf, "standard error"),
                ("count_se", missing, "standard error"),
            )
            for (field, value, reason) in cases
                response = R"modifyList($base_response, setNames(list($value), $field))"
                fit_function = R"local({ answer <- $response; function(...) answer })"
                with_pscl_fixture(fit_function) do
                    res = fit(counts, groups, sizes)
                    @test only(res.fits) isa TaxonFailed
                    @test occursin(reason, only(failed(res)).message)
                    @test isempty(pvalues(res))
                    @test isempty(adjusted_pvalues(res))
                    @test (res.provenance.n_tested, res.provenance.n_failed) == (0, 1)
                end
            end
            for message in ("singular covariance from R", "")
                fit_function = R"""
                local({ msg <- $message; function(...) list(ok = FALSE, message = msg) })
                """
                with_pscl_fixture(fit_function) do
                    res = fit(counts, groups, sizes)
                    expected =
                        isempty(message) ? "pscl returned no fit and no message" : message
                    @test only(failed(res)).message == expected
                end
            end
            with_pscl_fixture(R"function(...) stop('controlled bridge error')") do
                res = fit(counts, groups, sizes)
                @test occursin("RCall:", only(failed(res)).message)
                @test occursin("controlled bridge error", only(failed(res)).message)
            end
        end

        @testset "$kind clamps negative LR and retains theta boundary fits" begin
            for theta in (1.0e-8, 1.0e7)
                response =
                    R"modifyList($base_response, list(loglik_full = -14, theta = $theta))"
                fit_function = R"local({ answer <- $response; function(...) answer })"
                with_pscl_fixture(fit_function) do
                    res = fit(counts, groups, sizes)
                    f = only(fitted(res))
                    @test f.lr_statistic == 0.0
                    @test f.pvalue == f.p_adjusted == 1.0
                    @test f.theta == theta
                    @test f.theta_note !== nothing
                    @test f.warnings[1:2] == ["first warning", "second warning"]
                    @test length(f.warnings) == 3
                    @test occursin("reduced model fitted better", last(f.warnings))
                    @test occursin(
                        "taxon-1: the reduced model",
                        last(res.provenance.warnings),
                    )
                end
            end
        end
    end

    @testset "ZINB collapse, empty probabilities, and the exact floor" begin
        for (probabilities, expected_type) in (
            (Float64[], TaxonFailed),
            ([prevfloat(1.0e-6)], TaxonRefused),
            ([0.0, 1.0e-6], TaxonFit),
        )
            response = R"modifyList($base_response, list(zero_part_prob = $probabilities))"
            fit_function = R"local({ answer <- $response; function(...) answer })"
            with_pscl_fixture(fit_function) do
                res = zinb(counts, groups, sizes)
                @test only(res.fits) isa expected_type
                if expected_type === TaxonRefused
                    @test only(refused(res)).reason ==
                          "no zero inflation is detectable; the NB GLM result applies"
                    @test isempty(pvalues(res))
                elseif expected_type === TaxonFailed
                    @test occursin(
                        "no fitted zero-part probability",
                        only(failed(res)).message,
                    )
                    @test isempty(pvalues(res))
                end
                # Hurdle fits do not apply the ZINB collapse rule.
                @test length(fitted(hurdle_nb(counts, groups, sizes))) == 1
            end
        end
    end

    @testset "missing pscl and validation precedence" begin
        with_pscl_fixture(R"function(...) stop('must not fit')"; available = false) do
            for fit in (hurdle_nb, zinb)
                @test_throws r"pscl.*not installed" fit(counts, groups, sizes)
                @test_throws ArgumentError fit(Float64.(counts), groups, sizes)
            end
        end
    end

    @testset "only successful taxa enter the BH family" begin
        fit_function = R"""
        local({
            answer <- $base_response
            function(y, ...) {
                if (y[1] == 0) stop("prefit refusal must not reach R")
                if (y[1] == 2) return(list(ok = FALSE, message = "failed middle taxon"))
                p <- if (y[1] == 1) 0.04 else 0.01
                answer[["loglik_full"]] <- -10
                answer[["loglik_null"]] <- -10 + log(p)
                answer[["warnings"]] <- "fixture warning"
                answer
            }
        })
        """
        with_pscl_fixture(fit_function) do
            mixed = hcat(
                [1, 0, 2, 3, 0, 4, 5, 6],
                zeros(Int, 8),
                [2, 0, 2, 3, 0, 4, 5, 6],
                [3, 0, 2, 3, 0, 4, 5, 6],
            )
            taxa = ["larger-p", "refused", "failed", "smaller-p"]
            for fit in (hurdle_nb, zinb)
                res = fit(mixed, groups, sizes; taxa = taxa)
                @test [r.taxon for r in res.fits] == taxa
                @test [r.taxon for r in fitted(res)] == ["larger-p", "smaller-p"]
                @test only(refused(res)).taxon == "refused"
                @test only(failed(res)).message == "failed middle taxon"
                @test pvalues(res) ≈ [0.04, 0.01]
                @test adjusted_pvalues(res) ≈ [0.04, 0.02]
                p = res.provenance
                @test (p.n_tested, p.n_refused, p.n_failed) == (2, 1, 1)
                @test p.warnings ==
                      ["larger-p: fixture warning", "smaller-p: fixture warning"]

                all_refused = fit(zeros(Int, 8, 2), groups, sizes)
                @test length(refused(all_refused)) == 2
                @test isempty(pvalues(all_refused))
                @test isempty(adjusted_pvalues(all_refused))
                @test isempty(all_refused.provenance.warnings)
                @test all_refused.provenance.n_tested == 0
            end
        end
    end
end

@testset "R fixtures restore the engine after an exception" begin
    saved_fit = R"zic_fit"
    saved_available = R"zic_has_pscl"
    saved_versions = R"zic_versions"
    @test_throws ErrorException with_pscl_fixture(R"function(...) NULL") do
        error("intentional fixture failure")
    end
    @test rcopy(R"identical(zic_fit, $saved_fit)")
    @test rcopy(R"identical(zic_has_pscl, $saved_available)")
    @test rcopy(R"identical(zic_versions, $saved_versions)")
end
