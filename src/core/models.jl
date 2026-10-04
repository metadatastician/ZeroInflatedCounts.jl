# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# The two models, and the one test they report.
#
# `hurdle_nb` and `zinb` differ only in the R function they call and in the
# refusal floors; everything else — validation, the per-taxon test, the
# multiple-testing rule, the provenance — is shared, which is the point of
# keeping them in one file.
#
# The per-taxon test is a likelihood-ratio test of the group effect in both
# parts together, against the same model with the group term dropped from both
# parts. It has 2 degrees of freedom, so the p-value is exp(-LR/2), the exact
# chi-squared(2) upper tail; no approximation and no table is involved.

# The NB size parameter's bounds, as in the NB GLM the consumer already runs.
const THETA_UPPER = 1.0e7
const THETA_LOWER = 1.0e-8

# The zero-part weight below which a ZINB's mixture is indistinguishable from
# an ordinary NB. docs/method-conditions/: "zero part collapsed".
const ZERO_WEIGHT_FLOOR = 1.0e-6

"""
    bh_adjust(p) -> Vector{Float64}

Benjamini–Hochberg step-up adjustment, computed directly from the definition:
with `p` sorted ascending, the adjusted value of the k-th is
`min_{j ≥ k} (n/j) p_(j)`, clamped at 1. Written out rather than called from a
package so that the rule the method-conditions document names is the rule that
runs.
"""
function bh_adjust(p::AbstractVector{<:Real})
    n = length(p)
    n == 0 && return Float64[]
    ord = sortperm(p)
    adj = Vector{Float64}(undef, n)
    running = 1.0
    for k in n:-1:1
        i = ord[k]
        running = min(running, Float64(p[i]) * n / k)
        adj[i] = running
    end
    return adj
end

"""
    chi2_2_sf(x) -> Float64

Upper tail of the chi-squared distribution with 2 degrees of freedom,
`exp(-x/2)`. The likelihood-ratio test in this package always has 2 degrees of
freedom and this is its exact tail, not an approximation.
"""
chi2_2_sf(x::Real) = exp(-x / 2)

"""
    theta_note_for(theta) -> Union{Nothing,String}

The note the NB GLM already gives when the size parameter sits on a bound, or
`nothing` when it is interior. `theta` is kept either way: this is a note, not a
refusal.
"""
function theta_note_for(theta::Real)
    if theta >= THETA_UPPER
        return "theta at the upper bound (>= 1e7): the count part is effectively " *
               "Poisson for this taxon"
    elseif theta <= THETA_LOWER
        return "theta at the lower bound (<= 1e-8): the count part is effectively a " *
               "point mass at the mean for this taxon"
    end
    return nothing
end

"""
    zero_inflation_detectable(zero_part_prob) -> Bool

`true` when at least one fitted zero weight reaches the floor. A ZINB whose
zero weight is below the floor in *every* sample is indistinguishable from an
ordinary NB: the test is refused rather than reported as a null result. An
empty vector (the fit returned no probabilities at all) is not detectable, so
the taxon cannot be reported as tested.
"""
zero_inflation_detectable(zero_part_prob::AbstractVector) =
    !isempty(zero_part_prob) && !all(w -> w < ZERO_WEIGHT_FLOOR, zero_part_prob)

"""
    prefit_refusal(kind, inputs, j) -> Union{Nothing,String}

The conditions the method-conditions document requires to be detected *from the
counts, before fitting*, in the order they are checked here:

1. the sample floors (ZINB: more than 6 samples; hurdle: at least 4 samples);
2. the zero part is separated in either group — no zeros at all, or nothing but
   zeros;
3. the hurdle count part's floor of 3 positive counts per group;
4. the caller's `min_prevalence`, when one was configured.

Returns the reason as the document words it, or `nothing` when the taxon may be
fitted.
"""
function prefit_refusal(kind::Symbol, inputs::ZIInputs, j::Int)
    y = view(inputs.counts, :, j)
    n = length(y)

    if kind === :zinb && n <= 6
        return "a zinb is fitted jointly with 6 parameters: a run with 6 or fewer " *
               "samples is refused (n = $n)"
    end
    if kind === :hurdle_nb && n < 4
        return "the hurdle zero part has 3 parameters fitted on all samples: at " *
               "least 4 samples are required (n = $n)"
    end

    positive = Dict{String,Int}()
    for grp in (inputs.ref, inputs.contrast)
        sel = inputs.groups .== grp
        positive[grp] = count(>(0), y[sel])
        if positive[grp] == 0 || positive[grp] == count(sel)
            return "the zero part is separated in group $grp " *
                   "($(positive[grp]) of $(count(sel)) samples are positive)"
        end
    end

    if kind === :hurdle_nb
        fewest = min(positive[inputs.ref], positive[inputs.contrast])
        if fewest < 3
            return "the hurdle count part needs at least 3 positive counts in each " *
                   "group (reference $(inputs.ref): $(positive[inputs.ref]), " *
                   "contrast $(inputs.contrast): $(positive[inputs.contrast]))"
        end
    end

    if inputs.min_prevalence > 0
        prev = count(>(0), y) / n
        if prev < inputs.min_prevalence
            return "prevalence $(round(prev; digits = 4)) is below the configured " *
                   "minimum $(inputs.min_prevalence)"
        end
    end

    return nothing
end

"""
    fit_taxon(kind, inputs, j) -> TaxonResult

Fit and test one taxon. Never throws for a fitting problem: a taxon that cannot
be fitted is a [`TaxonFailed`](@ref) carrying R's own message.
"""
function fit_taxon(kind::Symbol, inputs::ZIInputs, j::Int)
    taxon = inputs.taxa[j]

    refusal = prefit_refusal(kind, inputs, j)
    refusal === nothing || return TaxonRefused(taxon, refusal)

    y = Float64.(view(inputs.counts, :, j))
    g = inputs.groups
    s = inputs.size_factors
    kind_string = kind === :hurdle_nb ? "hurdle" : "zeroinfl"
    levels = [inputs.ref, inputs.contrast]

    raw = try
        rcopy(R"zic_fit($y, $g, $s, $kind_string, $levels)")
    catch err
        return TaxonFailed(taxon, "RCall: " * sprint(showerror, err))
    end

    warnings = String[_as_string(w) for w in _as_strings(raw, "warnings")]

    if !_as_bool(raw, "ok", false)
        message = _as_string(_rfield(raw, "message"))
        isempty(message) && (message = "pscl returned no fit and no message")
        return TaxonFailed(taxon, message)
    end

    loglik_full = _as_float(raw, "loglik_full")
    loglik_null = _as_float(raw, "loglik_null")
    theta = _as_float(raw, "theta")
    se_count = _as_float(raw, "count_se")
    se_zero = _as_float(raw, "zero_se")

    if !(isfinite(loglik_full) && isfinite(loglik_null) && isfinite(theta))
        return TaxonFailed(taxon, "the fit returned a non-finite log-likelihood or " *
                                  "negative-binomial size parameter")
    end

    if !_as_bool(raw, "converged", true)
        return TaxonFailed(taxon, "pscl reported non-convergence (converged = FALSE)")
    end

    if !isfinite(se_count) || !isfinite(se_zero)
        return TaxonFailed(taxon, "the covariance matrix is not positive definite: " *
                                  "the standard error of the group coefficient is not finite")
    end

    theta_note = theta_note_for(theta)

    if kind === :zinb
        zero_prob = _as_floats(raw, "zero_part_prob")
        if isempty(zero_prob)
            return TaxonFailed(taxon, "the pscl result carries no fitted zero-part probability")
        end
        if !zero_inflation_detectable(zero_prob)
            return TaxonRefused(taxon, "no zero inflation is detectable; the NB GLM result applies")
        end
    end

    lr = 2 * (loglik_full - loglik_null)
    if lr < 0
        push!(warnings, "the reduced model fitted better than the full model " *
                        "(2*delta = $(lr)); the likelihood-ratio statistic is reported as 0")
        lr = 0.0
    end

    return TaxonFit(taxon, _as_float(raw, "count_coef"), se_count,
        _as_float(raw, "zero_coef"), se_zero, theta, loglik_full, loglik_null, lr,
        2, chi2_2_sf(lr), NaN, theta_note, warnings)
end

"""
    fit_table(kind, counts, groups, size_factors; kwargs...) -> ZIFit

The shared body of [`hurdle_nb`](@ref) and [`zinb`](@ref): validate once, fit
taxon by taxon, then apply Benjamini–Hochberg across exactly the taxa that
produced a p-value.
"""
function fit_table(kind::Symbol, counts::AbstractMatrix, groups::AbstractVector,
    size_factors::AbstractVector; taxa = nothing, ref = nothing,
    min_prevalence::Real = 0.0)
    inputs = check_inputs(counts, groups, size_factors; taxa = taxa, ref = ref,
        min_prevalence = min_prevalence)

    if !pscl_available()
        error("the R package `pscl` is not installed for this RCall session: the " *
              "$(kind) method is unavailable. Install it in the R that RCall uses " *
              "with install.packages(\"pscl\").")
    end

    versions = pscl_versions()

    fits = TaxonResult[fit_taxon(kind, inputs, j) for j in 1:size(inputs.counts, 2)]

    tested = findall(r -> r isa TaxonFit, fits)
    if !isempty(tested)
        adjusted = bh_adjust([fits[i].pvalue for i in tested])
        for (k, i) in enumerate(tested)
            fits[i] = set_adjusted_p(fits[i], adjusted[k])
        end
    end

    n_tested = length(tested)
    n_refused = count(r -> r isa TaxonRefused, fits)
    n_failed = count(r -> r isa TaxonFailed, fits)

    all_warnings = String[]
    for r in fits
        r isa TaxonFit || continue
        for w in r.warnings
            push!(all_warnings, r.taxon * ": " * w)
        end
    end

    count_formula = "count ~ group + offset(log(size_factor))"
    zero_formula = "~ group + log(size_factor)"
    offset_method = "count part: log(size_factor) as an offset; " *
                    "zero part: log(size_factor) as a covariate"

    provenance = Provenance(kind, count_formula, zero_formula,
        "count ~ offset(log(size_factor))", "~ log(size_factor)", offset_method,
        versions.r_version, versions.pscl_version, versions.mass_version,
        n_tested, n_refused, n_failed, all_warnings)

    return ZIFit(kind, inputs.taxa, fits, provenance)
end
