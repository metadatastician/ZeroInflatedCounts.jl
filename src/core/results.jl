# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# Result types for the two models.
#
# One taxon is tested at a time (see docs/method-conditions/). A taxon that
# cannot be fitted is never given a stand-in p-value: it is reported as
# `TaxonRefused` (a condition visible in the counts) or `TaxonFailed` (the
# fit was attempted and R reported a failure).

"""
    TaxonFit

A single taxon that was fitted and tested.

Fields are the descriptive per-part estimates plus the likelihood-ratio test of
the group effect in both parts together:

* `log_fold_change`, `se_log_fold_change` — count part, log scale, contrast vs
  reference (the `g` coefficient of the negative-binomial part).
* `zero_log_odds`, `se_zero_log_odds` — the zero part's `g` coefficient, contrast
  vs reference, on the logit scale. In pscl's parameterisation the zero part
  models P(count > 0) for `hurdle_nb` and P(structural zero) for `zinb`, so the
  same field reads as the log-odds that *any* reads were observed (hurdle) or of
  a structural zero (ZINB). That convention is measured against pscl in
  `test/cases/40-pscl-cross-check.jl`, not assumed here.
* `theta` — fitted negative-binomial size (dispersion) parameter.
* `loglik_full`, `loglik_null`, `lr_statistic`, `df`, `pvalue` — the test.
* `p_adjusted` — Benjamini–Hochberg over the taxa that produced a p-value.
* `warnings` — every warning R raised for this fit, verbatim.
* `notes` — interpretation notes (for example a θ at a bound).
"""
struct TaxonFit
    taxon::String
    log_fold_change::Float64
    se_log_fold_change::Float64
    zero_log_odds::Float64
    se_zero_log_odds::Float64
    theta::Float64
    loglik_full::Float64
    loglik_null::Float64
    lr_statistic::Float64
    df::Int
    pvalue::Float64
    p_adjusted::Float64
    theta_note::Union{Nothing,String}
    warnings::Vector{String}
end

"""
    TaxonRefused

A taxon refused before fitting, with the reason and no p-value.
"""
struct TaxonRefused
    taxon::String
    reason::String
end

"""
    TaxonFailed

A taxon whose fit was attempted and failed (non-convergence, a covariance
matrix that is not positive definite, or any other error from `pscl`), with
R's message verbatim and no p-value.
"""
struct TaxonFailed
    taxon::String
    message::String
end

const TaxonResult = Union{TaxonFit,TaxonRefused,TaxonFailed}

"""
    Provenance

What was run, with what, recorded for every call (docs/method-conditions/:
"Provenance written for every run").
"""
struct Provenance
    method::Symbol
    count_formula::String
    zero_formula::String
    null_count_formula::String
    null_zero_formula::String
    offset_method::String
    r_version::String
    pscl_version::String
    mass_version::String
    n_tested::Int
    n_refused::Int
    n_failed::Int
    warnings::Vector{String}
end

"""
    ZIFit

The result of a whole-table call to [`hurdle_nb`](@ref) or [`zinb`](@ref).

`fits` holds one entry per taxon, in the column order of `counts`, each of
which is a [`TaxonFit`](@ref), a [`TaxonRefused`](@ref) or a
[`TaxonFailed`](@ref). `p_adjusted` values in the `TaxonFit` entries are
Benjamini–Hochberg across exactly the taxa that produced a p-value; refused and
failed taxa are outside that family.
"""
struct ZIFit
    method::Symbol
    taxa::Vector{String}
    fits::Vector{TaxonResult}
    provenance::Provenance
end

"""
    fitted(res::ZIFit) -> Vector{TaxonFit}

The taxa that were fitted and tested, in input order.
"""
fitted(res::ZIFit) = TaxonFit[r for r in res.fits if r isa TaxonFit]

"""
    refused(res::ZIFit) -> Vector{TaxonRefused}

The taxa refused before fitting, with their reasons.
"""
refused(res::ZIFit) = TaxonRefused[r for r in res.fits if r isa TaxonRefused]

"""
    failed(res::ZIFit) -> Vector{TaxonFailed}

The taxa whose fit failed, with R's message.
"""
failed(res::ZIFit) = TaxonFailed[r for r in res.fits if r isa TaxonFailed]

"""
    pvalues(res::ZIFit) -> Vector{Float64}

The unadjusted likelihood-ratio p-values of the fitted taxa, in input order.
"""
pvalues(res::ZIFit) = [r.pvalue for r in fitted(res)]

"""
    adjusted_pvalues(res::ZIFit) -> Vector{Float64}

The Benjamini–Hochberg adjusted p-values of the fitted taxa, in input order.
"""
adjusted_pvalues(res::ZIFit) = [r.p_adjusted for r in fitted(res)]

function Base.show(io::IO, res::ZIFit)
    p = res.provenance
    print(io, "ZIFit(", res.method, ": ", p.n_tested, " tested, ", p.n_refused,
        " refused, ", p.n_failed, " failed)")
end

function Base.show(io::IO, fit::TaxonFit)
    print(io, "TaxonFit(", fit.taxon, ": LFC=", round(fit.log_fold_change; digits = 4),
        " (SE ", round(fit.se_log_fold_change; digits = 4), "), zero-log-odds=",
        round(fit.zero_log_odds; digits = 4), ", p=", fit.pvalue,
        ", p_adj=", fit.p_adjusted, ")")
end

Base.show(io::IO, r::TaxonRefused) = print(io, "TaxonRefused(", r.taxon, ": ", r.reason, ")")
Base.show(io::IO, r::TaxonFailed) = print(io, "TaxonFailed(", r.taxon, ": ", r.message, ")")

"""
    set_adjusted_p(fit::TaxonFit, p::Float64) -> TaxonFit

Return `fit` with `p_adjusted` replaced. Used once, after the whole table has
been fitted and Benjamini–Hochberg has been applied across the family.
"""
function set_adjusted_p(fit::TaxonFit, p::Float64)
    return TaxonFit(fit.taxon, fit.log_fold_change, fit.se_log_fold_change, fit.zero_log_odds,
        fit.se_zero_log_odds, fit.theta, fit.loglik_full, fit.loglik_null,
        fit.lr_statistic, fit.df, fit.pvalue, p, fit.theta_note, fit.warnings)
end

"""
    summary_table(res::ZIFit) -> String

A plain-text table, one row per taxon, for logs and for a caller that has no
table renderer of its own.
"""
function summary_table(res::ZIFit)
    out = IOBuffer()
    println(out, "method: ", res.method)
    println(out, res.method === :hurdle_nb ?
                "zero part reads as P(count > 0) on the logit scale (pscl's hurdle convention)" :
                "zero part reads as P(structural zero) on the logit scale (pscl's zeroinfl convention)")
    println(out, rpad("taxon", 24), rpad("LFC", 12), rpad("SE", 12), rpad("zero-log-odds", 16),
        rpad("theta", 12), rpad("p", 12), rpad("p_adj", 12), "state")
    for r in res.fits
        if r isa TaxonFit
            println(out, rpad(r.taxon, 24), rpad(string(round(r.log_fold_change; digits = 4)), 12),
                rpad(string(round(r.se_log_fold_change; digits = 4)), 12),
                rpad(string(round(r.zero_log_odds; digits = 4)), 16),
                rpad(string(round(r.theta; digits = 4)), 12),
                rpad(string(round(r.pvalue; digits = 6)), 12),
                rpad(string(round(r.p_adjusted; digits = 6)), 12), "tested")
        elseif r isa TaxonRefused
            println(out, rpad(r.taxon, 24), rpad("-", 12), rpad("-", 12), rpad("-", 16),
                rpad("-", 12), rpad("-", 12), rpad("-", 12), "refused: ", r.reason)
        else
            println(out, rpad(r.taxon, 24), rpad("-", 12), rpad("-", 12), rpad("-", 16),
                rpad("-", 12), rpad("-", 12), rpad("-", 12), "failed: ", r.message)
        end
    end
    print(out, "R ", res.provenance.r_version, "; pscl ", res.provenance.pscl_version,
        "; MASS ", res.provenance.mass_version)
    return String(take!(out))
end
