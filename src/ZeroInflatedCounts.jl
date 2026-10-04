# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# ZeroInflatedCounts — hurdle and zero-inflated negative binomial models for
# sequencing count data.
#
# Production fits are R `pscl` fits reached through RCall. This module's job is
# to make them answerable: it validates the data contract, decides per taxon
# whether a fit is computable at all, runs the same two-part formula for both
# models, reports one likelihood-ratio test per taxon with Benjamini–Hochberg
# across taxa, and records the provenance of every run.
#
# The specification this code is checked against is
# docs/method-conditions/zero-inflation-and-hurdle.adoc. Read it before
# changing a refusal, a formula or the test: every branch here is a line in
# that document.

module ZeroInflatedCounts

using RCall

include("core/results.jl")
include("core/pscl.jl")
include("core/inputs.jl")
include("core/models.jl")

export hurdle_nb, zinb, fitted, refused, failed, pvalues, adjusted_pvalues,
    summary_table, TaxonFit, TaxonRefused, TaxonFailed, ZIFit, Provenance

"""
    hurdle_nb(counts, groups, size_factors; taxa=nothing, ref=nothing,
              min_prevalence=0.0) -> ZIFit

Hurdle negative binomial for one count table.

`counts` is `samples × taxa`, non-negative integers, no pseudocounts.
`groups` labels each sample with one of exactly two groups; `size_factors` gives
each sample its sequencing-depth factor. Both are used in both parts of the
model:

* count part — `count ~ group + offset(log(size_factor))`, a zero-truncated
  negative binomial fitted on the positive counts only;
* zero part — `~ group + log(size_factor)`, a binomial GLM on whether *any*
  reads were seen, where `log(size_factor)` is a covariate because an offset
  has no meaning on the logit scale.

Both parts are fitted by `pscl::hurdle(dist = "negbin", zero.dist = "binomial")`
through RCall. No model is chosen automatically and nothing here replaces the
customer's NB GLM: this is a second fit to show beside it.

One p-value per taxon, from a likelihood-ratio test of the group effect in both
parts together (2 degrees of freedom), with Benjamini–Hochberg across the taxa
that produced one. Per-part estimates are descriptive and get no p-value of
their own.

A taxon is returned as a [`TaxonRefused`](@ref) when a condition visible in the
counts makes the fit impossible or undefined (a group with no zeros or nothing
but zeros; fewer than 3 positive counts per group; a prevalence below
`min_prevalence`), and as a [`TaxonFailed`](@ref) when the fit was attempted
and `pscl` reported a failure. Neither carries a p-value.

# Examples

```julia
counts = [0 12 3; 4 0 7; 9 21 0; 0 3 5]
groups = ["control", "control", "treated", "treated"]
sizes = [1.0, 1.1, 0.9, 1.0]
res = hurdle_nb(counts, groups, sizes; taxa = ["a", "b", "c"])
summary_table(res)
```
"""
function hurdle_nb(counts::AbstractMatrix, groups::AbstractVector,
    size_factors::AbstractVector; taxa = nothing, ref = nothing,
    min_prevalence::Real = 0.0)
    return fit_table(:hurdle_nb, counts, groups, size_factors; taxa = taxa, ref = ref,
        min_prevalence = min_prevalence)
end

"""
    zinb(counts, groups, size_factors; taxa=nothing, ref=nothing,
         min_prevalence=0.0) -> ZIFit

Zero-inflated negative binomial for one count table.

Same data contract as [`hurdle_nb`](@ref), and the same formula:

* count part — `count ~ group + offset(log(size_factor))`, negative binomial;
* zero part — `~ group + log(size_factor)`, the mixture weight for structural
  zeros, fitted jointly with the count part by
  `pscl::zeroinfl(dist = "negbin")` through RCall.

The zero weight is a mixture weight, separated from the negative binomial's own
zeros only by the shape of that distribution; it is never reported as the
probability that a taxon is absent.

Fitted jointly with 6 parameters, so a table with 6 or fewer samples is
refused. A taxon whose fitted zero weight is below 1e-6 in every sample is
refused with "no zero inflation is detectable; the NB GLM result applies" — the
test is not reported rather than reported as a null result.

One p-value per taxon, from a likelihood-ratio test of the group effect in both
parts together (2 degrees of freedom), Benjamini–Hochberg across taxa.
"""
function zinb(counts::AbstractMatrix, groups::AbstractVector,
    size_factors::AbstractVector; taxa = nothing, ref = nothing,
    min_prevalence::Real = 0.0)
    return fit_table(:zinb, counts, groups, size_factors; taxa = taxa, ref = ref,
        min_prevalence = min_prevalence)
end

end # module ZeroInflatedCounts
