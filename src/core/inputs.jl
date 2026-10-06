# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# The data contract, checked before anything is fitted.
#
# docs/method-conditions/ zero-inflation-and-hurdle.adoc, "Data accepted":
# non-negative integer read counts, no pseudocounts, two groups, one size
# factor per sample. Anything else is refused, not rounded.

"""
    ZIInputs

Validated inputs: the count table transposed to samples × taxa with taxa in
columns, the group labels as strings with `ref` and `contrast` fixed, and the
size factors. Constructed by [`check_inputs`](@ref); not exported.
"""
struct ZIInputs
    counts::Matrix{Int}
    groups::Vector{String}
    size_factors::Vector{Float64}
    taxa::Vector{String}
    ref::String
    contrast::String
    min_prevalence::Float64
end

"""
    check_inputs(counts, groups, size_factors; taxa=nothing, ref=nothing,
                 min_prevalence=0.0) -> ZIInputs

Validate and normalise a call to [`hurdle_nb`](@ref) or [`zinb`](@ref).

Refused here (an `ArgumentError`, naming what is wrong) rather than rounded or
silently repaired:

* counts that are not non-negative integers;
* a row count that does not match `groups` and `size_factors`;
* anything other than exactly two groups;
* size factors that are not finite and strictly positive;
* a `min_prevalence` outside `[0, 1]`.

`ref` chooses the reference group. Left at `nothing`, the reference is the
first group in sorted order, so the same table always gives the same contrast.
"""
function check_inputs(
    counts::AbstractMatrix,
    groups::AbstractVector,
    size_factors::AbstractVector;
    taxa = nothing,
    ref = nothing,
    min_prevalence::Real = 0.0,
)
    n, p = size(counts)

    if !(eltype(counts) <: Integer)
        throw(
            ArgumentError(
                "counts must be an integer matrix (got eltype " *
                "$(eltype(counts))); non-integer counts are refused, not rounded",
            ),
        )
    end
    if any(counts .< 0)
        throw(ArgumentError("counts must be non-negative; negative entries are refused"))
    end
    if n != length(groups) || n != length(size_factors)
        throw(
            ArgumentError(
                "counts has $n rows but groups has $(length(groups)) " *
                "and size_factors has $(length(size_factors)) entries",
            ),
        )
    end
    if p == 0 || n == 0
        throw(ArgumentError("counts must have at least one sample and one taxon"))
    end

    if length(size_factors) > 0 && any(!isfinite(s) for s in size_factors)
        throw(ArgumentError("size factors must be finite"))
    end
    if any(s -> s <= 0, size_factors)
        throw(ArgumentError("size factors must be strictly positive"))
    end

    levels = sort(unique(string.(groups)))
    if length(levels) != 2
        throw(
            ArgumentError(
                "exactly two groups are required (got " *
                "$(length(levels)): $(join(levels, ", ")))",
            ),
        )
    end
    if ref === nothing
        reference = levels[1]
        contrast = levels[2]
    else
        reference = string(ref)
        reference in levels || throw(
            ArgumentError(
                "ref = $(reference) is not one of the groups " * "($(join(levels, ", ")))",
            ),
        )
        contrast = levels[levels .!= reference][1]
    end

    if !(0 <= min_prevalence <= 1)
        throw(ArgumentError("min_prevalence must lie in [0, 1] (got $min_prevalence)"))
    end

    if taxa === nothing
        taxon_names = ["taxon-" * string(j) for j = 1:p]
    else
        length(taxa) == p || throw(
            ArgumentError("taxa has $(length(taxa)) entries but counts has $p columns"),
        )
        taxon_names = string.(taxa)
    end

    return ZIInputs(
        Matrix{Int}(counts),
        string.(groups),
        Float64.(size_factors),
        taxon_names,
        reference,
        contrast,
        Float64(min_prevalence),
    )
end
