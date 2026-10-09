# SPDX-License-Identifier: MPL-2.0
# Copyright (c) 2026 Jonathan D.A. Jewell (hyperpolymath) <j.d.a.jewell@open.ac.uk>
#
# The R side.
#
# Production fits are `pscl` fits, reached through RCall. The R functions below
# are the whole of the interface: they take plain vectors, return plain scalars
# and vectors, and never return anything a caller has to interpret. In
# particular they do NOT compute a test statistic — the likelihood-ratio
# statistic, its degrees of freedom and the p-value are computed on the Julia
# side, so that the same arithmetic is used for both models and so that the
# test in test/ can recompute them independently.
#
# `pscl::hurdle` is called with `dist = "negbin", zero.dist = "binomial"` and
# `pscl::zeroinfl` with `dist = "negbin"`, exactly as docs/method-conditions/
# specifies. Both take the two-part formula
#
#     count ~ group + offset(log(size_factor)) | group + log(size_factor)
#
# and the null model drops the group term from BOTH parts, as specified.

const _R_HELPERS = raw"""
zic_fit <- function(y, g, s, kind, levels) {
    # `levels` is [reference, contrast], passed from Julia. Without it R would
    # sort the levels itself and the "g" coefficient would silently change sign
    # whenever the caller's reference group is not the alphabetically first one.
    d <- data.frame(y = y, g = factor(g, levels = levels), s = s)
    w <- character(0)
    run <- function(expr) withCallingHandlers(
        tryCatch(expr, error = function(e) e),
        warning = function(cond) {
            w <<- c(w, conditionMessage(cond))
            invokeRestart("muffleWarning")
        })
    full <- if (kind == "hurdle") {
        run(pscl::hurdle(y ~ g + offset(log(s)) | g + log(s), data = d,
                         dist = "negbin", zero.dist = "binomial"))
    } else {
        run(pscl::zeroinfl(y ~ g + offset(log(s)) | g + log(s), data = d,
                           dist = "negbin"))
    }
    if (inherits(full, "error")) {
        return(list(ok = FALSE, message = conditionMessage(full),
                    warnings = w, stage = "full"))
    }
    null <- if (kind == "hurdle") {
        run(pscl::hurdle(y ~ offset(log(s)) | log(s), data = d,
                         dist = "negbin", zero.dist = "binomial"))
    } else {
        run(pscl::zeroinfl(y ~ offset(log(s)) | log(s), data = d,
                           dist = "negbin"))
    }
    if (inherits(null, "error")) {
        return(list(ok = FALSE, message = conditionMessage(null),
                    warnings = w, stage = "null"))
    }
    sumc <- summary(full)$coefficients
    cn <- names(coef(full, model = "count"))
    zn <- names(coef(full, model = "zero"))
    cig <- which(startsWith(cn, "g"))
    zig <- which(startsWith(zn, "g"))
    if (length(cig) != 1 || length(zig) != 1) {
        return(list(ok = FALSE,
                    message = paste0("pscl returned ", length(cig),
                                     " count-part and ", length(zig),
                                     " zero-part group coefficients"),
                    warnings = w, stage = "extract"))
    }
    Xz <- model.matrix(~ g + log(s), data = d)
    # The zero part's fitted probability, in pscl's own parameterisation:
    # for `hurdle` this is P(count > 0) (the spec's "whether any reads were
    # observed"), for `zeroinfl` it is P(structural zero) (the mixture weight).
    zero_part_prob <- as.numeric(plogis(Xz %*% coef(full, model = "zero")))
    list(ok = TRUE,
         count_coef = as.numeric(coef(full, model = "count")[cig]),
         count_se = as.numeric(sumc$count[cig, 2]),
         zero_coef = as.numeric(coef(full, model = "zero")[zig]),
         zero_se = as.numeric(sumc$zero[zig, 2]),
         theta = as.numeric(full$theta),
         loglik_full = as.numeric(logLik(full)),
         loglik_null = as.numeric(logLik(null)),
         zero_part_prob = zero_part_prob,
         converged = (if (is.null(full$converged)) TRUE else isTRUE(full$converged)) &&
                     (if (is.null(null$converged)) TRUE else isTRUE(null$converged)),
         warnings = w,
         stage = "")
}

zic_versions <- function() {
    list(r_version = R.version.string,
         pscl_version = if (requireNamespace("pscl", quietly = TRUE))
             as.character(utils::packageVersion("pscl")) else NA_character_,
         mass_version = if (requireNamespace("MASS", quietly = TRUE))
             as.character(utils::packageVersion("MASS")) else NA_character_)
}

zic_has_pscl <- function() requireNamespace("pscl", quietly = TRUE)
"""

# The helpers are (re)defined once per session, when the package is loaded. The
# call is cheap and idempotent; doing it here rather than at top level keeps
# `Pkg.precompile` from needing R at build time... it still needs R at load
# time, because RCall is a hard dependency: a machine without R cannot load
# this package, which is the honest outcome for a package whose production path
# IS R.
function __init__()
    R"eval(parse(text = $(_R_HELPERS)))"
    return nothing
end

"""
    pscl_versions() -> NamedTuple

The R, `pscl` and `MASS` versions of the R this session is talking to. Called
once per table fit, and recorded in the [`Provenance`](@ref).
"""
function pscl_versions()
    raw = rcopy(R"zic_versions()")
    return (
        r_version = _as_string(_rfield(raw, "r_version")),
        pscl_version = _as_string(_rfield(raw, "pscl_version")),
        mass_version = _as_string(_rfield(raw, "mass_version")),
    )
end

"""
    pscl_available() -> Bool

Whether R's `pscl` package is installed and loadable. The method-conditions
document requires this to be reported with the package named, not guessed at.
"""
pscl_available() = rcopy(R"zic_has_pscl()")::Bool

# --- reading an R named list --------------------------------------------------
# RCall turns an R named list into a Dict whose values are vectors or scalars,
# depending on length. Everything below is tolerant of the spellings RCall has
# used for that conversion, and of `missing` from an R `NA`: a field that is
# absent or unusable becomes `nothing`, `NaN` or "", and every caller checks
# before it acts. It never becomes a plausible-looking number.

function _as_dict(res)
    if res isa AbstractDict
        return Dict{String,Any}(string(k) => v for (k, v) in res)
    elseif res isa NamedTuple
        return Dict{String,Any}(string(k) => v for (k, v) in pairs(res))
    elseif res isa AbstractVector
        out = Dict{String,Any}()
        for entry in res
            if entry isa Pair
                out[string(entry.first)] = entry.second
            end
        end
        return out
    end
    return Dict{String,Any}()
end

function _rfield(res, name::String)
    d = _as_dict(res)
    return get(d, name, nothing)
end

function _as_number(x)
    if x === nothing || x isa Missing
        return NaN
    elseif x isa Number
        return Float64(x)
    elseif x isa AbstractVector && length(x) == 1
        return _as_number(x[1])
    end
    return NaN
end

function _as_float(res, name::String)
    return _as_number(_rfield(res, name))
end

function _as_bool(res, name::String, default::Bool)
    v = _rfield(res, name)
    if v isa AbstractVector && length(v) == 1
        v = v[1]
    end
    return v isa Bool ? v : default
end

function _as_floats(res, name::String)
    v = _rfield(res, name)
    v === nothing && return Float64[]
    v isa Missing && return Float64[]
    v isa Number && return [Float64(v)]
    if v isa AbstractVector
        return [x isa Number ? Float64(x) : NaN for x in v]
    end
    return Float64[]
end

function _as_strings(res, name::String)
    v = _rfield(res, name)
    v === nothing && return String[]
    v isa Missing && return String[]
    v isa AbstractString && return [String(v)]
    if v isa AbstractVector
        return [_as_string(x) for x in v]
    end
    return String[]
end

function _as_string(x)
    if x === nothing || x isa Missing
        return ""
    elseif x isa AbstractString
        return String(x)
    elseif x isa AbstractVector
        return join(string.(x), ", ")
    end
    return string(x)
end
