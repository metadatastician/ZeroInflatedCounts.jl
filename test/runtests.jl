# ZeroInflatedCounts — test entry point (julia-library archetype)
#
# Load order matters:
#   1. Aqua package-shape gate — runs FIRST: a package that does not load
#      cleanly, whose deps are not compat-closed, or whose exports are
#      ambiguous is rejected before any behaviour test can pass.
#   2. The oracle and the test helpers — both are included into this scope
#      (Main) rather than imported, so a direct `julia test/runtests.jl` and
#      `Pkg.test()` load exactly the same thing.
#   3. Behaviour — every test/cases/*.jl is included, one file per concern.

using Test
using ZeroInflatedCounts

@testset "ZeroInflatedCounts — Aqua package shape" begin
    using Aqua
    Aqua.test_all(ZeroInflatedCounts)
end

# The independent likelihood oracle (pure Julia, no R) and the simulation
# helpers. `Likelihoods` is a module defined in this scope, which is why the
# helpers can name it directly.
include(joinpath(@__DIR__, "oracle", "Likelihoods.jl"))
include(joinpath(@__DIR__, "helpers", "simulation.jl"))

@testset "ZeroInflatedCounts — behaviour" begin
    # NOTE: no `continue`/`break` inside a `for` inside @testset — Test's
    # loop-form detection re-parses the loop body at top level, where
    # `continue` is a ParseError. Filter with a guarded call instead.
    # (Measured trap, 2026-09-19 — see the Julia testing guide.)
    cases = joinpath(@__DIR__, "cases")
    if isdir(cases)
        for path in sort(readdir(cases, join = true))
            endswith(path, ".jl") && include(path)
        end
    end
end
