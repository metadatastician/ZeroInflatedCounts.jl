# SPDX-License-Identifier: MPL-2.0
# RCall's scalar/vector and named-list representations share one contract.

@testset "R named-list conversion" begin
    zic = ZeroInflatedCounts
    for raw in (Dict("theta" => 2), Dict(:theta => 2), (theta = 2,), [:theta => 2])
        @test zic._rfield(raw, "theta") == 2
        @test zic._as_float(raw, "theta") === 2.0
        @test zic._rfield(raw, "absent") === nothing
    end
    @test zic._as_dict(Any[:theta=>2, "ignored"]) == Dict("theta" => 2)
    @test isempty(zic._as_dict(nothing))
    @test zic._as_float((theta = [2.5],), "theta") == 2.5
    for unusable in (nothing, missing, "2.5", [], [1, 2], [missing])
        @test isnan(zic._as_float((theta = unusable,), "theta"))
    end
    for value in (true, false), spelling in (value, [value])
        @test zic._as_bool((ok = spelling,), "ok", !value) === value
    end
    for default in (true, false),
        unusable in (nothing, missing, 1, "true", [], [true, false])

        @test zic._as_bool((ok = unusable,), "ok", default) === default
        @test zic._as_bool((;), "absent", default) === default
    end

    @test zic._as_floats((p = 0.25,), "p") == [0.25]
    @test zic._as_floats((p = [0, 1],), "p") == [0.0, 1.0]
    @test isequal(zic._as_floats((p = Any[0.25, missing, "bad"],), "p"), [0.25, NaN, NaN])
    @test isempty(zic._as_floats((;), "absent"))
    @test isempty(zic._as_strings((;), "absent"))
    for unusable in (nothing, missing, "bad")
        @test isempty(zic._as_floats((p = unusable,), "p"))
    end
    @test zic._as_strings((warnings = "one warning",), "warnings") == ["one warning"]
    @test zic._as_strings((warnings = ["first", "second"],), "warnings") ==
          ["first", "second"]
    @test zic._as_strings((warnings = missing,), "warnings") == String[]
    @test zic._as_string(missing) == ""
    @test zic._as_string(nothing) == ""
    @test zic._as_string(["first", "second"]) == "first, second"
end
