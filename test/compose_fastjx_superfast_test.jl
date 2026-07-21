@testitem "2wayCoupling" begin
    using GasChem, EarthSciMLBase
    using OrdinaryDiffEqRosenbrock
    using ModelingToolkit
    sol_middle = 39.999599729358856
    sf = couple(SuperFast(), FastJX(0.0))
    sys = convert(System, sf)
    tspan = (0.0, 3600 * 24)
    prob = ODEProblem(sys, [], tspan; build_initializeprob = false)
    sol = solve(prob, Rosenbrock23(), saveat = 10.0)
    @test sol[sys.SuperFast₊O3][4320] ≈ sol_middle rtol = 1.0e-4
end

@testitem "2wayCoupling interpolation" begin
    using GasChem, EarthSciMLBase
    using OrdinaryDiffEqRosenbrock
    using ModelingToolkit
    sol_middle = 39.999599729358856

    sf = couple(SuperFast(), FastJX_interpolation_troposphere(0.0))
    sys = convert(System, sf)
    tspan = (0.0, 3600 * 24)
    prob = ODEProblem(sys, [], tspan; build_initializeprob = false)
    sol = solve(prob, Rosenbrock23(), saveat = 10.0)
    @test sol[sys.SuperFast₊O3][4320] ≈ sol_middle rtol = 1.0e-4
end

# The two composition goldens above both sample index 4320 of the 24 h solve, which at the
# operators' default coordinates (lat 40N, long -97W, t_ref = 0) is cosSZA = -0.35 — the middle
# of the night, where every J is identically zero. They therefore cannot see a change in the
# actinic-flux table at all: the bin-17/18 correction moved surface J_NO2 by 2.77x and left
# `sol_middle` unchanged to 1.7e-12. Nor is O3 a usable probe here — over the whole 24 h it
# only ranges 39.9996 to 40.0, so day and night differ by 4e-6 relative, well inside the 1e-4
# tolerance. Assert the J values themselves, at an index that is actually in daylight.
@testitem "FastJX_interpolation_troposphere daytime J values" begin
    using GasChem, EarthSciMLBase
    using OrdinaryDiffEqRosenbrock
    using ModelingToolkit

    night, day = 4320, 6660          # t = 43190 s (cosSZA -0.35) and 66590 s (cosSZA +0.45)

    sf = couple(SuperFast(), FastJX_interpolation_troposphere(0.0))
    sys = convert(System, sf)
    prob = ODEProblem(sys, [], (0.0, 3600 * 24); build_initializeprob = false)
    sol = solve(prob, Rosenbrock23(), saveat = 10.0)

    getj(name) = sol[only(filter(eq -> string(eq.lhs) == name, observed(sys))).lhs]

    j_NO2 = getj("FastJX₊j_NO2(t)")
    j_H2O2 = getj("FastJX₊j_H2O2(t)")

    # Night: the table's geometric zero must survive any regeneration of the flux tables.
    @test j_NO2[night] == 0.0
    @test j_H2O2[night] == 0.0

    # Day. j_NO2 is dominated by bin 17, so it guards the TOA-flux/bin-definition pairing;
    # j_H2O2 lives in the UV bins, so it guards the bins 1-16 half of the table independently.
    @test j_NO2[day]≈0.004144177136858269 rtol=1.0e-4
    @test j_H2O2[day]≈1.1366858279977815e-6 rtol=1.0e-4
end
