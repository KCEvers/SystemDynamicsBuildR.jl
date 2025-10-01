using Test
using Unitful
using SystemDynamicsBuildR.sdbuildR_units

@testset "sdbuildR_units tests" begin

    @testset "Time units - Common year based" begin
        @testset "common_yr" begin
            @test 1u"common_yr" == 365u"d"
            @test 2u"common_yr" == 730u"d"
            @test uconvert(u"d", 1u"common_yr") == 365u"d"
        end

        @testset "common_quarter" begin
            @test 4u"common_quarter" ≈ 365u"d"
            @test uconvert(u"d", 1u"common_quarter") ≈ 91.25u"d"
        end

        @testset "common_month" begin
            @test 12u"common_month" ≈ 365u"d"
            @test uconvert(u"d", 1u"common_month") ≈ 30.416666666666668u"d" atol=1e-10u"d"
        end
    end

    @testset "Time units - Julian year based" begin
        @testset "quarter" begin
            @test 4u"quarter" == 1u"yr"
            @test uconvert(u"d", 1u"quarter") ≈ 91.3125u"d"
        end

        @testset "month" begin
            @test 12u"month" == 1u"yr"
            @test uconvert(u"d", 1u"month") ≈ 30.4375u"d"
        end

        @testset "Difference between common and Julian" begin
            # Common year (365 days) vs Julian year (365.25 days)
            @test uconvert(u"d", 1u"common_yr") < uconvert(u"d", 1u"yr")
            @test uconvert(u"d", 1u"common_quarter") < uconvert(u"d", 1u"quarter")
            @test uconvert(u"d", 1u"common_month") < uconvert(u"d", 1u"month")
        end
    end

    @testset "Volume units" begin
        @testset "quart" begin
            @test 1u"quart" == 946.35u"cm^3"
            @test uconvert(u"mL", 1u"quart") ≈ 946.35u"mL"
        end

        @testset "US_gal" begin
            @test 1u"US_gal" == 0.003785411784u"m^3"
            @test uconvert(u"L", 1u"US_gal") ≈ 3.785411784u"L"
            
            # 1 gallon = 4 quarts
            @test uconvert(u"quart", 1u"US_gal") ≈ 4u"quart" atol=0.001u"quart"
        end

        @testset "fluidOunce" begin
            @test 1u"fluidOunce" ≈ 29.5735295625u"mL"
            
            # 128 fluid ounces = 1 gallon
            @test uconvert(u"US_gal", 128u"fluidOunce") ≈ 1u"US_gal" atol=0.001u"US_gal"
        end
    end

    @testset "Mass units" begin
        @testset "tonne (metric ton)" begin
            @test 1u"tonne" == 1000u"kg"
            @test uconvert(u"kg", 2.5u"tonne") == 2500u"kg"
        end

        @testset "ton (US short ton)" begin
            @test 1u"ton" == 907.18474u"kg"
            @test uconvert(u"kg", 1u"ton") ≈ 907.18474u"kg"
        end

        @testset "Difference between ton and tonne" begin
            # Metric tonne is heavier than US ton
            @test uconvert(u"kg", 1u"tonne") > uconvert(u"kg", 1u"ton")
            @test uconvert(u"ton", 1u"tonne") ≈ 1.10231u"ton" atol=0.001u"ton"
        end
    end

    @testset "Amount units (molecular)" begin
        @testset "atom" begin
            avogadro = 6.02214076e23
            @test 1u"atom" ≈ 1/avogadro * u"mol"
            @test avogadro * u"atom" ≈ 1u"mol" rtol=1e-10
        end

        @testset "molecule" begin
            avogadro = 6.02214076e23
            @test 1u"molecule" ≈ 1/avogadro * u"mol"
            @test avogadro * u"molecule" ≈ 1u"mol" rtol=1e-10
        end

        @testset "atom vs molecule" begin
            # They should be equivalent
            @test 1u"atom" == 1u"molecule"
        end

        @testset "Conversions to moles" begin
            @test uconvert(u"mol", 1e23u"atom") ≈ 0.166u"mol" atol=0.001u"mol"
            @test uconvert(u"atom", 1u"mol") ≈ 6.02214076e23u"atom" rtol=1e-10
        end
    end

    @testset "Currency units" begin
        @testset "EUR" begin
            @test 100u"EUR" isa Quantity
            @test dimension(1u"EUR") == Unitful.NoDims
        end

        @testset "USD" begin
            @test 50u"USD" isa Quantity
            @test dimension(1u"USD") == Unitful.NoDims
        end

        @testset "GBP" begin
            @test 75u"GBP" isa Quantity
            @test dimension(1u"GBP") == Unitful.NoDims
        end

        @testset "Currency arithmetic" begin
            # Can add same currencies
            @test 100u"USD" + 50u"USD" == 150u"USD"
            @test 200u"EUR" - 50u"EUR" == 150u"EUR"
            
            # Can multiply by scalars
            @test 2 * 50u"GBP" == 100u"GBP"
        end
    end

    @testset "Angular units" begin
        @testset "deg_ (degree)" begin
            @test 1u"deg_" ≈ π/180 atol=1e-10
            @test 180u"deg_" ≈ π rtol=1e-10
            @test 360u"deg_" ≈ 2π rtol=1e-10
        end

        @testset "Degree conversions" begin
            @test uconvert(u"rad", 90u"deg_") ≈ π/2 * u"rad" rtol=1e-10
            @test uconvert(u"deg_", π*u"rad") ≈ 180u"deg_" rtol=1e-10
        end
    end

    @testset "Electromagnetic units" begin
        @testset "ohm_" begin
            @test 1u"ohm_" == 1u"V"/u"A"
            @test dimension(1u"ohm_") == dimension(1u"Ω")
        end

        @testset "Ohm's law" begin
            V = 10u"V"
            I = 2u"A"
            R = V / I
            
            @test uconvert(u"ohm_", R) ≈ 5u"ohm_"
        end
    end

    @testset "Physical constants as units" begin
        @testset "reduced_Planck_constant" begin
            @test 1u"reduced_Planck_constant" ≈ u"h"/(2π)
            @test dimension(1u"reduced_Planck_constant") == dimension(1u"h")
        end

        @testset "anghertz" begin
            @test 1u"anghertz" ≈ 2π/u"s"
            @test dimension(1u"anghertz") == dimension(1u"Hz")
        end

        @testset "Physical constant units are valid" begin
            # Just check they can be created without error
            @test 1u"superconducting_magnetic_flux_quantum" isa Quantity
            @test 1u"Stefan_Boltzmann_constant" isa Quantity
            @test 1u"Bohr_magneton" isa Quantity
            @test 1u"Rydberg_constant" isa Quantity
            @test 1u"magnetic_constant" isa Quantity
            @test 1u"electric_constant" isa Quantity
        end
    end

    @testset "Temperature units" begin
        @testset "degC (Celsius)" begin
            @test 1u"degC" isa Quantity
            # Celsius is affine, can't directly compare equality
            @test dimension(1u"degC") == dimension(1u"K")
        end

        @testset "degF (Fahrenheit)" begin
            @test 1u"degF" isa Quantity
            @test dimension(1u"degF") == dimension(1u"K")
        end
    end

    @testset "Unit arithmetic and conversions" begin
        @testset "Time period calculations" begin
            # Project duration
            duration = 2.5u"common_yr"
            @test uconvert(u"d", duration) == 912.5u"d"
            
            # Quarterly reporting
            quarters = 4u"quarter"
            @test uconvert(u"yr", quarters) == 1u"yr"
        end

        @testset "Volume conversions" begin
            # Recipe scaling
            @test uconvert(u"mL", 2u"fluidOunce") ≈ 59.147u"mL" atol=0.01u"mL"
            
            # Fuel tank
            @test uconvert(u"L", 15u"US_gal") ≈ 56.781u"L" atol=0.01u"L"
        end

        @testset "Mass conversions" begin
            # Shipping weight
            @test uconvert(u"ton", 5u"tonne") ≈ 5.512u"ton" atol=0.001u"ton"
        end

        @testset "Molecular calculations" begin
            # Number of molecules in a sample
            sample = 0.5u"mol"
            n_molecules = uconvert(u"molecule", sample)
            @test n_molecules ≈ 3.01107e23u"molecule" rtol=1e-4
        end

        @testset "Mixed unit operations" begin
            # Rate calculation
            distance = 100u"km"
            time = 2u"quarter"
            speed = distance / time
            
            @test dimension(speed) == dimension(1u"m/s")
        end
    end

    @testset "Edge cases and special scenarios" begin
        @testset "Zero values" begin
            @test 0u"common_yr" == 0u"d"
            @test 0u"USD" isa Quantity
            @test 0u"tonne" == 0u"kg"
        end

        @testset "Fractional units" begin
            @test 0.5u"common_yr" == 182.5u"d"
            @test 0.25u"US_gal" ≈ 1u"quart" atol=0.01u"quart"
        end

        @testset "Large values" begin
            @test 1000u"common_yr" == 365000u"d"
            @test 1e6u"atom" isa Quantity
        end

        @testset "Unit combination" begin
            # Area calculation
            length = 10u"m"
            width = 5u"m"
            area = length * width
            @test area == 50u"m^2"
            
            # Density
            mass = 1u"tonne"
            volume = 1u"m^3"
            density = mass / volume
            @test uconvert(u"kg/m^3", density) == 1000u"kg/m^3"
        end

        @testset "Inverse units" begin
            rate = 0.05u"1/common_yr"
            @test rate isa Quantity
            @test dimension(rate) == dimension(1u"Hz")
            
            # Interest rate application
            principal = 1000u"USD"
            time = 1u"common_yr"
            interest = principal * rate * time
            @test ustrip(interest) ≈ 50.0
        end
    end

    @testset "Real-world application examples" begin
        @testset "Financial modeling" begin
            # Quarterly interest
            annual_rate = 0.04u"1/common_yr"
            quarterly_rate = annual_rate / 4
            principal = 10000u"USD"
            quarters = 8u"common_quarter"
            
            # Simple interest
            interest = principal * quarterly_rate * quarters
            @test ustrip(interest) ≈ 800.0
        end

        @testset "Chemical calculations" begin
            # Moles to molecules
            amount = 2u"mol"
            molecules = uconvert(u"molecule", amount)
            @test molecules ≈ 1.204428152e24u"molecule" rtol=1e-4
        end

        @testset "Population dynamics" begin
            # Birth rate per common year
            birth_rate = 0.02u"1/common_yr"
            population = 1000
            time_span = 5u"common_yr"
            
            # Simple growth (births per year)
            births_per_year = population * birth_rate
            total_births = ustrip(births_per_year) * ustrip(uconvert(u"common_yr", time_span))
            @test total_births ≈ 100.0
        end

        @testset "Recipe conversion" begin
            # Convert recipe from US to metric
            butter = 8u"fluidOunce"
            butter_ml = uconvert(u"mL", butter)
            @test butter_ml ≈ 236.588u"mL" atol=0.01u"mL"
        end
    end

end