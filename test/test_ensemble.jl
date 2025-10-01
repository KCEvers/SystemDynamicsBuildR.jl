using SystemDynamicsBuildR.ensemble
using DataFrames
using Unitful
using Statistics

@testset "ensemble tests" begin

    @testset "Parameter combination generation" begin
        @testset "Crossed design - basic" begin
            param_ranges = Dict(
                :alpha => [0.1, 0.5, 1.0],
                :beta => [2.0, 5.0]
            )
            
            combinations, total = generate_param_combinations(
                param_ranges; crossed=true, n_replicates=10
            )
            
            # Should have 3 × 2 = 6 combinations
            @test length(combinations) == 6
            @test total == 60  # 6 × 10
            
            # Check that all combinations are unique
            @test length(unique(combinations)) == 6
        end

        @testset "Crossed design - three parameters" begin
            param_ranges = Dict(
                :a => [1, 2],
                :b => [10, 20],
                :c => [100, 200]
            )
            
            combinations, total = generate_param_combinations(
                param_ranges; crossed=true, n_replicates=5
            )
            
            @test length(combinations) == 8  # 2 × 2 × 2
            @test total == 40  # 8 × 5
        end

        @testset "Non-crossed design - paired" begin
            param_ranges = Dict(
                :alpha => [0.1, 0.5, 1.0],
                :beta => [2.0, 5.0, 10.0]
            )
            
            combinations, total = generate_param_combinations(
                param_ranges; crossed=false, n_replicates=10
            )
            
            # Should have 3 combinations (paired)
            @test length(combinations) == 3
            @test total == 30  # 3 × 10
            
            # Check pairing is correct (sorted keys: alpha, beta)
            @test combinations[1] == [0.1, 2.0]
            @test combinations[2] == [0.5, 5.0]
            @test combinations[3] == [1.0, 10.0]
        end

        @testset "Non-crossed design - error on mismatched lengths" begin
            param_ranges = Dict(
                :alpha => [0.1, 0.5],
                :beta => [2.0, 5.0, 10.0]  # Different length
            )
            
            @test_throws ArgumentError generate_param_combinations(
                param_ranges; crossed=false, n_replicates=10
            )
        end

        @testset "Single parameter" begin
            param_ranges = Dict(:alpha => [0.1, 0.5, 1.0])
            
            combinations, total = generate_param_combinations(
                param_ranges; crossed=true, n_replicates=5
            )
            
            @test length(combinations) == 3
            @test total == 15
        end

        @testset "Different replicate counts" begin
            param_ranges = Dict(:a => [1, 2], :b => [10, 20])
            
            _, total_10 = generate_param_combinations(
                param_ranges; crossed=true, n_replicates=10
            )
            _, total_100 = generate_param_combinations(
                param_ranges; crossed=true, n_replicates=100
            )
            
            @test total_10 == 40   # 4 × 10
            @test total_100 == 400 # 4 × 100
        end
    end

    @testset "Transform intermediaries" begin
        @testset "Basic transformation" begin
            # Create mock intermediate results
            intermediaries = [
                (t = [0.0, 1.0, 2.0], saveval = [10.0, 20.0, 30.0]),
                (t = [0.0, 1.0, 2.0], saveval = [15.0, 25.0, 35.0])
            ]
            
            transformed = transform_intermediaries(intermediaries, [:x])
            
            @test length(transformed) == 2
            @test transformed[1].t == [0.0, 1.0, 2.0]
            @test transformed[1].u == [10.0, 20.0, 30.0]
            @test isnothing(transformed[1].p)
        end

        @testset "Empty intermediaries" begin
            intermediaries = [
                (t = Float64[], saveval = Float64[])
            ]
            
            transformed = transform_intermediaries(intermediaries)
            
            @test length(transformed) == 1
            @test isempty(transformed[1].t)
            @test isempty(transformed[1].u)
        end

        @testset "Nothing intermediaries" begin
            intermediaries = [nothing, nothing]
            
            transformed = transform_intermediaries(intermediaries)
            
            @test length(transformed) == 2
            @test isempty(transformed[1].t)
            @test isempty(transformed[2].t)
        end
    end

    @testset "ensemble_to_df - basic functionality" begin
        # Create mock solution data
        @testset "Single variable, single trajectory" begin
            solve_out = [
                (
                    t = [0.0, 1.0, 2.0],
                    u = [10.0, 11.0, 12.0],
                    u0 = 10.0,
                    p = (alpha = 0.5, beta = 2.0)
                )
            ]
            
            ts_df, param_df, init_df = ensemble_to_df(
                solve_out, [:x], nothing, nothing, 1
            )
            
            # Check time series
            @test nrow(ts_df) == 3
            @test ts_df.j == [1, 1, 1]
            @test ts_df.i == [1, 1, 1]
            @test ts_df.time == [0.0, 1.0, 2.0]
            @test ts_df.variable == ["x", "x", "x"]
            @test ts_df.value == [10.0, 11.0, 12.0]
            
            # Check parameters
            @test nrow(param_df) == 2
            @test "alpha" in param_df.variable
            @test "beta" in param_df.variable
            
            # Check initial values
            @test nrow(init_df) == 1
            @test init_df.variable == ["x"]
            @test init_df.value == [10.0]
        end

        @testset "Multiple variables" begin
            solve_out = [
                (
                    t = [0.0, 1.0],
                    u = [[10.0, 20.0], [11.0, 21.0]],
                    u0 = [10.0, 20.0],
                    p = (alpha = 0.5,)
                )
            ]
            
            ts_df, param_df, init_df = ensemble_to_df(
                solve_out, [:x, :y], nothing, nothing, 1
            )
            
            @test nrow(ts_df) == 4  # 2 time points × 2 variables
            @test "x" in ts_df.variable
            @test "y" in ts_df.variable
            
            # Check values for x
            x_data = subset(ts_df, :variable => ByRow(==("x")))
            @test x_data.value == [10.0, 11.0]
            
            # Check values for y
            y_data = subset(ts_df, :variable => ByRow(==("y")))
            @test y_data.value == [20.0, 21.0]
        end

        @testset "Multiple trajectories with ensemble_n" begin
            solve_out = [
                (t = [0.0, 1.0], u = [10.0, 11.0], u0 = 10.0, p = (a = 1.0,)),
                (t = [0.0, 1.0], u = [10.5, 11.5], u0 = 10.5, p = (a = 1.0,)),
                (t = [0.0, 1.0], u = [20.0, 21.0], u0 = 20.0, p = (a = 2.0,)),
                (t = [0.0, 1.0], u = [20.5, 21.5], u0 = 20.5, p = (a = 2.0,))
            ]
            
            ts_df, param_df, init_df = ensemble_to_df(
                solve_out, [:x], nothing, nothing, 2  # 2 replicates per condition
            )
            
            # Check j indices (parameter combination)
            @test 1 in ts_df.j
            @test 2 in ts_df.j
            
            # Check i indices (replicate)
            @test all(i -> i in [1, 2], ts_df.i)
            
            # Trajectories 1-2 should have j=1, trajectories 3-4 should have j=2
            traj1_data = subset(ts_df, [:j, :i] => ByRow((j, i) -> j == 1 && i == 1))
            @test all(traj1_data.value .∈ Ref([10.0, 11.0]))
        end

        @testset "With Unitful quantities" begin
            solve_out = [
                (
                    t = [0.0u"s", 1.0u"s", 2.0u"s"],
                    u = [10.0u"m", 11.0u"m", 12.0u"m"],
                    u0 = 10.0u"m",
                    p = (alpha = 0.5u"m/s",)
                )
            ]
            
            ts_df, param_df, init_df = ensemble_to_df(
                solve_out, [:x], nothing, nothing, 1
            )
            
            # Units should be stripped
            @test ts_df.time == [0.0, 1.0, 2.0]
            @test ts_df.value == [10.0, 11.0, 12.0]
            @test param_df.value == [0.5]
        end

        @testset "With intermediaries" begin
            solve_out = [
                (t = [0.0, 1.0], u = [10.0, 11.0], u0 = 10.0, p = nothing)
            ]
            
            intermediaries = [
                (t = [0.0, 1.0], saveval = [100.0, 110.0])
            ]
            
            ts_df, _, _ = ensemble_to_df(
                solve_out, [:x], intermediaries, [:y], 1
            )
            
            # Should have both main and intermediate variables
            @test "x" in ts_df.variable
            @test "y" in ts_df.variable
            
            x_data = subset(ts_df, :variable => ByRow(==("x")))
            y_data = subset(ts_df, :variable => ByRow(==("y")))
            
            @test x_data.value == [10.0, 11.0]
            @test y_data.value == [100.0, 110.0]
        end

        @testset "Vector parameters" begin
            solve_out = [
                (t = [0.0, 1.0], u = [10.0, 11.0], u0 = 10.0, p = [0.5, 2.0, 3.0])
            ]
            
            _, param_df, _ = ensemble_to_df(
                solve_out, [:x], nothing, nothing, 1
            )
            
            @test nrow(param_df) == 3
            @test "p1" in param_df.variable
            @test "p2" in param_df.variable
            @test "p3" in param_df.variable
        end
    end

    @testset "ensemble_summ - statistical summaries" begin
        @testset "Basic statistics" begin
            # Create sample data
            timeseries_df = DataFrame(
                j = [1, 1, 1, 1, 1, 1],
                i = [1, 1, 2, 2, 3, 3],
                time = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0],
                variable = ["x", "x", "x", "x", "x", "x"],
                value = [10.0, 11.0, 12.0, 13.0, 14.0, 15.0]
            )
            
            stats = ensemble_summ(timeseries_df, [0.025, 0.975])
            
            @test nrow(stats) == 2  # 2 time points
            
            # Check time 0.0
            t0_stats = subset(stats, :time => ByRow(==(0.0)))
            @test t0_stats.mean[1] ≈ 12.0  # mean of [10, 12, 14]
            @test t0_stats.median[1] ≈ 12.0
            
            # Check time 1.0
            t1_stats = subset(stats, :time => ByRow(==(1.0)))
            @test t1_stats.mean[1] ≈ 13.0  # mean of [11, 13, 15]
        end

        @testset "Handling NaN and missing" begin
            timeseries_df = DataFrame(
                j = [1, 1, 1, 1],
                i = [1, 2, 3, 4],
                time = [0.0, 0.0, 0.0, 0.0],
                variable = ["x", "x", "x", "x"],
                value = [10.0, NaN, 12.0, missing]
            )
            
            stats = ensemble_summ(timeseries_df)
            
            @test nrow(stats) == 1
            @test stats.mean[1] ≈ 11.0  # mean of [10.0, 12.0]
            @test stats.missing_count[1] == 2
        end

        @testset "All NaN values" begin
            timeseries_df = DataFrame(
                j = [1, 1],
                i = [1, 2],
                time = [0.0, 0.0],
                variable = ["x", "x"],
                value = [NaN, NaN]
            )
            
            stats = ensemble_summ(timeseries_df)
            
            @test isnan(stats.mean[1])
            @test isnan(stats.median[1])
            @test isnan(stats.variance[1])
        end

        @testset "Multiple variables" begin
            timeseries_df = DataFrame(
                j = [1, 1, 1, 1],
                i = [1, 1, 2, 2],
                time = [0.0, 0.0, 0.0, 0.0],
                variable = ["x", "y", "x", "y"],
                value = [10.0, 100.0, 12.0, 120.0]
            )
            
            stats = ensemble_summ(timeseries_df)
            
            @test nrow(stats) == 2  # One row per variable
            
            x_stats = subset(stats, :variable => ByRow(==("x")))
            y_stats = subset(stats, :variable => ByRow(==("y")))
            
            @test x_stats.mean[1] ≈ 11.0
            @test y_stats.mean[1] ≈ 110.0
        end

        @testset "Custom quantiles" begin
            timeseries_df = DataFrame(
                j = [1, 1, 1, 1, 1],
                i = [1, 2, 3, 4, 5],
                time = [0.0, 0.0, 0.0, 0.0, 0.0],
                variable = ["x", "x", "x", "x", "x"],
                value = [1.0, 2.0, 3.0, 4.0, 5.0]
            )
            
            stats = ensemble_summ(timeseries_df, [0.1, 0.5, 0.9])
            
            @test hasproperty(stats, :q1)   # 0.1 → q1
            @test hasproperty(stats, :q5)   # 0.5 → q5
            @test hasproperty(stats, :q9)   # 0.9 → q9
            
            @test stats.q5[1] ≈ 3.0  # median
        end

        @testset "Multiple parameter combinations" begin
            timeseries_df = DataFrame(
                j = [1, 1, 2, 2],
                i = [1, 2, 1, 2],
                time = [0.0, 0.0, 0.0, 0.0],
                variable = ["x", "x", "x", "x"],
                value = [10.0, 12.0, 20.0, 22.0]
            )
            
            stats = ensemble_summ(timeseries_df)
            
            @test nrow(stats) == 2  # One row per j
            
            j1_stats = subset(stats, :j => ByRow(==(1)))
            j2_stats = subset(stats, :j => ByRow(==(2)))
            
            @test j1_stats.mean[1] ≈ 11.0
            @test j2_stats.mean[1] ≈ 21.0
        end
    end

    @testset "Threading variants" begin
        @testset "ensemble_to_df_threaded produces same results" begin
            # Create test data
            solve_out = [
                (t = [0.0, 1.0], u = [[10.0, 20.0], [11.0, 21.0]], 
                 u0 = [10.0, 20.0], p = (a = 1.0,)),
                (t = [0.0, 1.0], u = [[15.0, 25.0], [16.0, 26.0]], 
                 u0 = [15.0, 25.0], p = (a = 2.0,))
            ]
            
            ts1, p1, i1 = ensemble_to_df(solve_out, [:x, :y], nothing, nothing, 1)
            ts2, p2, i2 = ensemble_to_df_threaded(solve_out, [:x, :y], nothing, nothing, 1)
            
            # Results should be identical
            @test ts1.j == ts2.j
            @test ts1.i == ts2.i
            @test ts1.time == ts2.time
            @test ts1.variable == ts2.variable
            @test ts1.value == ts2.value
            
            @test p1 == p2
            @test i1 == i2
        end

        @testset "ensemble_summ_threaded produces same results" begin
            timeseries_df = DataFrame(
                j = repeat([1], 100),
                i = 1:100,
                time = repeat([0.0], 100),
                variable = repeat(["x"], 100),
                value = randn(100)
            )
            
            stats1 = ensemble_summ(timeseries_df)
            stats2 = ensemble_summ_threaded(timeseries_df)
            
            # Results should be very similar (small numerical differences possible)
            @test stats1.mean ≈ stats2.mean
            @test stats1.median ≈ stats2.median
            @test stats1.variance ≈ stats2.variance
        end
    end

    @testset "Integration tests" begin
        @testset "Full workflow" begin
            # Generate parameters
            param_ranges = Dict(:alpha => [0.1, 0.5], :beta => [1.0, 2.0])
            combinations, total = generate_param_combinations(
                param_ranges; crossed=true, n_replicates=2
            )
            
            @test length(combinations) == 4
            @test total == 8
            
            # Create mock ensemble results
            solve_out = []
            for (j, combo) in enumerate(combinations)
                for i in 1:2
                    push!(solve_out, (
                        t = [0.0, 1.0],
                        u = [combo[1] * 10, combo[1] * 10 + 1],
                        u0 = combo[1] * 10,
                        p = (alpha = combo[1], beta = combo[2])
                    ))
                end
            end
            
            # Convert to DataFrame
            ts_df, param_df, init_df = ensemble_to_df(
                solve_out, [:x], nothing, nothing, 2
            )
            
            @test nrow(ts_df) > 0
            @test maximum(ts_df.j) == 4  # 4 parameter combinations
            @test maximum(ts_df.i) == 2  # 2 replicates
            
            # Compute summaries
            stats = ensemble_summ(ts_df)
            
            @test nrow(stats) == 8  # 4 combinations × 2 time points
            @test all(stats.missing_count .== 0)
        end
    end

end