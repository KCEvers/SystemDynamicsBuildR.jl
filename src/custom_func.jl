"""
    custom_func

Custom utility functions for system dynamics modeling, including:
- Signal generation (ramps, steps, pulses, seasonal waves)
- Interpolation and extrapolation
- Mathematical utilities (logistic, logit, expit)
- Random sampling functions
- String/array utilities

This module is designed to work seamlessly with Unitful quantities.
"""
module custom_func

using Unitful
using DataInterpolations
using Distributions
using ..unit_func: convert_u

export itp, make_ramp, make_step, make_pulse, make_seasonal
export round_IM, logit, expit, logistic
export nonnegative, rbool, rdist
export indexof, contains_IM, round_
export is_function_or_interp, ⊕

# ============================================================================
# Type Checking Utilities
# ============================================================================

"""
    is_function_or_interp(x)

Check if `x` is a Function or an AbstractInterpolation object.

# Examples
```julia
julia> is_function_or_interp(sin)
true

julia> is_function_or_interp(itp([1, 2], [3, 4]))
true

julia> is_function_or_interp(5)
false
```
"""
is_function_or_interp(x) = isa(x, Function) || isa(x, DataInterpolations.AbstractInterpolation)

# ============================================================================
# Interpolation Functions
# ============================================================================

"""
    itp(x, y; method="linear", extrapolation="nearest")

Create an interpolation function from vectors `x` and `y`.

# Arguments
- `x::AbstractVector`: Independent variable values (will be sorted)
- `y::AbstractVector`: Dependent variable values
- `method::String="linear"`: Interpolation method ("linear" or "constant")
- `extrapolation::String="nearest"`: Extrapolation behavior ("nearest" or "NA")

# Returns
- `AbstractInterpolation`: Interpolation object that can be called as a function

# Examples
```julia
julia> f = itp([1, 2, 3], [10, 20, 30])
julia> f(1.5)
15.0

julia> f = itp([1, 3, 2], [10, 30, 20])  # Automatically sorted
julia> f(2.5)
25.0
```
"""
function itp(x, y; method="linear", extrapolation="nearest")
    # Ensure y is sorted along x
    idx = sortperm(x)
    x = x[idx]
    y = y[idx]

    # Extrapolation rule: What happens outside of defined values?
    # Rule "NA": return NaN; Rule "nearest": return closest value
    rule_method = if extrapolation == "NA"
        DataInterpolations.ExtrapolationType.None
    elseif extrapolation == "nearest"
        DataInterpolations.ExtrapolationType.Constant
    else
        extrapolation
    end

    if method == "constant"
        func = DataInterpolations.ConstantInterpolation(y, x; extrapolation=rule_method)
    elseif method == "linear"
        func = DataInterpolations.LinearInterpolation(y, x; extrapolation=rule_method)
    else
        throw(ArgumentError("Method must be 'constant' or 'linear', got: $method"))
    end

    return func
end

# ============================================================================
# Signal Generation Functions
# ============================================================================

"""
    make_ramp(time_units, times, start, finish, height=1.0)

Create a ramp signal that linearly increases from 0 to `height` between `start` and `finish` times.

The ramp starts at height 0 at time `start`, increases linearly, and reaches `height` at time `finish`.
Outside this range, the value is constant (0 before start, height after finish).

# Arguments
- `time_units`: Units for time (e.g., u"yr", u"d")
- `times`: Time vector or range (start, end)
- `start`: Start time of ramp
- `finish`: End time of ramp
- `height=1.0`: Maximum height of ramp (can be negative for decreasing ramp)

# Returns
- Interpolation function that can be evaluated at any time

# Examples
```julia
julia> r = make_ramp(u"yr", [0.0, 10.0], 2.0, 5.0, 10.0)
julia> r(3.5)  # Halfway through ramp
5.0
```
"""
function make_ramp(time_units, times, start, finish, height=1.0)
    @assert start < finish "The finish time of the ramp cannot be before the start time. To specify a decreasing ramp, set the height to a negative value."

    # Normalize units between times and ramp parameters
    start, finish = _normalize_time_units(times, time_units, start, finish)
    
    # Initialize ramp height
    start_h_ramp = 0.0
    add_y = 0.0
    
    # Match height units
    if eltype(height) <: Unitful.Quantity
        start_h_ramp = convert_u(start_h_ramp, Unitful.unit(height))
        add_y = convert_u(0.0, Unitful.unit(height))
    elseif !(eltype(height) <: Unitful.Quantity)
        height = convert_u(height, Unitful.unit(start_h_ramp))
        add_y = convert_u(0.0, Unitful.unit(start_h_ramp))
    end

    x = [start, finish]
    y = [start_h_ramp, height]

    # If the ramp is after the start time, add a zero at the start
    if start > first(times)
        x = [first(times); x]
        y = [add_y; y]
    end

    func = itp(x, y, method="linear", extrapolation="nearest")
    return func
end

"""
    make_step(time_units, times, start, height=1.0)

Create a step signal that jumps from 0 to `height` at time `start`.

# Arguments
- `time_units`: Units for time
- `times`: Time vector or range
- `start`: Time when step occurs
- `height=1.0`: Height of step

# Returns
- Interpolation function representing the step signal

# Examples
```julia
julia> s = make_step(u"s", [0.0, 10.0], 5.0, 2.0)
julia> s(4.9)  # Before step
0.0
julia> s(5.1)  # After step
2.0
```
"""
function make_step(time_units, times, start, height=1.0)
    # Normalize units
    start = _normalize_single_time(times, time_units, start)
    
    add_y = eltype(height) <: Unitful.Quantity ? convert_u(0.0, Unitful.unit(height)) : 0.0

    x = [start, times[2]]
    y = [height, height]

    # If the step is after the start time, add a zero at the start
    if start > first(times)
        x = [first(times); x]
        y = [add_y; y]
    end

    func = itp(x, y, method="constant", extrapolation="nearest")
    return func
end

"""
    make_pulse(time_units, times, start, height=1.0, width=1.0*time_units, repeat_interval=nothing)

Create a pulse signal with specified width and optional repetition.

# Arguments
- `time_units`: Units for time
- `times`: Time vector or range
- `start`: Start time of first pulse
- `height=1.0`: Height of pulse
- `width=1.0*time_units`: Duration of each pulse
- `repeat_interval=nothing`: Time between pulse starts (nothing = single pulse)

# Returns
- Interpolation function representing the pulse train

# Examples
```julia
julia> p = make_pulse(u"s", [0.0, 20.0], 5.0, 1.0, 2.0, 10.0)  # Pulse every 10s
julia> p(6.0)  # During first pulse
1.0
julia> p(8.0)  # Between pulses
0.0
```
"""
function make_pulse(time_units, times, start, height=1.0, width=1.0 * time_units, repeat_interval=nothing)
    # Validate width
    width_value = eltype(width) <: Unitful.Quantity ? Unitful.ustrip(convert_u(width, time_units)) : width
    if width_value <= 0.0
        throw(ArgumentError("The width of the pulse cannot be equal to or less than 0; to indicate an 'instantaneous' pulse, specify the simulation step size (dt)."))
    end

    # Normalize units
    start, width, repeat_interval = _normalize_pulse_units(times, time_units, start, width, repeat_interval)

    # Define start and end times of pulses
    last_time = last(times)
    step_size = isnothing(repeat_interval) ? last_time * 2 : repeat_interval
    start_ts = collect(start:step_size:last_time)
    end_ts = start_ts .+ width

    # Build signal as vectors of times and y-values
    signal_times = [start_ts; end_ts]
    signal_y = [fill(height, length(start_ts)); fill(0, length(end_ts))]

    add_y = eltype(height) <: Unitful.Quantity ? convert_u(0.0, Unitful.unit(height)) : 0.0

    # Add zeros at boundaries if needed
    if minimum(start_ts) > first(times)
        signal_times = [first(times); signal_times]
        signal_y = [add_y; signal_y]
    end

    if maximum(end_ts) < last_time
        signal_times = [signal_times; last_time]
        signal_y = [signal_y; add_y]
    end

    # Sort by time
    perm = sortperm(signal_times)
    x = signal_times[perm]
    y = signal_y[perm]
    func = itp(x, y, method="constant", extrapolation="nearest")

    return func
end

"""
    make_seasonal(times, dt, period=u"1yr", shift=u"0yr")

Create a seasonal cosine wave with specified period and phase shift.

The wave oscillates between -1 and 1 with the formula: cos(2π(t - shift)/period)

# Arguments
- `dt`: Time step for sampling
- `times`: Time range [start, end]
- `period=u"1yr"`: Period of oscillation
- `shift=u"0yr"`: Phase shift (positive = delay)

# Returns
- Interpolation function representing the seasonal pattern

# Examples
```julia
julia> wave = make_seasonal(0.1u"yr", [0.0u"yr", 2.0u"yr"], 1.0u"yr")
julia> wave(0.0u"yr")  # Peak of cosine
1.0
julia> wave(0.5u"yr")  # Trough
-1.0
```
"""
function make_seasonal(dt, times, period=u"1yr", shift=u"0yr")
    @assert Unitful.ustrip(period) > 0 "The period of the seasonal wave must be greater than 0."

    time_vec = times[1]:dt:times[2]
    phase = 2 * pi .* (time_vec .- shift) ./ period
    y = cos.(phase)
    func = itp(time_vec, y, method="linear", extrapolation="nearest")

    return func
end

# ============================================================================
# Helper Functions for Unit Conversion
# ============================================================================

"""
    _normalize_time_units(times, time_units, start, finish)

Internal helper to normalize time units between simulation times and signal parameters.
"""
function _normalize_time_units(times, time_units, start, finish)
    if eltype(times) <: Unitful.Quantity
        # Times have units, ensure start/finish match
        if !(eltype(start) <: Unitful.Quantity)
            start = convert_u(start, time_units)
        end
        if !(eltype(finish) <: Unitful.Quantity)
            finish = convert_u(finish, time_units)
        end
    else
        # Times are unitless, convert start/finish if they have units
        if eltype(start) <: Unitful.Quantity
            start = Unitful.ustrip(convert_u(start, time_units))
        end
        if eltype(finish) <: Unitful.Quantity
            finish = Unitful.ustrip(convert_u(finish, time_units))
        end
    end
    return start, finish
end

"""
    _normalize_single_time(times, time_units, start)

Internal helper to normalize a single time value.
"""
function _normalize_single_time(times, time_units, start)
    if eltype(times) <: Unitful.Quantity
        if !(eltype(start) <: Unitful.Quantity)
            start = convert_u(start, time_units)
        end
    else
        if eltype(start) <: Unitful.Quantity
            start = Unitful.ustrip(convert_u(start, time_units))
        end
    end
    return start
end

"""
    _normalize_pulse_units(times, time_units, start, width, repeat_interval)

Internal helper to normalize time units for pulse signals.
"""
function _normalize_pulse_units(times, time_units, start, width, repeat_interval)
    if eltype(times) <: Unitful.Quantity
        if !(eltype(start) <: Unitful.Quantity)
            start = convert_u(start, time_units)
        end
        if !(eltype(width) <: Unitful.Quantity)
            width = convert_u(width, time_units)
        end
        if !isnothing(repeat_interval) && !(eltype(repeat_interval) <: Unitful.Quantity)
            repeat_interval = convert_u(repeat_interval, time_units)
        end
    else
        if eltype(start) <: Unitful.Quantity
            start = Unitful.ustrip(convert_u(start, time_units))
        end
        if eltype(width) <: Unitful.Quantity
            width = Unitful.ustrip(convert_u(width, time_units))
        end
        if !isnothing(repeat_interval) && eltype(repeat_interval) <: Unitful.Quantity
            repeat_interval = Unitful.ustrip(convert_u(repeat_interval, time_units))
        end
    end
    return start, width, repeat_interval
end

# ============================================================================
# Mathematical Functions
# ============================================================================

"""
    round_IM(x::Real, digits::Int=0)

Round a number using Insight Maker's convention where 0.5 rounds up.

Note: Julia's default `round()` uses banker's rounding where 0.5 rounds to nearest even.
This function always rounds 0.5 up to match Insight Maker behavior.

# Examples
```julia
julia> round_IM(0.5)
1.0
julia> round_IM(1.5)
2.0
julia> round_IM(2.5)
3.0
```
"""
function round_IM(x::Real, digits::Int=0)
    scaled_x = x * 10.0^digits
    frac = scaled_x % 1
    
    # Check if fractional part is exactly ±0.5
    if abs(frac) == 0.5
        return ceil(scaled_x) / 10.0^digits
    else
        return round(scaled_x) / 10.0^digits
    end
end

"""
    logit(p)

Compute the logit (log-odds) function: log(p / (1 - p))

# Examples
```julia
julia> logit(0.5)
0.0
julia> logit(0.75)
1.0986122886681098
```
"""
logit(p) = log(p / (1 - p))

"""
    expit(x)

Compute the expit (inverse logit) function: 1 / (1 + exp(-x))

Also known as the logistic sigmoid function.

# Examples
```julia
julia> expit(0.0)
0.5
julia> expit(10.0)
0.9999546021312976
```
"""
expit(x) = 1 / (1 + exp(-x))

"""
    logistic(x, slope=1.0, midpoint=0.0, upper=1.0)

Compute a generalized logistic function with adjustable slope, midpoint, and upper bound.

Formula: upper / (1 + exp(-slope * (x - midpoint)))

# Arguments
- `x`: Input value
- `slope=1.0`: Steepness of the curve
- `midpoint=0.0`: x-value at the inflection point
- `upper=1.0`: Maximum asymptotic value

# Examples
```julia
julia> logistic(0.0, 1.0, 0.0, 1.0)  # Standard logistic at midpoint
0.5
julia> logistic(5.0, 2.0, 5.0, 10.0)  # Steeper curve, shifted
5.0
```
"""
function logistic(x, slope=1.0, midpoint=0.0, upper=1.0)
    @assert isfinite(Unitful.ustrip(slope)) && isfinite(Unitful.ustrip(midpoint)) && isfinite(Unitful.ustrip(upper)) "slope, midpoint, and upper must be finite numeric values"
    return upper / (1 + exp(-slope * (x - midpoint)))
end

"""
    nonnegative(x)

Ensure value(s) are non-negative by returning max(0, x).

Works with scalars, arrays, and Unitful quantities.

# Examples
```julia
julia> nonnegative(-5)
0.0
julia> nonnegative(3)
3.0
julia> nonnegative([-1, 2, -3])
3-element Vector{Float64}: [0.0, 2.0, 0.0]
```
"""
nonnegative(x::Real) = max(0.0, x)
nonnegative(x::Unitful.Quantity) = max(0.0, Unitful.ustrip(x)) * Unitful.unit(x)
nonnegative(x::AbstractArray{<:Real}) = max.(0.0, x)
nonnegative(x::AbstractArray{<:Unitful.Quantity}) = max.(0.0, Unitful.ustrip.(x)) .* Unitful.unit.(x)

# ============================================================================
# Rounding Utilities
# ============================================================================

"""
    round_(x, digits=0)

Flexible rounding function that handles both regular numbers and Unitful quantities.

# Examples
```julia
julia> round_(3.14159, digits=2)
3.14
julia> round_(5.6u"m", digits=0)
6.0 m
```
"""
round_(x, digits::Real) = round(x, digits=round(Int, digits))
round_(x; digits::Real=0) = round(x, digits=round(Int, digits))
round_(x::Unitful.Quantity, digits::Real) = round(Unitful.ustrip(x), digits=round(Int, digits)) * Unitful.unit(x)
round_(x::Unitful.Quantity; digits::Real=0) = round(Unitful.ustrip(x), digits=round(Int, digits)) * Unitful.unit(x)

# ============================================================================
# Random Sampling Functions
# ============================================================================

"""
    rbool(p)

Generate a random boolean value with probability `p` of being true.

Equivalent to Insight Maker's RandBoolean() function.

# Examples
```julia
julia> rbool(0.7)  # 70% chance of true
true
julia> rbool(0.0)  # Always false
false
```
"""
rbool(p) = rand() < p

"""
    rdist(a::Vector, b::Vector{<:Real})

Sample randomly from vector `a` with probabilities given by vector `b`.

Probabilities are automatically normalized to sum to 1.

# Arguments
- `a`: Vector of values to sample from
- `b`: Vector of probabilities (will be normalized, must all be non-negative)

# Examples
```julia
julia> rdist(["red", "green", "blue"], [0.5, 0.3, 0.2])
"red"  # (with 50% probability)
```
"""
function rdist(a::Vector{T}, b::Vector{<:Real}) where T
    if length(a) != length(b)
        throw(ArgumentError("Length of a and b must match"))
    end
    
    if any(x -> x < 0, b)
        throw(ArgumentError("All probabilities must be non-negative"))
    end
    
    b_sum = sum(b)
    if b_sum <= 0
        throw(ArgumentError("Sum of probabilities must be positive"))
    end
    
    b_normalized = b / b_sum
    return a[rand(Distributions.Categorical(b_normalized))]
end

# ============================================================================
# String and Array Utilities
# ============================================================================

"""
    indexof(haystack, needle)

Find the index of `needle` in `haystack`.

Works with both strings and arrays. Returns 0 if not found (Insight Maker convention).

# Examples
```julia
julia> indexof("hello", "ll")
3
julia> indexof([1, 2, 3, 4], 3)
3
julia> indexof("hello", "x")
0
```
"""
function indexof(haystack, needle)
    if isa(haystack, AbstractString) && isa(needle, AbstractString)
        pos = findfirst(needle, haystack)
        return isnothing(pos) ? 0 : first(pos)
    else
        pos = findfirst(==(needle), haystack)
        return isnothing(pos) ? 0 : pos
    end
end

"""
    contains_IM(haystack, needle)

Check if `haystack` contains `needle`.

Works with both strings and arrays.

# Examples
```julia
julia> contains_IM("hello world", "world")
true
julia> contains_IM([1, 2, 3], 2)
true
julia> contains_IM("hello", "x")
false
```
"""
function contains_IM(haystack, needle)
    if isa(haystack, AbstractString) && isa(needle, AbstractString)
        return occursin(needle, haystack)
    else
        return needle in haystack
    end
end

# ============================================================================
# Operators
# ============================================================================

"""
    ⊕(x, y)

Modulus operator (x mod y).

Unicode alternative to `mod(x, y)`.

# Examples
```julia
julia> 7 ⊕ 3
1
julia> 10 ⊕ 5
0
```
"""
⊕(x, y) = mod(x, y)

end # module
