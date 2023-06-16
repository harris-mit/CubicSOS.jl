# This file is a driver for SpherePacking.jl, the main second order implementations
# of the sphere packing problems. It is not use for any testing. We rather compute the
# minimum number of samples required for both the SOCP method and LP method to return
# the true Delsarte bound.

using BenchmarkTools

# Solve the Delsarte sphere packing in a sphere.
# Define a dictionary of the Delsarte bound for dimension n.
true_delsarte = Dict(
    2 => 6,
    3 => 13,
    4 => 25,
    8 => 240,
    24 => 196560
)
# This only involves polynomial evaluation, so we don't deal with large coefficients.
# The method sets the Gegenbauer coefficients to be the decision variables.
# We use a low degree cubic spline to certify the nonpositivity condition.

"""
    get_num_required_points_socp(d)
Find the number of required points in the discretization for the SOCP
for the Delsarte bound to return the true value.
"""
function get_num_required_points_socp(d)
    Amax = cos(pi / 3)
    kmax = 20
    # First find an upper bound on what's necessary
    estimated = 1
    density = 0
    while density != true_delsarte[d]
        estimated = 2 * estimated
        xs = range(-1, Amax, length = estimated)
        fk = solve_delsarte_socp(d, kmax, xs)
        density = floor(compute_from_expansion(d, fk, 1))
    end
    # now binary search in between
    large_enough = estimated
    too_small = 1
    while large_enough - too_small > 1
        in_between = Int(floor((too_small + large_enough)/2))
        xs = range(-1, Amax, length = in_between)
        fk = solve_delsarte_socp(d, kmax, xs)
        density = floor(compute_from_expansion(d, fk, 1))
        if density == true_delsarte[d]
            large_enough = in_between
        else
            too_small = in_between
        end
    end
    return large_enough
end

"""
    get_num_required_points_lp(d)
Find the number of required points in the discretization for the LP
for the Delsarte bound to return the true value.
"""
function get_num_required_points_lp(d)
    Amax = cos(pi / 3)
    kmax = 20
    # First find an upper bound on what's necessary
    estimated = 1
    density = 0
    while density != true_delsarte[d]
        estimated = 2 * estimated
        if estimated > 100000
            return Inf
        end
        xs = range(-1, Amax, length = estimated)
        fk = solve_delsarte_lp(d, kmax, xs)
        density = floor(compute_from_expansion(d, fk, 1))
    end
    # now binary search in between
    large_enough = estimated
    too_small = 1
    while large_enough - too_small > 1
        in_between = Int(floor((too_small + large_enough)/2))
        xs = range(-1, Amax, length = in_between)
        fk = solve_delsarte_lp(d, kmax, xs)
        density = floor(compute_from_expansion(d, fk, 1))
        if density == true_delsarte[d]
            large_enough = in_between
        else
            too_small = in_between
        end
    end
    return large_enough
end

ds = [2,3,4,8,24]
socp_ns = zeros(length(ds))
lp_ns = zeros(length(ds))
for di in 1:length(ds)
    socp_ns[di] = get_num_required_points_socp(ds[di])
    lp_ns[di] = get_num_required_points_lp(ds[di])
end

# rerun with times
socp_times = zeros(length(ds))
lp_times = zeros(length(ds))
kmax = 20
for di = 1:length(ds)
    d = ds[di]
    xs = range(-1, Amax, length = Int(socp_ns[di]))
    socp_times[di] = @belapsed solve_delsarte_socp(d, kmax, xs);
    xs = range(-1, Amax, length = Int(lp_ns[di]))
    lp_times[di] = @belapsed solve_delsarte_lp(d, kmax, xs);
end

# Testing our SOCP certificate:
gktrunc = value.(fk)[1:14] # Truncating off the ones that are basically 0 (and possibly slightly negative.)
gktrunc[1] = 1
all(gktrunc .> -1e-10 )
compute_from_expansion(d, gktrunc, 1)
xx = range(-1, Amax, length =1000)
yy = compute_from_expansion.(Ref(d), Ref(gktrunc), xx)
maximum(yy)
