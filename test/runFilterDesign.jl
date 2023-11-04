# Run the filter design code


N = 4;
M = 60; # so number of nodes is 61
omegas = range(0, .5, length = M+1);

# AR(1)
rk = [.95^k for k = 0:(2*N-1)];
x = get_optimal_filter(rk, omegas, N);
@show compute_coding_gain(x, rk)

# AR(2)
rho = .975;
theta = pi/3;
rk = zeros(2*N);
rk[1] = 1;
rk[2] = 2 * rho  * cos(theta) / (1 + rho^2);
for i = 3:length(rk)
    rk[i] = 2 * rho * cos(theta) * rk[i-1] - rho^2 * rk[i-2];
end
x = get_optimal_filter(rk, omegas, N);
@show compute_coding_gain(x, rk)

# box-spec
fs = .225;
rk = [sin(2 * pi * fs * k) / (2*pi*fs*k) for k = 1:(2*N-1)];
rk = [1; rk];
x = get_optimal_filter(rk, omegas, N);
@show compute_coding_gain(x, rk)



# See scaling of our method with 
Ms = [10, 20, 50, 90]
ar1gains = zeros(lastindex(Ms))
ar2gains = zeros(lastindex(Ms))
boxgains = zeros(lastindex(Ms))
for Mi = 1:lastindex(Ms)
    M = Ms[Mi]
    omegas = range(0, .5, length = M+1);
    # AR(1)
    rk = [.95^k for k = 0:(2*N-1)];
    x = get_optimal_filter(rk, omegas, N);
    ar1gains[Mi] = compute_coding_gain(x, rk)

    # AR(2)
    rho = .975;
    theta = pi/3;
    rk = zeros(2*N);
    rk[1] = 1;
    rk[2] = 2 * rho  * cos(theta) / (1 + rho^2);
    for i = 3:length(rk)
        rk[i] = 2 * rho * cos(theta) * rk[i-1] - rho^2 * rk[i-2];
    end
    x = get_optimal_filter(rk, omegas, N);
    ar2gains[Mi] = compute_coding_gain(x, rk)

    # box-spec
    fs = .225;
    rk = [sin(2 * pi * fs * k) / (2*pi*fs*k) for k = 1:(2*N-1)];
    rk = [1; rk];
    x = get_optimal_filter(rk, omegas, N);
    boxgains[Mi] = compute_coding_gain(x, rk)
end
