# Run just the Delsarte sphere packing bounds for tests vs truth
# Compared to Pfender, Florian. "Improved Delsarte bounds for spherical codes in small dimensions."
# Journal of Combinatorial Theory, Series A 114.6 (2007): 1133-1147.
# They use n for dimension, we use d.
using Test

kmax = 30 # maximum number of Gegenbauer coefficients to include

# Upper bound of N for d = 3, angle = 44.3 degrees. Should have N = 24.
d = 3
Amax = cos(44.3 * pi / 180)
xs = range(-1, Amax, length =20)
fk = solve_delsarte_socp2(d, kmax, xs)
@test(floor(compute_from_expansion(d, fk, 1)) == 24)

# Upper bound of N for d = 5, angle = 85.39 degrees. Should have N = 11.
d = 5
Amax = cos(85.39 * pi / 180)
xs = range(-1, Amax, length =10)
fk = solve_delsarte_socp2(d, kmax, xs)
@test(floor(compute_from_expansion(d, fk, 1)) == 11)

# Kissing number in d = 3
Amax = cos(pi / 3)
d = 3
xs = range(-1, Amax, length =8)
fk = solve_delsarte_socp2(d, kmax, xs)
@test(floor(compute_from_expansion(d, fk, 1)) == 13)

# Kissing number in d = 8
d = 8
xs = range(-1, Amax, length =50)
fk = solve_delsarte_socp2(d, kmax, xs)
@test(floor(compute_from_expansion(d, fk, 1)) == 240)

 # Kissing number in d = 9
xs = range(-1, Amax, length =75)
d = 9
fk = solve_delsarte_socp2(d, kmax, xs)
@test(floor(compute_from_expansion(d, fk, 1)) == 380)

# Kissing number in d = 16
d = 16
xs = range(-1, Amax, length =100)
fk = solve_delsarte_socp2(d, kmax, xs)
@test(floor(compute_from_expansion(d, fk, 1)) >=  8313)
# Will need more refined xs to meet bound more closely
