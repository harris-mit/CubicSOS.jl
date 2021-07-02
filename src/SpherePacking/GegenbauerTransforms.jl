using PyCall

export populate_gegenbauer_transform!
export populate_gegenbauer_transform_analytic!, populate_gegenbauer_transform_d3!
export compute_from_expansion

"""
    scaled_gegenbauer(d, k)
Scale the C_k (with specialpolynomials) to Q_k
"""
function scaled_gegenbauer(d, k)
    v = zeros(k+1)
    v[k+1] = 1
    dfactor = 1
    p = Gegenbauer{0}([0])
    if d >= 3
        dfactor = (d + 2 * k - 2) / (d - 2)
        # get a vector that is an indicator for C_k
        p =  Gegenbauer{(d-2)/2}(v)
    else #d = 2
        dfactor = 2
        if k == 0
            dfactor = 1 #factor is 1 if k = 0
        end
        p = Chebyshev(v)
    end
    return dfactor * p
end

"""
    get_hermite_basis_poly(x, j, xi, xip1)
Version of get_hermite_basis that returns a polynomial
"""
function get_hermite_basis_poly(j, xi, xip1)
    delta = xip1 - xi
    p = Polynomial([1])
    if j == 1
        p = Polynomial([1 - 2 * xi^3 / delta^3 - 3 * xi^2 / delta^2,
                        6*xi^2 / delta^3 + 6 * xi / delta^2,
                        -6 * xi / delta^3 - 3/delta^2, 2/delta^3])
    elseif j == 2
        p = Polynomial([-xi - xi^3/delta^2 - 2 * xi^2 / delta,
                        1 + 3 * xi^2 / delta^2 + 4*xi/delta,
                        -3 * xi / delta^2 - 2 / delta,
                        1/delta^2])
    elseif j == 3
        p = Polynomial([2 * xi^3 / delta^3 + 3 * xi^2 / delta^2,
                        -6 * xi^2 / delta^3 - 6 * xi / delta^2,
                        3 / delta^2 + 6 * xi / delta^3,
                        -2 / delta^3])
    elseif j == 4
        p = Polynomial([-xi^3 / delta^2 - xi^2 / delta,
                        2 * xi / delta + 3 * xi^2 / delta^2,
                        -3 * xi / delta^2 - 1/delta,
                        1/delta^2])
    else
        error("Invalid basis number")
    end
    return p
end


"""
    populate_gegenbauer_transform!(Gcoefs, d, xs)
Compute the integrals
C_{ikj} = int_{I_i} Q_k(x) b_j(x) (1-x^2)^((d-3)/2)dx
for each subinterval i, polynomial k, and basis element j.
Fixing k and taking the dot product of Gcoefs with the spline coefficients for each
subinterval will give a_d g_k Q_k(1) by the orthogonality of Q_k.

Note that matrix index k refers to polynomial index k - 1.
(i.e. the first entry corresponds to the constant polynomial)
The computation method uses adaptive quadrature unless the error does not get better
than tol and then the integral is computed symbolically.
"""
function populate_gegenbauer_transform!(Gcoefs, d, xs, tol = 1e-7)
    ad = get_ad(d)
    w = x -> (1 - x^2)^((d-3)/2)# weight function
    x = pyimport("sympy").symbols("x")
    wpoly = (1 - x^2)^((d-3)/2)# weight function
    for k = 1:size(Gcoefs, 1) # coefficients
        print("Finished k = ", k,"\n")
        for i = 1:size(Gcoefs, 2) # x interval
            for j = 1:size(Gcoefs, 3) # 4
                # Define function to integrate
                xi = xs[i]
                xip1 = xs[i + 1]
                Qk = convert(Polynomial, scaled_gegenbauer(d, k-1)) # start at 0
                f = x -> Qk(x) * get_hermite_basis(x, j, xi, xip1) * w(x)
                (v, err) = hquadrature(f, xi, xip1, maxevals = 10^7, rtol=1e-4)
                if err > tol
                    Qkpoly = poly2sym(Qk)
                    print("The adaptive quadrature premature terminated")
                    @show (k,i,j,err,v)
                    b = get_hermite_basis_poly(j, xi, xip1)
                    b = poly2sym(b)
                    f = Qkpoly * b * wpoly
                    v = pyimport("sympy").N(pyimport("sympy").integrate(f, (x, xi, xip1)))
                    @show "Analytic v", v
                end
                Gcoefs[k,i,j] = v / Qk(1) / ad
            end
        end
    end
end


"""
    populate_gegenbauer_transform_d3!(Gcoefs, d, xs)
Computes the same integrals as populate_gegenbauer_transform! but assumes that
d = 3, so the integrand is a polynomial, and therefore all integrals can be computed
analytically with Polynomials.jl
"""
function populate_gegenbauer_transform_d3!(Gcoefs, d, xs)
    if d != 3
        error("This function assumes d = 3, so w = 1")
    end
    ad = get_ad(d)
    # weight function is 1 when d = 3
    for k = 1:size(Gcoefs, 1) # coefficients
        print("Started k = ", k,"\n")
        for i = 1:size(Gcoefs, 2) # x interval
            for j = 1:size(Gcoefs, 3) # 4
                # Define function to integrate
                xi = xs[i]
                xip1 = xs[i + 1]
                Qk = scaled_gegenbauer(d, k-1)# start at 0
                f = Qk * get_hermite_basis_poly(j, xi, xip1)
                v = integrate(f, xi, xip1)
                Gcoefs[k,i,j] = v / Qk(1) / ad
            end
        end
    end
end


"""
    populate_gegenbauer_transform_analytic!(Gcoefs_analytic, d, xs)
Computes the same integrals as populate_gegenbauer_transform! with symbolic
integration. This is much slower.
"""
function populate_gegenbauer_transform_analytic!(Gcoefs_analytic, d, xs)
    # use analytic integrals
    ad = get_ad(d)
    x = pyimport("sympy").symbols("x")
    w = (1 - x^2)^((d-3)/2)# weight function
    for k = 1:size(Gcoefs_analytic, 1) # coefficients
        print("Finished k = ", k)
        Qk = convert(Polynomial, scaled_gegenbauer(d, k-1))
        Qk1 = Qk(1) # start at 0
        Qk = poly2sym(Qk)
        for i = 1:size(Gcoefs_analytic, 2) # x interval
            for j = 1:size(Gcoefs_analytic, 3) # 4
                # Define function to integrate
                xi = xs[i]
                xip1 = xs[i + 1]
                b = get_hermite_basis_poly(j, xi, xip1)
                b = poly2sym(b)
                f = Qk * b * w
                v = pyimport("sympy").N(pyimport("sympy").integrate(f, (x, xi, xip1)))
                Gcoefs_analytic[k,i,j] = v / Qk1 / ad
            end
        end
    end
end

"""
    get_ad(d)
The normalizing factor that depends only on the dimension.
"""
function get_ad(d)
    w = x -> (1 - x^2)^((d-3)/2)# weight function
    Qk = convert(Polynomial, scaled_gegenbauer(d, 0)) # any k will do, here use 0
    f = x -> Qk(x)^2 * w(x) # since Q_0 = 1, this is just w(x)
    (v, err) = hquadrature(f, -1, 1, maxevals = 10^6)
    if err > 1e-7
        error("The adaptive quadrature premature terminated")
    end
    return v / Qk(1) # This is the same as just v, since Q_0(1) = 1
end

"""
    poly2sym(p)
Convert a polynomial to a sympy expression.
"""
function poly2sym(p)
    x = pyimport("sympy").symbols("x")
    c = p.coeffs
    result = 0
    for i = 1:length(c)
        result += c[i] * x^(i-1)
    end
    return result
end

"""
    compute_from_expansion(d, gk, x)
If a function is defined by the Gegenbauer coefficients in gk, then compute
the evaluation of the function at x.
"""
function compute_from_expansion(d, gk, x)
    val = 0
    for k = 1:length(gk)
        val += gk[k] * scaled_gegenbauer(d, k-1)(x)
    end
    return val
end

"""
    compute_deriv_from_expansion(d, gk, x)
If a function is defined by the Gegenbauer coefficients in gk, then compute
the evaluation of the derivative of the function at x.
"""
function compute_deriv_from_expansion(d, gk, x)
    val = 0
    for k = 1:length(gk)
        val += gk[k] * derivative(scaled_gegenbauer(d, k-1))(x)
    end
    return val
end


"""
    compute_4th_div_bound(d, gk)
Computes a bound on the 4th derivative of a function defined by the Gegenbauer
coefficients gk by computing the l1 norm of the coefficients of each term.
This assumes that gk >= 0.
This bound is rather weak; consider improving.
"""
function compute_4th_div_bound(d, gk)
    val = 0
    for k = 1:length(gk)
        fthdiv = derivative(scaled_gegenbauer(d, k-1), 4)
        p = convert(Polynomial, fthdiv)
        val += gk[k] * sum(abs, p.coeffs)
    end
    return val
end
