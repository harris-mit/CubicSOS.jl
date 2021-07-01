# Compute the Fourier transforms necessary to set up the SOCP.
# The analytic form of the integral
# \int_0^\pi e^{-2 \pi i s r \cos(\theta)} \sin(\theta)^{n-2} d\theta
# was computed with Mathematica.


"""
We've integrated out the theta dependence
r = |x|, the argument of the model F(r)
d = the order of the derivative we need
s = the frequency argument of the fourier transform
n = the dimension of the problem
j = which of the four bases to use
xi = the left hand endpoint of the interval (used for Hermite basis)
xip1 = x_{i+1}, or the right hand endpoint of the interval
Returns the integrand of the Fourier integral.
TODO: Rather than computing the bessel functions on each evaluation,
the thetaint integral can be computed one time elsewhere.
"""
function get_fourier_integrand(r, d, s, n, j, xi, xip1, omega)
    if d != 0
        error("Fourier integrand requires d = 0")
    end
    thetaint = 0
    basis_factor = get_hermite_basis(r, j, xi, xip1)
    if n == 8
        tpirs = 2 * pi * r * s
        thetaint = 15*(-tpirs * besselj(0, tpirs) + (2 - pi^2 * r^2 * s^2) * besselj(1, tpirs)) / (8 * pi^4 * r^5 * s^5)
        # now evaluate the limit in the case this is NaN
        if r == 0 || s == 0
            thetaint = 5 * pi / 16
        end
    elseif n == 5
        tpirs = 2 * pi * r * s
        thetaint = (-tpirs * cos(tpirs) + sin(tpirs))/(2 * pi^3 * r^3 * s^3)
        if r == 0 || s == 0
            thetaint = 4/3
        end
    elseif n == 4
        thetaint = besselj(1, 2 * pi * r * s) / (2 * r * s)
        if r == 0 || s == 0
            thetaint = pi / 2
        end
    elseif n == 3
        thetaint = sin(2*pi*r*s)/(pi*r*s)
        if r == 0 || s == 0
            thetaint = 2
        end
    elseif n == 2
        thetaint = pi * besselj(0, 2 * pi * r * s)
        if r == 0 || s == 0
            thetaint = pi
        end
    else
        error("Unimplemented value of n.")
    end
    return omega * basis_factor * exp(-pi * r^2) * r^(n-1) * thetaint
end

"""
We've integrated out the theta dependence
r = |x|, the argument of the model F(r)
d = the order of the derivative we need
s = the frequency argument of the fourier transform
n = the dimension of the problem
j = which of the four bases to use
xi = the left hand endpoint of the interval (used for Hermite basis)
xip1 = x_{i+1}, or the right hand endpoint of the interval
Returns the integrand of the Fourier integral.
"""
function get_fourier_deriv_integrand(r, d, s, n, j, xi, xip1, omega)
    if d != 1
        error("Fourier integrand derivative requires d = 1")
    end
    thetaint = 0
    basis_factor = get_hermite_basis(r, j, xi, xip1)
    if n == 8
        tpirs = 2 * pi * r * s
        thetaint = 15 * 1im * (3 * pi * r * s * besselj(1, tpirs) + (-6 + pi^2 * r^2 * s^2) * besselj(2,tpirs)) / (8 * pi^4 * r^5 * s^5)
        if s == 0 || r == 0
            thetaint = 0 # in that case the limit is 0
        end
    elseif n == 5
        thetaint = 1im * (6 * pi * r * s * cos(2 * pi * r * s) + (-3 + 4 * pi^2 * r^2 * s^2) * sin(2 * pi * r * s))/(4 * pi^4 * r^4 * s^4)
        if s == 0 || r == 0
            thetaint = 0
        end
    elseif n == 4
        thetaint = -1im * besselj(2, 2 * pi * r * s) / (2 * r * s)
        if s == 0 || r == 0
            thetaint = 0
        end
    elseif n == 3
        thetaint = -1im * (-2 * pi * r * s * cos(2 * pi * r * s) + sin(2 * pi * r * s))/(2 * pi^2 * r^2 * s^2)
        if s == 0 || r == 0
           thetaint = 0
        end
    elseif n == 2
        thetaint = -1im * pi * besselj(1, 2 * pi * r * s)
        if r == 0 || s == 0
            thetaint = 0
        end
    else
        error("Unimplemented value of n.")
    end
    return -2 * pi * 1im * omega * basis_factor * exp(-pi * r^2) * r^n * thetaint
end

"""
This is a fourth derivative bound on the Fourier transform
Additional constant factors are included in the SpherePacking run script.
"""
function fourth_div_bound(r, j, xi, xip1, n)
    basis_factor = get_hermite_basis(r, j, xi, xip1)
    return r^(n+3) * exp(-pi * r^2) * abs(basis_factor) # The integral is of the absolute value
end

"""
Compute the integrals of the above form.
n is the ambient dimension
Fcoefs, Fpcoefs and their integral errors indexed by yindex, xinterval, basis_index
F4bnd indexed by the x interval and the basis index
"""
function populate_fourier_integrals!(n, xs, ys, omega, Fcoefs, Fpcoefs, F4bnd, Fcoefserr, Fpcoefserr)
    for yinterval = 1:size(Fcoefs)[1]
        for xinterval = 1:size(Fcoefs)[2]
            for basis_ind = 1:4
                d = 0
                xi = xs[xinterval]
                xip1 = xs[xinterval + 1]
                f = x -> get_fourier_integrand(x, d, ys[yinterval], n, basis_ind, xi, xip1, omega)
                (v, err) = hquadrature(f, xi, xip1, maxevals = 10^6)
                if err > 1e-8
                    error("The adaptive quadrature prematurely terminated for Fourier transform computation.")
                end
                Fcoefs[yinterval, xinterval, basis_ind] = v
                Fcoefserr[yinterval, xinterval, basis_ind] = err
                d = 1
                fp = x -> get_fourier_deriv_integrand(x, d, ys[yinterval], n, basis_ind, xi, xip1, omega)
                (v, err) = hquadrature(fp, xi, xip1, maxevals = 10^6)
                if err > 1e-8
                    error("The adaptive quadrature prematurely terminated for Fourier transform computation.")
                end
                Fpcoefs[yinterval, xinterval, basis_ind] = v
                Fpcoefserr[yinterval, xinterval, basis_ind] = err
                if yinterval == 1 # bound is independent of what y (or s) interval you are in
                    fp4 = x -> fourth_div_bound(x, basis_ind, xi, xip1, n)
                    (v, err) = hquadrature(fp4, xi, xip1, maxevals = 10^6)
                    if err > 1e-8
                        error("The adaptive quadrature prematurely terminated for Fourier transform computation.")
                    end
                    F4bnd[xinterval, basis_ind] = v
                    # All other constant factors are included outside of this function
                end
            end
        end
    end
end
