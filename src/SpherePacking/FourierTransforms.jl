# Compute the Fourier transforms necessary to set up the SOCP.

function get_basis_factor(r, j, xi, xip1)
    # r =  the value to evaluate at ("x")
    # j = which of the four basis elements to use
    # xi = the left hand endpoint of the interval
    # xip1 = x_{i+1}, or the right hand endpoint of the interval
    # Switch to monomial basis for jump!!!!!
    if j == 1
        basis = (r - xi)^3 / (xip1 - xi)^3
    elseif j == 2
        basis = (r - xi)^2 * (r - xip1) / (xip1 - xi)^3
    elseif j == 3
        basis = (r - xi) * (r - xip1)^2 / (xip1 - xi)^3
    else
        basis = (r - xip1)^3 / (xip1 - xi)^3
    end
    return basis
end

function get_fourier_integrand(r, d, s, n, j, xi, xip1)
    # We've integrated out the theta dependence
    # r = |x|, the argument of the model F(r)
    # d = the order of the derivative we need
    # s = the frequency argument of the fourier transform
    # n = the dimension of the problem
    # j = which of the four bases to use
    # xi = the left hand endpoint of the interval (used for Lagrange basis)
    # xip1 = x_{i+1}, or the right hand endpoint of the interval
    # Use the following to pick out which basis to use
    if d != 0
        error("Fourier integrand requires d = 0")
    end
    thetaint = 0
    basis_factor = get_basis_factor(r, j, xi, xip1)
    if n == 8
        tpirs = 2 * pi * r * s
        thetaint = 15*(-tpirs * besselj(0, tpirs) + (2 - pi^2 * r^2 * s^2) * besselj(1, tpirs)) / (8 * pi^4 * r^5 * s^5)
        # now evaluate the limit in the case this is NaN
        if r == 0 || s == 0
            thetaint = 5 * pi / 16
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

function get_fourier_deriv_integrand(r, d, s, n, j, xi, xip1)
    # We've integrated out the theta dependence
    # r = |x|, the argument of the model F(r)
    # d = the order of the derivative we need
    # s = the frequency argument of the fourier transform
    # n = the dimension of the problem
    # j = which of the four bases to use
    # xi = the left hand endpoint of the interval (used for Lagrange basis)
    # xip1 = x_{i+1}, or the right hand endpoint of the interval
    # Use the following to pick out which basis to use
    if d != 1
        error("Fourier integrand derivative requires d = 1")
    end
    thetaint = 0
    basis_factor = get_basis_factor(r, j, xi, xip1)
    if n == 8
        tpirs = 2 * pi * r * s
        thetaint = 15 * 1im * (3 * pi * r * s * besselj(1, tpirs) + (-6 + pi^2 * r^2 * s^2) * besselj(2,tpirs)) / (8 * pi^4 * r^5 * s^5)
        if s == 0 || r == 0
            thetaint = 0 # in that case the limit is 0
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

function fourth_div_bound(r, j, xi, xip1, n)
    # This is a fourth derivative bound on the Fourier transform
    # Additional constant factors are included in the SpherePacking run script.
    basis_factor = get_basis_factor(r, j, xi, xip1)
    return r^(n+3) * exp(-pi * r^2) * abs(basis_factor) # The integral is of the absolute value
end

function populate_fourier_integrals!(n, Fcoefs, Fpcoefs, F4bnd, Fcoefserr, Fpcoefserr)
    # Compute the integrals of the above form.
    # n is the ambient dimension
    # Fcoefs, Fpcoefs and their integral errors indexed by yindex, xinterval, basis_index
    # F4bnd indexed by the x interval and the basis index
    for yinterval = 1:num_ys
        for xinterval = 1:(num_xs - 1)
            for basis_ind = 1:4
                d = 0
                xi = xs[xinterval]
                xip1 = xs[xinterval + 1]
                f = x -> get_fourier_integrand(x, d, ys[yinterval], n, basis_ind, xi, xip1)
                (v, err) = hquadrature(f, xi, xip1, maxevals = 10^6)
                if err > 1e-8
                    error("The adaptive quadrature prematurely terminated for Fourier transform computation.")
                end
                Fcoefs[yinterval, xinterval, basis_ind] = v
                Fcoefserr[yinterval, xinterval, basis_ind] = err
                d = 1
                fp = x -> get_fourier_deriv_integrand(x, d, ys[yinterval], n, basis_ind, xi, xip1)
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
