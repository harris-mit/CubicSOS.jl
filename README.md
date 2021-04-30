#CubicSOS.jl

CubicSOS.jl serves two primary purposes:
<ol>
<li> defines a CubicSpline data type
<li> provides a constraint that the spline is nonnegative over an interval
</ol>

## The CubicSpline type
Suppose we want to interpolate the function $f$. A cubic spline can be instantiated as follows

> `p = CubicSpline(x_vals, y_vals, deriv_vals)`

Here, `x_vals` is a list of the interpolation nodes, `y_vals = f.(x_vals)`, and
`deriv_vals = f'.(x_vals)`. The constructor works with `y_vals` and `deriv_vals` real numbers, but
more frequently we will want to define JuMP variables

> `model = Model()`<br>
> `@variable(model, y_vals[1:length(x_vals)])`<br>
> `@variable(model, deriv_vals[1:length(x_vals)])`

and instantiate `CubicSpline` with those variables. Then, we can use `p` in our optimization problem.

### Basic operations with splines

Suppose we have a real number `a` and two splines `p1` and `p2`. The following operations are defined in the
natural way:
<ol>
<li> `p1 + p2` returns a spline whose values and derivatives have been added
<li> `p1 - p2` returns a spline whose values and derivatives have been subtracted
<li> `p1 + a` returns a spline whose values only have been incremented by `a`.
<li> `a * p1` return a spline where the value and derivatives have been multiplied by `a`.
</ol>

TODO: Evaluation, coefficients

## Nonnegativity constraint of spline
The spline determines a cubic polynomial over each subinterval. If `p` is a spline, then
the number of intervals is `length(p) - 1`. We can constrain that `p` is nonnegative
on the `i`th interval by adding the constraint
> constrain_spline_nonnegative!(model, p, i)

This can be repeated for every interval on which we want to enforce nonnegativity.
It is possible and in practice useful to combine this nonnegativity operation with the
basic operations available. For instance we could impose
| Property      | Code |
| ----------- | ----------- |
| $p_1 \leq 0$  | `constrain_spline_nonnegative!(model, -p1, i)`       |
| $p_1 \geq a$  | `constrain_spline_nonnegative!(model, p1 - a, i) `        |
| $p_1 \geq p_2$      | `constrain_spline_nonnegative!(model, p1 - p2, i) `      |
In this case, it may be useful for $a$ to be a constant or a variable, and one or both of the
splines may have values that are decision variables.
