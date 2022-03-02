# CubicSOS.jl

CubicSOS.jl serves two primary purposes:

1. defines a CubicInterpolant data type
2. provides a constraint that the interpolant is nonnegative over an interval


## The CubicInterpolant type
Suppose we want to interpolate the function `f`. A cubic interpolant can be instantiated as follows

> `p = CubicInterpolant(x_vals, y_vals, deriv_vals)`

Here, `x_vals` is a list of the interpolation nodes, `y_vals = f.(x_vals)`, and
`deriv_vals = f'.(x_vals)`. The constructor works with `y_vals` and `deriv_vals` real numbers, but
more frequently we will want to define JuMP variables

> `model = Model()`<br>
> `@variable(model, y_vals[1:length(x_vals)])`<br>
> `@variable(model, deriv_vals[1:length(x_vals)])`

and instantiate `CubicInterpolant` with those variables. Then, we can use `p` in our optimization problem.

### Evaluation of interpolant

The most basic operation one can do with an interpolant is to evaluate it at a point on its domain. The function

> `evaluate_cubic(p, x)`,

or more simply

> `p(x)`,

calculates the interpolated value of the interpolant `p` at the point `x`.

### Basic operations with interpolants

Suppose we have a real number `a` and two interpolants `p1` and `p2`. The following operations are defined in the
natural way:

- `p1 + p2` returns an interpolant whose values and derivatives have been added
- `p1 - p2` returns an interpolant whose values and derivatives have been subtracted
- `p1 + a` returns an interpolant whose values only have been incremented by `a`.
- `a * p1` return an interpolant where the value and derivatives have been multiplied by `a`.

If we are adding or subtracting two interpolants that are defined on the same interpolation points, the resulting interpolant is also computed on the same points. The two interpolants must always have the same endpoints. If they are based on different interior points, the result is based on the union of those grids.

## Nonnegativity constraint of interpolant
The interpolant determines a cubic polynomial over each subinterval. If `p` is a interpolant, then
the number of intervals is `length(p) - 1`. We can constrain that `p` is nonnegative
on the `i`th interval by adding the constraint
> `constrain_interpolant_nonnegative!(model, p, i)`

If `i` is an AbstractArray of interval indices, such constraints are added for each interval. Alternatively,

> `constrain_interpolant_nonnegative!(model, p)`

enforces the constraint on the entire interval on which `p` is defined.

It is possible and in practice useful to combine the nonnegativity condition with the
basic operations available. For instance we could impose
| Property      | Code |
| ----------- | ----------- |
| `p1 <= 0` | `constrain_interpolant_nonnegative!(model, -p1, i)`       |
| `p1 >= a`  | `constrain_interpolant_nonnegative!(model, p1 - a, i) `        |
| `p1 >= p2`      | `constrain_interpolant_nonnegative!(model, p1 - p2, i) `      |

In this case, it may be useful for `a` to be a constant or a variable, and one or both of the
interpolants may have values that are decision variables.
