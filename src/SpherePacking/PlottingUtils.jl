# Utils for plotting what the functions look like given the output

function eval_F(x, Q1s, Q2s, xs)
    # evaluate at x
    # mats Q1s, Q2s,
    # list of intervals given by xs
    if x > maximum(xs)
        return 0
    end
    idx = x2idx(x, xs)
    xi = xs[idx]
    xip1 = xs[idx+1]
    Q1 = Q1s[idx].value
    Q2 = Q2s[idx].value
    v = [(x-xi)/(xip1 - xi); (x-xip1)/(xi-xip1)]
    return -((x-xi)/(xip1 - xi)* v' * Q1 * v + (xip1 - x)/(xip1 - xi) * v' * Q2 * v) * exp(-pi * x^2)
end

function eval_Fhat(y, X1s, X2s, ys)
    # evaluate at x
    # mats Q1s, Q2s,
    # list of intervals given by xs
    idx = x2idx(y, ys)
    yi = ys[idx]
    yip1 = ys[idx+1]
    X1 = X1s[idx].value
    X2 = X2s[idx].value
    v = [(y-yi)/(yip1 - yi); (y-yip1)/(yi-yip1)]
    return (y-yi)/(yip1 - yi)* v' * X1 * v + (yip1 - y)/(yip1 - yi) * v' * X2 * v
end
