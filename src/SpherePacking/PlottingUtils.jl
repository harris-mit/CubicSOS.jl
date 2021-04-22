# Utils for plotting what the functions look like given the output

"""
Find which bucket x is in
"""
function x2idx(x, xs)
    for j = 1:length(xs)
        if x >= xs[j] && j == (length(xs) - 1)
            return j
        elseif x >= xs[j] && x < xs[j+1]
            return j
        end
    end
end


"""
Evaluate at x
Mats Q1s, Q2s,
List of intervals given by xs
"""
function eval_F(x, Q1s, Q2s, xs)
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


"""
Evaluate at x
Mats Q1s, Q2s,
List of intervals given by xs
"""
function eval_Fhat(y, X1s, X2s, ys)
    idx = x2idx(y, ys)
    yi = ys[idx]
    yip1 = ys[idx+1]
    X1 = X1s[idx].value
    X2 = X2s[idx].value
    v = [(y-yi)/(yip1 - yi); (y-yip1)/(yi-yip1)]
    return (y-yi)/(yip1 - yi)* v' * X1 * v + (yip1 - y)/(yip1 - yi) * v' * X2 * v
end

"""
Plotting the functions F and Fhat that are output.
"""
function get_function_plots(Q1s, Q2s, X1s, X2s, xs, ys, num_plot_pts = 100)
    xx = 0:maximum(xs)/num_plot_pts:maximum(xs)
    yy = 0:maximum(ys)/num_plot_pts:maximum(ys)
    F = x -> eval_F(x, Q1s, Q2s, xs)
    Fhat = y -> eval_Fhat(y, X1s, X2s, ys)

    fig, axs = subplots(2,2);
    fig.suptitle("", fontsize=20);
    axs[1,1].plot(xx, F.(xx));
    axs[1,1].set_title("F");
    axs[1,1].plot([rad, rad], axs[1,1].get_ylim(), linestyle = "--");
    axs[2,1].plot(yy, Fhat.(yy));
    axs[2,1].set_title("Fhat");
    xx = rad:(maximum(xs) - rad)/num_plot_pts:maximum(xs);
    axs[1,2].plot(xx, F.(xx));
    axs[1,2].set_title("F on subdomain");
    yy = 9:(maximum(ys) - 10)/num_plot_pts:maximum(ys);
    axs[2,2].plot(yy, Fhat.(yy));
    axs[2,2].set_title("Fhat on subdomain")
    return fig, axs
end
