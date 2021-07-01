# Utils for plotting what the functions look like given the output
#ENV["MPLBACKEND"]="qt5agg"
#using PyPlot

"""
Plotting the splines F and Fhat that are output.
"""
function get_function_plots(f::CubicSpline, fhat::CubicSpline, xs, ys, rad, num_plot_pts = 200)
    xx = 0:maximum(xs)/num_plot_pts:maximum(xs)
    yy = 0:maximum(ys)/num_plot_pts:maximum(ys)

    fig, axs = subplots(2,2);
    fig.suptitle("", fontsize=20);
    feval = x -> value(evaluate_cubic(f, x)) * exp(-pi * x^2)
    fhateval = x -> value(evaluate_cubic(fhat, x))
    axs[1,1].plot(xx, feval.(xx) );
    axs[1,1].set_title("F");
    axs[1,1].plot([rad, rad], axs[1,1].get_ylim(), linestyle = "--");
    axs[2,1].plot(yy, fhateval.(yy));
    axs[2,1].set_title("Fhat");
    xx = rad:(maximum(xs) - rad)/num_plot_pts:maximum(xs);
    axs[1,2].plot(xx, feval.(xx));
    axs[1,2].set_title("F on subdomain");
    yy = 9:(maximum(ys) - 10)/num_plot_pts:maximum(ys);
    axs[2,2].plot(yy, fhateval.(yy));
    axs[2,2].set_title("Fhat on subdomain")
    return fig, axs
end
