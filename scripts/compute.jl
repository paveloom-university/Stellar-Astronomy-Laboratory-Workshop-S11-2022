# This script computes rotation curves from galaxy potentials
# with known parameters, toys with the parameters a bit to see
# the effects, then tries to infer these parameters back

# Masses are in units of M_0 = 2.325 * 10^7 M_Sun, other parameters
# are in kpc. Initial values are taken from Bajkova, Bobylev (2020),
# ranges are taken from Granados et al. (2021)

"Parse the string, taking more arguments if it's quoted"
function parse_string(i)::String
    # Start from the first argument after the flag
    j = i
    # If the arguments starts with an apostrophe,
    s = if startswith(ARGS[i], "'")
        # Find the index of the argument
        # which ends with an apostrophe
        while !endswith(ARGS[j], "'")
            j += 1
        end
        # Join the arguments in one string
        # and remove the apostrophes
        chop(join(ARGS[i:j], ' '), head=1, tail=1)
    else
        # Return the next argument
        ARGS[i]
    end
    return s
end

# Define default values for optional arguments
POSTFIX = ""

# Parse the options
for i in eachindex(ARGS)
    # A postfix for the names of output files
    if ARGS[i] == "--postfix"
        try
            global POSTFIX = " ($(parse_string(i+1)))"
        catch
            println("Couldn't parse the value of the `--postfix` argument.")
            exit(1)
        end
    end
end

# Prepare color codes
RESET = "\e[0m"
GREEN = "\e[32m"
YELLOW = "\e[33m"

# Check for required arguments
if "--help" in ARGS
    println("""
        $(YELLOW)USAGE:$(RESET)
        { julia --project=. | ./julia.bash } scripts/compute.jl [--postfix <POSTFIX>]

        $(YELLOW)OPTIONS:$(RESET)
            $(GREEN)--postfix <POSTFIX>$(RESET)    A postfix for the names of output files"""
    )
    exit(1)
end

"Padding in the output"
pad = " "^4

"Floating point type used across the script"
F = Float64

"Integer type used across the script"
I = UInt64

# Define the paths
CURRENT_DIR = @__DIR__
ROOT_DIR = basename(CURRENT_DIR) == "scripts" ? dirname(CURRENT_DIR) : CURRENT_DIR
PLOTS_DIR = joinpath(ROOT_DIR, "plots")
TRACES_DIR = joinpath(ROOT_DIR, "traces")

# Make sure the needed directories exist
mkpath(PLOTS_DIR)
mkpath(TRACES_DIR)

println('\n', " "^4, "> Loading the packages...")

using LaTeXStrings
using Optim
using Plots

# Use the PGFPlotsX backend for plots
pgfplotsx()

# Change some of the default parameters for plots
default(fontfamily="Computer Modern", dpi=300)

# Define the parameters of the Galaxy potential model

# Bulge (Plummer potential, which is a reduced version of the Miyamoto & Nagai potential)
M_B = 443.0
B_B = 0.2672

# Disk (Miyamoto & Nagai potential)
M_D = 2798.0
A_D = 4.40
B_D = 0.3084

# Halo (Navarro-Frenk-White potential)
M_H = 12474.0
A_H = 7.7

# Define the bounds of the parameters

# Bulge
M_B_bounds = [43.0, 2580.0]
B_B_bounds = [0.0, 1.5]

# Disk
M_D_bounds = [2150.0, 12903.0]
A_D_bounds = [1.0, 10.0]
B_D_bounds = [0.1, 15.0]

# Halo
M_H_bounds = [43.0, 12903.0]
A_H_bounds = [0.1, 30.0]

bounds = [M_B_bounds, B_B_bounds, M_D_bounds, A_D_bounds, B_D_bounds, M_H_bounds, A_H_bounds]

# Get lower bounds
lower_bounds = map(v -> first(v), bounds)

# Get upper bounds
upper_bounds = map(v -> last(v), bounds)

# Define the range of galactocentric distances
r_range = 1:0.1:20

"Compute the value of the R derivative of the Miyamoto & Nagai potential [100 km/s²]"
function miyamoto_nagai_phi_dr(r, z, m, a, b)::F
    return m * r / (r^2 + (a + √(z^2 + b^2))^2)^(1.5)
end

"Compute the value of the R derivative of the Navarro-Frenk-White potential [100 km/s²]"
function navarro_frenk_white_phi_dr(r, z, m, a)::F
    sq_sum = r^2 + z^2
    k_1 = √(sq_sum) / a + 1.0
    return m * r * log(k_1) / sq_sum^(1.5) - m * r / a / sq_sum / k_1
end

"Compute the value of ∂Φ(R, Z)/∂R [100 km/s²]"
function phi_dr(r, z, m_b, b_b, m_d)::F
    return miyamoto_nagai_phi_dr(r, z, m_b, 0.0, b_b) + # Bulge
           miyamoto_nagai_phi_dr(r, z, m_d, A_D, B_D) + # Disk
           navarro_frenk_white_phi_dr(r, z, M_H, A_H)   # Halo
end

"Compute the value of circular velocity according to v_c(R) = √(R * ∂Φ(R, 0)/∂R) [km/s]"
v_c(
    r;
    m_b=M_B,
    b_b=B_B,
    m_d=M_D
) = √(r * phi_dr(r, 0, m_b, b_b, m_d) * 100)

println(pad, "> Plotting the rotation curves...")

"Reset the global plot"
function reset()
    Plots.CURRENT_PLOT.nullableplot = nothing
end

"Plot the rotation curve with the specified parameters"
function plot_with(
    label="Initial";
    m_b::F=M_B,
    b_b::F=B_B,
    m_d::F=M_D
)
    plot!(
        (r) -> v_c(r; m_b, b_b, m_d), r_range;
        label,
        xlabel=L"R \;\, [\mathrm{kpc}]",
        ylabel=L"v_c \; [\mathrm{km \, s^{-1}}]",
        xminorticks=5
    )
end

# Place the legend under the rotation curves
default(legend=(0.25, 0.4))

# Plot rotation curves with different variations of masses
reset()
plot_with()
plot_with(L"M_b + 10%", m_b=M_B + 0.1 * M_B)
plot_with(L"M_b - 10%", m_b=M_B - 0.1 * M_B)
plot_with(L"M_d + 10%", m_d=M_D + 0.1 * M_D)
plot_with(L"M_d - 10%", m_d=M_D - 0.1 * M_D)
plot_with(L"M_b + 10%, M_d - 10%", m_b=M_B + 0.1 * M_B, m_d=M_D - 0.1 * M_D)
plot_with(L"M_b - 10%, M_d + 10%", m_b=M_B - 0.1 * M_B, m_d=M_D + 0.1 * M_D)

# Save the figure
savefig(joinpath(PLOTS_DIR, "Rotation curves$(POSTFIX).pdf"))

println(pad, "> Fitting rotation curves...")

# Reduce the mass of the disk by 10%. Try to fit the rotation
# curve to the initial one by varying the mass of the bulge

# Get the rotation curve with initial parameters
v₀ = v_c.(r_range)

"Compute the difference between the rotation curve with the specified
parameters and the rotation curve with initial parameters"
function residuals_with(m_b::F; b_b::F=B_B)
    # Compute the rotation curve
    v = v_c.(r_range; m_b, b_b, m_d=M_D - 0.1 * M_D)
    # Return the sum of squared differences
    sum((v - v₀) .^ 2)
end

println(pad * pad, "... by varying the mass of the bulge")

# Find the optimal parameter via the least squares method
res = optimize(
    θ -> residuals_with(θ[1]),
    [lower_bounds[1]],
    [upper_bounds[1]],
    [M_B],
    Fminbox(LBFGS()),
    Optim.Options(
        extended_trace=true,
        store_trace=true,
    )
)

# Save the trace to the file
open(joinpath(TRACES_DIR, "Fit via M_b$(POSTFIX).trace"), "w") do io
    println(io, summary(res))
    println(io, res.trace)
end

# Get the minimizer
m_b_optimal = res.minimizer[1]

println(pad * pad, "... by varying the mass of the bulge and its scale parameter")

# Do the same, but this time vary the scale parameter of the bulge, too
res = optimize(
    θ -> residuals_with(θ[1], b_b=θ[2]),
    lower_bounds[1:2],
    upper_bounds[1:2],
    [M_B, B_B],
    Fminbox(LBFGS()),
    Optim.Options(
        extended_trace=true,
        store_trace=true,
    )
)

# Save the trace to the file
open(joinpath(TRACES_DIR, "Fit via M_b and b_b$(POSTFIX).trace"), "w") do io
    println(io, summary(res))
    println(io, res.trace)
end

# Get the minimizer
θ_b_optimal = res.minimizer

# Place the legend in the top right corner
default(legend=:topright)

# Plot the results of the optimizations
reset()
plot_with()
plot_with(L"Fit via $M_b$", m_b=m_b_optimal)
plot_with(L"Fit via $M_b$ and $b_b$", m_b=θ_b_optimal[1], b_b=θ_b_optimal[2])

# Save the figure
savefig(joinpath(PLOTS_DIR, "Fit rotation curves$(POSTFIX).pdf"))

println()
