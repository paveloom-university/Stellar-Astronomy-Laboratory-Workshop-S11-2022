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

# Collect them in a vector
θ₀ = [M_B, B_B, M_D, A_D, B_D, M_H, A_H]

# Define the bounds of the parameters

# Bulge
M_Bᵣ = [43.0, 2580.0]
B_Bᵣ = [0.0, 1.5]

# Disk
M_Dᵣ = [2150.0, 12903.0]
A_Dᵣ = [1.0, 10.0]
B_Dᵣ = [0.1, 15.0]

# Halo
M_Hᵣ = [43.0, 12903.0]
A_Hᵣ = [0.1, 30.0]

# Collect them in a vector
θᵣ = [M_Bᵣ, B_Bᵣ, M_Dᵣ, A_Dᵣ, B_Dᵣ, M_Hᵣ, A_Hᵣ]

# Get lower bounds
θₗ = map(v -> first(v), θᵣ)

# Get upper bounds
θᵤ = map(v -> last(v), θᵣ)

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
function phi_dr(r, z, m_b, b_b, m_d, a_d, b_d, m_h, a_h)::F
    return miyamoto_nagai_phi_dr(r, z, m_b, 0.0, b_b) + # Bulge
           miyamoto_nagai_phi_dr(r, z, m_d, a_d, b_d) + # Disk
           navarro_frenk_white_phi_dr(r, z, m_h, a_h)   # Halo
end

"Compute the value of circular velocity according to v_c(R) = √(R * ∂Φ(R, 0)/∂R) [km/s]"
v_c(r; m_b=M_B, b_b=B_B, m_d=M_D, a_d=A_D, b_d=B_D, m_h=M_H, a_h=A_H) =
    √(r * phi_dr(r, 0, m_b, b_b, m_d, a_d, b_d, m_h, a_h) * 100)

println(pad, "> Plotting the rotation curves...")

"Reset the global plot"
function reset()
    Plots.CURRENT_PLOT.nullableplot = nothing
end

"Plot the rotation curve with the specified parameters"
function plot_with(label="Initial"; m_b=M_B, b_b=B_B, m_d=M_D, a_d=A_D, b_d=B_D, m_h=M_H, a_h=A_H, dashed=false)
    plot!(
        (r) -> v_c(r; m_b, b_b, m_d, a_d, b_d, m_h, a_h), r_range;
        label,
        xlabel=L"R \;\, [\mathrm{kpc}]",
        ylabel=L"v_c \; [\mathrm{km \, s^{-1}}]",
        xminorticks=5,
        linestyle=dashed ? :dash : :solid
    )
end

# Place the legend under the rotation curves
default(legend=(0.25, 0.4))

# Plot rotation curves with different variations of masses
reset()
plot_with()
plot_with(L"$ M_b $ +10% off", m_b=M_B + 0.1 * M_B)
plot_with(L"$ M_b $ --10% off", m_b=M_B - 0.1 * M_B)
plot_with(L"$ M_d $ +10% off", m_d=M_D + 0.1 * M_D)
plot_with(L"$ M_d $ --10% off", m_d=M_D - 0.1 * M_D)
plot_with(L"$ M_b $ +10% off, $ M_d $ --10% off", m_b=M_B + 0.1 * M_B, m_d=M_D - 0.1 * M_D)
plot_with(L"$ M_b $ --10% off, $ M_d $ +10% off", m_b=M_B - 0.1 * M_B, m_d=M_D + 0.1 * M_D)

# Save the figure
savefig(joinpath(PLOTS_DIR, "Rotation curves$(POSTFIX).pdf"))

println(pad, "> Fitting to the initial rotation curve...")

# Reduce the mass of the disk by 10%. Try to fit the rotation
# curve to the initial one by varying the mass of the bulge

# Get the rotation curve with initial parameters
v₀ = v_c.(r_range)

"Compute the difference between the rotation curve with the specified
parameters and the rotation curve with initial parameters"
function residuals_with(; m_b=M_B, b_b=B_B, m_d=M_D, a_d=A_D, b_d=B_D, m_h=M_H, a_h=A_H)
    # Compute the rotation curve
    v = v_c.(r_range; m_b, b_b, m_d, a_d, b_d, m_h, a_h)
    # Return the sum of squared differences
    sum((v - v₀) .^ 2)
end

println(pad, pad, "... by varying the mass of the bulge (with minus 10% of the disk mass)")

# Find the optimal parameter via the least squares method
res = optimize(
    θ -> residuals_with(m_b=θ[1], m_d=M_D - 0.1 * M_D),
    [θₗ[1]],
    [θᵤ[1]],
    [θ₀[1]],
    Fminbox(NelderMead()),
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

println(pad, pad, "... by varying the mass of the bulge and its scale parameter (with minus 10% of the disk mass)")

# Do the same, but this time vary the scale parameter of the bulge, too
res = optimize(
    θ -> residuals_with(m_b=θ[1], b_b=θ[2], m_d=M_D - 0.1 * M_D),
    θₗ[1:2],
    θᵤ[1:2],
    θ₀[1:2],
    Fminbox(NelderMead()),
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

println(pad, pad, "... by varying all parameters (off by -10%)")

# Let's reverse it. Try to get the initial parameters
# of the model from the initial rotation curve
res = optimize(
    θ -> residuals_with(m_b=θ[1], b_b=θ[2], m_d=θ[3], a_d=θ[4], b_d=θ[5], m_h=θ[6], a_h=θ[7]),
    θₗ,
    θᵤ,
    map(x -> x - 0.1 * x, θ₀),
    Fminbox(NelderMead()),
    Optim.Options(
        extended_trace=true,
        store_trace=true,
    )
)

# Save the trace to the file
open(joinpath(TRACES_DIR, "Fit via all$(POSTFIX).trace"), "w") do io
    println(io, summary(res))
    println(io, res.trace)
end

# Get the minimizer
θ_optimal = res.minimizer

# Place the legend under the rotation curves
default(legend=(0.25, 0.3))

# Plot the results of the optimizations
reset()
plot_with()
plot_with(L"Fit via $M_b$ with $M_d$ --10% off", m_b=m_b_optimal)
plot_with(L"Fit via $M_b$ and $b_b$ with $M_d$ --10% off", m_b=θ_b_optimal[1], b_b=θ_b_optimal[2])
plot_with(
    "Fit via all, starting --10% off",
    m_b=θ_optimal[1],
    b_b=θ_optimal[2],
    m_d=θ_optimal[3],
    a_d=θ_optimal[4],
    b_d=θ_optimal[5],
    m_h=θ_optimal[6],
    a_h=θ_optimal[7],
    dashed=true,
)

# Save the figure
savefig(joinpath(PLOTS_DIR, "Fit rotation curves$(POSTFIX).pdf"))

println()
