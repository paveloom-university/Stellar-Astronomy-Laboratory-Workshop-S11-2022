# This script computes rotation curves from galaxy potentials
# with known parameters, toys with the parameters a bit to see
# the effects, then tries to infer these parameters back

# Masses are in units of M_0 = 2.325 * 10^7 M_Sun, other parameters
# are in kpc. Initial values are taken from Bajkova, Bobylev (2020).

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

# Make sure the needed directories exist
mkpath(PLOTS_DIR)

println('\n', " "^4, "> Loading the packages...")

using LaTeXStrings
using Plots

# Use the PGFPlotsX backend for plots
pgfplotsx()

# Change some of the default parameters for plots
default(fontfamily="Computer Modern", dpi=300, legend=nothing)

# Define the parameters of the Galaxy potential model

# Bulge (Plummer potential, which is a reduced version of the Miyamoto & Nagai potential)
M_B = 443.0
A_B = 0.0
B_B = 0.2672

# Disk (Miyamoto & Nagai potential)
M_D = 2798.0
A_D = 4.40
B_D = 0.3084

# Halo (Navarro-Frenk-White potential)
M_H = 12474.0
A_H = 7.7

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
function phi_dr(r, z)::F
    return miyamoto_nagai_phi_dr(r, z, M_B, A_B, B_B) + # Bulge
           miyamoto_nagai_phi_dr(r, z, M_D, A_D, B_D) + # Disk
           navarro_frenk_white_phi_dr(r, z, M_H, A_H)   # Halo
end

"Compute the value of circular velocity according to v_c(R) = √(R * ∂Φ(R, 0)/∂R) [km/s]"
v_c(r) = √(r * phi_dr(r, 0) * 100)

println(" "^4, "> Plotting the rotation curves...")

# Plot the rotation curve
plot(
    v_c, 1, 20;
    xlabel=L"R \;\, [\mathrm{kpc}]",
    ylabel=L"v_c \; [\mathrm{km \, s^{-1}}]",
    xminorticks=5
)

# Save the figure
savefig(joinpath(PLOTS_DIR, "Rotation curve$(POSTFIX).pdf"))
savefig(joinpath(PLOTS_DIR, "Rotation curve$(POSTFIX).png"))

println()
