# #----------------------------------------------------------------------------------------------------------------------
# #----------------------------------------------------------------------------------------------------------------------
# This file contains functions and other snippets of code that are used in various calculations for Mimi-FAIRv2.0.
# #----------------------------------------------------------------------------------------------------------------------
# #----------------------------------------------------------------------------------------------------------------------


#######################################################################################################################
# CALCULATE CLIMATE & THERMAL RESPONSE PARAMETERS
#######################################################################################################################
# Description: From original Python FAIR code, "This function returns the FaIRv2.0.0-alpha default climate parameters.
#              Use the kwargs to specify pre-determined climate sensitivities. In both cases, the response timescales d1-3
#              (and the shortest-timescale coefficient, q1) are set to the central estimate of a CMIP6 inferred distribution
#              constrained with observational warming. The constraint does not significantly affect the central estimates of
#              the prior (ie. raw CMIP6 inference) distribution."
#
# Function Arguments:
#
#       TCR: Transient climate response (K).
#       RWF: Realized warming fraction (ratio of TCR/ECS).
#       F2x: Radiative forcing from a doubling of carbon dioxide concentrations (Wm⁻²).
#----------------------------------------------------------------------------------------------------------------------

function get_thermal_parameter_defaults(;TCR::Float64=1.79, RWF::Float64=0.552, F2x::Float64=3.759)

    # Set default values for d1, d2, d3, and q1 parameters.
    d1 = 0.903
    d2 = 7.92
    d3 = 355
    q1 = 0.180

    # Calculate equilibrium climate sensitivity based on user-specified TCR and RWF values.
    ECS = TCR/RWF

    # Intermediate calculations.
    v1 = (1-(d1/69.66) * (1-exp(-69.66/d1)))
    v2 = (1-(d2/69.66) * (1-exp(-69.66/d2)))
    v3 = (1-(d3/69.66) * (1-exp(-69.66/d3)))

    # Calculate equilibrium response of second (q2) and third (q3) thermal boxes.
    q3 = (((TCR/F2x) - q1*(v1-v2) - (ECS/F2x)*v2) / (v3-v2))
    q2 = (ECS/F2x - q1 -  q3)

    # Return thermal response parameters as a dataframe.
    df = DataFrame(d1=d1, d2=d2, d3=d3, q1=q1, q2=q2, q3=q3, v1=v1, v2=v2, v3=v3, tcr=TCR, rwf=RWF, F2x=F2x, ECS=ECS)

    return df
end



#######################################################################################################################
# CALCULATE CONSTANTS TO APPROXIMATE STATE-DEPENDENT TIMESCALE ADJUSTMENT FACTOR
#######################################################################################################################
# Description: This function calculates constants for estimating the state-dependent adjustment coefficient of a
#              reservior's lifetime (α) such that it approximates the Millar et al. (2017) numerical solution for a
#              iIRF100 carbon cycle parameterization at α=1.
#
# Function Arguments:
#
#       a: Fraction of emissions entering iᵗʰ atmospheric pool (by default, FAIR has four pools).
#       τ: Atmospheric lifetime of gas in iᵗʰ pool.
#----------------------------------------------------------------------------------------------------------------------

# Function version for a single gas (where a and τ are vectors).
function calculate_g0_g1(a::Array{Float64,1}, τ::Array{Float64,1})
    g1 = sum(a .* τ .* (1.0 .- (1.0 .+ 100.0 ./ τ) .* exp.(-100.0 ./ τ)))
    g0 = exp(-1.0 * sum(a .* τ .* (1.0 .- exp.(-100.0 ./ τ))) / g1)
    return g0, g1
end

# Function version for multiple gases (where a and τ are 2-d arrays).
function calculate_g0_g1(a::Array{Float64,2}, τ::Array{Float64,2})
    g1 = vec(sum(a .* τ .* (1.0 .- (1.0 .+ 100.0 ./ τ) .* exp.(-100.0 ./ τ)), dims=2))
    g0 = vec(exp.(-1.0 .* sum(a .* τ .* (1.0 .- exp.(-100.0 ./ τ)), dims=2) ./ g1))
    return g0, g1
end



#######################################################################################################################
# CALCULATE RADIATIVE FORCING
#######################################################################################################################
# Description: This function calculates the effective radiative forcing based on the atmospheric concentration of a 
#              greenhouse gas or other forcing agent. It follows Equation (6) in Leach et al. (2021).
#
# Function Arguments:
#
#       f:    A vector of three concentration-forcing coefficients (1 = logarithmic, 2 = linear, 3 = square root).
#       C:    Concentration of forcing agent in atmosphere.
#       C_pi: Pre-industrial concentration of forcing agent in atmosphere.
#----------------------------------------------------------------------------------------------------------------------

function calculate_RF(f, C, C_pi)

    if C <= 0.0
        # If concentration is negative, set log and square-root concentrations to 0.0. This follows the original Python FAIRv2.0 code.
        forcing = f[1] * 0.0 + f[2] * (C - C_pi) + f[3] * (0.0 - sqrt(C_pi))
    else
        # Otherwise, calculate forcing as the sum of logarithmic, linear, and square-root terms.
        forcing = f[1] * log(C / C_pi) + f[2] * (C - C_pi) + f[3] * (sqrt(C) - sqrt(C_pi))
    end

    return forcing
end
