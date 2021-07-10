using DataFrames
using CSVFiles
using Mimi

include(joinpath(@__DIR__, "..", "src/MimiFAIRv2.jl"))

#Load initial conditions (just need parameter names) and extract gas names.
init_gas_vals     = DataFrame(load(joinpath(@__DIR__, "..", "data", "fair_initial_gas_cycle_conditions_1750.csv"), skiplines_begin=7))
montreal_init     = filter(:gas_group => ==("montreal"), init_gas_vals)
flourinated_init  = filter(:gas_group => ==("flourinated"), init_gas_vals)
aerosol_plus_init = filter(:gas_group => ==("aerosol_plus"), init_gas_vals)

# Sort arrays of initial conditions for multiple gases so they are listed alphabetically.
sort!(montreal_init,     :gas_name)
sort!(flourinated_init,  :gas_name)
sort!(aerosol_plus_init, :gas_name)

# Set SSP (options = ['ssp119','ssp126','ssp245','ssp370','ssp585'])
ssp = "ssp370"

# Load replication emissions and forcing from Python.
python_emiss   = DataFrame(load(joinpath(@__DIR__, "..", "data", "python_replication_data", ssp*"_emissions.csv")))
python_exog_RF = DataFrame(load(joinpath(@__DIR__, "..", "data", "python_replication_data", ssp*"_forcing.csv")))[:,:External]

# Extract emissions arrays for multi-gas groupings.
montreal_emissions     = python_emiss[:, Symbol.(montreal_init.gas_name)]
flourinated_emissions  = python_emiss[:, Symbol.(flourinated_init.gas_name)]
aerosol_plus_emissions = python_emiss[:, Symbol.(aerosol_plus_init.gas_name)]

# Get a model instance for that SSP.
m = MimiFAIRv2.get_model(emissions_forcing_scenario=ssp, start_year=1750, end_year=2100)

# Update exogenous forcing and emissions data.
update_param!(m, :radiative_forcing, :exogenous_RF, python_exog_RF)
update_param!(m, :co2_cycle, :E_co2, python_emiss.carbon_dioxide)
update_param!(m, :ch4_cycle, :E_ch4, python_emiss.methane)
update_param!(m, :n2o_cycle, :E_n2o, python_emiss.nitrous_oxide)
update_param!(m, :flourinated_cycles, :E_flourinated, Array(flourinated_emissions))
update_param!(m, :montreal_cycles, :E_montreal, Array(montreal_emissions))
update_param!(m, :aerosol_plus_cycles, :E_aerosol_plus, Array(aerosol_plus_emissions))

run(m)

# Put CO₂ concentration, total radiative forcing, and global temperature anomaly into a dataframe.
results = DataFrame(co2 = m[:co2_cycle, :co2], rf = m[:radiative_forcing, :total_RF], temp = m[:temperature, :T])

# Save results.
save("Mimi_FAIR_replication_"*ssp*".csv", results)



#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------
# PLOTTING CODE IN R
#------------------------------------------------------
#------------------------------------------------------
#------------------------------------------------------

library(ggplot2)

setwd("C:/Users/fce/Desktop/fair2.0/MimiFAIRv2.0.jl")

# Load Mimi results
mimi_ssp119 = read.csv("C:/Users/fce/Desktop/fair2.0/MimiFAIRv2.0.jl/DONT PUT ON GIT - FAIR REPLICATION OF PYTHON RESULTS/Mimi_FAIR_replication_ssp119.csv")
mimi_ssp126 = read.csv("C:/Users/fce/Desktop/fair2.0/MimiFAIRv2.0.jl/DONT PUT ON GIT - FAIR REPLICATION OF PYTHON RESULTS/Mimi_FAIR_replication_ssp126.csv")
mimi_ssp245 = read.csv("C:/Users/fce/Desktop/fair2.0/MimiFAIRv2.0.jl/DONT PUT ON GIT - FAIR REPLICATION OF PYTHON RESULTS/Mimi_FAIR_replication_ssp245.csv")
mimi_ssp370 = read.csv("C:/Users/fce/Desktop/fair2.0/MimiFAIRv2.0.jl/DONT PUT ON GIT - FAIR REPLICATION OF PYTHON RESULTS/Mimi_FAIR_replication_ssp370.csv")
mimi_ssp585 = read.csv("C:/Users/fce/Desktop/fair2.0/MimiFAIRv2.0.jl/DONT PUT ON GIT - FAIR REPLICATION OF PYTHON RESULTS/Mimi_FAIR_replication_ssp585.csv")

# Load Python temperature results
python_ssp119 = read.csv("C:/Users/fce/Desktop/fair2.0/MimiFAIRv2.0.jl/data/python_replication_data/ssp119_temperature.csv")
python_ssp126 = read.csv("C:/Users/fce/Desktop/fair2.0/MimiFAIRv2.0.jl/data/python_replication_data/ssp126_temperature.csv")
python_ssp245 = read.csv("C:/Users/fce/Desktop/fair2.0/MimiFAIRv2.0.jl/data/python_replication_data/ssp245_temperature.csv")
python_ssp370 = read.csv("C:/Users/fce/Desktop/fair2.0/MimiFAIRv2.0.jl/data/python_replication_data/ssp370_temperature.csv")
python_ssp585 = read.csv("C:/Users/fce/Desktop/fair2.0/MimiFAIRv2.0.jl/data/python_replication_data/ssp585_temperature.csv")

# Put data into a dataframe.
mimi_data   = data.frame(year=1750:2100, ssp119=mimi_ssp119[,3], ssp126=mimi_ssp126[,3], ssp245=mimi_ssp245[,3], ssp370=mimi_ssp370[,3], ssp585=mimi_ssp585[,3])
python_data = data.frame(year=1750:2100, ssp119=python_ssp119[,2], ssp126=python_ssp126[,2], ssp245=python_ssp245[,2], ssp370=python_ssp370[,2], ssp585=python_ssp585[,2])

# Crop it to all start in the year 2000 (then do historical in black).
index_2000 = which(1750:2100 == 2000)

short_mimi_data = mimi_data[(index_2000+1):351, ]
short_python_data = python_data[(index_2000+1):351, ]

# Get historicla data
historic_mimi_data = mimi_data[1:(index_2000+1), ]
historic_python_data = python_data[1:(index_2000+1), ]

# Thin python data to every 5-years and to shortened time frame
thin_indices = seq(1,length(short_python_data[,1]), by=3)
thin_indices_historic = seq(1,length(historic_python_data[,1]), by=1)

python_data_thinned = short_python_data[thin_indices, ]
historic_python_data_thinned = historic_python_data[thin_indices_historic, ]

ssp119_col = "#6ace85"
ssp126_col = "#7f9be7"
ssp245_col = "#9c4ed2"
ssp370_col = "darkorange"
#ssp370_col = "#9c4ed2"
ssp585_col = "#ea3675"


p = ggplot()

p = p + geom_hline(yintercept=0, color="red", linetype="22")

p = p + geom_line(data=historic_mimi_data, aes(x=year, y=ssp119), color="black")
p = p + geom_point(data=historic_python_data_thinned, aes(x=year, y=ssp119), shape=21, size=1)

p = p + geom_line(data=short_mimi_data, aes(x=year, y=ssp119), color=ssp119_col)
p = p + geom_point(data=python_data_thinned, aes(x=year, y=ssp119), fill=ssp119_col, shape=21, size=1)

p = p + geom_line(data=short_mimi_data, aes(x=year, y=ssp126), color=ssp126_col)
p = p + geom_point(data=python_data_thinned, aes(x=year, y=ssp126), fill=ssp126_col, shape=21, size=1)

p = p + geom_line(data=short_mimi_data, aes(x=year, y=ssp245), color=ssp245_col)
p = p + geom_point(data=python_data_thinned, aes(x=year, y=ssp245), fill=ssp245_col, shape=21, size=1)

p = p + geom_line(data=short_mimi_data, aes(x=year, y=ssp370), color=ssp370_col)
p = p + geom_point(data=python_data_thinned, aes(x=year, y=ssp370), fill=ssp370_col, shape=21, size=1)

p = p + geom_line(data=short_mimi_data, aes(x=year, y=ssp585), color=ssp585_col)
p = p + geom_point(data=python_data_thinned, aes(x=year, y=ssp585), fill=ssp585_col, shape=21, size=1)

p = p + xlab("Year")
p = p + ylab("Global Surface Temperature Anomaly (K)")
#p = p + labs(title = "Python vs. Julia-Mimi versions of FaIR v2.0")
p = p + labs(title = "Python vs. Julia-Mimi versions of FaIR v2.0", subtitle = "RCMIP Emission & Forcing Scenarios (1750-2100)")

p = p + theme(panel.background = element_rect(fill = "transparent"),
                  plot.background = element_rect(fill = "transparent", color="NA"),
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_line(color="gray90", size=0.25),
                  #axis.line.y = element_blank(),
                  axis.line = element_line(color="black", size=0.25),
                  axis.ticks = element_line(color="black", size=0.25),
                  axis.text = element_text(size=9, colour="black"),
                  axis.title = element_text(size=9, colour="black"),
                  legend.position="none",
                  #plot.title = element_blank(),
                  #axis.title.y = element_blank(),
                  axis.ticks.length=unit(.1, "cm"))
                  #axis.text.y = element_blank(),
                  #axis.ticks.y = element_blank())

ggsave(p, file="Python_Mimi_FAIR2_temperature_comparison.jpg", device="jpeg", type="cairo", width=200, height=130, unit="mm", dpi=200)
ggsave(p, file="Python_Mimi_FAIR2_temperature_comparison.pdf", device="pdf", width=200, height=130, unit="mm", dpi=200, useDingbats=FALSE)

