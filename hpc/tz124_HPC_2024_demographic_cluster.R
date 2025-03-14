# CMEE 2024 HPC exercises R code pro forma
# For stochastic demographic model cluster run

rm(list=ls()) # good practice 
graphics.off()  

# Load `Demographic.R` and make sure you can call `stochastic_simulation`
source("Demographic.R")  
source("tz124_HPC_2024_main.R")
# Read the iter variable (provided by the HPC cluster)
iter <- as.numeric(Sys.getenv("PBS_ARRAY_INDEX")) 
#Set the random seed to ensure different tasks have different randomness
set.seed(iter)  
# Define the growth matrix (growth_matrix) and the reproduction matrix (reproduction_matrix)
growth_matrix <- matrix(c(0.1, 0.0, 0.0, 0.0,
                          0.5, 0.4, 0.0, 0.0,
                          0.0, 0.4, 0.7, 0.0,
                          0.0, 0.0, 0.25, 0.4), 
                        nrow=4, ncol=4, byrow=TRUE)

reproduction_matrix <- matrix(c(0.0, 0.0, 0.0, 2.6,
                                0.0, 0.0, 0.0, 0.0,
                                0.0, 0.0, 0.0, 0.0,
                                0.0, 0.0, 0.0, 0.0), 
                              nrow=4, ncol=4, byrow=TRUE)

# Computing `projection_matrix`
projection_matrix <- reproduction_matrix + growth_matrix

# Setting `clutch_distribution`
clutch_distribution <- c(0.06, 0.08, 0.13, 0.15, 0.16, 0.18, 0.15, 0.06, 0.03)

# Set `simulation_length` 120 (total number of simulation steps)
simulation_length <- 120

# ***Choose initial conditions** (iter ensures that the 4 cases are evenly distributed)
initial_conditions <- list(
  list(type = "adult_large",  state = state_initialise_adult(4, 100)),  # 100 individuals, adult stage
  list(type = "adult_small",  state = state_initialise_adult(4, 10)),   # 10 individuals, adult stage
  list(type = "spread_large", state = state_initialise_spread(4, 100)), # 100 individuals, evenly distributed
  list(type = "spread_small", state = state_initialise_spread(4, 10))   # 10 individuals, evenly distributed
)

# **Ensure that the 4 initial conditions are evenly distributed**
condition_index <- ((iter - 1) %% 4) + 1  
initial_state <- initial_conditions[[condition_index]]$state  
condition_type <- initial_conditions[[condition_index]]$type  

# **Create a list to store 150 simulations**
simulation_results <- vector("list",150)  

# **Run `stochastic_simulation` 150 times**
for (i in 1:150) {
  simulation_results[[i]] <- stochastic_simulation(
    initial_state, growth_matrix, reproduction_matrix, clutch_distribution, simulation_length
  )
}

# ***Store the results, the file name contains iter**
save(simulation_results, file = paste0("simulation_results_", iter, "_", condition_type, ".rda"))  
