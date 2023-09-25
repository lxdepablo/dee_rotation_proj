# Author: Luis X. de Pablo
# E-mail: luis.depablo@colorado.edu
# Last updated: September 25, 2023

# load packages and set working directory ----

library(tidyverse)
library(tictoc)
library(deSolve)
library(network)
library(igraph)
library(intergraph)

# set wd to active script location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load helper functions and ecosystem service functions
source("helper_functions.R")
source("ecosystem_services.R")


# run a simulation on empirical food webs ----

# read in some food webs from node and edge lists
food_webs_raw <- read_food_webs("data/food_webs/historical")

# run the simulation on every food web in list
full_sims <- lapply(1:length(food_webs_raw), function(i){
  sim_disturbance(nodes = food_webs_raw[[i]][[1]], 
                  edges = food_webs_raw[[i]][[2]], 
                  pre_dist_gens = 10, 
                  post_dist_gens = 10, 
                  dist_fun = rand_extinct)
})
names(full_sims) <- names(food_webs_raw)


# analysis ----

## plot biomasses over time ----
# add time column and site column to simulation results
full_sims_plotting <- full_sims
full_sims_plotting <- lapply(1:length(full_sims_plotting), function(i) {
  full_sims_plotting[[i]]$t <- as.numeric(rownames(full_sims_plotting[[i]]))
  full_sims_plotting[[i]]$site <- names(full_sims_plotting)[i]
  return(full_sims_plotting[[i]])
})

# make results long form for plotting
simulation_results_long <- lapply(full_sims_plotting, function(i) {
  # take every n row of array
  #n <- 1
  i %>%
    #slice(which(row_number() %% n == 0)) %>%
    pivot_longer(-c(t, site), names_to = "species_ID", values_to = "biomass")
})

# concatenate all simulation data into one dataframe
all_data <- bind_rows(simulation_results_long)

# draw plots
ggplot(data = all_data, aes(x = t, y = biomass, col = species_ID)) +
  facet_wrap(~site) +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_line()


## calculate relative ascendancies ----
# iterate over node and edge lists and generate adjacency matrices
adjacency_matrices <- lapply(1:length(food_webs_raw), function(i) {
  # get node and edge lists for current food web
  curr_nodes <- food_webs_raw[[i]][[1]]
  curr_edges <- food_webs_raw[[i]][[2]]
  
  # generate adjacency matrix
  curr_adj_mat <- make_adj(curr_nodes, curr_edges)
})
# get site names for adjacency matrices
names(adjacency_matrices) <- names(food_webs_raw)

# calculate relative ascendancies from adjacency matrices
relative_ascendancies <- lapply(1:length(adjacency_matrices), function(i) {
  # get current adjacency matrix
  curr_adj_mat <- adjacency_matrices[[i]]
  # create network object
  curr_network <- network(curr_adj_mat)
  #calculate rel ascendancy
  enaAscendency(curr_network)
})


## calculate ecosystem services ----

# use lapply to loop over all food webs
filtration_scores <- lapply(1:length(full_sims), function(i) {
  # apply again to loop over every row of time series dataframe
  filtration_t_series <- lapply(1:nrow(full_sims[[i]]), function(j) {
    # get current row from dataframe
    curr_row <- as.numeric(full_sims[[i]][j,])
    # get names from column headers
    names(curr_row) <- colnames(full_sims[[i]])
    
    # calculate score for that row
    curr_filtration <- water_filtration(curr_row, food_webs_raw[[i]][[1]])
  })
})

# make into time series dataframes
filtration_time_series <- lapply(1:length(filtration_scores), function(i) {
  data.frame(t = 1:length(filtration_scores[[i]]), site = names(full_sims)[i], score = unlist(filtration_scores[[i]]))
})

# plot filtration scores for all sites over time
all_filtration_data <- bind_rows(filtration_time_series)

ggplot(data = all_filtration_data, aes(x = t, y = score)) +
  geom_line() +
  facet_wrap(~site)





