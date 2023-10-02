# Author: Luis X. de Pablo
# E-mail: luis.depablo@colorado.edu
# Last updated: October 2, 2023

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

# run simulation on all sites 100 times
all_sims <- lapply(1:10, function(j) {
  # run the simulation on every food web in list
  pre_post_states <- lapply(1:length(food_webs_raw), function(i){
    # define number of generations to run simulation before and after disturbance
    pre_d_g <- 200
    post_d_g <- 200
    
    curr_sim <- sim_disturbance(nodes = food_webs_raw[[i]][[1]], 
                                edges = food_webs_raw[[i]][[2]], 
                                pre_dist_gens = pre_d_g, 
                                post_dist_gens = post_d_g, 
                                dist_fun = rand_extinct)
    # extract pre disturbance state and final state
    pre_dist_state <- curr_sim[pre_d_g,]
    final_state <- curr_sim[nrow(curr_sim),]
    rbind(pre_dist_state, final_state)
  })
  names(pre_post_states) <- names(food_webs_raw)
  pre_post_states
})

# loop over sim results and write them to CSV's
lapply(1:length(all_sims), function(i) {
  lapply(1:length(all_sims[[i]]), function(j) {
    # create new directory to save files to
    dir.create(file.path("data", "sim_results", names(all_sims[[i]][j])), showWarnings = FALSE)
    # set working directory
    setwd(file.path("data", "sim_results", names(all_sims[[i]][j])))
    # build filename
    fn <- paste0(names(all_sims[[i]][j]), "_", i, ".csv")
    # write sim results to csv files
    write_csv(all_sims[[i]][[j]], fn)
    # reset working directory
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  })
})



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
names(relative_ascendancies) <- names(food_webs_raw)

## calculate ecosystem services ----

# reshape sim data to be nicer to handle
clean_sims <- lapply(1:length(all_sims[[1]]), function(i) {
  curr_site <- names(all_sims[[1]][i])
  lapply(1:length(all_sims), function(j) {
    all_sims[[j]][[curr_site]]
  })
})
names(clean_sims) <- names(food_webs_raw)

# get ES data
# make vector with ecosystem service functions we want to run
es_functions <- c(total_carbon, water_filtration, fisheries, foundation_score)

# apply over list of ES functions
all_delta_eses <- lapply(es_functions, function(curr_fun) {
  curr_es_data <- lapply(1:length(clean_sims), function(j) {
    # calculate change in ecosystem service for all trials for a site
    delta_eses <- lapply(1:length(clean_sims[[j]]), function(i) {
      
      # get pre and post disturbance states
      pd_state <- clean_sims[[j]][[i]][1,]
      f_state <- clean_sims[[j]][[i]][2,]
      # get node list
      node_list <- food_webs_raw[[j]][[1]]
      # calculate change in ecosystem service
      delta_es <- calc_delta_es(pd_state, f_state, curr_fun, node_list)
    })
  })
  names(curr_es_data) <- names(food_webs_raw)
  curr_es_data
})
names(all_delta_eses) <- c("total_carbon", "water_filtration", "fisheries", "foundation_score")


## model change in es as a function of rel ascendancy ----

# pull es data out of nested lists and combine into one dataframe
es_dataframe <- lapply(1:length(all_delta_eses), function(h) {
  curr_service_name <- names(all_delta_eses)[h]
  curr_service_df <- lapply(1:length(all_delta_eses[[h]]), function(i) {
    # get name of current site
    curr_site_name <- names(all_delta_eses[[h]])[i]
    # apply over all trials for current site
    curr_site_values <- lapply(1:length(all_delta_eses[[h]][[i]]), function(j) {
      trial_num <- j
      # get delta es value for current trial
      curr_trial_df <- data.frame(delta_es = all_delta_eses[[h]][[i]][[j]])
      # add columns for site and trial number
      curr_trial_df$trial <- trial_num
      curr_trial_df$site <- curr_site_name
      curr_trial_df$service <- curr_service_name
      curr_trial_df
    })
    # combine into one dataframe
    bind_rows(curr_site_values)
  })
  bind_rows(curr_service_df)
})
es_dataframe <- bind_rows(es_dataframe)

# pull relative ascendancy data into a df as well
rel_asc_df <- data.frame(site = names(relative_ascendancies),
                         rel_asc = unlist(relative_ascendancies))
# combine es data and rel asc into one dataframe
model_data <- es_dataframe %>%
  left_join(rel_asc_df, by = "site")

ggplot(data = model_data, aes(x = rel_asc, y = delta_es, col = site)) +
  geom_point() +
  facet_wrap(~service, labeller = labeller(service = c(fisheries = "Fisheries",
                                                   foundation_score = "Foundation Species",
                                                   total_carbon = "Total Carbon",
                                                   water_filtration = "Water Filtration"))) +
  labs(x = "Relative Ascendancy", y = paste0(expression(Delta), " ES"))








