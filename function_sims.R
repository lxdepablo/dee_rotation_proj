# Author: Luis X. de Pablo
# E-mail: luis.depablo@colorado.edu
# Last updated: September 11, 2023

# load packages and set working directory to script location ----
library(tidyverse)
library(tictoc)
library(deSolve)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# define some useful functions ----
# function to extract the site name from the filename
extract_site_name <- function(fn) {
  # Extract substring between the first two periods
  site_name <- sub("^.*?\\.(.*?)\\..*$", "\\1", fn)
  # replace spaces with underscores
  site_name <- gsub(" ", "_", site_name)
}

# function to read in all node and edge lists in a directory and return a list of pairs for each site
read_food_webs <- function(directory){
  # List all files in the folder ending with "NODES.csv" or "EDGES.csv" and get site names
  node_file_list <- list.files(path = directory, pattern = "NODES.csv$", full.names = TRUE)
  node_site_list <- lapply(node_file_list, function(file) {
    extract_site_name(file)
  })
  edge_file_list <- list.files(path = directory, pattern = "EDGES.csv$", full.names = TRUE)
  edge_site_list <- lapply(edge_file_list, function(file) {
    extract_site_name(file)
  })
  
  # Use lapply to read all files and store them in a list of data frames
  nodes_lists_raw <- lapply(node_file_list, function(file) {
    # Read the CSV file and return it as a data frame
    read.csv(file)
  })
  edges_lists_raw <- lapply(edge_file_list, function(file) {
    # Read the CSV file and return it as a data frame
    read.csv(file)
  })
  # add site names to node and edge lists
  names(nodes_lists_raw) <- node_site_list
  names(edges_lists_raw) <- edge_site_list
  
  # match up node and edge list for each site
  node_edge_pairs <- lapply(1:length(nodes_lists_raw), function(i) {
    curr_site <- names(nodes_lists_raw)[i]
    list(nodes_lists_raw[[curr_site]], edges_lists_raw[[curr_site]])
  })
  # add names to matched list
  names(node_edge_pairs) <- names(nodes_lists_raw)
  
  return(node_edge_pairs)
}

# function to generate adjacency matrix from node and edge lists
make_adj <- function(nodes, edges) {
  # make matrix of zeroes with length and width equal to number of species in node list
  num_species <- nrow(nodes)
  adj_mat <- matrix(0, nrow = num_species, ncol = num_species)
  
  # name each row and column in adjacency matrix with node ID's from node list
  rownames(adj_mat) <- as.character(nodes$Node.ID)
  colnames(adj_mat) <- as.character(nodes$Node.ID)
  
  # use lapply to set matrix cell to 1 for each interaction in edge list
  for(i in 1:nrow(edges)) {
    # get node ID's for sp 1 and 2
    sp1_ID <- as.character(edges$species_1_Node.ID[i])
    sp2_ID <- as.character(edges$species_2_Node.ID[i])
    
    # set matrix[sp1, sp2] to 1
    adj_mat[sp1_ID, sp2_ID] <- 1
  }
  
  return(adj_mat)
}

# function to generate the initial state vector weighted by trophic level
make_state <- function(nodes, producer_biomass) {
  # for each step up in trophic level, biomass should divide by 10
  # use lapply to assign biomasses to each species in node list based on TL
  biomasses <- producer_biomass/(2^(nodes$TL-1))
  names(biomasses) <- nodes$Node.ID
  
  return(biomasses)
}

# placeholder function to calculate the relative ascendancy of a node in a food web
calculate_relative_ascendancy <- function(adjacency_matrix, node_ID) {
   # Extract the specified row and column of the node in the adjacency matrix
  row_of_node <- adjacency_matrix[node_ID, ]
  col_of_node <- adjacency_matrix[, node_ID]
  
  # Calculate in-degree (predators) and out-degree (prey) for the specified node
  in_degree <- sum(row_of_node)
  out_degree <- sum(col_of_node)
  
  # Calculate the relative ascendancy for the node
  relative_ascendancy <- in_degree / (in_degree + out_degree)
  
  return(relative_ascendancy)
}

# function from Henry's ATN script to get change in biomass over time using ATN
get_dB <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    S <- length(state) # get number of species in web
    dBdt <- rep(0,S) # create vector to store dBi for all i
    
    # zero out all species whose biomasses went negative after the previous step
    for (i in 1:S) {
      if (state[i] < 0) {
        state[i] <- 0
      }
    }
    
    # calculate dBi/dt for all i and store in dBdt
    for (i in 1:S) {
      # determine whether species i is basal and set intrinsic growth rate r accordingly
      if (sum(adjacency_matrix[i,]) == 0) {
        r <- 1 # set r to 1 if species i is basal
      }
      else {r <- 0} # otherwise set r to 0
      
      # get resource vector, the species that i eats
      resources <- which(adjacency_matrix[i,] == 1)
      
      # get consumer vector, the species that eat i
      consumers <- which(adjacency_matrix[,i] == 1)
      
      # calculate primary production
      primary_production <- r*state[i]*(1-(state[i]/K))
      
      # calculate metabolic loss
      metabolic_loss <- x*state[i]
      
      # calculate resource_gain
      resource_gain <- 0 # start by setting resource_gain to 0
      for (j in resources) {
        # calculate Fij
        b_sum <- 0
        for (k in resources) { # sum Bk^(1+q) for every k that i eats
          b_sum <- b_sum + (state[k])^(1+q)
        }
        Fij <- ((state[j])^(1+q))/(b_sum + B0^(1+q))
        
        # calculate gain from resource j
        resource_gain_j <- x*y*Fij*state[i]
        
        # add resource_gain_j to resource_gain
        resource_gain <- resource_gain + resource_gain_j
      }
      
      # calculate consumer_loss
      consumer_loss <- 0 # start by setting consumer_loss to 0
      for (j in consumers) {
        # calculate Fji
        j_resources <- which(adjacency_matrix[j,] == 1)
        b_sum <- 0
        for (k in j_resources) { #sum Bk^(1+q) for every k that j eats
          b_sum <- b_sum + (state[k])^(1+q)
        }
        Fji <- ((state[i])^(1+q))/(b_sum + B0^(1+q))
        
        # calculate loss to consumer j
        consumer_loss_j <- x*y*Fji*state[j]/e
        
        # add consumer_loss_j to consumer_loss
        consumer_loss <- consumer_loss + consumer_loss_j
      }
      
      # sum the above components and store in dBdt[i]
      dBdt[i] <- primary_production - metabolic_loss + resource_gain - consumer_loss
    }
    return(list(dBdt))
  })
}

# function to run a simulation for a given site
# input: node list, edge list, initial state vector
# output: dataframe of biomasses for each species at each time step, cols are species, rows are time
run_sim <- function(nodes, edges, state) {
  tic()
  # generate adjacency matrix from node and edge list
  adjacency_matrix <<- make_adj(nodes, edges)
  # set parameters (use values from Romanuk et al. (2009))
  pars <- c(K=1,x=0.5,y=6,e=1,B0=0.5,q=0.2)
  # set time vector times = (1,...,2000)
  times <- seq(1,100,1)
  
  # run simulation on web and store output as dataframe
  simulation_df <- as.data.frame(ode(state, times, get_dB, pars, method = "ode45"))
  simulation_df <- simulation_df[, 2:ncol(simulation_df)]
  # get species ID's from adjacency matrix
  species_IDs <- colnames(adjacency_matrix)
  colnames(simulation_df) <- species_IDs
  toc()
  
  return(simulation_df)
}


# run a simulation on empirical food webs ----

# read in some food webs from node and edge lists
food_webs_raw <- read_food_webs("data/food_webs")

# run the simulation on every food web in list
simulation_results <- lapply(1:length(food_webs_raw), function(i){
  # get node and edge lists for current food web
  curr_nodes <- food_webs_raw[[i]][[1]]
  curr_edges <- food_webs_raw[[i]][[2]]
  # generate initial state vector
  curr_state <- make_state(curr_nodes, producer_biomass = 10)
  # run simulation
  run_sim(nodes = curr_nodes, edges = curr_edges, curr_state)
})

# get site names for food webs
names(simulation_results) <- names(food_webs_raw)

# plot biomasses over time
# add time column to simulation results
simulation_results <- lapply(1:length(simulation_results), function(i) {
  simulation_results[[i]]$t <- as.numeric(rownames(simulation_results[[i]]))
  simulation_results[[i]]$site <- names(simulation_results)[i]
  return(simulation_results[[i]])
})

# make results long form for plotting
simulation_results_long <- lapply(simulation_results, function(i) {
  n <- 1
  i %>%
    slice(which(row_number() %% n == 0)) %>%
    pivot_longer(-c(t, site), names_to = "species_ID", values_to = "biomass")
})

# concatenate all simulation data into one dataframe
all_data <- bind_rows(simulation_results_long)

# draw plots
ggplot(data = all_data, aes(x = t, y = biomass, col = species_ID)) +
  facet_wrap(~site) +
  theme(legend.position = "none") +
  geom_line()


# define some ecosystem services ----

# example ecosystem service function to calculate fisheries score
# input: named vector of biomasses at a time point where names are Node ID's, node and edge lists
#   for food web of interest
# output: numeric score (total fish biomass)
evaluate_fisheries <- function(biomasses, nodes) {
  # get all the fish in the food web from node list
  fish <- nodes %>%
    filter(Class == "Actinopteri")
  # subset biomass list to include only fish
  biomasses <- biomasses[names(biomasses) %in% fish$Node.ID]
  # fisheries score is the sum of all fish biomasses (total fish biomass)
  score <- sum(biomasses)
}

# get biomasses at last time step
final_biomasses <- c(unlist(slice(simulation_df, nrow(simulation_df))))
# evaluate fisheries score on final biomasses
fisheries_score <- evaluate_fisheries(final_biomasses, nodes = food_webs_raw[[1]][[1]])

# use lapply to calculate fisheries score for every time step
fisheries_t_series <- data.frame(
  score = apply(simulation_df, 1, function(i) {
    evaluate_fisheries(i, nodes = food_webs_raw[[1]][[1]])}),
  t = 1:100
  )

# plot fisheries over time
ggplot(data = fisheries_t_series, aes(x = t, y = score)) +
  geom_line()





