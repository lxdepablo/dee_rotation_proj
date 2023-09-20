# Author: Luis X. de Pablo
# E-mail: luis.depablo@colorado.edu
# Last updated: September 11, 2023

# load packages and set working directory ----

library(tidyverse)
library(tictoc)
library(deSolve)
library(network)
library(igraph)
library(intergraph)
# set wd to active script location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# define some useful functions ----

# function to extract the site name from the filename
extract_site_name <- function(fn) {
  # Extract substring between the first two periods
  site_name <- sub("^[^.]+\\.([^.]+\\.[^.]+)\\..*", "\\1", fn)
  #site_name <- sub("(?:[^_]*_[^_]*_+([^.]+\\.[^.]+)\\.*)", "\\1", fn)
  
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

# helper function to measure the relative ascendancy of a network
as.extended <- function(x,zero.na=TRUE){
  #Check for network class object
  if (class(x) != "network"){warning('x is not a network class object')}
  #unpack the data from the network object
  flow <- as.matrix(x,attrname="flow")
  input <- x%v%'input'
  respiration <- x%v%'respiration'
  export <- x%v%'export'
  #recombine into the extended format
  import <- c(input,0,0,0)
  x <- cbind(flow,export,respiration,rep(0,nrow(flow)))
  x <- rbind(x,rep(0,length(import)),rep(0,length(import)),import)
  #make NA values zero if zero.na == TRUE
  if (zero.na){
    x[is.na(x)] = 0
  }
  return(x)
}

# function to measure the relative ascendancy of a network
enaAscendency <- function(x='network object'){
  if (class(x) != 'network'){warning('x is not a network class object')}
  
  
  #    if (any(is.na(x%v%'export'))){
  #        warning('Export data is absent from the model.')
  #    }
  #    if(any(is.na(x%v%'respiration'))){
  #           warning('Respiration data is absent from the model.')
  #   }
  
  ####### set initial conditions for calculations #########
  T.ulan <- as.extended(x)
  N <- ncol(T.ulan) # set up N
  r.td <- c.ld <- t.ulan <- ami <- mat.or.vec(N,N) # initialize ascendency matrix
  oh <- mat.or.vec(N,N) # initialize overhead matrix
  cap <- mat.or.vec(N,N) # initialize capacity matrix
  #calculate total system throughPUT
  TSTp <- sum(T.ulan)
  
  ## calculate H & CAPACITY  #######################################
  ## H = Total Flow Diversity
  
  h <- T.ulan/sum(T.ulan)
  h2 <- log2(h)
  h2[!is.finite(h2)] <- 0
  H = - sum(h * h2)   # Total Flow Diversity
  
  CAP <- H * TSTp     # Capactity
  
  ################### calculate AMI  #######################
  ## AMI = Average Mutual Informaiton
  
  # loop through T.ulan to calculate AMI
  for (i in 1:N){
    for (j in 1:N){
      if (T.ulan[i,j] == 0){
        ami[i,j] <- 0
      }else{
        ami[i,j] <- T.ulan[i,j]/TSTp * log2((T.ulan[i,j]*TSTp)/(sum(T.ulan[,j])*sum(T.ulan[i,])))
      }
    }
  }
  
  AMI <- sum(ami)
  
  ################ calculate ascendency and ascendency/capacity ###############
  
  ASC <- TSTp * AMI
  ASC.CAP <- ASC/CAP
  
  return(ASC.CAP)
}

# function from Henry's ATN script to get change in biomass over time using ATN
vectorized_get_dB <- function(t, state, parameters) {
  # extract parameters from list of model parameters
  r <- parameters$r
  K <- parameters$K
  x <- parameters$x
  y <- parameters$y
  e <- parameters$e
  B0 <- parameters$B0
  q <- parameters$q
  
  # get number of species in web
  S <- nrow(adjacency_matrix)
  
  # if a species' biomass fell below the extinction threshold, set its biomass to 0
  threshold <- 1*10^(-10)
  state <- ifelse(state < threshold, 0, state)
  
  # calculate primary production vector
  primary_production <- r*state*(1-state/K)
  
  # Calculate metabolic loss vector
  metabolic_loss <- x*state
  
  # calculate resource gain and consumer loss vectors
  resource_gain <- rep(0, S) # initiate a vector to store resource gain for each species
  consumer_loss <- rep(0, S) # initiate a vector to store consumer loss for each species
  
  for (i in 1:S) {
    # calculate resource_gain
    b_sum_i <- as.vector(state^(1+q) %*% adjacency_matrix[i,])
    F_ij <- (state^(1+q))/(b_sum_i + B0[i,]^(1+q))
    resource_gain[i] <- (x[i]*y[i,]*F_ij*state[i]) %*% adjacency_matrix[i,]
    
    # calculate consumer_loss
    b_sum_j <- t(adjacency_matrix %*% state^(1+q))
    F_ji <- state[i]^(1+q)/(b_sum_j + B0[,i]^(1+q))
    consumer_loss[i] <- (x*y[,i]*F_ji*state/e[,i]) %*% adjacency_matrix[,i]
  }
  
  # sum the above components and store in dBdt
  dBdt <- primary_production - metabolic_loss + resource_gain - consumer_loss
  
  return(list(dBdt))
}

# function to run a simulation for a given site
# input: node list, edge list, number of generations to run, initial state vector
# output: dataframe of biomasses for each species at each time step, cols are species, rows are time
vectorized_run_sim <- function(nodes, edges, num_generations, state) {
  tic()
  # nodes <- food_webs_raw[[1]][[1]]
  # edges <- food_webs_raw[[1]][[2]]
  # state <- make_state(nodes, 10)
  # num_generations <- 100
  
  # generate adjacency matrix from node and edge list
  adjacency_matrix <<- make_adj(nodes, edges)
  S <- nrow(adjacency_matrix)
  # set parameters
  r <- ifelse(rowSums(adjacency_matrix) == 0, 1, 0)
  K <- rep(1, S)
  x <- rep(0.5, S)
  y <- matrix(rep(6, S^2), nrow = S, ncol = S)
  e <- matrix(rep(1, S^2), nrow = S, ncol = S)
  B0 <- matrix(rep(0.5, S^2), nrow = S, ncol = S)
  q <- 0.2
  
  pars <- list(r = r, K = K, x = x, y = y, e = e, B0 = B0, q = q)
  # set time vector times = (1,...,2000)
  times <- seq(1,num_generations,1)
  
  # run simulation on web and store output as dataframe
  simulation_df <- as.data.frame(ode(state, times, vectorized_get_dB, pars, method = "ode45"))
  simulation_df <- simulation_df[, 2:ncol(simulation_df)]
  # get species ID's from adjacency matrix
  species_IDs <- colnames(adjacency_matrix)
  colnames(simulation_df) <- species_IDs
  toc()
  
  return(simulation_df)
}

# function to use a monte carlo simulation to optimize initial conditions for fewest extinctions
monte_carlo_state <- function(num_iterations, nodes, edges) {
  # get number of species in food web
  num_species <- nrow(nodes)
  
  # use lapply to generate random initial states and test each one on 25 generations of the ATN
  possible_states <- lapply(1:num_iterations, function(i) {
    # generate a random state
    rand_state <- runif(num_species)
    # test that state on the ATN
    curr_results <- vectorized_run_sim(nodes, edges, 100, rand_state)
    # get final biomasses
    final_state <- curr_results[nrow(curr_results), ]
    # check how many species went extinct
    num_extinct <- length(final_state[final_state < 0.0000000001])
    # return final state and store its score in the last position of the vector
    rand_state <- append(rand_state, num_extinct)
  })
  
  # return only the best state
  best_state <- NULL
  min_extinctions <- Inf
  
  for(i in 1:length(possible_states)) {
    curr_state <- possible_states[[i]]
    curr_score <- curr_state[length(curr_state)]
    if(curr_score < min_extinctions) {
      best_state <- curr_state[-length(curr_state)]
      min_extinctions <- curr_score
    }
  }
  print("FINISHED MONTE CARLO!")
  return(best_state)
}

# function to simulate a random extinction
# input: vector of biomasses
# output: new vector of biomasses with a random sp extinct
rand_extinct <- function(biomasses) {
  rand_sp <- runif(1, min = 1, max = length(biomasses) + 1)
  new_biomasses <- biomasses
  new_biomasses[rand_sp] <- 0
  return(new_biomasses)
}


# define some ecosystem services ----

# function to calculate fisheries score
# input: named vector of biomasses at a time point where names are Node ID's, node and edge lists
#   for food web of interest
# output: numeric score (total fish biomass)
fisheries <- function(biomasses, nodes) {
  # get all the fish in the food web from node list
  fish <- nodes %>%
    filter(Class == "Actinopteri")
  # subset biomass list to include only fish
  biomasses <- as.numeric(biomasses[names(biomasses) %in% fish$Node.ID])
  # fisheries score is the sum of all fish biomasses (total fish biomass)
  score <- sum(biomasses)
}

# function to calculate total carbon
total_carbon <- function(biomasses) {
  score <- sum(biomasses)
}

# function to calculate water filtration score
water_filtration <- function(biomasses, nodes) {
  # get all the filter feeding species
  filter_feeders <- nodes %>%
    filter(Class %in% c("Bivalvia", "Thecostraca"))
  
  # subset biomass list to include only filter feeders
  biomasses <- as.numeric(biomasses[names(biomasses) %in% filter_feeders$Node.ID])
  # fisheries score is the sum of all fish biomasses (total fish biomass)
  score <- sum(biomasses)
}

# function to calculate foundation species score
foundation_score <- function(biomasses, nodes) {
  # get all the foundation species
  foundation_species <- nodes %>%
    filter(Class %in% c("Thecostraca") | Order %in% c("Mytilida") | Genus %in% c("Fucus"))
  
  # subset biomass list to include only filter feeders
  biomasses <- as.numeric(biomasses[names(biomasses) %in% foundation_species$Node.ID])
  # fisheries score is the sum of all fish biomasses (total fish biomass)
  score <- sum(biomasses)
}


# run a simulation on empirical food webs ----

# read in some food webs from node and edge lists
food_webs_raw <- read_food_webs("data/food_webs/historical")

# run the simulation on every food web in list for 1000 generations
first_1000 <- lapply(1:length(food_webs_raw), function(i){
  # get node and edge lists for current food web
  curr_nodes <- food_webs_raw[[i]][[1]]
  curr_edges <- food_webs_raw[[i]][[2]]
  
  # generate initial state vector
  #curr_state <- make_state(curr_nodes, producer_biomass = 10)
  curr_state <- monte_carlo_state(1, curr_nodes, curr_edges)
  
  # run simulation
  vectorized_run_sim(nodes = curr_nodes, edges = curr_edges, num_generations = 50, curr_state)
})
# get site names for food webs
names(first_1000) <- names(food_webs_raw)

# get final states after 1000 generations
post_disturbance_states <- lapply(1:length(first_1000), function(i) {
  curr_food_web <- first_1000[[i]]
  final_state <- curr_food_web[nrow(curr_food_web),]
  new_state <- rand_extinct(final_state)
})

# run the simulation on every food web in list for 1000 more generations
second_1000 <- lapply(1:length(food_webs_raw), function(i){
  # get node and edge lists for current food web
  curr_nodes <- food_webs_raw[[i]][[1]]
  curr_edges <- food_webs_raw[[i]][[2]]
  
  # generate initial state vector
  #curr_state <- make_state(curr_nodes, producer_biomass = 10)
  curr_state <- unlist(post_disturbance_states[i])
  
  # run simulation
  vectorized_run_sim(nodes = curr_nodes, edges = curr_edges, num_generations = 50, curr_state)
})
# get site names for food webs
names(second_1000) <- names(food_webs_raw)

# combine together first 1000 and second 1000 generations
full_sims <- lapply(1:length(first_1000), function(i) {
  # get first 1000 and second 1000 generations for current site
  curr_first_1000 <- first_1000[[i]]
  curr_second_1000 <- second_1000[[i]]
  
  # combine them together into full sim
  full_sim <- rbind(curr_first_1000, curr_second_1000)
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





