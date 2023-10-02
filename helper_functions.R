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

# define some useful functions ----

# function to extract the site name from the filename
extract_site_name <- function(fn) {
  # Extract substring between the first and third periods
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
  
  return(simulation_df)
}

# function to use a monte carlo simulation to optimize initial conditions for fewest extinctions
monte_carlo_state <- function(num_iterations, nodes, edges) {
  # get number of species in food web
  num_species <- nrow(nodes)
  
  # use lapply to generate random initial states and test each one on 100 generations of the ATN
  possible_states <- lapply(1:num_iterations, function(i) {
    # generate a random state
    rand_state <- runif(num_species)
    # test that state on the ATN
    curr_results <- vectorized_run_sim(nodes, edges, 100, rand_state)
    # get final biomasses
    final_state <- curr_results[nrow(curr_results), ]
    # check how many species went extinct
    num_extinct <- length(final_state[final_state < 1e-10])
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

  return(best_state)
}

# function to generate a random initial state for a simulation
rand_state <- function(nodes) {
  # get number of species in food web
  num_species <- nrow(nodes)
  # generate state
  state <- runif(num_species)
}

# function to simulate a random extinction
# input: vector of biomasses
# output: new vector of biomasses with a random sp extinct
rand_extinct <- function(biomasses) {
  while(TRUE) {
    rand_sp <- runif(1, min = 1, max = length(biomasses) + 1)
    new_biomasses <- biomasses
    if(!(new_biomasses[rand_sp] < 1e-10)) {
      new_biomasses[rand_sp] <- 0
      return(new_biomasses)
    }
  }
  
}

# function to run a simulation with a disturbance part way through
# input: node/edge lists, num generations before disturbance, num gens post
#   disturbance, function to simulate disturbance
# output: dataframe or biomasses where columns are species and rows are time steps
sim_disturbance <- function(nodes, edges, pre_dist_gens, post_dist_gens, dist_fun, ...) {
  tic()
  # generate initial state vector
  init_state <- rand_state(nodes)
  
  # run first half of simulation
  first_half <- vectorized_run_sim(nodes, edges, pre_dist_gens, init_state)
  
  # simulate disturbance using disturbance function
  # get state at end of first half of sim
  pre_dist_state <- unlist(first_half[nrow(first_half),])
  # generate post disturbance state using dist_fun
  post_dist_state <- dist_fun(pre_dist_state, ...)
  
  # run second half of sim
  second_half <- vectorized_run_sim(nodes, edges, post_dist_gens, post_dist_state)
  
  # combine two halves of sim
  full_sim <- rbind(first_half, second_half)
  
  print("Finished simulation")
  toc()
  
  # Return the result
  return(full_sim)
}

# input: vectors of biomasses before disturbance and at end of simulation,
#   function to calculate ecosystem service
# output: numeric change in service between two time points
# function to calculate change in a service
calc_delta_es <- function(pre_dist_state, final_state, es_fun, ...) {
  # get service score before and after disturbance
  pd_es <- es_fun(pre_dist_state, ...)
  f_es <- es_fun(final_state, ...)
  # calculate change in service
  delta_es <- (f_es - pd_es)
}







