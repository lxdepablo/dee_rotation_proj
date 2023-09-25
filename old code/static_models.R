# load packages
library(igraph)
library(doBy)
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#-----------------------------------------------------------------------------#
# define helper functions

# create_species(): function to create and parameterize a species by providing beta for Gamma distribution
# input: positive real number beta
# output: vector of length 3
create_species <- function(beta) {
  new_species <- numeric(length = 3) # initialize a vector to describe a new species
  new_species[1] <- runif(1, min = 0, max = 1) # draw niche value n_i from Unif[0,1] and store in vector
  y <- runif(1, min = 0, max = 1) # set up to draw fundamental generality range r_i
  x <- 1 - (1 - y)^(1 / beta) # set up to draw fundamental generality range r_i
  new_species[2] <- x * new_species[1] # draw r_i, see Williams and Martinez (2000) for derivation
  new_species[3] <- runif(1, min = new_species[2] / 2, max = new_species[1]) # draw range center c_i from Unif[r_i/2, n_i] and store in vector
  return(new_species)
}

# make_Par(): function to construct a parameter matrix of n species
# input: positive integer S, positive real number beta
# output: n x 3 matrix
make_Par <- function(S, beta) {
  parameter_matrix <- matrix(NA, nrow = S, ncol = 3) # initialize a matrix to store species parameters
  colnames(parameter_matrix) <- c("n_i", "r_i", "c_i") # name matrix columns
  for (i in 1:S) { # populate matrix by looping through each species
    parameter_matrix[i,] <- t(create_species(beta)) # for each row create a new species
  }
  return(parameter_matrix)
}

# make_Adj(): function to construct an adjacency matrix from a parameter matrix
# input: S x 3 matrix for positive integer S
# output: S x S matrix
make_Adj <- function(parameter_matrix, B) {
  S <- nrow(parameter_matrix) # store the number of species
  adjacency_matrix <- matrix(NA, nrow = S, ncol = S) # initialize a matrix to store network edges
  for (i in 1:S) { # calculate feeding ranges for each species
    i_lower <- parameter_matrix[i,3] - parameter_matrix[i,2] / 2 # define lower boundary of a species' feeding range
    i_upper <- parameter_matrix[i,3] + parameter_matrix[i,2] / 2 # define upper boundary of a species' feeding range
    i_eats <- integer(length = S) # initialize a vector to store network edges for a species
    for (j in 1:S) { # populate vector by looping through each species
      if ((parameter_matrix[j,1] > i_lower) & (parameter_matrix[j, 1] < i_upper)) { # determine whether a species' n_i falls within the feeding range
        i_eats[j] <- 1
      }
    }
    adjacency_matrix[i,] <- i_eats # populate matrix by inserting vector for each species
  }
  
  # assign B basal species to to adjacency matrix
  basal_species <- which.minn(parameter_matrix[,1],B) # get a vector of the B species with lowest n_i values
  adjacency_matrix[basal_species, ] <- rep(0, S) # zero out rows for all basal species in the adjacency matrix
  
  return(adjacency_matrix)
}

# make_Edg(): function to construct an edge list from an adjacency matrix
# input: S x S matrix for positive integer S
# output: a x 2 matrix for positive integer a
make_Edg <- function(adjacency_matrix) {
  n <- sum(adjacency_matrix) # store the number of edges in an adjacency matrix
  edge_list <- matrix(NA, nrow = n, ncol = 2) # initialize a matrix to store network edges
  colnames(edge_list) <- c("consumer", "resource") # name matrix columns
  row_index <- 1 # initialize a counter to store the ith interaction to the ith row of the edge list
  for (i in 1:nrow(adjacency_matrix)) { # loop through all rows of adjacency matrix
    for (j in 1:ncol(adjacency_matrix)) { # loop through all columns of adjacency matrix
      if (adjacency_matrix[i,j] == 1) { # add i,j to the edge list if species i eats species j
        edge_list[row_index, 1] <- i
        edge_list[row_index, 2] <- j
        row_index <- row_index + 1 # increment row_index to the next row of edge list
      }
    }
  }
  return(edge_list)
}

# isolates(): function to detect identity of single isolated species in adjacency matrix
# input: S x S matrix for positive integer S
# output: vector of length a for non-negative integer a
isolates <- function(adjacency_matrix) {
  isolated_species <- numeric(length = 0) # initialize a vector to store indices of isolated species
  for (i in 1:nrow(adjacency_matrix)) { # loop through all species to identify isolated species
    if (sum(adjacency_matrix[i,]) == 0 & sum(adjacency_matrix[,i]) == 0) { # add species to vector if it eats no one and no one eats it
      isolated_species <- append(isolated_species, i)
    }
  }
  return(isolated_species)
}

# identicals(): function to detect identity of single trophically identical species in adjacency matrix
# input: S x S matrix for positive integer S
# output: vector of length a for non-negative integer a
identicals <- function(adjacency_matrix) {
  identical_species <- numeric(length = 0) # initialize a vector to store indices of identical species
  for (i in 1:nrow(adjacency_matrix)) { # loop through all species-pairs to identify identical species
    i_eats <- which(adjacency_matrix[i,] == 1) # initialize a vector to store all species that species i eats
    i_eaten <- which(adjacency_matrix[,i] == 1) # initialize a vector to store all species that eat species i
    for (j in 1:nrow(adjacency_matrix)) {
      if (j != i){
        j_eats <- which(adjacency_matrix[j,] == 1) # initialize a vector to store all species that species j eats
        j_eaten <- which(adjacency_matrix[,j] == 1) # initialize a vector to store all species that eats species j
        if (identical(i_eats, j_eats) & identical(i_eaten, j_eaten)) { # add species j to vector if it has the same eats and eaten vectors as species i
          identical_species <- append(identical_species, j)
        }
      }
    }
  }
  return(identical_species)
}

# lower_triangular(): function to make an acyclic adjacency matrix lower triangular
# input: S x S matrix for positive integer S
# output: list(S x S matrix, vector of length S)
lower_triangular <- function(adjacency_matrix) {
  graph <- graph_from_adjacency_matrix(adjacency_matrix) # convert adjacency matrix into igraph object
  graph <- topo_sort(graph, mode = "in") # sort adjacency matrix based on incoming edges
  return(list(adjacency_matrix[graph, graph], as.vector(graph)))
}

# DFS(): function to decycle an adjacency matrix
# input: S x S matrix for positive integer S
# output: S x S matrix
# taken from threat_functions.R by Gyorgy Barabas
DFS <- function(adjacency_matrix) {
  DFSCOLOR <- numeric(0)
  DFSBACKEDGE <- numeric(0)
  ORDERVERTICES <- numeric(0)
  DFSVisit <- function(adjacency_matrix, i) {
    DFSCOLOR[i] <<- 1
    for (j in 1:nrow(adjacency_matrix)) {
      if(adjacency_matrix[i,j] != 0) {
        if (DFSCOLOR[j] == 0) {
          DFSVisit(adjacency_matrix,j)
        } else {
          if(DFSCOLOR[j] == 1) {
            DFSBACKEDGE[i,j] <<- 1 # It's a back edge: list for removal
          }
        }
      }
    }
    DFSCOLOR[i] <<- 2
    ORDERVERTICES <<- c(i, ORDERVERTICES)
  }
  run_DFS <- function(adjacency_matrix) {
    S <- nrow(adjacency_matrix)
    DFSCOLOR <<- rep(0, S)
    DFSBACKEDGE <<- matrix(0, S, S)
    ORDERVERTICES <<- numeric(0)
    for (i in 1:S) if (DFSCOLOR[i] == 0) DFSVisit(adjacency_matrix, i)
    return(adjacency_matrix - DFSBACKEDGE)
  }
  return(run_DFS(adjacency_matrix))
}

# trophic_levels(): function to calculate the prey averaged trophic level for each species from an adjacency matrix
# input: S x S matrix for positive integer S
# output: vector of length S
# adapted from threat_functions.R by Gyorgy Barabas
trophic_levels <- function(adjacency_matrix) {
  A <- adjacency_matrix
  A <- A / rowSums(A)
  A[is.nan(A)] <- 0
  S <- nrow(A)
  tl <- solve(diag(S) - A, rep(1, S))
  return(tl)
}

# make_Att(): function to construct a matrix of trophic level and In Degree from an adjacency matrix
# input: S x S matrix for positive integer S
# output: S x 2 matrix
make_Att <- function(adjacency_matrix) {
  S <- nrow(adjacency_matrix) # store the number of species
  attribute_matrix <- matrix(NA, nrow = S, ncol = 2) # initialize a matrix to store attributes
  colnames(attribute_matrix) <- c("trophic_level", "in_degree") # name matrix columns
  attribute_matrix[, 1] <- trophic_levels(adjacency_matrix) # calculate and store trophic level for all species
  for (i in 1:S) { # store in degree for all species
    attribute_matrix[i, 2] <- sum(adjacency_matrix[i,]) # calculate the in degree for species i
  }
  return(attribute_matrix)
}

# remove_links(): function to remove links from existing network
# input: S x 5 matrix, S x S matrix and p, for positive integer S and p in [0,1]
# output: list(vector of length a, S x S matrix, S x 3 matrix) for positive integer a
remove_links <- function(parameter_matrix, adjacency_matrix, p) {
  links <- sum(adjacency_matrix) # get number of links in original web
  number_of_removals <- round(links * p, digits = 0) # get number of links to remove
  links_to_remove <- sample(which(adjacency_matrix == 1), number_of_removals) # get a vector of links to remove
  new_adjacency_matrix <- adjacency_matrix # initiate the new adjacency matrix
  new_adjacency_matrix[links_to_remove] <- 0 # set all links in removal vector to 0 in new adjacency matrix
  
  problem_species <- c(isolates(new_adjacency_matrix), identicals(new_adjacency_matrix)) # get vector of problem webs in new adjacency matrix
  new_matrix_graph <- graph_from_adjacency_matrix(new_adjacency_matrix, mode = "directed") # convert new adjacency matrix into igraph object to test for connectedness
  while(length(problem_species) > 0 | !is_connected(new_matrix_graph, mode = "weak")) { # keep redoing link removal process until the proposed new adjacency matrix satisfies conditions
    links_to_remove <- sample(which(adjacency_matrix == 1), number_of_removals)
    new_adjacency_matrix <- adjacency_matrix
    new_adjacency_matrix[links_to_remove] <- 0
    problem_species <- c(isolates(new_adjacency_matrix), identicals(new_adjacency_matrix))
    new_matrix_graph <- graph_from_adjacency_matrix(new_adjacency_matrix, mode = "directed")
  }
  
  lower_triangulate <- lower_triangular(new_adjacency_matrix) # lower triangulate the new adjacency matrix and store the species sequence
  new_adjacency_matrix <- lower_triangulate[[1]] # make the new adjacency matrix lower triangular
  new_parameter_matrix <- parameter_matrix[lower_triangulate[[2]], ] # reorder the new parameter matrix according to the species sequence
  new_edge_list <- make_Edg(new_adjacency_matrix) # make the new edge list based on the new adjacency matrix
  
  return(list(new_edge_list, new_adjacency_matrix, new_parameter_matrix))
}

#-----------------------------------------------------------------------------#
# Stratified Model
make_strat_web <- function(N, n, mixing_parameter) {
  mixing_matrix <- matrix(0, nrow = N, ncol = N) # initialize mixing matrix with all 0s
  for (i in 1:N-1) {
    mixing_matrix[i+1,i] <- mixing_parameter # populate mixing matrix
  }
  
  z <- c() # initialize group label vector
  for (i in 1:N) { # add all groups to label vector
    new_additions <- rep(i, n) # n species are added as trophic level i species
    z <- c(z, new_additions) # concatenate new species to existing label vector
  }
  
  adjacency_matrix <- matrix(0, nrow = length(z), ncol = length(z)) # initialize adjacency matrix with all 0s
  for (i in 1:length(z)) { # loop over all species in adjacency matrix
    i_group <- z[i] # get the group of species i
    for (j in 1:length(z)) { # loop over all possible prey of species i in adjacency matrix
      j_group <- z[j] # get the group of species j
      if (runif(1, min = 0, max = 1) < mixing_matrix[i_group, j_group]) {
        adjacency_matrix[i,j] <- 1 # i eats j with probability M_ij
      }
    }
  }
  problem_species <- unique(c(isolates(adjacency_matrix), identicals(adjacency_matrix))) # get vector of problem species in candidate web
  web_as_graph <- graph_from_adjacency_matrix(adjacency_matrix, mode = "directed") # convert candidate web into igraph object to check connectedness
  
  while (length(problem_species) > 0 | !is_connected(web_as_graph, mode = "weak")) { # keep generating candidate webs as long as problem species are present or web is not connected
    mixing_matrix <- matrix(0, nrow = N, ncol = N) # initialize mixing matrix with all 0s
    for (i in 1:N-1) {
      mixing_matrix[i+1,i] <- mixing_parameter # populate mixing matrix
    }
    
    z <- c() # initialize group label vector
    for (i in 1:N) { # add all groups to label vector
      new_additions <- rep(i, n) # n species are added as trophic level i species
      z <- c(z, new_additions) # concatenate new species to existing label vector
    }
    
    adjacency_matrix <- matrix(0, nrow = length(z), ncol = length(z)) # initialize adjacency matrix with all 0s
    for (i in 1:length(z)) { # loop over all species in adjacency matrix
      i_group <- z[i] # get the group of species i
      for (j in 1:length(z)) { # loop over all possible prey of species i in adjacency matrix
        j_group <- z[j] # get the group of species j
        if (runif(1, min = 0, max = 1) < mixing_matrix[i_group, j_group]) {
          adjacency_matrix[i,j] <- 1 # i eats j with probability M_ij
        }
      }
    }
    problem_species <- unique(c(isolates(adjacency_matrix), identicals(adjacency_matrix))) # reevaluate vector of problem species in new candidate web
    web_as_graph <- graph_from_adjacency_matrix(adjacency_matrix, mode = "directed") # convert new candidate web into igraph object to check connectedness
  }
  
  adjacency_matrix <- lower_triangular(adjacency_matrix)[[1]] # make the adjacency matrix lower triangular
  edge_list <- make_Edg(adjacency_matrix) # make the edge list based on the adjacency matrix
  
  return(list(edge_list, adjacency_matrix))
}

#-----------------------------------------------------------------------------#
# Stochastic Block Model
make_sbm_web <- function(N, n, C) {
  alpha_denom <- 0 # initialize the denominator of alpha
  for (i in 1:N-1) { # sum together n_1*n_2 +...+ n_(N-1)*n_N
    x <- n[i]*n[i+1]
    alpha_denom <- alpha_denom + x
  }
  alpha <- (C*(sum(n))^2) / alpha_denom
  
  mixing_matrix <- matrix(0, nrow = N, ncol = N) # initialize mixing matrix with all 0s
  for (i in 1:N-1) {
    mixing_matrix[i+1,i] <- alpha # populate mixing matrix
  }
  
  z <- c() # initialize group label vector
  for (i in 1:N) { # add all groups to label vector
    new_additions <- rep(i, n[i]) # n[i] species are added as trophic level i species
    z <- c(z, new_additions) # concatenate new species to existing label vector
  }
  
  adjacency_matrix <- matrix(0, nrow = length(z), ncol = length(z)) # initialize adjacency matrix with all 0s
  for (i in 1:length(z)) { # loop over all species in adjacency matrix
    i_group <- z[i] # get the group of species i
    for (j in 1:length(z)) { # loop over all possible prey of species i in adjacency matrix
      j_group <- z[j] # get the group of species j
      if (runif(1, min = 0, max = 1) < mixing_matrix[i_group, j_group]) {
        adjacency_matrix[i,j] <- 1 # i eats j with probability M_ij
      }
    }
  }
  adjacency_matrix <- lower_triangular(adjacency_matrix) # make the adjacency matrix lower triangular
  edge_list <- make_Edg(adjacency_matrix) # make the edge list based on the adjacency matrix
  
  return(list(edge_list, adjacency_matrix))
}

#-----------------------------------------------------------------------------#
# Cascade Model
# returns a model-generated food web with species S and connectance C
# see Cohen and Newman (1985) for model
# Cohen and Newman (1985): 10.1098/rspb.1985.0042

# make_cascade_web(): master function that generates cascade model network
# input: Number of species S, connectance C for positive integer S and C in [0,1]
# output: list(a x 2 matrix, S x S matrix, S x 2 matrix) for positive integer a
make_cascade_web <- function(S, C) {
  p <- (2*C*S)/ (S - 1) # store probability parameter for Bernoulli trials
  nsims <- S*(S-1) / 2 # store the number of Bernoulli trials to run
  sim_vector <- rbinom(nsims, 1, p) # get a vector of Bernoulli trial outcomes
  adjacency_matrix <- matrix(0, nrow = S, ncol = S) # initialize an empty adjacency matrix
  adjacency_matrix[lower.tri(adjacency_matrix, diag = FALSE)] <- sim_vector # store vector of outcomes in upper triangular elements of adjacency matrix
  adjacency_matrix <- DFS(adjacency_matrix) # remove all cycles from the adjacency matrix
  problem_species <- c(isolates(adjacency_matrix), identicals(adjacency_matrix)) # initialize a vector of all isolated or identical species
  adjacency_matrix_graph <- graph_from_adjacency_matrix(adjacency_matrix, mode = "directed") # convert adjacency matrix into igraph object to test for connectedness
  
  while (length(problem_species) > 0 | !is_connected(adjacency_matrix_graph, mode = "weak")) { # keep regenerating adjacency matrices until the proposed new adjacency matrix satisfies conditions
    sim_vector <- rbinom(nsims, 1, p)
    adjacency_matrix <- matrix(0, nrow = S, ncol = S)
    adjacency_matrix[lower.tri(adjacency_matrix, diag = FALSE)] <- sim_vector
    adjacency_matrix <- DFS(adjacency_matrix)
    problem_species <- c(isolates(adjacency_matrix), identicals(adjacency_matrix))
    adjacency_matrix_graph <- graph_from_adjacency_matrix(adjacency_matrix, mode = "directed")
  }
  
  lower_triangulated <- lower_triangular(adjacency_matrix) # lower triangulate the adjacency matrix
  adjacency_matrix <- lower_triangulated[[1]]
  edge_list <- make_Edg(adjacency_matrix)
  
  return(list(edge_list, adjacency_matrix))
}

#-----------------------------------------------------------------------------#
# Niche Model
# returns a model-generated food web with species S and connectance C
# see Williams and Martinez (2000) for model
# Williams and Martinez (2000): 10.1038/35004572

# make_niche_web(): master function that generates niche model network
# input: Number of species S, number of basal species B, connectance C for positive integers S and B, and C in [0,1]
# output: list(a x 2 matrix, S x S matrix, S x 5 matrix) for positive integer a
make_niche_web <- function(S, B, C) {
  beta <- (1 - 2*C) / (2*C) # store beta parameter of Gamma distribution
  parameter_matrix <- make_Par(S, beta) # make a parameter matrix with S species
  adjacency_matrix <- make_Adj(parameter_matrix, B) # make an adjacency matrix using the parameter matrix with B basal species
  adjacency_matrix <- DFS(adjacency_matrix) # remove all cycles from the adjacency matrix
  
  problem_species <- c(isolates(adjacency_matrix), identicals(adjacency_matrix)) # initialize a vector of all isolated or identical species
  adjacency_matrix_graph <- graph_from_adjacency_matrix(adjacency_matrix, mode = "directed") # convert adjacency matrix into igraph object to test for connectedness
  while (length(problem_species) > 0 | !is_connected(adjacency_matrix_graph, mode = "weak")) { # keep regenerating parameter and adjacency matrices until the proposed new adjacency matrix satisfies conditions
    parameter_matrix <- make_Par(S, beta)
    adjacency_matrix <- make_Adj(parameter_matrix, B)
    adjacency_matrix <- DFS(adjacency_matrix)
    problem_species <- c(isolates(adjacency_matrix), identicals(adjacency_matrix))
    adjacency_matrix_graph <- graph_from_adjacency_matrix(adjacency_matrix, mode = "directed")
  }
  
  lower_triangulated <- lower_triangular(adjacency_matrix) # lower triangulate the adjacency matrix and store the species sequence
  adjacency_matrix <- lower_triangulated[[1]] # make the adjacency matrix lower triangular
  edge_list <- make_Edg(adjacency_matrix) # make the edge list based on the adjacency matrix
  parameter_matrix <- cbind(parameter_matrix[lower_triangulated[[2]], ], make_Att(adjacency_matrix)) # reorder the parameter matrix according to the species sequence and add the trophic level and indegree of each species
  
  return(list(edge_list, adjacency_matrix, parameter_matrix))
}