# measure structure metrics for all generated null webs and all empirical subwebs
# Henry Li
# July 26, 2023

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(FRACTION)
library(NetSwan)
library(igraph)
library(cheddar)
library(network)
library(intergraph)
library(NetIndices)
library(rnetcarto)
library(vegan)
library(plyr)
library(tidyverse)

#-----------------------------------------------------------------------------#
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

# function to measure trophic coherency of a network
trophic.coherence <- function(network) {
  adj_matrix <- as.matrix(as_adjacency_matrix(network))
  TL_vector <- V(network)$TL
  interactions <- which(adj_matrix == 1, arr.ind = TRUE)
  TC_vector <- rep(NA, nrow(interactions))
  
  for (i in 1:length(TC_vector)) {
    interaction <- unname(interactions[i, ])
    resource <- interaction[1]
    consumer <- interaction[2]
    
    resource_TL <- TL_vector[resource]
    consumer_TL <- TL_vector[consumer]
    trophic_distance <- consumer_TL - resource_TL
    TC_vector[i] <- trophic_distance
  }
  
  TC <- sqrt(mean(TC_vector^2) - (mean(TC_vector)^2))
  return(TC)
}

# function to measure degree heterogeneity of a network
degree.heterogeneity <- function(network) {
  adj_matrix <- as.matrix(as_adjacency_matrix(network))
  degree_vector <- rep(NA, nrow(adj_matrix))
  
  for (i in 1:length(degree_vector)) {
    degree <- sum(adj_matrix[i, ]) + sum(adj_matrix[, i]) - adj_matrix[i, i]
    degree_vector[i] <- degree
  }
  
  DH <- mean(degree_vector^2) / (mean(degree_vector))^2
  return(DH)
}

#-----------------------------------------------------------------------------#
# master function to calculate structure metrics for all null webs of a given subweb and export dataframe
calculate_structure_metrics_null <- function(subweb_name, assembly, resolution) {
  # load in subweb and its null webs
  subweb <- list.webs[[which(SY == subweb_name)]]
  load(paste("C:/Users/henry/OneDrive/Boulder/Projects/GOM community data/generated_subwebs/null_webs/", assembly, ".", resolution, "/", subweb_name, ".null.webs.RData", sep = ""))
  web_list <- c(null_web_list, list(subweb))
  
  # initiate dataframe to store structure metrics
  metrics <- as.data.frame(matrix(NA, nrow = length(web_list), ncol = 17))
  colnames(metrics) <- c("subweb",
                         "type",
                         "assembly",
                         "resolution",
                         "time",
                         "site",
                         "tidalzone",
                         "size",
                         "links",
                         "transitivity",
                         "undirected.modularity",
                         "relative.ascendancy",
                         "motif.food.chain",
                         "motif.omnivory",
                         "motif.loop",
                         "motif.apparent.competition",
                         "motif.direct.competition")
  
  # calculate structure metrics for all webs
  for (i in 1:length(web_list)) {
    # read in web
    web <- web_list[[i]]
    if (class(web) != "igraph") next
    
    # get rid of multiple edges between the same two nodes
    if (sum(which_multiple(web)) > 0) {
      web <- delete_edges(web, which(which_multiple(web) == 1))
    }
    
    # get motif count in web
    motif_vector <- triad_census(web)
    
    metrics$subweb[i] <- subweb_name
    metrics$type[i] <- ifelse(i == length(web_list), "empirical", "null")
    metrics$assembly[i] <- assembly
    metrics$resolution[i] <- resolution
    metrics$time[i] <- ifelse(substr(subweb_name, 1, 4) %in% c("2013", "2014"), "contemporary", "historical")
    metrics$site[i] <- word(subweb_name, 2, sep = "\\.")
    metrics$tidalzone[i] <- word(subweb_name, 3, sep = "\\.")
    metrics$size[i] <- vcount(web)
    metrics$links[i] <- ecount(web)
    metrics$transitivity[i] <- transitivity(web)
    metrics$undirected.modularity[i] <- modularity(cluster_louvain(as.undirected(web, mode = "collapse")))
    metrics$relative.ascendancy[i] <- enaAscendency(asNetwork(web))
    metrics$motif.food.chain[i] <- motif_vector[6] / sum(motif_vector)
    metrics$motif.omnivory[i] <- motif_vector[9] / sum(motif_vector)
    metrics$motif.loop[i] <- motif_vector[10] / sum(motif_vector)
    metrics$motif.apparent.competition[i] <- motif_vector[5] / sum(motif_vector)
    metrics$motif.direct.competition[i] <- motif_vector[4] / sum(motif_vector)
  }
  metrics <- na.omit(metrics)
  
  # export structure metrics as a .csv file
  write.csv(metrics, paste("C:/Users/henry/OneDrive/Boulder/Projects/GOM community data/henry_ch1/metrics/null_webs/", assembly, ".", resolution, "/", subweb_name, ".metrics.csv", sep = ""), row.names = F)
}

# master function to calculate structure metrics for all empirical subwebs and export dataframe
calculate_structure_metrics_empirical <- function(assembly, resolution) {
  # load in subwebs
  load(paste("C:/Users/henry/OneDrive/Boulder/Projects/GOM community data/generated_subwebs/empirical_subwebs/", assembly, ".", resolution, "/web.lists.tidalzone.", assembly, ".", resolution, ".historical.RData", sep = ""))
  
  # initiate dataframe to store structure metrics
  metrics <- as.data.frame(matrix(NA, nrow = length(list.webs), ncol = 16))
  colnames(metrics) <- c("subweb",
                         "type",
                         "assembly",
                         "resolution",
                         "time",
                         "site",
                         "tidalzone",
                         "size",
                         "links",
                         "mean.degree",
                         "connectance",
                         "mean.tl",
                         "fraction.basal",
                         "fraction.intermediate",
                         "fraction.top",
                         "fraction.omnivorous")
  
  # calculate structure metrics for all webs
  for (i in 1:length(list.webs)) {
    # read in web
    web <- list.webs[[i]]
    web_cheddar <- cheddar.webs[[i]]
    
    if (class(web) != "igraph") next
    
    # get rid of multiple edges between the same two nodes
    if (sum(which_multiple(web)) > 0) {
      web <- delete_edges(web, which(which_multiple(web) == 1))
    }
    
    # list web as an adjacency matrix
    adjacency_matrix <- as.matrix(as_adjacency_matrix(web))# list of networks as adjacency matrices
    
    metrics$subweb[i] <- SY[i]
    metrics$type[i] <- "empirical"
    metrics$assembly[i] <- assembly
    metrics$resolution[i] <- resolution
    metrics$time[i] <- ifelse(substr(SY[i], 1, 4) %in% c("2013", "2014"), "contemporary", "historical")
    metrics$site[i] <- word(SY[i], 2, sep = "\\.")
    metrics$tidalzone[i] <- word(SY[i], 2, sep = "\\.")
    metrics$size[i] <- vcount(web)
    metrics$links[i] <- ecount(web)
    metrics$mean.degree[i] <- 2 * ecount(web) / vcount(web)
    metrics$connectance[i] <- ecount(web) / (vcount(web))^2
    metrics$mean.tl[i] <- mean(TrophInd(adjacency_matrix)$TL)
    metrics$fraction.basal[i] <- FractionBasalNodes(web_cheddar)
    metrics$fraction.intermediate[i] <- FractionIntermediateNodes(web_cheddar) 
    metrics$fraction.top[i] <- FractionTopLevelNodes(web_cheddar)
    metrics$fraction.omnivorous[i] <- FractionOmnivorous(web_cheddar)
  }
  metrics <- na.omit(metrics)
  
  # export structure metrics as a .csv file
  write.csv(metrics, paste("C:/Users/henry/OneDrive/Boulder/Projects/GOM community data/henry_ch1/metrics/empirical_webs/", assembly, ".", resolution, "/", assembly, ".", resolution, ".metrics.csv", sep = ""), row.names = F)
}

#-----------------------------------------------------------------------------#
assembly <- c("simple", "full.free.living", "full.parasitic")
resolution <- c("aggregated", "disaggregated")
config_df <- expand_grid(assembly, resolution)

# measure structure metrics for all empirical subwebs for all subweb configurations
for (i in 1:nrow(config_df)) {
  calculate_structure_metrics_empirical(config_df$assembly[i], config_df$resolution[i])
}

df1 <- read.csv("C:/Users/henry/OneDrive/Boulder/Projects/GOM community data/henry_ch1/metrics/empirical_webs/full.free.living.aggregated/full.free.living.aggregated.metrics.csv")
df2 <- read.csv("C:/Users/henry/OneDrive/Boulder/Projects/GOM community data/henry_ch1/metrics/empirical_webs/full.free.living.disaggregated/full.free.living.disaggregated.metrics.csv")
df3 <- read.csv("C:/Users/henry/OneDrive/Boulder/Projects/GOM community data/henry_ch1/metrics/empirical_webs/full.parasitic.aggregated/full.parasitic.aggregated.metrics.csv")
df4 <- read.csv("C:/Users/henry/OneDrive/Boulder/Projects/GOM community data/henry_ch1/metrics/empirical_webs/full.parasitic.disaggregated/full.parasitic.disaggregated.metrics.csv")
df5 <- read.csv("C:/Users/henry/OneDrive/Boulder/Projects/GOM community data/henry_ch1/metrics/empirical_webs/simple.aggregated/simple.aggregated.metrics.csv")
df6 <- read.csv("C:/Users/henry/OneDrive/Boulder/Projects/GOM community data/henry_ch1/metrics/empirical_webs/simple.disaggregated/simple.disaggregated.metrics.csv")
df <- rbind(df1, df2, df3, df4, df5, df6)
write.csv(df, "C:/Users/henry/OneDrive/Boulder/Projects/GOM community data/henry_ch1/metrics/empirical_webs/all/all.metrics.csv", row.names = F)

# measure structure metrics for all null webs of all subwebs
for (i in 1:nrow(config_df)) {
  load(paste("C:/Users/henry/OneDrive/Boulder/Projects/GOM community data/generated_subwebs/empirical_subwebs/", config_df$assembly[i], ".", config_df$resolution[i], "/web.lists.tidalzone.", config_df$assembly[i], ".", config_df$resolution[i], ".historical.RData", sep = ""))
  
  for (j in 1:length(list.webs)) {
    if (j %% 10 == 0) {print(j)}
    if (class(list.webs[[j]]) != "igraph") next
    calculate_structure_metrics_null(SY[j], config_df$assembly[i], config_df$resolution[i])
  }
  
  # rbind all .csvs for a given subweb configuration into one .csv file
  file_list <- Sys.glob(paste("C:/Users/henry/OneDrive/Boulder/Projects/GOM community data/henry_ch1/metrics/null_webs/", config_df$assembly[i], ".", config_df$resolution[i], "/*.csv", sep = ""))
  df <- ldply(file_list, read.csv, header = T)
  write.csv(df, paste("C:/Users/henry/OneDrive/Boulder/Projects/GOM community data/henry_ch1/metrics/null_webs/", config_df$assembly[i], ".", config_df$resolution[i], "/", config_df$assembly[i], ".", config_df$resolution[i], ".metrics.csv", sep = ""), row.names = F)
}

df1 <- read.csv("C:/Users/henry/OneDrive/Boulder/Projects/GOM community data/henry_ch1/metrics/null_webs/full.free.living.aggregated/full.free.living.aggregated.metrics.csv")
df2 <- read.csv("C:/Users/henry/OneDrive/Boulder/Projects/GOM community data/henry_ch1/metrics/null_webs/full.free.living.disaggregated/full.free.living.disaggregated.metrics.csv")
df3 <- read.csv("C:/Users/henry/OneDrive/Boulder/Projects/GOM community data/henry_ch1/metrics/null_webs/full.parasitic.aggregated/full.parasitic.aggregated.metrics.csv")
df4 <- read.csv("C:/Users/henry/OneDrive/Boulder/Projects/GOM community data/henry_ch1/metrics/null_webs/full.parasitic.disaggregated/full.parasitic.disaggregated.metrics.csv")
df5 <- read.csv("C:/Users/henry/OneDrive/Boulder/Projects/GOM community data/henry_ch1/metrics/null_webs/simple.aggregated/simple.aggregated.metrics.csv")
df6 <- read.csv("C:/Users/henry/OneDrive/Boulder/Projects/GOM community data/henry_ch1/metrics/null_webs/simple.disaggregated/simple.disaggregated.metrics.csv")
df <- rbind(df1, df2, df3, df4, df5, df6)
write.csv(df, "C:/Users/henry/OneDrive/Boulder/Projects/GOM community data/henry_ch1/metrics/null_webs/all/all.metrics.csv", row.names = F)