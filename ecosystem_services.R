# Author: Luis X. de Pablo
# E-mail: luis.depablo@colorado.edu
# Last updated: September 11, 2023

# load packages and set working directory ----

library(tidyverse)


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