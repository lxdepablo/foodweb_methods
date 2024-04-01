# load libraries
library(tidyverse)
library(kableExtra)
library(igraph)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

# read in raw data
data_raw <- read_csv("data/gateway-cleaned.csv")

# explore
colnames(data_raw)
unique(data_raw$link.methodology)

citations <- unique(data_raw$link.citation)

kbl(citations, col.names = "citation")

# wrangle data
data_clean <- data_raw %>%
  select(c(link.citation, interaction.type, con.taxonomy, con.mass.mean.g., res.taxonomy, res.mass.mean.g., longitude, latitude, ecosystem.type, foodweb.name, study.site))

# build foodwebs from edge list
food_webs <- lapply(unique(data_clean$foodweb.name), function(fw){
  curr_site <- data_clean %>%
    filter(foodweb.name == fw) %>%
    select(c(con.taxonomy, res.taxonomy))
  curr_edgelist <- as.matrix(curr_site)
  curr_graph <- graph_from_edgelist(curr_edgelist, directed = T)
})

# extract structure metrics
mean_degrees <- lapply(food_webs)
  
  
  
  
  
  
