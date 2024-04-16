# load libraries
library(tidyverse)
library(kableExtra)
library(igraph)
library(lme4)
library(report)
library(car)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

# read in raw data
data_raw <- read_csv("data/gateway-cleaned.csv")
metadata_raw <- read_csv("data/metadata.csv")

# explore
colnames(data_raw)
unique(data_raw$link.methodology)

citations <- unique(data_raw$link.citation)

kbl(citations, col.names = "citation")

# wrangle data
data_clean <- data_raw %>%
  select(c(link.citation, interaction.type, con.taxonomy, con.mass.mean.g., res.taxonomy, res.mass.mean.g., longitude, latitude, ecosystem.type, foodweb.name, study.site)) %>%
  rename(foodweb_name = foodweb.name)

data_summary <- data_clean %>%
  group_by(link.citation, foodweb_name, study.site, ecosystem.type, latitude, longitude) %>%
  summarize(mean_consumer_mass_g = mean(con.mass.mean.g.), mean_resource_mass_g = mean(res.mass.mean.g.))

# build foodwebs from edge list
food_webs <- lapply(unique(data_clean$foodweb_name), function(fw){
  curr_site <- data_clean %>%
    filter(foodweb_name == fw) %>%
    select(c(con.taxonomy, res.taxonomy))
  curr_edgelist <- as.matrix(curr_site)
  curr_graph <- list(graph_from_edgelist(curr_edgelist, directed = T))
  # add food web name back in
  names(curr_graph) <- fw
  curr_graph
})


# extract structure metrics
food_web_stats <- do.call(rbind, lapply(food_webs, function(i){
  # extract current food web and its name
  curr_name <- names(i)[1]
  curr_web <- i[[1]]
  # calculate some stats
  curr_all_degrees <- degree(curr_web, mode = "total")
  curr_in_degrees <- degree(curr_web, mode = "in")
  curr_out_degrees <- degree(curr_web, mode = "out")
  in_degree_dist <- list(curr_in_degrees)
  out_degree_dist <- list(curr_out_degrees)
  
  stats <- tibble(foodweb_name = curr_name,
             n = vcount(curr_web),
             mean_degree = mean(curr_all_degrees),
             diameter = diameter(curr_web),
             in_degree_distribution = in_degree_dist,
             out_degree_distribution = out_degree_dist,
             mean_distance = mean_distance(curr_web),
             connectance = edge_density(curr_web))
}))

# standardize metadata columns
metadata_clean <- metadata_raw %>%
  mutate(github_citation = tolower(github_citation),
         study.site = tolower(study.site),
         study.site = case_when(
           study.site == "grand caricale" ~ "grand caricaie",
           TRUE ~ study.site
         )) %>%
  rename("link.citation" = "github_citation")

# join in metadata
stats_and_metadata <- food_web_stats %>%
  left_join(select(data_summary, c(link.citation, foodweb_name, study.site, ecosystem.type, longitude, latitude)), by = "foodweb_name") %>%
  # standardize study sites to metadata sheet
  mutate(study.site = case_when(
    study.site == "intertidal study sites in chile: el quisco" ~ "intertidal study sites in chile",
    study.site == "intertidal study sites in chile: los molles" ~ "intertidal study sites in chile",
    study.site == "intertidal study sites in chile : curaumilla" ~ "intertidal study sites in chile",
    study.site == "intertidal study sites in chile: las cruces (ecim)" ~ "intertidal study sites in chile",
    study.site == "bylot" ~ "bylot island",
    study.site == "ythan estuary; tidal estuary of river ythan; forvie nature reserve" ~ "ythan estuary",
    study.site == "herschel" ~ "herschel island",
    TRUE ~ study.site
  )) %>%
  left_join(metadata_clean, by = c("link.citation", "study.site")) %>%
  # drop food webs we couldn't find data for
  filter(include == "yes",
         !is.na(gut_content),
         !is.na(expert),
         !is.na(literature),
         !is.na(field_obs),
         !is.na(feeding_trial),
         !is.na(natural_history),
         !is.na(stable_isotope),
         !is.na(pred_analysis),
         !is.na(molecular),
         !is.na(presence_sampling),
         !is.na(fatty_acid)) %>%
  # recode yes/no to 1/0
  mutate(expert = ifelse(expert == "yes", 1, 0),
         gut_content = ifelse(gut_content == "yes", 1, 0),
         literature = ifelse(literature == "yes", 1, 0),
         field_obs = ifelse(field_obs == "yes", 1, 0),
         feeding_trial = ifelse(feeding_trial == "yes", 1, 0),
         natural_history = ifelse(natural_history == "yes", 1, 0),
         stable_isotope = ifelse(stable_isotope == "yes", 1, 0),
         pred_analysis = ifelse(pred_analysis == "yes", 1, 0),
         molecular = ifelse(molecular == "yes", 1, 0),
         presence_sampling = ifelse(presence_sampling == "yes", 1, 0),
         fatty_acid = ifelse(fatty_acid == "yes", 1, 0))

# make some models
linear_model <- lm(connectance ~ gut_content + 
                     expert + 
                     literature + 
                     field_obs + 
                     feeding_trial + 
                     natural_history + 
                     pred_analysis + 
                     molecular + 
                     presence_sampling +
                     paper_1_year, 
                   data = stats_and_metadata,
                   na.action = na.fail)
car::vif(linear_model)
summary(linear_model)
MuMIn::dredge(linear_model)

reduced_model_1 <- lm(n ~ feeding_trial +
                      field_obs +
                      gut_content +
                      literature +
                      molecular +
                      paper_1_year +
                      pred_analysis,
                    data = stats_and_metadata)

reduced_model_2 <- lm(connectance ~ literature +
                        molecular +
                        natural_history +
                        paper_1_year,
                      data = stats_and_metadata)

car::vif(reduced_model_2)

summary(reduced_model_1)
summary(reduced_model_2)


mixed_model <- lmer(connectance ~ feeding_trial +
                      field_obs + 
                      molecular +
                      natural_history + 
                      paper_1_year +
                      (1|link.citation), 
                   data = stats_and_metadata)
mixed_model_n <- lmer(n ~ feeding_trial +
                        field_obs +
                        gut_content +
                        literature +
                        molecular +
                        paper_1_year +
                        pred_analysis +
                        (1|link.citation), 
                      data = stats_and_metadata)

summary(mixed_model)
report(mixed_model)
report(mixed_model_n)

AIC(linear_model)
AIC(mixed_model)
AIC(reduced_model_2)

# make some visualizations
ggplot(data = stats_and_metadata, aes(x = connectance)) +
  geom_histogram()

ggplot(data = stats_and_metadata, aes(x = n)) +
  geom_histogram()

ggplot(data = stats_and_metadata, aes(x = mean_degree)) +
  geom_histogram()

ggplot(data = stats_and_metadata, aes(x = (diameter/n))) +
  geom_histogram()

ggplot(data = stats_and_metadata, aes(x = mean_distance)) +
  geom_histogram()

