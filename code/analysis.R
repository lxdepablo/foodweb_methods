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

stats_and_metadata <- read_csv("data/stats_and_metadata.csv")

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


# build null models
# loop over every food web and build n null models for each food web
null_models <- lapply(1:nrow(stats_and_metadata), function(i){
  lapply(1:10, function(n){
    # curr_in_degrees <- stats_and_metadata$in_degree_distribution[[i]]
    # curr_out_degrees <- stats_and_metadata$out_degree_distribution[[i]]
    
    n_vertices <- stats_and_metadata$n[i]
    edge_prob <- stats_and_metadata$connectance[i]
    
    # make null model
    #graph <- sample_degseq(out.deg=curr_out_degrees, in.deg=curr_in_degrees, method="simple.no.multiple")
    graph <- sample_gnp(n_vertices, edge_prob, directed = TRUE)
  })
})

# sample null models for stats of interest
null_model_stats <- do.call(rbind, lapply(null_models, function(i){
  curr_stats <- do.call(rbind, lapply(i, function(n){
    data.frame(diameter = diameter(n),
               mean_distance = mean_distance(n))
  }))
  data.frame(null_diameter = mean(curr_stats$diameter),
             null_mean_distance = mean(curr_stats$mean_distance))
}))

# join back to other data
stats_and_metadata <- cbind(stats_and_metadata, null_model_stats) %>%
  mutate(diff_diameter = diameter - null_diameter,
         diff_distance = mean_distance - null_mean_distance)

ggplot(data = stats_and_metadata, aes(x = diff_diameter)) +
  geom_histogram()

# make some models
global_model_diameter <- lm(diff_diameter ~ gut_content + 
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

global_model_distance <- lm(diff_distance ~ gut_content + 
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

# do model selection
MuMIn::dredge(global_model_diameter)
MuMIn::dredge(global_model_distance)

# build reduced models
reduced_model_diameter <- lm(diff_diameter ~ literature + 
                               natural_history + 
                               pred_analysis + 
                               molecular + 
                               paper_1_year, 
                             data = stats_and_metadata,
                             na.action = na.fail)

reduced_model_distance <- lm(diff_distance ~ literature + 
                               natural_history + 
                               pred_analysis + 
                               molecular + 
                               paper_1_year, 
                             data = stats_and_metadata,
                             na.action = na.fail)
car::vif(reduced_model_diameter)
summary(reduced_model_diameter)
summary(reduced_model_distance)

# # null models
# stats_and_metadata$n_null <- NA
# stats_and_metadata$diameter_null <- NA
# stats_and_metadata$mean_distance_null <- NA
# 
# for (i in 1:nrow(stats_and_metadata)){
#   # i = 6
#   print(i)
#   in_degree <- unname(stats_and_metadata$in_degree_distribution[[i]])
#   out_degree <- unname(stats_and_metadata$out_degree_distribution[[i]])
#   graph <- sample_degseq(out.deg=out_degree, in.deg=in_degree, method="simple.no.multiple")
#   n_list <- list()
#   diameter_list <- list()
#   dist_list <- list()
#   # loop through and generate 10 random graphs each
#   for (i in 1:10){
#     g <- sample_degseq(out.deg=out_degree, in.deg=in_degree, method="simple.no.multiple")
#     n_list[[length(n_list)+1]] = vcount(g)
#     diameter_list[[length(diameter_list)+1]] = diameter(g)
#     dist_list[[length(dist_list)+1]] = mean_distance(g)
#   }
#   avg_n <- mean(unlist(n_list))
#   avg_diam <- mean(unlist(diameter_list))
#   avg_dist <- mean(unlist(dist_list))
#   
#   # stats_and_metadata$null_model[i] <- graph
#   stats_and_metadata$n_null[i] <- avg_n
#   stats_and_metadata$diameter_null[i] <- avg_diam
#   stats_and_metadata$mean_distance_null[i] <- avg_dist
#   
# }
# 
# # testing
# in_degree = unname(stats_and_metadata$in_degree_distribution[[1]])
# out_degree = unname(stats_and_metadata$out_degree_distribution[[1]])
# 
# test_graph <- sample_degseq(out.deg=out_degree, in.deg=in_degree, method="simple.no.multiple")