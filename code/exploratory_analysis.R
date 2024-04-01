# load libraries
library(tidyverse)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../")

# read in raw data
data_raw <- read_csv("data/gateway-cleaned.csv")

