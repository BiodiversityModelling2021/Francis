#### Biodiversity modeling course 2021

## Day 2 - Probabilities 

## Problem to solve 

library(tidyverse)

# Load and clean data
mass <- readr::read_csv("data/invasion_NZ.txt")

mass_clean <- mass %>% filter(BodyMass >0, !is.na(BodyMass))

# Plot body mass of invaders and non invaders 
mass_clean %>%
  ggplot(aes(x = log(BodyMass), fill = as.factor(Invader_NZ))) +
  geom_histogram() +
  facet_wrap(~Invader_NZ, scales="free_y", ncol=1)

# Describe data
summary(mass_clean$BodyMass %>% log)

bell_curve <- function(x) exp(-3*(x-4)^2)

bell_curve(mass_clean$BodyMass) %>% summary
