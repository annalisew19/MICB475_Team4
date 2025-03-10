#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

#### Load data ####
# Include this R file in your zipped project file
# Otherwise include the whole script for regenerating this (better practice)
load("pdmouse.RData")

# Convert to relative abundance
mouse_RA <- transform_sample_counts(pdmouse, fun=function(x) x/sum(x))

# Subset dataset into treatment and control groups
healthy_stat <- subset_samples(mouse_RA, `donor_status`=="Healthy")
pd_stat <- subset_samples(mouse_RA, `donor_status`=="PD")

# Set a prevalence threshold and abundance threshold. Be prepared to justify
healthy_list <- core_members(healthy_stat, detection=0.001, prevalence = 0.10)
pd_list <- core_members(pd_stat, detection=0.001, prevalence = 0.10)

list_full <- list(Parkinson = pd_list, Healthy = healthy_list)

# Create a Venn diagram using all the ASVs shared and unique to antibiotic users and non users
mouse_venn <- ggVennDiagram(x = list_full)

ggsave("venn_mouse_status.png", mouse_venn)
