#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(indicspecies)

#### Load data ####
# Include this R file in your zipped project file
# Otherwise include the whole script for regenerating this (better practice)
load("pdmouse.RData")

#### Indicator Species/Taxa Analysis ####
# glom to Genus
mouse_genus <- tax_glom(pdmouse, "Genus", NArm = FALSE)
mouse_genus_RA <- transform_sample_counts(mouse_genus, fun=function(x) x/sum(x))

#ISA
isa_mouse <- multipatt(t(otu_table(mouse_genus_RA)), cluster = sample_data(mouse_genus_RA)$`donor_status`)
summary(isa_mouse)
taxtable <- tax_table(pdmouse) %>% as.data.frame() %>% rownames_to_column(var="ASV")

res <- isa_mouse$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05)

View(res)
