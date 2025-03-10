#### Install ####
install.packages("indicspecies")

library(tidyverse)
library(phyloseq)
library(indicspecies)

#### Load data ####
load("ivf.RData")


#### Indicator Species/Taxa Analysis ####
# glom to Genus
ivf_genus <- tax_glom(ivf, "Genus", NArm = FALSE)
ivf_genus_RA <- transform_sample_counts(ivf_genus, fun=function(x) x/sum(x))

# Remove samples with NA in pregnancy_success
ivf_genus_RA_clean <- subset_samples(ivf_genus_RA, !is.na(pregnancy_success))

#ISA
isa_ivf <- multipatt(t(otu_table(ivf_genus_RA_clean)), cluster = sample_data(ivf_genus_RA_clean)$`pregnancy_success`)
summary(isa_ivf)
taxtable <- tax_table(ivf) %>% as.data.frame() %>% rownames_to_column(var="ASV")

res <- isa_ivf$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05)

View(res)
