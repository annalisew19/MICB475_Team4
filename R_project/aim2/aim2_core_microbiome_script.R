#load
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

#### Load data ####
load("ivf_phyloseq.RData") #unrarafied object, analysis take vary seq depth into consideration already

#### "core" microbiome ####

# Convert to relative abundance
phyloseq_RA <- transform_sample_counts(ivf_phyloseq, fun=function(x) x/sum(x))

# Filter dataset by outcome
phyloseq_successful_outcome <- subset_samples(phyloseq_RA, `outcome`=="successful")
phyloseq_unsuccesful_outcome <- subset_samples(phyloseq_RA, `outcome`=="unsuccessful")

# What ASVs are found in more than 70% of samples in each outcome category?
# trying changing the prevalence to see what happens
successful_ASVs <- core_members(phyloseq_successful_outcome, detection=0, prevalence = 0.4)
unsuccessful_ASVs <- core_members(phyloseq_unsuccesful_outcome, detection=0, prevalence = 0.4)

# What are these ASVs? you can code it in two different ways to see the same things
succ_table <- tax_table(prune_taxa(successful_ASVs, ivf_phyloseq))
unsucc_table <- tax_table(prune_taxa(unsuccessful_ASVs, ivf_phyloseq))


# can plot those ASVs' relative abundance
prune_taxa(unsuccessful_ASVs,phyloseq_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`outcome`, scales ="free") #x-axis don't depend on each other
#less groups from antibiotic group

# Notice that in this dataset, there are very few CORE microbiome members. This is common
### What if we wanted to make a Venn diagram of all the ASVs that showed up in each treatment?
successful_list <- core_members(phyloseq_successful_outcome, detection=0.001, prevalence = 0.10)
unsuccessful_list <- core_members(phyloseq_unsuccesful_outcome, detection=0.001, prevalence = 0.10)

outcome_list_full <- list(Successful = successful_list, Unsuccessful = unsuccessful_list)


# Create a Venn diagram using all the ASVs shared and unique to antibiotic users and non users
first_venn <- ggVennDiagram(x = outcome_list_full)

ggsave("venn_antibiotic.png", first_venn)
