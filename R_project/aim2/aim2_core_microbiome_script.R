#load
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

#### Load data ####
load("ivf_phyloseq.RData") #unrarafied object, analysis take vary seq depth into consideration already

#### "core" microbiome ####

# Convert raw read counts to relative abundance. counts -> percentage, so more comparable
phyloseq_RA <- transform_sample_counts(ivf_phyloseq, fun=function(x) x/sum(x))

# Filter dataset by outcome
phyloseq_successful_outcome <- subset_samples(phyloseq_RA, `outcome`=="successful")
phyloseq_unsuccesful_outcome <- subset_samples(phyloseq_RA, `outcome`=="unsuccessful")

# core microbiome calculation
  # detection: ASV must have relative abundance greater than 0; must be present to be considered
  # prevalence: ASV must be present in at least xx% of the sample
  #tried that the higehst prevalence is 0.4
successful_ASVs <- core_members(phyloseq_successful_outcome, detection=0, prevalence = 0.4)
unsuccessful_ASVs <- core_members(phyloseq_unsuccesful_outcome, detection=0, prevalence = 0.4)

# retrieve taxonomic classifications of the ASVs identified as core members
tax_table(prune_taxa(successful_ASVs, ivf_phyloseq))
tax_table(prune_taxa(unsuccessful_ASVs, ivf_phyloseq))


# can plot those ASVs' relative abundance
prune_taxa(unsuccessful_ASVs,phyloseq_RA) %>% 
  plot_bar(fill="Genus")
#  facet_wrap(.~`outcome`, scales ="free") #x-axis don't depend on each other
#less groups from antibiotic group

# Notice that in this dataset, there are very few CORE microbiome members. This is common
### What if we wanted to make a Venn diagram of all the ASVs that showed up in each treatment?
 #detection more tha 0.001 (0.1%)
successful_list <- core_members(phyloseq_successful_outcome, detection=0.001, prevalence = 0.10)
unsuccessful_list <- core_members(phyloseq_unsuccesful_outcome, detection=0.001, prevalence = 0.10)

outcome_list_full <- list(Successful = successful_list, Unsuccessful = unsuccessful_list)


# Create a Venn diagram using all the ASVs shared and unique to antibiotic users and non users
first_venn <- ggVennDiagram(x = outcome_list_full)

ggsave("venn_outcome", first_venn)


# Get the taxonomy table from your original phyloseq object
taxonomy <- as.data.frame(tax_table(ivf_phyloseq))

# Add ASV names as rownames (makes subsetting easier)
taxonomy$ASV <- rownames(taxonomy)

# Get unique ASVs for successful outcome
unique_successful_ASVs <- setdiff(successful_list, unsuccessful_list)

# Get unique ASVs for unsuccessful outcome
unique_unsuccessful_ASVs <- setdiff(unsuccessful_list, successful_list)

# Create tables of unique ASVs and their taxonomy
successful_taxa <- taxonomy[taxonomy$ASV %in% unique_successful_ASVs, ]
unsuccessful_taxa <- taxonomy[taxonomy$ASV %in% unique_unsuccessful_ASVs, ]
