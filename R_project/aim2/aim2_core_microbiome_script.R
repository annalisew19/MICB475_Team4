#load
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

#### Load data ####
load("ivf_phyloseq.RData") #unrarafied object, analysis take vary seq depth into consideration already

# Convert raw read counts to relative abundance. counts -> percentage, so more comparable
phyloseq_RA <- transform_sample_counts(ivf_phyloseq, fun=function(x) x/sum(x))

#### "core" microbiome of outcome (successful or unsuccesful) ####

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

outcome_list_full <- list(S= successful_list, U= unsuccessful_list)


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


#### "core" microbiome of disease ####

# Filter dataset by disease
phyloseq_no_pregnancy <- subset_samples(phyloseq_RA, `disease`=="No pregnancy")
phyloseq_ongoing_preg <- subset_samples(phyloseq_RA, `disease`=="Ongoing pregnancy")
phyloseq_live_birth <- subset_samples(phyloseq_RA, `disease`=="Live birth")
phyloseq_biochem_preg <- subset_samples(phyloseq_RA, `disease`=="Biochemical pregnancy")
phyloseq_clinical_miscarriage <- subset_samples(phyloseq_RA, `disease`=="Clinical miscarriage")

# core microbiome calculation
# detection: ASV must have relative abundance greater than 0; must be present to be considered
# prevalence: ASV must be present in at least xx% of the sample
#tried that the higehst prevalence is 0.4
no_pregnancy_ASVs <- core_members(phyloseq_no_pregnancy, detection=0, prevalence = 0.4)
ongoing_preg_ASVs <- core_members(phyloseq_ongoing_preg, detection=0, prevalence = 0.4)
live_birth_ASVs <- core_members(phyloseq_live_birth, detection=0, prevalence = 0.4)
biochem_preg_ASVs <- core_members(phyloseq_biochem_preg, detection=0, prevalence = 0.4)
clinical_miscarriage_ASVs <- core_members(phyloseq_clinical_miscarriage, detection=0, prevalence = 0.4)


# retrieve taxonomic classifications of the ASVs identified as core members
tax_table(prune_taxa(no_pregnancy_ASVs, ivf_phyloseq))
tax_table(prune_taxa(ongoing_preg_ASVs, ivf_phyloseq))
tax_table(prune_taxa(live_birth_ASVs, ivf_phyloseq))
tax_table(prune_taxa(biochem_preg_ASVs, ivf_phyloseq))
tax_table(prune_taxa(clinical_miscarriage_ASVs, ivf_phyloseq))


# Notice that in this dataset, there are very few CORE microbiome members. This is common
### What if we wanted to make a Venn diagram of all the ASVs that showed up in each treatment?
#detection more tha 0.001 (0.1%)
no_pregnancy_list <- core_members(phyloseq_no_pregnancy, detection=0.001, prevalence = 0.1)
on_going_list <- core_members(phyloseq_ongoing_preg, detection=0.001, prevalence = 0.1)
live_birth_list <- core_members(phyloseq_live_birth, detection=0.001, prevalence = 0.1)
biochem_preg_list <- core_members(phyloseq_biochem_preg, detection=0.001, prevalence = 0.1)
clinical_miscarriage_list <- core_members(phyloseq_clinical_miscarriage, detection=0.001, prevalence = 0.1)

outcome_list_full <- list(no_preg = no_pregnancy_list, 
                          on_going = on_going_list,
                          live_birth = live_birth_list,
                          biochem_preg = biochem_preg_list,
                          clinical_miscarriage = clinical_miscarriage_list)


# Create a Venn diagram using all the ASVs shared and unique to antibiotic users and non users
first_venn <- ggVennDiagram(x = outcome_list_full)

ggsave("venn_outcome", first_venn)


# Get all core member lists
all_lists <- list(no_pregnancy_list, 
                  on_going_list, 
                  live_birth_list, 
                  biochem_preg_list, 
                  clinical_miscarriage_list)
names(all_lists) <- c("no_pregnancy", "ongoing_pregnancy", "live_birth", "biochem_pregnancy", "clinical_miscarriage")

# Create a list to store the unique ASVs for each group
unique_asvs_per_group <- list()

# Loop through each group
for (group_name in names(all_lists)) {
  # Get the core members for the current group
  current_group_list <- all_lists[[group_name]]
  
  # Get the core members for all OTHER groups
  other_groups_lists <- all_lists[names(all_lists) != group_name]
  other_groups_combined <- unlist(other_groups_lists)
  
  # Find the unique ASVs for the current group
  unique_asvs <- setdiff(current_group_list, other_groups_combined)
  
  # Store the unique ASVs in the list
  unique_asvs_per_group[[group_name]] <- unique_asvs
}

# Example: Print the unique ASVs for "No pregnancy"
print(unique_asvs_per_group$no_pregnancy)

# Get the taxonomy table from your original phyloseq object
taxonomy <- as.data.frame(tax_table(ivf_phyloseq))

# Add ASV names as rownames (makes subsetting easier)
taxonomy$ASV <- rownames(taxonomy)

# Loop through the groups and create the tables
taxa_tables_list <- list()

for (group_name in names(unique_asvs_per_group)) {
  unique_asvs <- unique_asvs_per_group[[group_name]]
  taxa_table <- taxonomy[taxonomy$ASV %in% unique_asvs, ]
  taxa_tables_list[[group_name]] <- taxa_table # Store the table in the list
  
}

# Access the tables from the list
View(taxa_tables_list$no_pregnancy)
View(taxa_tables_list$ongoing_pregnancy)
View(taxa_tables_list$live_birth)
View(taxa_tables_list$biochem_pregnancy)
View(taxa_tables_list$clinical_miscarriage)

