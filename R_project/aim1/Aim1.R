# Goal: Aim1 -  Categorize the endometrial microbiota and microbial diversity in relation to 
# age groups and their corresponding reproductive success

# Load libraries
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(vegan)
library(DESeq2)
library(ape)


#### Binning Reproductive outcome and age ####

# Load metadata
ivf_metaFP <- "../data_files/IVF_metadata.tsv"
ivf_meta <- read_delim(ivf_metaFP)

# Update metadate file to only include: sample-id, AGE, disease columns
ivf_meta_updated <- ivf_meta %>%
  select(`sample-id`, AGE, disease)

# Create new column named "outcome" 
# filter out NA in "disease" column
# successful if "disease" = Live birth, Ongoing pregnancy 
# unsuccessful if "disease" = Clinical miscarriage, Biochemical pregnancy, No pregnancy

ivf_meta_updated <- ivf_meta_updated %>%
  filter(disease != "NA: Not Applicable") %>%
  mutate(outcome = ifelse(disease %in% c("Live birth", 
                                         "Ongoing pregnancy"),
                          "successful",
                          ifelse(disease %in% c("Clinical miscarriage",
                                                "Biochemical pregnancy",
                                                "No pregnancy", 
                                                "Ectopic pregnancy"),
                                 "unsuccessful", NA)))

# In AGE column filter out AGE <= 25
# Bin samples into the following age groups: 26-30, 31-35, 36-40, 41-45, 46-50
ivf_meta_updated <- ivf_meta_updated %>%
  filter(AGE > 25) %>%
  mutate(age_group = cut(AGE,
                         breaks = c(25, 30, 35, 40, 45, 50),
                         labels = c("26-30", "31-35", "36-40", "41-45", "46-50"),
                         right = TRUE)) # specifies that the right age boundary should be included



#### Creating phyloseq object ####
# Load feature table, taxonomy, tree 
otufp <- "../data_files/ivf_export/ivf-table_export/feature-table.txt"
ivf_otu <- read_delim(file = otufp, delim = "\t", skip = 1)

taxfp <- "../data_files/ivf_export/ivf_taxonomy_export/taxonomy.tsv"
ivf_tax <- read_delim(taxfp)

phylotreefp <- "../data_files/ivf_export/ivf-rooted-tree_export/tree.nwk"
ivf_phylotree <- read.tree(phylotreefp)

# Format OTU table
# with rownames and colnames as OTUs and sampleIDs, respectively

