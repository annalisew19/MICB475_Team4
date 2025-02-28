######## Aim 2: Indicator Species Analysis

#load packages
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(vegan)
library(DESeq2)
library(ape)
library(indicspecies)

########  METADATA MANIPULATION

# Load metadata
IVFmeta <- read_delim("IVF_metadata.tsv", delim = "\t")

# select for coloumns, include only sample id, age, prg outcome (disease), sample name, and tissue
IVFmeta_select <- select(IVFmeta, `sample-id`,`AGE`, `disease`, `tissue`)

# filters rows that have NA in the column AGE and disease
IVFmeta_age_group <- filter(IVFmeta_select, 
                         AGE != "NA",
                         AGE > 25,
                         disease != "NA", 
                         disease != "NA: Not Applicable") %>%
   mutate(age_group = case_when(
    AGE >= 26 & AGE <= 30 ~ "26-30",
    AGE >= 31 & AGE <= 35 ~ "31-35",
    AGE >= 36 & AGE <= 40 ~ "36-40",
    AGE >= 41 & AGE <= 45 ~ "41-45",
    AGE >= 46 & AGE <= 50 ~ "46-50"))


# group by age group and summarize number of samples in each age group
IVFmeta_age_group %>%
  group_by(age_group) %>%  # group by 'age_group'
  summarize(count_samples = n())  

# group by age group and summarize number of samples in each age group
IVFmeta_age_group %>%
  group_by(age_group) %>%  # group by 'age_group'
  summarize(count_samples = n())  

# create new column called preg_outcome, and group "disease" into successful or unsuccessful
IVFmeta_age_and_outcome <- IVFmeta_age_group %>%
  mutate(preg_outcome = case_when(
    disease %in% c("Live birth", "Ongoing pregnancy") ~ "successful",
    disease %in% c("Clinical miscarriage", "Biochemical pregnancy", "No pregnancy", "Ectopic pregnancy") ~ "unsucessful"))

# group by disease and age group, summarize number of samples in each
IVFmeta_age_and_outcome %>%
  group_by(preg_outcome, age_group) %>%  # Group by both preg_outcome and age_group
  summarize(count_samples = n())


############ CREATING PHYLOSEQ OBJECT

# Load features table, taxonomy file, and tree
IVFotu <- read_delim("feature-table.txt", delim = "\t", skip=1)
IVFtax <- read_delim("taxonomy.tsv", delim="\t")
IVFphylotree <- read.tree("tree.nwk")

### format OTU table
# save everything except first column (OTU ID) into a matrix
ivf_otu_mat <- as.matrix(IVFotu[,-1]) #remove the first index column
# Make first column (#OTU ID) the rownames of the new matrix
rownames(ivf_otu_mat) <- IVFotu$'#OTU ID'
# Use the "otu_table" function to make an OTU table
ivf_OTU <- otu_table(ivf_otu_mat, taxa_are_rows = TRUE) 
class(ivf_OTU)

### Format sample metadata
# Save everything except sampleid as new data frame
ivf_samp_df <- as.data.frame(IVFmeta_age_and_outcome[,-1])
# Make sampleids the rownames
rownames(ivf_samp_df)<- IVFmeta_age_and_outcome$'sample-id'
# Make phyloseq sample data with sample_data() function
ivf_SAMP <- sample_data(ivf_samp_df)
class(ivf_SAMP)

### Format taxonomy
# Convert taxon strings to a table with separate taxa rank columns
ivf_tax_mat <- IVFtax %>% select(-Confidence)%>% #remove Confidence columm
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
ivf_tax_mat <- ivf_tax_mat[,-1]
# Make sampleids the rownames
rownames(ivf_tax_mat) <- IVFtax$'Feature ID'
# Make taxa table
ivf_TAX <- tax_table(ivf_tax_mat)
class(ivf_TAX)

### Create phyloseq object 
# Merge all into a phyloseq object
ivf_mpt <- phyloseq(ivf_OTU, ivf_SAMP, ivf_TAX, IVFphylotree)

############ INDICATOR SPECIES ANALYSIS

