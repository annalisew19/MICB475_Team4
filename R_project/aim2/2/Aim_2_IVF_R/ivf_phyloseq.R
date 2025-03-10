##### Phyloseq #####
#### Install Packages ####
install.packages("tidyverse")
install.packages("vegan")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

#### Load packages ####
library(phyloseq)
library(tidyverse)
library(ape)

#### Load data ####
metaFP <- "IVF_metadata_status_age.tsv"
meta <- read_delim(file = metaFP, delim = "\t")

otufp  <- "Files/feature-table.txt"
otu <- read_delim(file=otufp, delim = "\t", skip=1)

taxfp <- "Files/taxonomy.tsv"
tax <- read_delim(file = taxfp, delim="\t")

phylotreefp <- "Files/tree.nwk"
phylotree <- read.tree(phylotreefp)

#### Formatting ####
### OTU ###
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$`#OTU ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

### Metadata ###
# Create new data frame omitting sample_name
samp_df <- as.data.frame(meta[,-1])
# Make sample_name the rownames
rownames(samp_df)<- meta$'sample-id'
# Make phyloseq sample data
SAMP <- sample_data(samp_df)

### Taxonomy ###
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$`Feature ID`
# Make taxa table
TAX <- tax_table(tax_mat)

#### Create phyloseq object ####
# Merge all as a phyloseq object
ivf <- phyloseq(OTU, SAMP, TAX, phylotree)

#### Looking at phyloseq object ####
otu_table(ivf)
sample_data(ivf)
tax_table(ivf)
phy_tree(ivf)

#### Saving ####
save(ivf, file="ivf.RData")
