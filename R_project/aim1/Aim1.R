# Goal: Aim1 -  Categorize the endometrial microbiota and microbial diversity in relation to 
# age groups and their corresponding reproductive success

# Load libraries
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(vegan)
library(DESeq2)
library(ape)
library(picante)


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

## Format OTU table
# with rownames and colnames as OTUs and sampleIDs, respectively
ivf_otu_mat <- as.matrix(ivf_otu[,-1])
# Make the first column (#OTU ID) the rownames of the new matrix
rownames(ivf_otu_mat) <- ivf_otu$'#OTU ID'
OTU <- otu_table(ivf_otu_mat, taxa_are_rows = TRUE)
class(OTU)

## Format ivf_meta_updated
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(ivf_meta_updated[,-1])
# Make sampleids the rownames
rownames(samp_df) <- ivf_meta_updated$'sample-id'
# Make phyloseq sample data with sample_data() function
META <- sample_data(samp_df)
class(META)

## Format taxonomy
# Conver taxon strings to a table with separate taxa rank columns
tax_mat <- ivf_tax %>% select(-Confidence) %>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix()
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- ivf_tax$'Feature ID'
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)

## Create phyloseq object
ivf_phyloseq <- phyloseq(OTU, META, TAX, ivf_phylotree)

## Filtering and Rarefying phyloseq object
# Remove non-bacterial sequences, if any
ivf_phyloseq_filt <- subset_taxa(ivf_phyloseq, Domain == "d__Bacteria"
                                 & Class!="c__Chloroplast"
                                 & Family !="f_Mitochondria")
# Remove ASVs that have less than 5 counts total
ivf_phyloseq_nolow <- filter_taxa(ivf_phyloseq_filt, function(x) sum(x)>5, prune = TRUE)
# Remove samples with less than 100 reads
ivf_final <- prune_samples(sample_sums(ivf_phyloseq_nolow)>100, ivf_phyloseq_nolow)

# Adjust plot margins to make sure there is enough space for the plot
par(mar = c(4, 4, 2, 2))  # Adjust margins (bottom, left, top, right)

# Rarefy samples
rarecurve(t(as.data.frame(otu_table(ivf_final))), cex=0.1)
ivf_rare <- rarefy_even_depth(ivf_final, rngseed = 1, sample.size = 2500)

#### Alpha Diversity #### 
## Shannon's Diversity
gg_richness <- plot_richness(ivf_rare, x = "age_group", color = "outcome", measures = c("Shannon")) +
  facet_wrap(~ outcome) +
  geom_boxplot() +
  labs(title = "Shannon Diversity Across Age Groups by IVF Outcome",
       x = "Age Group",
       y = "Shannon Diversity")
gg_richness

ggsave(filename = "shannon_diversity.png",
       gg_richness,
       height = 4, width = 6)

## Faith's PD
# calculate Faith's phylogenetic diversity as PD
phylo_dist <- pd(t(otu_table(ivf_rare)), phy_tree(ivf_rare),
                 include.root = F)

# add PD to metadata table
sample_data(ivf_rare)$PD <- phylo_dist$PD

# plot metadata category against hte PD
plot.pd <- ggplot(sample_data(ivf_rare), aes(age_group, PD)) + 
  geom_boxplot() +
  facet_wrap(~ outcome) +
  xlab("Age Group") +
  ylab("Phylogenetic Diversity")

ggsave("faithpd_boxplot.png",
       height = 4,
       width = 6)

## Statistical Analysis



#### Beta Diversity ####







