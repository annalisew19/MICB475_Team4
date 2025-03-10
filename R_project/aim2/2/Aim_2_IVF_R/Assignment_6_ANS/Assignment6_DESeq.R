#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(DESeq2)


#### Load data ####
load("pdmouse.RData")

# Add one to the remove any zeros
mouse_plus1 <- transform_sample_counts(pdmouse, function(x) x+1)

# DESeq comparing donor status
mouse_deseq <- phyloseq_to_deseq2(mouse_plus1, ~`donor_status`)
DESEQ_mouse <- DESeq(mouse_deseq)
res <- results(DESEQ_mouse, tidy=TRUE, 
               #this will ensure that Healthy is your reference group
               contrast = c("donor_status","PD","Healthy"))
View(res)

# Create bar plot
# To get table of results
sigASVs <- res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

View(sigASVs)

# Get only asv names
sigASVs_vec <- sigASVs %>%
  pull(ASV)

# Prune phyloseq file
mouse_DESeq <- prune_taxa(sigASVs_vec,pdmouse)
#Add taxonomy onto DESeq results table
merged_results <- tax_table(mouse_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

# Make DESeq plot
ggplot(merged_results) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
