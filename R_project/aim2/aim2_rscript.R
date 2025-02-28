#load packages
library(tidyverse)

# Load metadata
# the read_delim is a version of the other read.delim command but from tidyverse
IVFmeta <- read_delim("IVF_metadata.tsv", delim = "\t")

# Load features table
IVFotu <- read_delim("feature-table.txt", delim = "\t", skip=1)

# Load taxonomy file
IVFtax <- read_delim("taxonomy.tsv", delim="\t")

#dplyr data manipulation

#select for coloumns, include only sample id, age, prg outcome (disease), sample name, and tissue
IVFmeta_select <- select(IVFmeta, `sample-id`,`AGE`, `disease`, `Sample Name`, `tissue`)

# Filters rows based on values  
IVFmeta_filter <- filter(IVFmeta_select, AGE != "NA", disease != "NA") #remove NA in age and disease
