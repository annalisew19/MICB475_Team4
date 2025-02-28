#load packages
library(tidyverse)

# Load metadata
# the read_delim is a version of the other read.delim command but from tidyverse
IVFmeta <- read_delim("IVF_metadata.tsv", delim = "\t")

# Load features table
IVFotu <- read_delim("feature-table.txt", delim = "\t", skip=1)

# Load taxonomy file
IVFtax <- read_delim("taxonomy.tsv", delim="\t")

### dplyr data manipulation

# select for coloumns, include only sample id, age, prg outcome (disease), sample name, and tissue
IVFmeta_select <- select(IVFmeta, `sample-id`,`AGE`, `disease`, `Sample Name`, `tissue`)

# filters rows that have NA in the column AGE and disease
IVFmeta_filter <- filter(IVFmeta_select, 
                         AGE != "NA",
                         AGE > 26,
                         disease != "NA", 
                         disease != "NA: Not Applicable") 

# create new column called age_group
IVFmeta_age_group <- IVFmeta_filter %>%
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
          
