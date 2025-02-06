#### Histogram of ages and sample sizes ####

# Goal: Load metadata into R, plot histogram of number of samples according to each age

# load libraries
library(tidyverse)
library(phyloseq)

# Set the metadata file path and read in metadata
metadatafp <- "../data_files/IVF_metadata.tsv"
meta <- read_delim(metadatafp)


# make a histogram with sample size of age
age_histo <- ggplot(meta, aes(x = AGE)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  scale_x_continuous(breaks = seq(20, 50, by = 5)) +
  scale_y_continuous(breaks = seq(0, 70, by = 10)) +
  labs(title = "Age Distribution of Samples",
       x = "Age",
       y = "Number of Samples")
age_histo

# save age_histo into age_samplesizes folder
ggsave("age_samplesizes.png",
       plot = age_histo)

#### Frequency of different pregnancy outcomes with each age ####
# Count the number of samples per AGE and disease
age_preg_outcome <- meta %>%
  group_by(AGE, disease) %>%
  summarise(Count = n(), .groups = "drop")

# Create a stacked barplot to visualize the counts of different pregnancy outcome
# for each age
age_outcome_count <- ggplot(age_preg_outcome, aes(x = AGE, y = Count, fill = disease)) +
  geom_col(position = "stack") +
  labs(title = "Pregnancy Outcomes by Age",
       x = "Age",
       y = "Number of Samples",
       fill = "Pregnancy Outcome") 
age_outcome_count

# save age_histo into age_samplesizes folder
ggsave("pregnancy_outcomes_by_age.png",
       plot = age_outcome_count)

#Modify metadata to use agegroups in the graph
#add new column age group and filter out all the NAs
meta <- meta %>%
  filter(!is.na(AGE)) %>%
  mutate(AGE_GROUP = case_when(
    AGE >= 20 & AGE <= 25 ~ "20-25",
    AGE >= 26 & AGE <= 30 ~ "26-30",
    AGE >= 31 & AGE <= 35 ~ "31-35",
    AGE >= 36 & AGE <= 40 ~ "36-40",
    AGE >= 41 & AGE <= 45 ~ "41-45",
    AGE >= 46 & AGE <= 50 ~ "46-50"
  ))

#updated histogram using age_group
agegroup_histo <- ggplot(meta, aes(x = AGE_GROUP)) +
  geom_bar(fill = "skyblue", color = "black") +
  labs(title = "Age Group Distribution of Samples",
       x = "Age Group",
       y = "Number of Samples") 
agegroup_histo

#save as png
ggsave("agegroup_samplesizes.png",
       plot = agegroup_histo)

# Create a stacked barplot to visualize the counts of different pregnancy outcome
# for each age group
agegroup_preg_outcome <- meta %>%
  group_by(AGE_GROUP, disease) %>%
  summarise(Count = n(), .groups = "drop")

agegroup_outcome_count <- ggplot(agegroup_preg_outcome, aes(x = AGE_GROUP, y = Count, fill = disease)) +
  geom_col(position = "stack") +
  labs(title = "Pregnancy Outcomes by Age Group",
       x = "Age",
       y = "Number of Samples",
       fill = "Pregnancy Outcome") 
agegroup_outcome_count

#save as png
ggsave("pregnancy_outcomes_by_agegroup.png",
       plot = agegroup_outcome_count)
