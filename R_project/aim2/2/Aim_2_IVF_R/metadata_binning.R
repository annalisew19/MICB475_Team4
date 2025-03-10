install.packages("dplyr")
library(dplyr)

metaFP <- "IVF_metadata.tsv"
meta <- read_delim(file = metaFP, delim = "\t")

unique_outcomes <- unique(meta$disease)
print(unique_outcomes)



IVF_metadata_status <- meta %>%
  mutate(pregnancy_success = case_when(
    disease %in% c("Live birth", "Ongoing pregnancy") ~ "Successful",
    disease %in% c("No pregnancy", "Biochemical pregnancy", "Clinical miscarriage","Ectopic pregnancy") ~ "Not Successful",
    TRUE ~ NA_character_  # Handle any other unexpected values
  ))

IVF_metadata_status_age <- IVF_metadata_status %>%
  mutate(age_group = case_when(
    AGE >= 26 & AGE <= 30 ~ "26-30",
    AGE >= 31 & AGE <= 35 ~ "31-35",
    AGE >= 36 & AGE <= 40 ~ "36-40",
    AGE >= 41 & AGE <= 45 ~ "41-45",
    AGE >= 46 & AGE <= 50 ~ "46-50",
    TRUE ~ "Other"  # Assign "Other" to ages outside the specified range
  ))

# Save the new metadata as a TSV file
write.table(IVF_metadata_status_age, file = "IVF_metadata_status_age.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
