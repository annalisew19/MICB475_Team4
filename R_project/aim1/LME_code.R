# Load necessary libraries
library(lme4)       # For linear mixed-effects models
library(ggplot2)    # For visualization
library(dplyr)      # For data manipulation
library(ggeffects)

# Example: Load your dataset (replace 'your_data.csv' with actual dataset)
# Data should have columns: 'age' (continuous), 'alpha_diversity' (continuous), 'cancer_stage' (categorical)
phylo_dist <- pd(t(otu_table(gc_rare)), phy_tree(gc_rare),
                 include.root=F) 

# add PD to metadata table
sample_data(gc_rare)$PD <- phylo_dist$PD

# Need to extract information
alphadiv <- estimate_richness(gc_rare)
samp_dat <- sample_data(gc_rare)
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)


data <- samp_dat_wdiv

# Convert 'cancer_stage' to a factor (if it's not already)
data$cancer_stage <- as.factor(data$Group)

# Fit the Linear Mixed Effects Model (LME)
lme_model <- lmer(PD ~ Age + (1 | cancer_stage), data = data)

# Create a new dataset with model predictions
data$predicted <- predict(lme_model, re.form = ~ (1 | cancer_stage))  # Include random effects


# Print model summary
summary(lme_model)

# Plot: Regression lines for each cancer stage
ggplot(data, aes(x = Age, y = PD, color = cancer_stage)) +
  geom_point(alpha = 0.6) +  # Scatter plot
  geom_line(aes(y = predicted), size = 1) +  # Regression lines
  labs(title = "PD Across Age for Different Cancer Stages",
       x = "Age",
       y = "PD",
       color = "Cancer Stage") +
  theme_minimal()



# Spearman's correlation for each cancer stage separately
cor_results <- data %>%
  group_by(cancer_stage) %>%
  summarise(spearman_cor = cor.test(Age, PD, method = "spearman")$estimate,
            p_value = cor.test(Age, PD, method = "spearman")$p.value)

# Print correlation results per cancer stage
print(cor_results)
