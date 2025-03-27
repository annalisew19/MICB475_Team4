# Random Forest Model 

# Based on this guide: https://www.r-bloggers.com/2021/04/random-forest-in-r/

# Load library
library(randomForest)
library(caret)

#### Make model using the 9 taxa ####
# Load in data that has the ASVs as columns with their counts and the metadata column of interest per sample but remove sample_id
data <- read.delim("kn_class_filt.csv", sep = ",")

# Ensure sample size is the same across categories ideally
data$ionizing_radiation <- as.factor(data$ionizing_radiation)
table(data$ionizing_radiation)

# Lets start with random seed so the outcome will be repeatable and store train and test data.
set.seed(222)
ind <- sample(2, nrow(data), replace = TRUE, prob = c(0.7, 0.3))
train <- data[ind==1,]
test <- data[ind==2,]

# Run the random forest code on the training data
rf <- randomForest(ionizing_radiation~., data=train, proximity=TRUE)
print(rf)

# you can plot the model
plot(rf) 
MDSplot(rf, train$ionizing_radiation)
varImpPlot(rf)

# prediction and confusion matrix
p1 <- predict(rf, train)
confusionMatrix(p1, train$ ionizing_radiation)

# run prediction and confusion on the test data now
p2 <- predict(rf, test)
confusionMatrix(p2, test$ ionizing_radiation)

#### Make model using all taxa ####
all_data <- read.delim("kn_class.csv", sep = ",")

all_data$ionizing_radiation <- as.factor(all_data$ionizing_radiation)
table(all_data$ionizing_radiation)

train_all <- all_data[ind==1,]
test_all <- all_data[ind==2,]

rf_all <- randomForest(ionizing_radiation~., data=train_all, proximity=TRUE)
print(rf_all)

# you can plot the model
plot(rf_all)
MDSplot(rf_all, train_all$ionizing_radiation)
varImpPlot(rf_all)

# prediction and confusion matrix
p1_all <- predict(rf_all, train_all)
confusionMatrix(p1_all, train_all$ ionizing_radiation)

# run prediction and confusion on the test data now
p2_all <- predict(rf_all, test_all)
confusionMatrix(p2_all, test_all$ ionizing_radiation)
