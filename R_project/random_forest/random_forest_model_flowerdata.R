# Random Forest Model 

# Based on this guide: https://www.r-bloggers.com/2021/04/random-forest-in-r/

# Load library
library(randomForest)
library(datasets)
library(caret)

# Using sample dataset iris
data<-iris
str(data)

# Ensure sample size is the same across categories
data$Species <- as.factor(data$Species)
table(data$Species)

# Lets start with random seed so the outcome will be repeatable and store train and test data.
set.seed(222)
ind <- sample(2, nrow(data), replace = TRUE, prob = c(0.7, 0.3))
train <- data[ind==1,]
test <- data[ind==2,]

# Run the random forest code on the training data
rf <- randomForest(Species~., data=train, proximity=TRUE)
print(rf)

# you can plot the model
plot(rf)
MDSplot(rf, train$Species)
varImpPlot(rf)

# prediction and confusion matrix
p1 <- predict(rf, train)
confusionMatrix(p1, train$ Species)

# run prediction and confusion on the test data now
p2 <- predict(rf, test)
confusionMatrix(p2, test$ Species)

