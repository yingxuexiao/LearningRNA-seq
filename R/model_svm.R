
library(e1071)

index <- sample(1:150, 50)
iris

training_set <- iris[-index, ]

test_set <- iris[index, ]

svm.model <- svm(Species ~ ., data = training_set)

unkown_set <- test_set[, -5]

predicted_label <- predict(svm.model, unkown_set)

true_label <- test_set[, 5]

table(predicted_label, true_label)

