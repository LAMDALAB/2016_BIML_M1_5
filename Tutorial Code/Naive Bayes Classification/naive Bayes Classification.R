# Naive Bayes Classification

## BioConductor - made4 데이터 사용

### 필요한 라이브러리 불러오기
library(made4)
library(e1071)

### 데이터 로드
data(khan)

### Train, Test data 생성
train_data <- khan$train
train_class <- khan$train.classes
test_data <- khan$test
test_class <- khan$test.classes

### 데이터 전처리
#### made4에서 주는 데이터는 행/열을 바꿔야함
train_data <- t(train_data)
test_data <- t(test_data)

### Naive Bayes Classifier model 구축
naive_bayes_model <- naiveBayes(x=train_data, y=train_class)

### test 데이터 predict
test_result_nb <- predict(naive_bayes_model, test_data)

### predict class 비교
table(test_result_nb, test_class) 

## plsgenomics 데이터 사용

### 필요한 라이브러리 불러오기
library(plsgenomics)
library(e1071)

### 데이터 로드
data(SRBCT)

### Train, Test data 생성
train_row <- sample(1:83,60)
train_data <- SRBCT$X[train_row,]
test_data <- SRBCT$X[-train_row,]
train_class <- factor(SRBCT$Y[train_row],levels = c("1","2","3","4"))
test_class <- factor(SRBCT$Y[-train_row],levels = c("1","2","3","4"))

### Naive Bayes Classifier model 구축
naive_bayes_model <- naiveBayes(x=train_data, y=train_class)

### test 데이터 predict
test_result_nb <- predict(naive_bayes_model, test_data)

### predict class 비교
table(test_result_nb, test_class) 

## Cross Validation
k=5
index <- sample(1:k,nrow(SRBCT$X),replace=TRUE)
folds <- 1:k
myRes=data.frame()
for (i in 1:k){
  # create training set
  training <- subset(SRBCT$X, index %in% folds[-i]) 
  # create training set label
  training_class <- factor(subset(SRBCT$Y, index %in% folds[-i]),levels = c("1","2","3","4")) 
  # create test set
  test <- subset(SRBCT$X, index %in% c(i))
  # create test set label
  test_class <- factor(subset(SRBCT$Y, index %in% c(i)),levels = c("1","2","3","4")) 
  # train model
  naive_bayes_model <- naiveBayes(x=training, y=training_class) 
  # run model on test set
  temp <- data.frame(predict(naive_bayes_model, test)) 
  colnames(temp)="Predicted"
  # create data.frame for results
  results <- data.frame(Predicted=temp, Actual=test_class) 
  # append results for each iteration
  myRes <- rbind(myRes, results)
}
table(myRes)


## Evaluation - Confusion Matrix, ROC curve
library(caret)
library(pROC)

### Confusion matrix
confusionMatrix(myRes$Predicted,myRes$Actual)

### ROC curve
roc_nb <- plot.roc(as.numeric(myRes$Predicted), as.numeric(myRes$Actual))
lines(roc_nb, col="black")
legend("bottomright", c("Naive Bayes"), fill = c("black"))

### AUC 
multiclass.roc(as.numeric(myRes$Predicted), as.numeric(myRes$Actual), percent=TRUE)

