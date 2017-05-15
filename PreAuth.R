# Remove environment variables from R cache
rm(list = ls(all = TRUE))

# Set the working directory
setwd("C:\\Users\\gopal\\Desktop\\Project Study")

# Read the data from csv file
data <- read.csv(file = "PriorAuth_Data.csv", header = T)

# Initial Summary of the dataset
summary(data)
# Structure of data
str(data)
attr <- colnames(data)

# Removing duplicate records
data <- unique(data.frame(data))

# Changing logical attribute to numeric
data <- within(data, {
  Target <- as.numeric(Target)
})

# Drop the attributes which are not required(Most probably)
drop_Attr <- c("UserID", "DoctorID")
attr <- setdiff(attr, drop_Attr)
data <- data[, attr]

# Removing unused references
rm(attr)
rm(drop_Attr)

# Tried removing duplicate records after initial preprocessing
# data <- unique(data.frame(data))

# Checking Summary of the unique data
summary(data)
str(data)

# Download required libraries
library(dplyr)
library(dummies)
library(car)

# Clustering of attributes steps
##########################################################

# Attribute "Drug" w.r.t Target
data_Drug <- select(data, Drug)
dummy_Drug <- dummy(data$Target)
data_Drug <- data.frame("Drug" = data$Drug, dummy_Drug)
data_Drug <-
  data_Drug %>% group_by(Drug) %>% summarise("False" = sum(Target0), "True" = sum(Target1))

# Clustering with mclust package
library(mclust)
mclust_Drug <- Mclust(data_Drug)
# plot results
plot(mclust_Drug, "Classification")
# display the best model
summary(mclust_Drug)


# clustering using Kmeans on Drug.
# kmeans returns an object of class "kmeans" which has a print and a fitted method.
# It is a list with at least following components:
# Total within-cluster sum of squares, i.e. sum(withinss)
# A vector of integers indicating the cluster
# A matrix of cluster centres.

# Random no.generator to restore same(constant result everytime)
set.seed(1234)
tot.wss <- 0
for (i in 1:15) {
  tot.wss[i] <- kmeans(data_Drug[, -c(1)], centers = i)$tot.withinss
}
# Plotting cluster plot to select K
plot(1:15,
     tot.wss,
     type = 'b',
     xlab = "Number of Drug Clusters",
     ylab = "Total Within Group Sum of Squares")
# Removing unused references
rm(i, tot.wss)


# Kmeans with 4 cluster
km_drug <- kmeans(data_Drug[, -c(1)], 4, nstart = 25)
km_drug$cluster

data_Drug <-
  data.frame(data_Drug, "Drug_Cluster" = as.factor(km_drug$cluster))

# Doing innerjoin with original data by data_drug
join_Drug <- inner_join(data, data_Drug, by = "Drug")
join_Drug <-
  subset(data.frame(join_Drug), select = -c(False, True, Drug))


# Removing unused references
rm(dummy_Drug)
rm(data_Drug)
rm(km_drug)

##########################################################

# Attribute "GPI" w.r.t Target

data_GPI <- select(data, GPI)
dummay_GPI <- dummy(data$Target)
data_GPI <- data.frame("GPI" = data$GPI, dummay_GPI)
data_GPI <-
  data_GPI %>% group_by(GPI) %>% summarise("False" = sum(Target0), "True" = sum(Target1))

# Random no.generator to restore same(constant result everytime)
set.seed(1234)
tot.wss <- 0
for (i in 1:15) {
  tot.wss[i] <- kmeans(data_GPI[, -c(1)], centers = i)$tot.withinss
}
# Plotting cluster plot to select K
plot(1:15,
     tot.wss,
     type = 'b',
     xlab = "Number of GPI Clusters",
     ylab = "Total Within Group Sum of Squares")
# Removing unused references
rm(i, tot.wss)


# Kmeans with 4 cluster
km_GPI <- kmeans(data_GPI[, -c(1)], 4, nstart = 25)
km_GPI$cluster

data_GPI <-
  data.frame(data_GPI, "GPI_Cluster" = as.factor(km_GPI$cluster))

# Doing innerjoin with original data by data_GPI
join_GPI <- inner_join(join_Drug, data_GPI, by = "GPI")
join_GPI <- subset(data.frame(join_GPI), select = -c(False, True, GPI))


# Removing unused references
rm(dummay_GPI)
rm(join_Drug)
rm(data_GPI)
rm(km_GPI)

##########################################################

# Attribute "NDC" w.r.t Target

data_NDC <- select(data, NDC)
dummy_NDC <- dummy(data$Target)
data_NDC <- data.frame("NDC" = data$NDC, dummy_NDC)
data_NDC <-
  data_NDC %>% group_by(NDC) %>% summarise("False" = sum(Target0), "True" = sum(Target1))

# Random no.generator to restore same(constant result everytime)
set.seed(1234)
tot.wss <- 0
for (i in 1:15) {
  tot.wss[i] <- kmeans(data_NDC[, -c(1)], centers = i)$tot.withinss
}
# Plotting cluster plot to select K
plot(1:15,
     tot.wss,
     type = 'b',
     xlab = "Number of NDC Clusters",
     ylab = "Total Within Group Sum of Squares")
# Removing unused references
rm(i, tot.wss)


# Kmeans with 4 cluster
km_NDC <- kmeans(data_NDC[, -c(1)], 4, nstart = 25)
km_NDC$cluster

data_NDC <-
  data.frame(data_NDC, "NDC_Cluster" = as.factor(km_NDC$cluster))

# Doing innerjoin with original data by data_GPI
join_NDC <- inner_join(join_GPI, data_NDC, by = "NDC")
join_NDC <- subset(data.frame(join_NDC), select = -c(False, True, NDC))


# Removing unused references
rm(dummy_NDC)
rm(join_GPI)
rm(data_NDC)
rm(km_NDC)

##########################################################

# Attribute RxGroupId w.r.t Target

data_RxGroupId <- select(data, RxGroupId)
dummy_RxGroupId <- dummy(data$Target)
data_RxGroupId <-
  data.frame("RxGroupId" = data$RxGroupId, dummy_RxGroupId)
data_RxGroupId <-
  data_RxGroupId %>% group_by(RxGroupId) %>% summarise("False" = sum(Target0), "True" = sum(Target1))

# Random no.generator to restore same(constant result everytime)
set.seed(1234)
tot.wss <- 0
for (i in 1:15) {
  tot.wss[i] <- kmeans(data_RxGroupId[, -c(1)], centers = i)$tot.withinss
}
# Plotting cluster plot to select K
plot(1:15,
     tot.wss,
     type = 'b',
     xlab = "Number of RxGroupId Clusters",
     ylab = "Total Within Group Sum of Squares")
# Removing unused references
rm(i, tot.wss)


# Kmeans with 4 cluster
km_RxGroupId <- kmeans(data_RxGroupId[, -c(1)], 4, nstart = 25)
km_RxGroupId$cluster

data_RxGroupId = data.frame(data_RxGroupId, "RxGroupId_Cluster" = as.factor(km_RxGroupId$cluster))

# Doing innerjoin with original data by data_GPI
join_RxGroupId <-
  inner_join(join_NDC, data_RxGroupId, by = "RxGroupId")
join_RxGroupId <-
  subset(data.frame(join_RxGroupId), select = -c(False, True, RxGroupId))


# Removing unused references
rm(dummy_RxGroupId)
rm(join_NDC)
rm(data_RxGroupId)
rm(km_RxGroupId)

##########################################################

# Attribute Bin w.r.t Target

data_Bin <- select(data, Bin)
dummy_Bin <- dummy(data$Target)
data_Bin <- data.frame("Bin" = data$Bin, dummy_Bin)
data_Bin <-
  data_Bin %>% group_by(Bin) %>% summarise("False" = sum(Target0), "True" = sum(Target1))

# Random no.generator to restore same(constant result everytime)
set.seed(1234)
tot.wss <- 0
for (i in 1:15) {
  tot.wss[i] <- kmeans(data_Bin[, -c(1)], centers = i)$tot.withinss
}
# Plotting cluster plot to select K
plot(1:15,
     tot.wss,
     type = 'b',
     xlab = "Number of Bin Clusters",
     ylab = "Total Within Group Sum of Squares")
# Removing unused references
rm(i, tot.wss)


# Kmeans with 4 cluster
km_Bin <- kmeans(data_Bin[, -c(1)], 4, nstart = 25)
km_Bin$cluster

data_Bin = data.frame(data_Bin, "Bin_Cluster" = as.factor(km_Bin$cluster))

# Doing innerjoin with original data by data_GPI
join_Bin <- inner_join(join_RxGroupId, data_Bin, by = "Bin")
join_Bin <- subset(data.frame(join_Bin), select = -c(False, True, Bin))


# Removing unused references
rm(dummy_Bin)
rm(join_RxGroupId)
rm(data_Bin)
rm(km_Bin)

##########################################################

# Attribute PCN w.r.t Target

data_PCN <- select(data, PCN)
dummy_PCN <- dummy(data$Target)
data_PCN <- data.frame("PCN" = data$PCN, dummy_PCN)
data_PCN <-
  data_PCN %>% group_by(PCN) %>% summarise("False" = sum(Target0), "True" = sum(Target1))

# Random no.generator to restore same(constant result everytime)
set.seed(1234)
tot.wss <- 0
for (i in 1:15) {
  tot.wss[i] <- kmeans(data_PCN[, -c(1)], centers = i)$tot.withinss
}
# Plotting cluster plot to select K
plot(1:15,
     tot.wss,
     type = 'b',
     xlab = "Number of PCN Clusters",
     ylab = "Total Within Group Sum of Squares")
# Removing unused references
rm(i, tot.wss)


# Kmeans with 4 cluster
km_PCN <- kmeans(data_PCN[, -c(1)], 4, nstart = 25)
km_PCN$cluster

data_PCN <-
  data.frame(data_PCN, "PCN_Cluster" = as.factor(km_PCN$cluster))

# Doing innerjoin with original data by data_GPI
join_PCN <- inner_join(join_Bin, data_PCN, by = "PCN")
join_PCN <- subset(data.frame(join_PCN), select = -c(False, True, PCN))


# Removing unused references
rm(dummy_PCN)
rm(join_Bin)
rm(data_PCN)
rm(km_PCN)

##########################################################

# Attribute State w.r.t Target

data_State <- select(data, State)
dummy_State <- dummy(data$Target)
data_State <- data.frame("State" = data$State, dummy_State)
data_State <-
  data_State %>% group_by(State) %>% summarise("False" = sum(Target0), "True" = sum(Target1))

# Random no.generator to restore same(constant result everytime)
set.seed(1234)
tot.wss <- 0
for (i in 1:15) {
  tot.wss[i] <- kmeans(data_State[, -c(1)], centers = i)$tot.withinss
}
# Plotting cluster plot to select K
plot(1:15,
     tot.wss,
     type = 'b',
     xlab = "Number of State Clusters",
     ylab = "Total Within Group Sum of Squares")
# Removing unused references
rm(i, tot.wss)


# Kmeans with 4 cluster
km_State <- kmeans(data_State[, -c(1)], 4, nstart = 25)
km_State$cluster

data_State <-
  data.frame(data_State, "State_Cluster" = as.factor(km_State$cluster))

# Doing innerjoin with original data by data_GPI
join_State <- inner_join(join_PCN, data_State, by = "State")
join_State <-
  subset(data.frame(join_State), select = -c(False, True, State))


# Removing unused references
rm(dummy_State)
rm(join_PCN)
rm(data_State)
rm(km_State)

##########################################################

# Attribute DrugClass w.r.t Target

data_DrugClass <- select(data, DrugClass)
dummy_DrugClass <- dummy(data$Target)
data_DrugClass <-
  data.frame("DrugClass" = data$DrugClass, dummy_DrugClass)
data_DrugClass <-
  data_DrugClass %>% group_by(DrugClass) %>% summarise("False" = sum(Target0), "True" = sum(Target1))

# Random no.generator to restore same(constant result everytime)
set.seed(1234)
tot.wss <- 0
for (i in 1:15) {
  tot.wss[i] <- kmeans(data_DrugClass[, -c(1)], centers = i)$tot.withinss
}
plot(1:15,
     tot.wss,
     type = 'b',
     xlab = "Number of DrugClass Clusters",
     ylab = "Total Within Group Sum of Squares")
# Removing unused references
rm(i, tot.wss)


# Kmeans with 4 cluster
km_DrugClass <- kmeans(data_DrugClass[, -c(1)], 4, nstart = 25)
km_DrugClass$cluster

data_DrugClass = data.frame(data_DrugClass, "DrugClass_Cluster" = as.factor(km_DrugClass$cluster))

# Doing innerjoin with original data by data_GPI
join_DrugClass <-
  inner_join(join_State, data_DrugClass, by = "DrugClass")
join_DrugClass <-
  subset(data.frame(join_DrugClass), select = -c(False, True, DrugClass))


# Removing unused references
rm(dummy_DrugClass)
rm(join_State)
rm(data_DrugClass)
rm(km_DrugClass)

##########################################################

# Attribute DrugSubClass w.r.t Target

data_DrugSubClass <- select(data, DrugSubClass)
dummy_DrugSubClass <- dummy(data$Target)
data_DrugSubClass <-
  data.frame("DrugSubClass" = data_DrugSubClass, dummy_DrugSubClass)
data_DrugSubClass <-
  data_DrugSubClass %>% group_by(DrugSubClass) %>% summarise("False" = sum(Target0), "True" = sum(Target1))

# Random no.generator to restore same(constant result everytime)
set.seed(1234)
tot.wss <- 0
for (i in 1:15) {
  tot.wss[i] <-
    kmeans(data_DrugSubClass[, -c(1)], centers = i)$tot.withinss
}
plot(1:15,
     tot.wss,
     type = 'b',
     xlab = "Number of DrugSubClass Clusters",
     ylab = "Total Within Group Sum of Squares")
# Removing unused references
rm(i, tot.wss)


# Kmeans with 4 cluster
km_DrugSubClass <- kmeans(data_DrugSubClass[, -c(1)], 4, nstart = 25)
km_DrugSubClass$cluster

data_DrugSubClass = data.frame(data_DrugSubClass,
                               "DrugSubClass_Cluster" = as.factor(km_DrugSubClass$cluster))

# Doing innerjoin with original data by data_GPI
join_DrugSubClass <-
  inner_join(join_DrugClass, data_DrugSubClass, by = "DrugSubClass")
join_DrugSubClass <-
  subset(data.frame(join_DrugSubClass),
         select = -c(False, True, DrugSubClass))


# Removing unused references
rm(dummy_DrugSubClass)
rm(join_DrugClass)
rm(data_DrugSubClass)
rm(km_DrugSubClass)

##########################################################

# Attribute DrugChemicalName w.r.t Target

data_DrugChemical <- select(data, Drug_Chemical_Name)
dummy_Drug_Chemical_Name <- dummy(data$Target)
data_DrugChemical <-
  data.frame("Drug_Chemical_Name" = data$Drug_Chemical_Name,
             dummy_Drug_Chemical_Name)
data_DrugChemical <-
  data_DrugChemical %>% group_by(Drug_Chemical_Name) %>% summarise("False" = sum(Target0), "True" = sum(Target1))

# Random no.generator to restore same(constant result everytime)
set.seed(1234)

tot.wss <- 0
for (i in 1:15) {
  tot.wss[i] <-
    kmeans(data_DrugChemical[, -c(1)], centers = i)$tot.withinss
}
# Plotting cluster plot to select K
plot(1:15,
     tot.wss,
     type = 'b',
     xlab = "Number of Drug_Chemical Clusters",
     ylab = "Total Within Group Sum of Squares")
# Removing unused references
rm(i, tot.wss)


# Kmeans with 4 cluster
km_Drug_Chemical <- kmeans(data_DrugChemical[, -c(1)], 4, nstart = 25)
km_Drug_Chemical$cluster

data_DrugChemical <-
  data.frame(data_DrugChemical,
             "Drug_Chemical_Cluster" = as.factor(km_Drug_Chemical$cluster))

# Doing innerjoin with original data by data_GPI
join_Drug_Chemical <-
  inner_join(join_DrugSubClass, data_DrugChemical, by = "Drug_Chemical_Name")
join_Drug_Chemical <-
  subset(data.frame(join_Drug_Chemical),
         select = -c(False, True, Drug_Chemical_Name))

# Removing unused references
rm(dummy_Drug_Chemical_Name)
rm(join_DrugSubClass)
rm(km_Drug_Chemical)
rm(data_DrugChemical)

##########################################################

# Attribute DrugGroup w.r.t Target

data_DrugGroup <- select(data, DrugGroup)
dummy_DrugGroup <- dummy(data$Target)
data_DrugGroup <-
  data.frame("DrugGroup" = data$DrugGroup, dummy_DrugGroup)
data_DrugGroup <-
  data_DrugGroup %>% group_by(DrugGroup) %>% summarise("False" = sum(Target0), "True" = sum(Target1))

# Random no.generator to restore same(constant result everytime)
set.seed(1234)

tot.wss <- 0
for (i in 1:15) {
  tot.wss[i] <- kmeans(data_DrugGroup[, -c(1)], centers = i)$tot.withinss
}
# Plotting cluster plot to select K
plot(1:15,
     tot.wss,
     type = 'b',
     xlab = "Number of DrugGroup Clusters",
     ylab = "Total Within Group Sum of Squares")
# Removing unused references
rm(i, tot.wss)


# Kmeans with 4 cluster
km_DrugGroup <- kmeans(data_DrugGroup[, -c(1)], 4, nstart = 25)
km_DrugGroup$cluster

data_DrugGroup <-
  data.frame(data_DrugGroup, "DrugGroup_Cluster" = as.factor(km_DrugGroup$cluster))

# Doing innerjoin with original data by data_GPI
join_DrugGroup <-
  inner_join(join_Drug_Chemical, data_DrugGroup, by = "DrugGroup")
cluster_Data <-
  subset(data.frame(join_DrugGroup), select = -c(False, True, DrugGroup))

# Removing unused references
rm(dummy_DrugGroup)
rm(join_DrugGroup)
rm(join_Drug_Chemical)
rm(km_DrugGroup)
rm(data_DrugGroup)

##########################################################

#Checking the structure of the clustered data
str(cluster_Data)
cluster_Data <- as.data.frame(unclass(cluster_Data))
str(cluster_Data)

# Download required library for date format
library(lubridate)

# Change Transdate format
cluster_Data$TransDate <- mdy(cluster_Data$TransDate)

# Random no.generator to restore same(constant result everytime)
set.seed(1234)
# Time-based sampling for Train(91%) and Test(9%) data (Tried to fit in 80:20)
#train_Data <- subset(cluster_Data,cluster_Data$TransDate>="2013-07-22" & cluster_Data$TransDate<="2013-08-06")
# Time-based sampling for Train(77%) and Test(23%) data
train_Data <-
  subset(
    cluster_Data,
    cluster_Data$TransDate >= "2013-07-22" &
      cluster_Data$TransDate <= "2013-08-04"
  )

# Removing Transdate and GPI_Cluster from Train and Test data
train_Data <-
  subset(data.frame(train_Data), select = -c(TransDate, GPI_Cluster))
train_Data$Target <- factor(train_Data$Target)
train_row <- row.names(train_Data)
test_row <- setdiff(row.names(cluster_Data), train_row)
test_Data <- cluster_Data[test_row, ]
test_Data <-
  subset(data.frame(test_Data), select = -c(TransDate, GPI_Cluster))

# Removing unused references
rm(cluster_Data, train_row, test_row)

# Train data without target
train_NoTarget <-
  subset(data.frame(train_Data), select = -c(Target))
# test data without target
test_NoTarget <- subset(data.frame(test_Data), select = -c(Target))

# Written to CSV to verify the records taken for Train and Test  data
write.csv(train_Data,
          "F:/INSOFE/Course Material/Project/train_data.csv")
write.csv(test_Data, "F:/INSOFE/Course Material/Project/test_data.csv")


# Stratified sampling
# train_RowIDs <- cluster_Data(1:nrow(dt), nrow(dt)*0.7)
# train_Data <- dt[train_RowIDs,]
# str(train_Data)
# test_Data <- dt[-train_RowIDs,]

##########################################################

###----------------Ensemble:Stacking--------------------

# Download required libraries for diferent models
library(randomForest)
library(rpart)
library(C50)

# Random no.generator to restore same(constant result everytime)
set.seed(1234)
# Build randomforest model on the training dataset
rf_Model <-
  randomForest(
    train_Data$Target ~ .,
    data = train_NoTarget,
    keep.forest = T,
    ntree = 50
  )
summary(rf_Model)

print(rf_Model)

rf_Model$importance
round(importance(rf_Model), 2)

# plot tree: Visualizing RF model by variable importance
varImpPlot(rf_Model)


# Build CART model on the training dataset
cart_Model <-
  rpart(train_Data$Target ~ ., train_NoTarget, method = "class")
summary(cart_Model)

# plot tree: Visualizing CART model decision tree
plot(cart_Model, uniform = TRUE,
     main = "Classification Tree for PreAuth")
text(cart_Model,
     use.n = TRUE,
     all = TRUE,
     cex = .8)


# Build C5.0 model on the training dataset
# Tree is decomposed into a rule-based model
#c50_Model <- C5.0(train_Data$Target ~ ., train_NoTarget, rules = T)

c50_Model <- C5.0(train_Data$Target ~ ., train_NoTarget, rules = F)
summary(c50_Model)

# plot tree: Visualizing C50 model decision tree
plot(c50_Model)


# Build Logistic regression on the training dataset
glm_Model <-
  glm(train_Data$Target ~ ., train_NoTarget, family = binomial)
summary(glm_Model)

# plot tree: Visualizing regression
plot(glm_Model)
abline(glm_Model)


##########################################################

#---------------Predict on Train Data(without target)--------------------
# Using randomforest Model predicting with the train dataset
rf_Train <-
  predict(rf_Model,
          train_NoTarget,
          type = "response",
          norm.votes = TRUE)

# Build confusion matrix to find accuracy and recall using randomforest model
cm_rfTrain  <-
  table("actual" = train_Data$Target, "predicted" = rf_Train)
cm_rfTrain
#accuracy <- sum(number of correctly classified instances per class)/number of instances
accu_rfTrain <- sum(diag(cm_rfTrain)) / sum(cm_rfTrain)
accu_rfTrain
#recall = number of correctly classified instances per class/number of instances per class
recall_rfTrain <-  cm_rfTrain[2, 2] / (sum(cm_rfTrain[2, ]))
recall_rfTrain

#--------------------
# Using CART Model predict on train data
cart_Train <- predict(cart_Model, train_NoTarget, type = "vector")
table(cart_Train)

# if we choose type=vector, then replace 1 with 0 and 2 with 1
cart_Train <- ifelse(cart_Train == 1, 0, 1)
table(cart_Train)

# Build confusion matrix to find accuracy and recall using CART model
cm_cartTrain <-
  table("actual" = train_Data$Target, "predicted" = cart_Train)
cm_cartTrain
# accuracy <- sum(number of correctly classified instances per class)/number of instances
accu_cartTrain <- sum(diag(cm_cartTrain)) / sum(cm_cartTrain)
accu_cartTrain
# recall = number of correctly classified instances per class/number of instances per class
recall_cartTrain <- cm_cartTrain[2, 2] / (sum(cm_cartTrain[2, ]))
recall_cartTrain

#--------------------
# Using C5.0 Model predicting with the train dataset
c50_Train <- predict(c50_Model, train_NoTarget, type = "class")
c50_Train <- as.vector(c50_Train)
table(c50_Train)

# Build confusion matrix to find accuracy and recall using C5.0 model
cm_c50Train <-
  table("actual" = train_Data$Target, "predicted" = c50_Train)
cm_c50Train
# accuracy <- sum(number of correctly classified instances per class)/number of instances
accu_c50Train <- sum(diag(cm_c50Train)) / sum(cm_c50Train)
accu_c50Train
# recall = number of correctly classified instances per class/number of instances per class
recall_c50Train <-  cm_c50Train[2, 2] / (sum(cm_c50Train[2, ]))
recall_c50Train

#--------------------
# Using GLM Model predicting on train dataset
glm_Train = predict(glm_Model, type = "response")
#it gives probabilities, so we #need to convert to 1's and 0's;
# if >0.5 show as 1 or else show as 0.
glm_Train <- ifelse(glm_Train > 0.5, 1, 0)
table(glm_Train)

# Build confusion matrix to find accuracy and recall using GLM Model
#cm_glmTrain <- table(glm_Train, train_Data$Target)
cm_glmTrain <-
  table("actual" = train_Data$Target, "predicted" = glm_Train)
cm_glmTrain
# accuracy <- sum(number of correctly classified instances per class)/number of instances
accu_glmTrain <- sum(diag(cm_glmTrain)) / sum(cm_glmTrain)
accu_glmTrain
# recall = number of correctly classified instances per class/number of instances per class
recall_glmTrain <- cm_glmTrain[2, 2] / (sum(cm_glmTrain[2, ]))
recall_glmTrain

# Removing unused references
rm(cm_rfTrain, cm_cartTrain, cm_c50Train, cm_glmTrain)


#----------------------------Ensemble with Train data-----------------#

# Combining only 3 algorithm for Train predictions: CART, C5.0 & Logistic Regression
PredAll_TrainModels <- data.frame(CART = cart_Train,
                                  C50 = c50_Train,
                                  RF = rf_Train)
# Doing factoring
PredAll_TrainModels <-
  data.frame(sapply(PredAll_TrainModels, as.factor))

# Summary and structure of all train data models together
str(PredAll_TrainModels)
summary(PredAll_TrainModels)

# Removing unused references
rm(cart_Train, rf_Train, c50_Train)

# Viewing the predictions of each modelc50_Train
table(PredAll_TrainModels$CART) #CART
table(PredAll_TrainModels$C50)  #C5.0
table(PredAll_TrainModels$RF)  #Logistic Regression
table(train_Data$Target) #Original Dataset DV

# Adding the original DV to the dataframe
PredAll_TrainModels <-
  cbind(PredAll_TrainModels, Target = train_Data$Target)

# Ensemble Model with GLM as Meta Learner
str(PredAll_TrainModels)
head(PredAll_TrainModels)
ensemble_Model_Train <-
  glm(train_Data$Target ~ ., PredAll_TrainModels, family = binomial)
summary(ensemble_Model_Train)

# Check the "ensemble_Model model" on the train data
ensemble_Train <-
  predict(ensemble_Model_Train, PredAll_TrainModels,
          type = "response")
ensemble_Train <- ifelse(ensemble_Train > 0.5, 1, 0)
table(ensemble_Train)

# Build confusion matrix to find accuracy and recall using glm_ensemble model : On Train data
cm_EnsembleTrain <-
  table("actual" = PredAll_TrainModels$Target, "predicted" = ensemble_Train)
cm_EnsembleTrain
# accuracy <- sum(number of correctly classified instances per class)/number of instances
accu_EnsembleTrain <-
  sum(diag(cm_EnsembleTrain)) / sum(cm_EnsembleTrain)
accu_EnsembleTrain
# recall = number of correctly classified instances per class/number of instances per class
recall_EnsembleTrain <-
  cm_EnsembleTrain[2, 2] / (sum(cm_EnsembleTrain[2, ]))
recall_EnsembleTrain


###########################################################

#---------------Predict on Test Data(without target)--------------------



# Using randomforest Model prediction on test dataset
rf_Test <-
  predict(rf_Model,
          test_NoTarget,
          type = "response",
          norm.votes = TRUE)

# Build confusion matrix to find accuracy and recall
cm_rfTest <-
  table("actual" = test_Data$Target, "predicted" = rf_Test)

cm_rfTest
#accuracy <- sum(number of correctly classified instances per class)/number of instances
accu_rfTest <- sum(diag(cm_rfTest)) / sum(cm_rfTest)
accu_rfTest
#recall = number of correctly classified instances per class/number of instances per class
recall_rfTest <-  cm_rfTest[2, 2] / (sum(cm_rfTest[2, ]))
recall_rfTest

#--------------------
# Using CART Model prediction on test dataset
cart_Test <- predict(cart_Model, test_NoTarget, type = "vector")
cart_Test <- ifelse(cart_Test == 1, 0, 1)

# Build confusion matrix to find accuracy and recall using CART model
cm_cartTest <-
  table("actual" = test_Data$Target, "predicted" = cart_Test)
cm_cartTest
# accuracy <- sum(number of correctly classified instances per class)/number of instances
accu_cartTest <- sum(diag(cm_cartTest)) / sum(cm_cartTest)
accu_cartTest
# recall = number of correctly classified instances per class/number of instances per class
recall_cartTest <-  cm_cartTest[2, 2] / (sum(cm_cartTest[2, ]))
recall_cartTest
#--------------------
# Using C50 Model prediction on test dataset
c50_Test <- predict(c50_Model, test_NoTarget, type = "class")
c50_Test <- as.vector(c50_Test)

# Build confusion matrix to find accuracy and recall using C50 model
#cm_C50Test <- table(c50_Test, test_Data$Target)
cm_C50Test <-
  table("actual" = test_Data$Target, "predicted" = c50_Test)
cm_C50Test
# accuracy <- sum(number of correctly classified instances per class)/number of instances
accu_c50Test <- sum(diag(cm_C50Test)) / sum(cm_C50Test)
accu_c50Test
# recall = number of correctly classified instances per class/number of instances per class
recall_c50Test <-  cm_C50Test[2, 2] / (sum(cm_C50Test[2, ]))
recall_c50Test


#--------------------
# Using GLM Model prediction on test dataset
glm_Test <- predict(glm_Model, test_NoTarget, type = "response")
glm_Test <- ifelse(glm_Test > 0.5, 1, 0)

# Build confusion matrix to find accuracy and recall using GLM Model
cm_glmTest <-
  table("actual" = test_Data$Target, "predicted" = glm_Test)
cm_glmTest
# accuracy <- sum(number of correctly classified instances per class)/number of instances
accu_glmTest <- sum(diag(cm_glmTest)) / sum(cm_glmTest)
accu_glmTest
# recall = number of correctly classified instances per class/number of instances per class
recall_glmtTest <-  cm_glmTest[2, 2] / (sum(cm_glmTest[2, ]))
recall_glmtTest

# Removing unused references
rm(cm_rfTest, cm_cartTest, cm_C50Test, cm_glmTest)
rm(rf_Model, cart_Model, c50_Model, glm_Model)

#----------------------------Ensemble with Test data-----------------#
# Combining test predictions of CART, C5.0 & Log Regression together
PredAll_TestModels <- data.frame(CART = cart_Test,
                                 C50 = c50_Test,
                                 RF = rf_Test)
# Removing unused references
rm(cart_Test, c50_Test, rf_Test)

# Viewing the predictions of each test model
table(PredAll_TestModels$CART) #CART
table(PredAll_TestModels$C50)  #C5.0
table(PredAll_TestModels$RF)  #Logistic Regression
table(test_Data$Target) #Original Dataset DV


# Doing factoring
PredAll_TestModels <-
  data.frame(sapply(PredAll_TestModels, as.factor))

# Adding the original DV to the dataframe
PredAll_TestModels <-
  cbind(PredAll_TestModels, Target = test_Data$Target)

# Ensemble Model with GLM as Meta Learner
str(PredAll_TestModels)
head(PredAll_TestModels)
ensemble_Model_Test <-
  glm(test_Data$Target ~ ., PredAll_TestModels, family = binomial)
summary(ensemble_Model_Test)

# Combiner "glm_ensemble model" on the test data
ensemble_Test <-
  predict(ensemble_Model_Test, PredAll_TestModels, type = "response")
ensemble_Test <- ifelse(ensemble_Test > 0.5, 1, 0)
table(ensemble_Test)

# Build confusion matrix to find accuracy and recall using glm_ensemble model : On Test data
cm_EnsembleTest <-
  table("actual" = test_Data$Target, "predicted" = ensemble_Test)
cm_EnsembleTest
# accuracy <- sum(number of correctly classified instances per class)/number of instances
accu_EnsembleTest <- sum(diag(cm_EnsembleTest)) / sum(cm_EnsembleTest)
accu_EnsembleTest
# recall = number of correctly classified instances per class/number of instances per class
recall_EnsembleTest <-
  cm_EnsembleTest[2, 2] / (sum(cm_EnsembleTest[2, ]))
recall_EnsembleTest

# Removing unused references
rm(train_Data, test_Data, train_NoTarget, test_NoTarget)

# Saving ensemble models
saveRDS(ensemble_Model_Train, file = "F:/INSOFE/Course Material/Project/ensemble_Train_save.rds")
saveRDS(ensemble_Model_Test, file = "F:/INSOFE/Course Material/Project/ensemble_Test_save.rds")
