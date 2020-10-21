# Paulina Kurowska

# -------------------------- Removing all and closes the windows ---------------------
rm(list = ls())
graphics.off()
# -------------------------------- Packages ------------------------------------------
library(TH.data)
library(mboost)
library(ggplot2)
library(l2boost)
library(caret)
library(rpart)
library(rpart.plot)
# ---------------------------------- Data loading -----------------------------------

data("bodyfat", package = "TH.data")     
head(bodyfat)

# data partition into trainset and testset
# trainset will be used to build=train the models
# testset will act as a set for performance validation (out of sample)
sample_size <- as.integer(0.8*nrow(bodyfat))
set.seed(2018) # reproducibility for sample()
train <- sample(1:nrow(bodyfat), size = sample_size, replace = F)
trainset <- bodyfat[train,]
testset <- bodyfat[-train,]

# ----------------------------------OLS------------------------------------------

# ols estimation
ols_model <- lm(DEXfat ~ hipcirc  + kneebreadth + anthro3a, data = trainset)
coef(ols_model)
ols_ssr <- sum(resid(ols_model)^2) # 664.5877
ols_mse <- (1/nrow(trainset))*ols_ssr
ols_rmse <- sqrt(ols_mse) # 3.44

# ----------------------------- Boosting applied --------------------------------
# with boosting
# the results below were used in the presentation
boost_model <-glmboost(DEXfat ~ hipcirc  + kneebreadth + anthro3a, data = trainset)
coef(boost_model)
sum(resid(boost_model)^2)
predict_boost <- predict(boost_model,trainset)
ssr_boost <- sum((predict_boost - trainset$DEXfat)^2) #664.588
boost_mse <- (1/nrow(trainset))*ssr_boost
boost_rmse <- sqrt(boost_mse)
#plot
plot(boost_model)
#plot
plot(boost_model, ylim=(range(coef(boost_model, 
                                   which = c("hipcirc","kneebreadth","anthro3a")))))
# -------------------------------- AIC --------------------------------------
AIC(ols_model)

aic<-AIC(boost_model)
mstop_aic <- mstop(aic)
plot(aic)
# par(mai = par("mai") * c(1, 1, 1, 1.8), mfrow=c(1,2))
par(mfrow=c(1,1))
plot(boost_model, ylim=(range(coef(boost_model, 
                                   which = c("hipcirc","kneebreadth","anthro3a")))))
abline(v=mstop_aic, col="green")
abline(v=mstop_aic,col="red")

boost_model_aic<-glmboost(DEXfat ~ hipcirc  + kneebreadth + anthro3a, data = trainset, control = boost_control(mstop=mstop_aic))
boost_model_aic
# calculating ssr, mse, rmse
aic_boost_ssr <- sum(resid(boost_model_aic)^2)
aic_boost_mse <- (1/nrow(trainset))*aic_boost_ssr
aic_boost_rmse <- sqrt(aic_boost_mse) #3.449346

set.seed(15091987)
cv<-cvrisk(boost_model )
mstop_cv <- mstop(cv)
boost_model_cv<-glmboost(DEXfat ~ hipcirc  + kneebreadth + anthro3a, data = trainset, control = boost_control(mstop = mstop_cv ))
coef(boost_model_cv)
sum(resid(boost_model_cv)^2)

# ---------------------------- ALL Variable for Report ----------------------------
# all variables 
boost_model_all <-glmboost(DEXfat ~ ., data = trainset)
boost_model_all
plot(boost_model_all)
plot(boost_model_all, ylim=(range(coef(boost_model_all, 
        which = names(bodyfat)[-2]))))
ols_model_all<-lm(DEXfat ~ ., data = trainset)

ols_model_all
summary(ols_model_all)
ssr_ols_train <- sum(resid(ols_model_all)^2)# 510.8855, smaller ssr
mse_ols_train <- (1/length(trainset$DEXfat))*ssr_ols_train # 9.122956
rmse_ols_train <- sqrt(mse_ols_train) # 3.020423

ssr_boost_train <- sum(resid(boost_model_all)^2) # 522.46
mse_boost_train <- (1/length(trainset$DEXfat))*ssr_boost_train # 9.329643
rmse_boost_train <- sqrt(mse_boost_train) # 3.054446

ols_model_corrected <- step(ols_model_all)
summary(ols_model_corrected)
ssr_ols_corrected <- sum(ols_model_corrected$residuals^2) # 534.6465
mse_ols_corrected <- (1/length(trainset$DEXfat))*ssr_ols_corrected
rmse_ols_corrected <- sqrt(mse_ols_corrected) #3.089864

# how does it perform on a testset?
predict_boost_all <- predict(boost_model_all,testset)
ssr_boost_all <- sum((predict_boost_all - testset$DEXfat)^2) #152.0508
mse_boost_all <- (1/length(testset$DEXfat))*ssr_boost_all #10.13672
rmse_boost_all <- sqrt(mse_boost_all)#3.183821

predict_ols_all <- predict(ols_model_all,testset)
ssr_ols_all <- sum((predict_ols_all - testset$DEXfat)^2) #155.8104
mse_ols_all <- (1/length(testset$DEXfat))*ssr_ols_all #10.38736
rmse_ols_all <- sqrt(mse_ols_all) # 3.222943

ssr_ols_corrected_test <- sum((predict(ols_model_corrected, testset)-testset$DEXfat)^2) # 150.6031
mse_ols_corrected_test <- (1/length(testset$DEXfat))*sum(ssr_ols_corrected_test) # 10.04021
rmse_ols_corrected_test <- sqrt(mse_ols_corrected_test) # 3.168629
# conclusion boosting is better on a testset, but has a slightly poorer performance on a trainset

summary(step(ols_model_all)) # the same as indicates the p-values

# AIC and CV optimalization
aic_boost_all <- AIC(boost_model_all)
aic_mstop_all <- mstop(aic_boost_all) # 48
set.seed(1509)
cv_setting <- cv(model.weights(boost_model_all), type= "kfold")
cv_boost_all <- cvrisk(boost_model_all, folds =cv_setting)
cv_mstop_all <- mstop(cv_boost_all) #54

variables <- names(trainset)
variables <- variables[variables!="DEXfat"]
par(mfrow=c(1,2))
plot(boost_model_all, ylim=(range(coef(boost_model_all, 
                                   which = c(variables)))), main = paste("Gradient Boosting ", "\n","for bodyfat"))
abline(v=aic_mstop_all, col="green")
abline(v=cv_mstop_all,col="red")
# Diagnostics: ssr, mse, rmse for train and for test
# mstop according to cv = 91, metrics for trainset
boost_model_cv_train <- glmboost(DEXfat ~ ., data = trainset, control = boost_control(mstop = cv_mstop_all)) 
ssr_boost_all_train_cv <- sum((boost_model_cv_train$resid())^2)
mse_boost_all_train_cv <- (1/length(trainset$DEXfat))*ssr_boost_all_train_cv#9.339361
rmse_boost_all_train_cv <- sqrt(mse_boost_all_train_cv)
# the same for testset
ssr_boost_all_test_cv <- sum((predict(boost_model_cv_train, testset)-testset$DEXfat)^2)
mse_boost_all_test_cv <- (1/length(testset$DEXfat))*ssr_boost_all_test_cv #10.08784
rmse_boost_all_test_cv <- sqrt(mse_boost_all_test_cv)#3.176137

# now for aic_mstop
# trainset, nu = 0.1
boost_model_aic_train <- glmboost(DEXfat ~ ., data = trainset, control = boost_control(mstop = aic_mstop_all)) 
ssr_boost_all_train_aic <- sum((boost_model_aic_train$resid())^2)
mse_boost_all_train_aic <- (1/length(trainset$DEXfat))*ssr_boost_all_train_aic#9.504008
rmse_boost_all_train_aic <- sqrt(mse_boost_all_train_aic)#3.082857
# testset
ssr_boost_all_test_aic <- sum((predict(boost_model_aic_train,testset)-testset$DEXfat)^2)
mse_boost_all_test_aic <- (1/length(testset$DEXfat))*ssr_boost_all_test_aic #5.882201
rmse_boost_all_test_aic <- sqrt(mse_boost_all_test_aic)#2.425325

# blackboost

blackboost_model<- blackboost(DEXfat ~ ., data = trainset)
blackboost_model2<- blackboost(DEXfat ~ ., data = trainset, control = boost_control(mstop = 25, nu=0.2))
predict
set.seed(15091987)
blackboost_cv <- cvrisk(blackboost_model)#48


# ---------------------- Regression tree and Random Forest (for comparison purposes)-------------
tree<-rpart(DEXfat ~ ., data = trainset)
par(mfrow=c(1,1))
rpart.plot(tree)
ssr_tree <- sum((residuals(tree))^2)# 1133.72
mse_tree <- (1/length(trainset$DEXfat))*ssr_tree # 20.245
rmse_tree <- sqrt(mse_tree) # 4.499445

tree_predict <- predict(tree, testset)
ssr_tree_test <- sum((tree_predict-testset$DEXfat)^2)# 381.9813
mse_tree_test <- (1/length(testset$DEXfat))*ssr_tree_test # 25.46542
rmse_tree_test <- sqrt(mse_tree_test) # 5.046327
# boosted decision tree with caret
library(caret)
# based on https://www.youtube.com/watch?v=GZafoSe5DvI
# settings for cross-validation : 5-fold cross validation, repeated 3 times
trainctrl <- trainControl(method= "repeatedcv", number = 5, repeats = 3)
# training tree in caret
set.seed(15091987)
caret_tree <- train( DEXfat~ ., data = trainset, method= "rpart", trControl =trainctrl)
predict_caret_tree <- predict(caret_tree, trainset)
rsq_caret_tree <- cor(predict_caret_tree, trainset$DEXfat)^2 # 0.6430407
ssr_caret_tree <- sum(residuals(caret_tree)^2)
mse_caret_tree <- (1/length(trainset$DEXfat))*ssr_caret_tree#45.10273
rmse_caret_tree <- sqrt(mse_caret_tree)# 6.715857
# random forest in caret
caret_randomforest <- train( DEXfat~ ., data = trainset, method= "rf", trControl =trainctrl)
compare <- resamples(list(rpart_tree= caret_tree, randomForest=caret_randomforest))
rmse_compare <- summary(compare)$statistics$RMSE[,4] #we are interested in mean RMSE
rsq_compare <- summary(compare)$statistics$Rsquared[,4] #we are interested in mean RSQ

#testset
predict_caret_tree_test <- predict(caret_tree, testset)
ssr_caret_tree_test <- sum((predict_caret_tree_test-testset$DEXfat)^2)# 352.6408
mse_caret_tree_test <- (1/length(trainset$DEXfat))*ssr_caret_tree_test# 6.297157
rmse_caret_tree_test <- sqrt(mse_caret_tree_test)# 2.509414

predict_caret_randomforest <- predict(caret_randomforest, testset)
ssr_caret_randomforest_test <- sum((predict_caret_randomforest-testset$DEXfat)^2)# 201.3924
mse_caret_randomforest_test <- (1/length(trainset$DEXfat))*ssr_caret_randomforest_test
rmse_caret_randomforest_test <- sqrt(mse_caret_randomforest_test)

# all in one table
results_all <-as.data.frame(matrix(NA, nrow=7, ncol = 8))
colnames(results_all) <- c("model", "ssr_train","mse_train","rmse_train","ssr_test","mse_test", "rmse_test", "RSQ")
results_all[,1] <- c("boost_CV", "boost_AIC", "tree", "tree_cv", "random_forest", "OLS", "OLS_tunned")
results_all[1,2] <- ssr_boost_all_train_cv
results_all[1,3] <- mse_boost_all_train_cv
results_all[1,4] <- rmse_boost_all_train_cv
results_all[1,5] <- ssr_boost_all_test_cv
results_all[1,6] <- mse_boost_all_test_cv
results_all[1,7] <- rmse_boost_all_test_cv
results_all[1,8] <- cor(boost_model_cv_train$fitted(), trainset$DEXfat)^2

results_all[2,2] <- ssr_boost_all_train_aic
results_all[2,3] <- mse_boost_all_train_aic
results_all[2,4] <- rmse_boost_all_train_aic
results_all[2,5] <- ssr_boost_all_test_aic
results_all[2,6] <- mse_boost_all_test_aic
results_all[2,7] <- rmse_boost_all_test_aic
results_all[2,8] <- cor(boost_model_aic_train$fitted(), trainset$DEXfat)^2

results_all[3,2] <- ssr_tree
results_all[3,3] <- mse_tree
results_all[3,4] <- rmse_tree
results_all[3,5] <- ssr_tree_test
results_all[3,6] <- mse_tree_test
results_all[3,7] <- rmse_tree_test
results_all[3,8] <- cor(predict(tree), trainset$DEXfat)^2

results_all[4,4] <- rmse_compare[1]
results_all[4,7] <- rmse_caret_tree_test
results_all[4,8] <- rsq_compare[1]


results_all[5,4] <- rmse_compare[2]
results_all[5,7] <- rmse_caret_randomforest_test
results_all[5,8] <- rsq_compare[2]

results_all[6,2] <- ssr_ols_train
results_all[6,3] <- mse_ols_train
results_all[6,4] <- rmse_ols_train
results_all[6,5] <- ssr_ols_all
results_all[6,6] <- mse_ols_all
results_all[6,7] <- rmse_ols_all
results_all[6,8] <- summary(ols_model_all)$adj.r.squared # adjusted

results_all[7,2] <- ssr_ols_corrected
results_all[7,3] <- mse_ols_corrected
results_all[7,4] <- rmse_ols_corrected
results_all[7,5] <- ssr_ols_corrected_test
results_all[7,6] <- mse_ols_corrected_test
results_all[7,7] <- rmse_ols_corrected_test
results_all[7,8] <- summary(ols_model_corrected)$adj.r.squared # adjusted

results_all[,2:8] <- round(results_all[,2:8],2)
print(results_all)
# ===================================================================================
# ---------------------------L2Boosting ---------------------------------------------
# ===================================================================================
# l2boosting

l2_boost<-l2boost(DEXfat ~ ., data = trainset,M=1000, nu=.2)
l2_predict<-predict(l2_boost, testset)
l2_predict_res<-l2_predict$yhat
ssr_l2<- sum((l2_predict_res-testset$DEXfat)^2)
cbind(l2_predict_res,testset$DEXfat)
par(mfrow=c(1,2))
plot(l2_boost)
plot(l2_boost, type="coef")

# -----------------------L2 loop for best nr of steps acc to cv -------------------
nu.list <- seq(from = 0.1, to= 1, by = 0.1)
m.list <- c(100)
results.l2 <-as.data.frame(matrix(NA, nrow = 10, ncol = 6))
results.l2[,1] <- nu.list
colnames(results.l2)<- c("nu_value", "M=100", "MSE_cv","MSE_test", "SSR_train","SSR_test")
for (n in 1:length(nu.list)){
  nu_value <- nu.list[n]
  for (m in 1:length(m.list)){
    m_value <- m.list[m]
    set.seed(15091987)#for reproducibility
    l2.object.cv.sim      <- cv.l2boost(trainset[,-c(2)],  trainset$DEXfat, M=m_value, nu=nu_value)
    optimal_step          <- l2.object.cv.sim$opt.step
    results.l2[n,2]     <- optimal_step
    results.l2[n,3]     <- l2.object.cv.sim$mse
    l2.object.opt         <- l2boost(trainset[,-c(2)],  trainset$DEXfat, M=optimal_step, nu=nu_value)
    l2.object.opt.pred    <- predict(l2.object.opt, testset[,-c(2)])
    l2.train.opt.pred     <- predict(l2.object.opt, trainset[,-c(2)])
    results.l2[n,4]     <- (1/length(testset$DEXfat))*(sum((l2.object.opt.pred$yhat-testset$DEXfat)^2))
    results.l2[n,5]      <- sum((l2.train.opt.pred$yhat-trainset$DEXfat)^2)
    results.l2[n,6]      <- sum((l2.object.opt.pred$yhat-testset$DEXfat)^2)
  }
}
results.l2[,3:6] <- round(results.l2[,3:6],2)
print(results.l2)
min(results.l2$MSE_cv)
best_nu <- results.l2$nu[results.l2$MSE_cv==min(results.l2$MSE_cv)]
best_m <- results.l2$`M=100`[results.l2$MSE_cv==min(results.l2$MSE_cv)]
best.l2 <- l2boost(trainset[,-c(2)],  trainset$DEXfat, M=best_m, nu=best_nu)
#png("~/Downloads/k4rtik-latex-project-report-template-4a5e948/l2_best.png")
par(mfrow=c(1,2))
plot(best.l2)
plot(best.l2, type = "coef")
#dev.off()
# ================================================================================================
# ----------------------------------- ADABOOST ---------------------------------------------------
# ================================================================================================
# example for calculating adapted from MIT tutorial available on: https://www.youtube.com/watch?v=gmok1h8wG-Q
x<- c(1,5,1,5,3)
y<- c(5,5,1,1,3)
letter<- c("A", "B", "D", "E", "C")
example <- as.data.frame(cbind(letter,x,y))

ggplot(example[1:4,], aes(x,y, color=as.factor(letter)))+
  geom_point(size=10, pch="+")+
  geom_text(data = example[1:4,],aes(label = letter), position=position_dodge(0.9))+
  geom_point(data=example[5,],aes(x,y, color=as.factor(letter)),size=10, pch="-")+
  geom_text(data = example[5,],aes(label = letter))+
  geom_vline(xintercept=c(1.5, 2.5, 3.5))+
  ggtitle("Visualization of AdaBoost example")

# ----------------------------------------- THE END ---------------------------------------------

