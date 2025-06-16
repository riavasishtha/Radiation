# estimate the causal effect of occupational radiation exposure (a continuous treatment) on survival outcomes, using data from nuclear industry workers 

library(tidyverse)
library(caret)
library(rsample)

data <- read.csv(file = "/Users/riavasishtha/Desktop/Radiation/Radiation/cleaned_data.csv", header = T)
glimpse(data)  # look at variables & types

# use classification model to predict a discrete category based on input features
# here, we're trying to find if it will be censored or not for the ipcw

data$class_uncensored <- ifelse(data$Death_event == 1, "yes", "no")  
data$class_uncensored <- as.factor(data$class_uncensored)

set.seed(123)
split <- initial_split(data, prop = 0.7, strata = class_uncensored)
train_data <- training(split)
test_data <- testing(split)

rf_grid <- expand.grid(
  .mtry = c(3, 5, 7),
  .splitrule = "gini",
  .min.node.size = c(5, 10)
)

fit_control <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

set.seed(123)
rf_ipcw_model <- train(
  class_uncensored ~ Birth_year + Total_dose_mSv + Facilities_worked + 
    Pay_code_1 + Pay_code_2 + Pay_code_3 + 
    Race_W + Race_B + Race_R + Year_of_hire + 
    Sex_numeric + Years_worked + Age_at_termination,
  data = train_data,
  method = "ranger",
  tuneGrid = rf_grid,
  trControl = fit_control,
  metric = "ROC",
  importance = 'impurity'
)

data$prob_uncensored <- predict(rf_ipcw_model, newdata = data, type = "prob")[, "yes"]
data$ipcw_weight <- 1 / data$prob_uncensored
summary(data$ipcw_weight)
hist(data$ipcw_weight,
     breaks = 30,
     main = "Histogram of IPC Weights",
     xlab = "IPCW Weight",
     col = "skyblue",
     border = "white")

ggplot(rf_ipcw_model, metric = "ROC") +
  labs(title = "ROC Performance Across Parameter Combinations")

# heat map
trellis.par.set(caretTheme())
plot(rf_ipcw_model,
     metric = "ROC",
     plotType = "level",
     scales = list(x = list(rot = 90)),
     main = "Hyperparameter Tuning (ROC)")

var_imp <- varImp(rf_ipcw_model)
plot(var_imp, top = 10, main = "Top Predictors of Uncensoring")
