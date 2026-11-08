############################################################
# STATISTICAL DATA ANALYSIS
# From Functional Connectivity to Macro-Structure: a multimodal analysis 
# of Language Lateralisation in the Human Brain
# Dataset: BIL&GIN: A neuroimaging, cognitive, behavioral, and genetic
# database for the study of human brain lateralization
############################################################

# ==========================================================
# 1. DATA IMPORT AND INITIAL INSPECTION
# ==========================================================

library(tidyverse)
library(ggplot2)

#install.packages("corrplot")
library(corrplot)

#install.packages("glmnet")
library(glmnet)

#install.packages("car")
library(car)

library(dplyr)

library(lmtest)

library(emmeans)

#install.packages("patchwork")
library(patchwork) 

getwd()

# Import the dataset 
data <- read_delim(
  "dataset/287_Labache_2020_eLife.txt",
  delim = "\t",
  locale = locale(decimal_mark = ",")
)

# Initial inspection of dataset structure and data quality
dim(data)
glimpse(data)
head(data)
tail(data)
summary(data)
colSums(is.na(data))

# ==========================================================
# 2. DATA PREPARATION
# ==========================================================

# Select variables relevant to the research question
data_project <- data %>%
  select(
    ID,
    Age,
    EducLevel,
    Sex,
    MP,
    ES,
    Rs_mIHHC,#resting-state interhemispheric connectivity
    Rs_DC_Asymmetry, #degree centrality
    sentCORE_PROD_Asym,
    sentCORE_LISN_Asym,
    sentCORE_READ_Asym,
    GMasym,
    WMasym,
    CCvol,
    ITV,
    Classif_MMHC,
    DISSOCIATIONS_2categories
  )

# Convert categorical variables to factors
data_project <- data_project %>%
  mutate(
    Sex = factor(Sex),
    MP = factor(MP),
    Classif_MMHC = factor(Classif_MMHC),
    DISSOCIATIONS_2categories = factor(DISSOCIATIONS_2categories)
  )

# Inspect the analysis dataset structure
glimpse(data_project)
#check missing values
colSums(is.na(data_project))

# ID=218 valori atipici

# ==========================================================
# 3. DESCRIPTIVE ANALYSIS
# ==========================================================

#summary for continuous variables
data_project %>%
  select(
    Age, EducLevel, ES, Rs_mIHHC, Rs_DC_Asymmetry,
    sentCORE_PROD_Asym, sentCORE_LISN_Asym, sentCORE_READ_Asym,
    GMasym, WMasym, CCvol, ITV
  ) %>%
  summary()

#summary for categorical variables
table(data_project$Sex); prop.table(table(data_project$Sex))
table(data_project$MP); prop.table(table(data_project$MP))
table(data_project$Classif_MMHC); prop.table(table(data_project$Classif_MMHC))
table(data_project$DISSOCIATIONS_2categories); prop.table(table(data_project$DISSOCIATIONS_2categories))

#generate histograms and boxplots side by side
continuous_vars <- c(
  "Age", "EducLevel", "ES", "Rs_mIHHC", "Rs_DC_Asymmetry",
  "sentCORE_PROD_Asym", "sentCORE_LISN_Asym", "sentCORE_READ_Asym",
  "GMasym", "WMasym", "CCvol", "ITV"
)

colors_list <- c('violet', 'lightblue', '#339999', '#9999ff', 'cyan', 
                 'yellow', 'orange', 'red', 'pink', 'lightgreen', 'gray', 'brown')

for (i in seq_along(continuous_vars)) {
  var_name <- continuous_vars[i]
  var_data <- data_project[[var_name]]
  col_choice <- colors_list[i]
  
  # automatic saving
  # png(filename = paste0("plot_side_", var_name, ".png"), width = 900, height = 450)
  
  # 1 row 2 cols setting
  par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
  
  # Histogram on the left
  hist(var_data, prob = TRUE, col = col_choice, 
       xlab = var_name, main = paste(var_name, "- Histogram"), breaks = 15)
  
  # Vertical boxplot on the right
  boxplot(var_data, col = col_choice, 
          ylab = var_name, main = paste(var_name, "- Boxplot"))
  
  # dev.off() 
}

#bring standard R setting back
par(mfrow = c(1, 1))

# ==========================================================
# 4. EXPLORATORY ANALYSIS
# ==========================================================

vars_to_correlate <- data_project[, c("sentCORE_PROD_Asym", "sentCORE_LISN_Asym",
                                      "sentCORE_READ_Asym", "EducLevel", "Rs_mIHHC", 
                                      "Rs_DC_Asymmetry", "ES", "Age", "GMasym", 
                                      "WMasym", "CCvol", "ITV")]

#calculate correlation matrix (screen for multicollinearity)
cor_matrix <- cor(vars_to_correlate, use = "complete.obs")

corrplot(cor_matrix, 
         method = "color", 
         type = "upper", 
         col = colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", 
                                  "#F7F7F7", "#92C5DE", "#4393C3", "#2166AC", "#084081"))(200),
         addCoef.col = "black",   # add numbers inside cells
         number.cex = 0.6,        # adjust numbers dimension
         tl.col = "black",        # color of variables' labels
         tl.srt = 45,             # labels' slope
         title = "Correlation matrix",
         tl.cex=0.9,
         mar = c(0,0,1,0))

tasks <- c("sentCORE_PROD_Asym", "sentCORE_LISN_Asym", "sentCORE_READ_Asym")

predictors <- c("Rs_mIHHC", "ES", "Rs_DC_Asymmetry", "GMasym", "WMasym", "CCvol")

for (predictor in predictors) {
  for (task in tasks) {
    
    p <- ggplot(
      data_project,
      aes(
        x = .data[[predictor]],
        y = .data[[task]]
      )
    ) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "blue") +
      geom_smooth(method = "loess", formula = y ~ x, se = FALSE,color = "red",
        linetype = "dashed") +
      labs(
        title = paste(task, "vs", predictor),
        x = predictor,
        y = task
      ) +
      theme_classic()
    
    print(p)
  }
}

# Part 1: Resting-State Connectivity (Rs_mIHHC) vs Language Lateralisation

for (task in tasks) {
  cat("\n--- Rs_mIHHC vs", task, "---\n")
  print(cor.test(data_project$Rs_mIHHC, data_project[[task]], method = "pearson"))
  print(cor.test(data_project$Rs_mIHHC, data_project[[task]], method = "spearman"))
}

# Part 2: Degree Centrality (Rs_DC_Asymmetry) vs Language Lateralisation

for (task in tasks) {
  cat("\n--- Rs_DC_Asymmetry vs", task, "---\n")
  print(cor.test(data_project$Rs_DC_Asymmetry, data_project[[task]], method = "pearson"))
  print(cor.test(data_project$Rs_DC_Asymmetry, data_project[[task]], method = "spearman"))
}

# Part 3: Manual Preference (ES) vs Language Lateralisation

for (task in tasks) {
  cat("\n--- Edinburgh Score (ES) vs", task, "---\n")
  print(cor.test(data_project$ES, data_project[[task]], method = "pearson"))
  print(cor.test(data_project$ES, data_project[[task]], method = "spearman"))
}

# Part 4: Structural Organization (GMasym, WMasym, CCvol) vs Language Lateralisation

structural_vars <- c("GMasym", "WMasym", "CCvol")

for (struct in structural_vars) {
  for (task in tasks) {
    cat("\n---", struct, "vs", task, "---\n")
    print(cor.test(data_project[[struct]], data_project[[task]], method = "pearson"))
    print(cor.test(data_project[[struct]], data_project[[task]], method = "spearman"))
  }
}

# ===============================================================================
# 5. STATISTICAL ANALYSIS 1 - Functional Connectivity and Language Lateralisation
# ===============================================================================

#MULTIPLE LINEAR REGRESSION MODELS

# Predictors and outcomes
predictors_A1 <- c("Rs_mIHHC", "Age", "Sex", "ES", "EducLevel")

outcomes <- c(
  PROD = "sentCORE_PROD_Asym",
  LISN = "sentCORE_LISN_Asym",
  READ = "sentCORE_READ_Asym"
)

analyze_lm_A1 <- function(data, task, outcome, predictors_A1) {
  
  # Prepare complete cases and preserve original row numbers
  d <- data
  d$.row <- seq_len(nrow(d))
  vars <- c(outcome, predictors_A1)
  d <- d[complete.cases(d[, vars, drop = FALSE]), ]
  
  # Original model
  f <- reformulate(predictors_A1, response = outcome)
  model <- lm(f, data = d)
  
  cat("\n\n==========", task, "- ORIGINAL MODEL ==========\n")
  print(summary(model))
  
  # Diagnostic plots
  par(mfrow = c(2, 2))
  plot(model, which = 1:4)
  par(mfrow = c(1, 1))
  
  # Cook's distance, correctly matched to IDs
  used_rows <- as.integer(rownames(model.frame(model)))
  used_data <- d[match(used_rows, d$.row), ]
  
  cook <- data.frame(
    ID = used_data$ID,
    row = used_data$.row,
    Cook = as.numeric(cooks.distance(model))
  )
  
  cook <- cook[order(cook$Cook, decreasing = TRUE), ]
  
  cat("\n--- Top 10 Cook's distances ---\n")
  print(head(cook, 10))
  
  # Select the three highest Cook's distances
  top3 <- head(cook, 3)
  
  cat("\n--- Participants selected for sensitivity analysis ---\n")
  print(data[top3$row, ])
  
  # Sensitivity model excluding those three observations
  d_sens <- d[!d$.row %in% top3$row, ]
  model_sens <- lm(f, data = d_sens)
  
  cat("\n--- SENSITIVITY MODEL ---\n")
  print(summary(model_sens))
  
  # Save results
  list(
    model = model,
    sensitivity = model_sens,
    cook = cook,
    top3 = top3
  )
}

# Run all three tasks
results_lm_A1 <- lapply(names(outcomes), function(task) {
  analyze_lm(
    data = data_project,
    task = task,
    outcome = outcomes[[task]],
    predictors_A1 = predictors_A1
  )
})

names(results_lm_A1) <- names(outcomes)

#Extract p-values only for Rs_mIHHC
p_mIHHC_A1 <- c(
  PROD = coef(summary(results_lm_A1$PROD$model))["Rs_mIHHC", "Pr(>|t|)"],
  LISN = coef(summary(results_lm_A1$LISN$model))["Rs_mIHHC", "Pr(>|t|)"],
  READ = coef(summary(results_lm_A1$READ$model))["Rs_mIHHC", "Pr(>|t|)"]
)

#multiple test correction on Rs_mIHHC
results_p_mIHHC_A1 <- data.frame(
  Outcome = names(p_mIHHC_A1),
  p_raw = unname(p_mIHHC_A1),
  p_Holm = unname(p.adjust(p_mIHHC_A1, method = "holm"))
)

cat("\n--- MULTIPLE-TESTING CORRECTION FOR Rs_mIHHC ONLY ---\n")
print(results_p_mIHHC_A1)

#Variance Inflation factor (multicollinearity)
cat("\n--- VIF ANALYSIS 1 ---\n")
print(vif(results_lm_A1$PROD$model))
print(vif(results_lm_A1$LISN$model))
print(vif(results_lm_A1$READ$model))

#PENALIZED REGRESSION MODELS: RIDGE AND LASSO

fit_ridge_lasso <- function(data, task, outcome, predictors_A1) {
  
  # Keep complete cases for this task
  vars <- c(outcome, predictors_A1)
  
  d <- data[
    complete.cases(data[, vars, drop = FALSE]),
    ,
    drop = FALSE
  ]
  
  # Design matrix: automatically encodes categorical predictors
  X <- model.matrix(
    reformulate(predictors_A1),
    data = d
  )[, -1, drop = FALSE]
  
  y <- d[[outcome]]
  
  # Same cross-validation folds for Ridge and Lasso
  set.seed(123)
  foldid <- sample(rep(1:10, length.out = nrow(X)))
  
  # Ridge: alpha = 0
  ridge <- cv.glmnet(
    x = X,
    y = y,
    alpha = 0,
    family = "gaussian",
    foldid = foldid,
    type.measure = "mse",
    standardize = TRUE
  )
  
  # Lasso: alpha = 1
  lasso <- cv.glmnet(
    x = X,
    y = y,
    alpha = 1,
    family = "gaussian",
    foldid = foldid,
    type.measure = "mse",
    standardize = TRUE
  )
  
  cat("\n==========", task, "==========\n")
  cat("N =", nrow(d), "\n")
  
  cat("\n--- RIDGE ---\n")
  cat("Lambda min:", ridge$lambda.min, "\n")
  cat("CV MSE:", min(ridge$cvm), "\n")
  
  cat("\n--- LASSO ---\n")
  cat("Lambda min:", lasso$lambda.min, "\n")
  cat("CV MSE:", min(lasso$cvm), "\n")
  
  return(list(
    ridge = ridge,
    lasso = lasso,
    n = nrow(d),
    X = X,
    y = y
  ))
}

results_penalized_A1 <- lapply(names(outcomes), function(task) {
  
  fit_ridge_lasso(
    data = data_project,
    task = task,
    outcome = outcomes[[task]],
    predictors_A1 = predictors_A1
  )
  
})

names(results_penalized_A1) <- names(outcomes)

for (task in names(results_penalized_A1)) {
  
  cat("\n==========", task, "==========\n")
  
  for (method in c("ridge", "lasso")) {
    
    model <- results_penalized_A1[[task]][[method]]
    
    cat("\n", toupper(method), "\n")
    cat("Lambda min:", model$lambda.min, "\n")
    cat("MSE (lambda min):",
        min(model$cvm), "\n")
    
    cat("Lambda 1se:", model$lambda.1se, "\n")
    
    cat("Coefficients:\n")
    print(coef(model, s = "lambda.min"))
  }
}

# ============================================================
# COMPARE THREE MODELS: LINEAR REGRESSION, RIDGE AND LASSO
# Repeated 10-fold cross-validation (5 repetitions)
# ============================================================

set.seed(123)

# Function to compare the three models
compare_predictive_A1 <- function(data, outcome, predictors_A1,
                               K = 10, repetition = 1) {
  
  # Keep complete cases for the outcome and predictors
  vars <- c(outcome, predictors_A1)
  d <- data[complete.cases(data[, vars, drop = FALSE]), ,
            drop = FALSE]
  
  # Ensure Sex is categorical
  d$Sex <- factor(d$Sex)
  
  n <- nrow(d)
  y <- d[[outcome]]
  
  # Generate outer folds
  set.seed(1000 + repetition)
  folds <- sample(rep(1:K, length.out = n))
  
  # Store out-of-fold predictions
  pred <- matrix(
    NA_real_,
    nrow = n,
    ncol = 3,
    dimnames = list(NULL, c("Linear", "Ridge", "Lasso"))
  )
  
  f <- reformulate(predictors_A1, response = outcome)
  x_formula <- reformulate(predictors_A1)
  
  # Outer cross-validation
  for (k in seq_len(K)) {
    
    train <- d[folds != k, , drop = FALSE]
    test  <- d[folds == k, , drop = FALSE]
    
    # 1. Multiple linear regression
    fit_lm <- lm(f, data = train)
    
    pred[folds == k, "Linear"] <-
      predict(fit_lm, newdata = test)
    
    # Design matrices for Ridge and Lasso
    X_train <- model.matrix(
      x_formula, data = train
    )[, -1, drop = FALSE]
    
    X_test <- model.matrix(
      x_formula, data = test
    )[, -1, drop = FALSE]
    
    y_train <- train[[outcome]]
    
    # Inner folds for selecting lambda
    set.seed(2000 + repetition * 100 + k)
    
    inner_folds <- sample(
      rep(1:10, length.out = nrow(train))
    )
    
    # 2. Ridge regression
    fit_ridge <- cv.glmnet(
      X_train,
      y_train,
      alpha = 0,
      family = "gaussian",
      foldid = inner_folds,
      type.measure = "mse"
    )
    
    pred[folds == k, "Ridge"] <- as.numeric(
      predict(
        fit_ridge,
        newx = X_test,
        s = "lambda.min"
      )
    )
    
    # 3. Lasso regression
    fit_lasso <- cv.glmnet(
      X_train,
      y_train,
      alpha = 1,
      family = "gaussian",
      foldid = inner_folds,
      type.measure = "mse"
    )
    
    pred[folds == k, "Lasso"] <- as.numeric(
      predict(
        fit_lasso,
        newx = X_test,
        s = "lambda.min"
      )
    )
  }
  
  # Calculate predictive performance
  metrics <- do.call(
    rbind,
    lapply(colnames(pred), function(method) {
      
      errors <- y - pred[, method]
      
      data.frame(
        Outcome = outcome,
        Model = method,
        RMSE = sqrt(mean(errors^2)),
        MAE = mean(abs(errors)),
        R2_predictive =
          1 - sum(errors^2) / sum((y - mean(y))^2)
      )
    })
  )
  
  metrics$Repetition <- repetition
  
  return(metrics)
}

# Run the comparison for all outcomes, 5 repetitions
results_repeated_A1 <- list()

for (r in 1:5) {
  
  results_repeated_A1[[r]] <- do.call(
    rbind,
    lapply(names(outcomes), function(task) {
      
      compare_predictive_A1(
        data = data_project,
        outcome = outcomes[[task]],
        predictors_A1 = predictors_A1,
        K = 10,
        repetition = r
      )
    })
  )
}

# Combine results from all repetitions
results_repeated_A1 <- do.call(rbind, results_repeated_A1)

# View the results of each repetition
#print(results_repeated_A1)

# Summarise mean and standard deviation across repetitions
summary_repeated_A1 <- results_repeated_A1 %>%
  group_by(Outcome, Model) %>%
  summarise(
    RMSE_mean = mean(RMSE),
    RMSE_SD = sd(RMSE),
    MAE_mean = mean(MAE),
    MAE_SD = sd(MAE),
    R2_mean = mean(R2_predictive),
    R2_SD = sd(R2_predictive),
    .groups = "drop"
  )

cat("\n--- PREDICTIVE PERFORMANCE SUMMARY ANALYSIS 1---\n")
print(summary_repeated_A1)

----------------------------------------------------------------------------------------------

#Rs_mIHHC values to calculate predictions
newdata <- data.frame(
  Rs_mIHHC = seq(
    min(data_project$Rs_mIHHC, na.rm = TRUE),
    max(data_project$Rs_mIHHC, na.rm = TRUE),
    length.out = 100
  ),
  Age = mean(data_project$Age, na.rm = TRUE),
  Sex = factor("F", levels = levels(data_project$Sex)),
  ES = mean(data_project$ES, na.rm = TRUE),
  EducLevel=mean(data_project$EducLevel, na.rm = TRUE)
)

#Predictions for each model
pred_prod <- predict(results_lm_A1$PROD$model, newdata, interval = "confidence")
pred_lisn <- predict(results_lm_A1$LISN$model, newdata, interval = "confidence")
pred_read <- predict(results_lm_A1$READ$model, newdata, interval = "confidence")

#Organizing results
plot_data <- rbind(
  data.frame(newdata, as.data.frame(pred_prod), Task = "Production"),
  data.frame(newdata, as.data.frame(pred_lisn), Task = "Listening"),
  data.frame(newdata, as.data.frame(pred_read), Task = "Reading")
)

ggplot(plot_data, aes(x = Rs_mIHHC, y = fit, color = Task, fill = Task)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c(
    "Production" = "#3399FF",
    "Listening" = "#00CC66",
    "Reading" = "#9933FF"
  )) +
  scale_fill_manual(values = c(
    "Production" = "#3399FF",
    "Listening" = "#00CC66",
    "Reading" = "#9933FF"
  )) +
  labs(
    title = "Functional interhemispheric connectivity and language lateralisation",
    subtitle = "Adjusted predictions for Age, Sex, ES and Educational level (95% confidence intervals)",
    x = "Functional connectivity (Rs_mIHHC)",
    y = "Functional language asymmetry",
    color = "Task",
    fill = "Task"
  ) +
  theme_minimal()+
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 12)
  )

 # ggsave(
 #   "analysis1_regression.png",
 #  width = 10,
 #  height = 7,
 #  dpi = 300
 # )

#MODEL COMPARISON GRAPH

# Design matrix per Ridge e Lasso
X_new <- model.matrix(
  reformulate(predictors_A1),
  data = newdata
)[, -1, drop = FALSE]

# Predictions of three methods for each outcome
compare_models <- function(task) {
  
  # Linear
  pred_lm <- predict(
    results_lm_A1[[task]]$model,
    newdata = newdata
  )
  
  # Ridge
  pred_ridge <- as.numeric(predict(
    results_penalized_A1[[task]]$ridge,
    newx = X_new,
    s = "lambda.min"
  ))
  
  # Lasso
  pred_lasso <- as.numeric(predict(
    results_penalized_A1[[task]]$lasso,
    newx = X_new,
    s = "lambda.min"
  ))
  
  data.frame(
    Rs_mIHHC = rep(newdata$Rs_mIHHC, 3),
    Prediction = c(pred_lm, pred_ridge, pred_lasso),
    Model = factor(
      rep(c("Linear regression", "Ridge", "Lasso"),
          each = nrow(newdata)),
      levels = c("Linear regression", "Ridge", "Lasso")
    ),
    Task = task
  )
}

# Unire i risultati dei tre outcome
comparison_data <- do.call(
  rbind,
  lapply(names(outcomes), compare_models)
)

# Etichette più leggibili
comparison_data$Task <- factor(
  comparison_data$Task,
  levels = c("PROD", "LISN", "READ"),
  labels = c("Production", "Listening", "Reading")
)

# Grafico comparativo
ggplot(
  comparison_data,
  aes(x = Rs_mIHHC, y = Prediction, color = Model)
) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Task, scales = "free_y") +
  scale_color_manual(values = c(
    "Linear regression" = "#00ccff",
    "Ridge" = "#66ff00",
    "Lasso" = "#ff3300"
  )) +
  labs(
    title = "Comparison of regression models",
    subtitle = "Predictions adjusted for Age, Sex, ES and Educational level",
    x = "Functional connectivity (Rs_mIHHC)",
    y = "Predicted language asymmetry",
    color = "Model"
  ) +
  theme_minimal()+
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 12)
  )

 ggsave(
   "analysis1_comparemodels.png",
  width = 10,
  height = 7,
   dpi = 300
 )

# =========================================================================
# 6. STATISTICAL ANALYSIS 2 - Degree centrality and language lateralisation
# =========================================================================

#Does Rs_DC_Asymmetry explain additional variance in language functional 
#beyond Rs_mIHHC lateralization?

#MULTIPLE LINEAR REGRESSION MODELS

predictors_A2 <- c("Rs_DC_Asymmetry", "Rs_mIHHC", "Age", "Sex", "ES", "EducLevel")

outcomes <- c(
  PROD = "sentCORE_PROD_Asym",
  LISN = "sentCORE_LISN_Asym",
  READ = "sentCORE_READ_Asym"
)

# Function for regression, diagnostics and sensitivity analysis
analyze_lm_A2 <- function(data, task, outcome, predictors_A2) {
  
  # Prepare complete cases and preserve original row numbers
  d <- data
  d$.row <- seq_len(nrow(d))
  vars <- c(outcome, predictors_A2)
  d <- d[complete.cases(d[, vars, drop = FALSE]), ]
  
  # Original model
  f <- reformulate(predictors_A2, response = outcome)
  model <- lm(f, data = d)
  
  cat("\n\n==========", task, "- ORIGINAL MODEL ==========\n")
  print(summary(model))
  
  # Diagnostic plots
  par(mfrow = c(2, 2))
  plot(model, which = 1:4)
  par(mfrow = c(1, 1))
  
  # Cook's distance, correctly matched to IDs
  used_rows <- as.integer(rownames(model.frame(model)))
  used_data <- d[match(used_rows, d$.row), ]
  
  cook <- data.frame(
    ID = used_data$ID,
    row = used_data$.row,
    Cook = as.numeric(cooks.distance(model))
  )
  
  cook <- cook[order(cook$Cook, decreasing = TRUE), ]
  
  cat("\n--- Top 10 Cook's distances ---\n")
  print(head(cook, 10))
  
  # Select the three highest Cook's distances
  top3 <- head(cook, 3)
  
  cat("\n--- Participants selected for sensitivity analysis ---\n")
  print(data[top3$row, ])
  
  # Sensitivity model excluding those three observations
  d_sens <- d[!d$.row %in% top3$row, ]
  model_sens <- lm(f, data = d_sens)
  
  cat("\n--- SENSITIVITY MODEL ---\n")
  print(summary(model_sens))
  
  list(
    model = model,
    sensitivity = model_sens,
    cook = cook,
    top3 = top3
  )
}

# Run all three tasks
results_lm_A2 <- lapply(names(outcomes), function(task) {
  analyze_lm_A2(
    data = data_project,
    task = task,
    outcome = outcomes[[task]],
    predictors_A2 = predictors_A2
  )
})

names(results_lm_A2) <- names(outcomes)

# Extract p-values for Rs_DC_Asymmetry only
p_DC_A2 <- c(
  PROD = coef(summary(results_lm_A2$PROD$model))[
    "Rs_DC_Asymmetry", "Pr(>|t|)"
  ],
  LISN = coef(summary(results_lm_A2$LISN$model))[
    "Rs_DC_Asymmetry", "Pr(>|t|)"
  ],
  READ = coef(summary(results_lm_A2$READ$model))[
    "Rs_DC_Asymmetry", "Pr(>|t|)"
  ]
)

# Multiple-testing correction across the three outcomes
results_p_A2 <- data.frame(
  Outcome = names(p_DC_A2),
  p_raw = unname(p_DC_A2),
  p_Holm = unname(p.adjust(p_DC_A2, method = "holm"))
)

print(results_p_A2)

print(vif(results_lm_A2$PROD$model))
print(vif(results_lm_A2$LISN$model))
print(vif(results_lm_A2$READ$model))

# COMPARE A1 and A2 for all three outcomes

compare_A1_A2 <- function(model_A1, model_A2) {
  
  list(
    ANOVA = anova(model_A1, model_A2),
    
    AIC_BIC = data.frame(
      Model = c("A1_base", "A2_extended"),
      AIC = c(AIC(model_A1), AIC(model_A2)),
      BIC = c(BIC(model_A1), BIC(model_A2))
    )
  )
}

comparison_A1_A2 <- list(
  PROD = compare_A1_A2(
    results_lm_A1$PROD$model,
    results_lm_A2$PROD$model
  ),
  LISN = compare_A1_A2(
    results_lm_A1$LISN$model,
    results_lm_A2$LISN$model
  ),
  READ = compare_A1_A2(
    results_lm_A1$READ$model,
    results_lm_A2$READ$model
  )
)

# Display results
comparison_A1_A2$PROD$ANOVA
comparison_A1_A2$PROD$AIC_BIC

comparison_A1_A2$LISN$ANOVA
comparison_A1_A2$LISN$AIC_BIC

comparison_A1_A2$READ$ANOVA
comparison_A1_A2$READ$AIC_BIC

#PENALIZED REGRESSION MODELS: RIDGE AND LASSO

fit_ridge_lasso <- function(data, task, outcome, predictors_A2) {
  
  # Keep complete cases for this task
  vars <- c(outcome, predictors_A2)
  
  d <- data[
    complete.cases(data[, vars, drop = FALSE]),
    ,
    drop = FALSE
  ]
  
  # Design matrix: automatically encodes categorical predictors
  X <- model.matrix(
    reformulate(predictors_A2),
    data = d
  )[, -1, drop = FALSE]
  
  y <- d[[outcome]]
  
  # Same cross-validation folds for Ridge and Lasso
  set.seed(123)
  foldid <- sample(rep(1:10, length.out = nrow(X)))
  
  # Ridge: alpha = 0
  ridge <- cv.glmnet(
    x = X,
    y = y,
    alpha = 0,
    family = "gaussian",
    foldid = foldid,
    type.measure = "mse",
    standardize = TRUE
  )
  
  # Lasso: alpha = 1
  lasso <- cv.glmnet(
    x = X,
    y = y,
    alpha = 1,
    family = "gaussian",
    foldid = foldid,
    type.measure = "mse",
    standardize = TRUE
  )
  
  cat("\n==========", task, "==========\n")
  cat("N =", nrow(d), "\n")
  
  cat("\n--- RIDGE ---\n")
  cat("Lambda min:", ridge$lambda.min, "\n")
  cat("CV MSE:", min(ridge$cvm), "\n")
  
  cat("\n--- LASSO ---\n")
  cat("Lambda min:", lasso$lambda.min, "\n")
  cat("CV MSE:", min(lasso$cvm), "\n")
  
  return(list(
    ridge = ridge,
    lasso = lasso,
    n = nrow(d),
    X = X,
    y = y
  ))
}

results_penalized_A2 <- lapply(names(outcomes), function(task) {
  
  fit_ridge_lasso(
    data = data_project,
    task = task,
    outcome = outcomes[[task]],
    predictors_A2 = predictors_A2
  )
  
})

names(results_penalized_A2) <- names(outcomes)

for (task in names(results_penalized_A2)) {
  
  cat("\n==========", task, "==========\n")
  
  for (method in c("ridge", "lasso")) {
    
    model <- results_penalized_A2[[task]][[method]]
    
    cat("\n", toupper(method), "\n")
    cat("Lambda min:", model$lambda.min, "\n")
    cat("MSE (lambda min):",
        min(model$cvm), "\n")
    
    cat("Lambda 1se:", model$lambda.1se, "\n")
    
    cat("Coefficients:\n")
    print(coef(model, s = "lambda.min"))
  }
}

# ============================================================
# COMPARE THREE MODELS: LINEAR REGRESSION, RIDGE AND LASSO
# Repeated 10-fold cross-validation (5 repetitions)
# ============================================================

set.seed(123)

compare_predictive_A2 <- function(data, outcome, predictors_A2,
                               K = 10, repetition = 1) {
  
  # Keep complete cases for the outcome and predictors
  vars <- c(outcome, predictors_A2)
  d <- data[complete.cases(data[, vars, drop = FALSE]), ,
            drop = FALSE]
  
  # Ensure Sex is categorical
  d$Sex <- factor(d$Sex)
  
  n <- nrow(d)
  y <- d[[outcome]]
  
  # Generate outer folds
  set.seed(1000 + repetition)
  folds <- sample(rep(1:K, length.out = n))
  
  # Store out-of-fold predictions
  pred <- matrix(
    NA_real_,
    nrow = n,
    ncol = 3,
    dimnames = list(NULL, c("Linear", "Ridge", "Lasso"))
  )
  
  f <- reformulate(predictors_A2, response = outcome)
  x_formula <- reformulate(predictors_A2)
  
  # Outer cross-validation
  for (k in seq_len(K)) {
    
    train <- d[folds != k, , drop = FALSE]
    test  <- d[folds == k, , drop = FALSE]
    
    # 1. Multiple linear regression
    fit_lm <- lm(f, data = train)
    
    pred[folds == k, "Linear"] <-
      predict(fit_lm, newdata = test)
    
    # Design matrices for Ridge and Lasso
    X_train <- model.matrix(
      x_formula, data = train
    )[, -1, drop = FALSE]
    
    X_test <- model.matrix(
      x_formula, data = test
    )[, -1, drop = FALSE]
    
    y_train <- train[[outcome]]
    
    # Inner folds for selecting lambda
    set.seed(2000 + repetition * 100 + k)
    
    inner_folds <- sample(
      rep(1:10, length.out = nrow(train))
    )
    
    # 2. Ridge regression
    fit_ridge <- cv.glmnet(
      X_train,
      y_train,
      alpha = 0,
      family = "gaussian",
      foldid = inner_folds,
      type.measure = "mse"
    )
    
    pred[folds == k, "Ridge"] <- as.numeric(
      predict(
        fit_ridge,
        newx = X_test,
        s = "lambda.min"
      )
    )
    
    # 3. Lasso regression
    fit_lasso <- cv.glmnet(
      X_train,
      y_train,
      alpha = 1,
      family = "gaussian",
      foldid = inner_folds,
      type.measure = "mse"
    )
    
    pred[folds == k, "Lasso"] <- as.numeric(
      predict(
        fit_lasso,
        newx = X_test,
        s = "lambda.min"
      )
    )
  }
  
  # Calculate predictive performance
  metrics <- do.call(
    rbind,
    lapply(colnames(pred), function(method) {
      
      errors <- y - pred[, method]
      
      data.frame(
        Outcome = outcome,
        Model = method,
        RMSE = sqrt(mean(errors^2)),
        MAE = mean(abs(errors)),
        R2_predictive =
          1 - sum(errors^2) / sum((y - mean(y))^2)
      )
    })
  )
  
  metrics$Repetition <- repetition
  
  return(metrics)
}

# Run the comparison for all outcomes, 5 repetitions
results_repeated_A2 <- list()

for (r in 1:5) {
  
  results_repeated_A2[[r]] <- do.call(
    rbind,
    lapply(names(outcomes), function(task) {
      
      compare_predictive_A2(
        data = data_project,
        outcome = outcomes[[task]],
        predictors_A2 = predictors_A2,
        K = 10,
        repetition = r
      )
    })
  )
}

# Combine results from all repetitions
results_repeated_A2 <- do.call(rbind, results_repeated_A2)

# View the results of each repetition
#print(results_repeated_A2)

# Summarise mean and standard deviation across repetitions
summary_repeated_A2 <- results_repeated_A2 %>%
  group_by(Outcome, Model) %>%
  summarise(
    RMSE_mean = mean(RMSE),
    RMSE_SD = sd(RMSE),
    MAE_mean = mean(MAE),
    MAE_SD = sd(MAE),
    R2_mean = mean(R2_predictive),
    R2_SD = sd(R2_predictive),
    .groups = "drop"
  )

cat("\n--- PREDICTIVE PERFORMANCE SUMMARY ANALYSIS 2---\n")
print(summary_repeated_A2)

-------------------------------------------------------------------------------------
  
#update newdata including the new variable Rs_DC_Asymmetry
newdata <- data.frame(
  Rs_mIHHC = seq(
  min(data_project$Rs_mIHHC, na.rm = TRUE),
  max(data_project$Rs_mIHHC, na.rm = TRUE),
  length.out = 100
),
  Rs_DC_Asymmetry = mean(data_project$Rs_DC_Asymmetry, na.rm = TRUE), 
  Age = mean(data_project$Age, na.rm = TRUE),
  Sex = factor("F", levels = levels(data_project$Sex)),
  ES = mean(data_project$ES, na.rm = TRUE),
  EducLevel= mean(data_project$EducLevel, na.rm = TRUE)
)

#new predictions
pred_prod_dc <- predict(results_lm_A2$PROD$model, newdata, interval = "confidence")
pred_lisn_dc <- predict(results_lm_A2$LISN$model, newdata, interval = "confidence")
pred_read_dc <- predict(results_lm_A2$READ$model, newdata, interval = "confidence")

plot_data_dc <- rbind(
  data.frame(newdata, as.data.frame(pred_prod_dc), Task = "Production"),
  data.frame(newdata, as.data.frame(pred_lisn_dc), Task = "Listening"),
  data.frame(newdata, as.data.frame(pred_read_dc), Task = "Reading")
)

ggplot(plot_data_dc, aes(x = Rs_mIHHC, y = fit, color = Task, fill = Task)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c(
    "Production" = "#3399FF",
    "Listening" = "#00CC66",
    "Reading" = "#9933FF"
  )) +
  scale_fill_manual(values = c(
    "Production" = "#3399FF",
    "Listening" = "#00CC66",
    "Reading" = "#9933FF"
  )) +
  labs(
    title = "Functional connectivity and language lateralisation (Extended Model)",
    subtitle = "Adjusted predictions for Rs_DC_Asymmetry, Age, Sex, ES and Educational level (95% confidence intervals)",
    x = "Functional connectivity (Rs_mIHHC)",
    y = "Functional language asymmetry",
    color = "Task",
    fill = "Task"
  ) +
  theme_minimal()+
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 12)
  )

 # ggsave(
 #   "analysis2_regression.png",
 #   width = 10,
 #   height = 7,
 #   dpi = 300
 #  )

# Design matrix per Ridge e Lasso
X_new <- model.matrix(
  reformulate(predictors_A2),
  data = newdata
)[, -1, drop = FALSE]

# Predizioni dei tre metodi per ciascun outcome
compare_models <- function(task) {
  
  # Regressione lineare multipla
  pred_lm <- predict(
    results_lm_A2[[task]]$model,
    newdata = newdata
  )
  
  # Ridge
  pred_ridge <- as.numeric(predict(
    results_penalized_A2[[task]]$ridge,
    newx = X_new,
    s = "lambda.min"
  ))
  
  # Lasso
  pred_lasso <- as.numeric(predict(
    results_penalized_A2[[task]]$lasso,
    newx = X_new,
    s = "lambda.min"
  ))
  
  data.frame(
    Rs_mIHHC = rep(newdata$Rs_mIHHC, 3),
    Prediction = c(pred_lm, pred_ridge, pred_lasso),
    Model = factor(
      rep(c("Linear regression", "Ridge", "Lasso"),
          each = nrow(newdata)),
      levels = c("Linear regression", "Ridge", "Lasso")
    ),
    Task = task
  )
}

# Unire i risultati dei tre outcome
comparison_data <- do.call(
  rbind,
  lapply(names(outcomes), compare_models)
)

# Etichette più leggibili
comparison_data$Task <- factor(
  comparison_data$Task,
  levels = c("PROD", "LISN", "READ"),
  labels = c("Production", "Listening", "Reading")
)

ggplot(
  comparison_data,
  aes(x = Rs_mIHHC, y = Prediction, color = Model)
) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Task, scales = "free_y") +
  scale_color_manual(values = c(
    "Linear regression" = "#00ccff",
    "Ridge" = "#66ff00",
    "Lasso" = "#ff3300"
  )) +
  labs(
    title = "Comparison of regression models",
    subtitle = "Predictions adjusted for Rs_DC_Asymmetry, Age, Sex, ES and Educational level",
    x = "Functional connectivity (Rs_mIHHC)",
    y = "Predicted language asymmetry",
    color = "Model"
  ) +
  theme_minimal()+
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 12)
  )

  # 
  # ggsave(
  #   "analysis2_comparemodels.png",
  #  width = 10,
  #  height = 7,
  #   dpi = 300
  # )

# =========================================================================
# 7. STATISTICAL ANALYSIS 3 - Manual preference and language lateralisation
# =========================================================================

#Is manual preference associated to functional language asymmetry?

# MULTIPLE LINEAR REGRESSION MODELS

predictors_A3_ES <- c("ES", "Age", "Sex")
predictors_A3_MP <- c("MP", "Age", "Sex")

outcomes <- c(
  PROD = "sentCORE_PROD_Asym",
  LISN = "sentCORE_LISN_Asym",
  READ = "sentCORE_READ_Asym"
)

analyze_lm_A3 <- function(data, task, outcome, predictors) {
  
  # Prepare complete cases and preserve original row numbers
  d <- data
  d$.row <- seq_len(nrow(d))
  vars <- c(outcome, predictors)
  d <- d[complete.cases(d[, vars, drop = FALSE]), ]
  
  # Original model
  f <- reformulate(predictors, response = outcome)
  model <- lm(f, data = d)
  
  cat("\n\n==========", task, "- ORIGINAL MODEL ==========\n")
  print(summary(model))
  
  # Diagnostic plots
  par(mfrow = c(2, 2))
  plot(model, which = 1:4)
  par(mfrow = c(1, 1))
  
  # Cook's distance, correctly matched to IDs
  used_rows <- as.integer(rownames(model.frame(model)))
  used_data <- d[match(used_rows, d$.row), ]
  
  cook <- data.frame(
    ID = used_data$ID,
    row = used_data$.row,
    Cook = as.numeric(cooks.distance(model))
  )
  
  cook <- cook[order(cook$Cook, decreasing = TRUE), ]
  
  cat("\n--- Top 10 Cook's distances ---\n")
  print(head(cook, 10))
  
  top3 <- head(cook, 3)
  
  cat("\n--- Participants selected for sensitivity analysis ---\n")
  print(data[top3$row, ])
  
  # Sensitivity model excluding those three observations
  d_sens <- d[!d$.row %in% top3$row, ]
  model_sens <- lm(f, data = d_sens)
  
  cat("\n--- SENSITIVITY MODEL ---\n")
  print(summary(model_sens))
  
  list(
    model = model,
    sensitivity = model_sens,
    cook = cook,
    top3 = top3
  )
}

#RUN LINEAR MODELS FOR BOTH ES AND MP

# Run all three tasks for ES
results_lm_A3_ES <- lapply(names(outcomes), function(task) {
  analyze_lm_A3(data_project, task, outcomes[[task]], predictors_A3_ES)
})
names(results_lm_A3_ES) <- names(outcomes)

# Run all three tasks for MP
results_lm_A3_MP <- lapply(names(outcomes), function(task) {
  analyze_lm_A3(data_project, task, outcomes[[task]], predictors_A3_MP)
})
names(results_lm_A3_MP) <- names(outcomes)

# Extract p-values for ES
p_ES_A3 <- c(
  PROD = coef(summary(results_lm_A3_ES$PROD$model))["ES", "Pr(>|t|)"],
  LISN = coef(summary(results_lm_A3_ES$LISN$model))["ES", "Pr(>|t|)"],
  READ = coef(summary(results_lm_A3_ES$READ$model))["ES", "Pr(>|t|)"]
)

results_p_ES_A3 <- data.frame(
  Outcome = names(p_ES_A3),
  p_raw = unname(p_ES_A3),
  p_Holm = unname(p.adjust(p_ES_A3, method = "holm"))
)

cat("\n--- Multiple-testing correction: ES ---\n")
print(results_p_ES_A3)

# Extract minimum p-value for categorical MP across its dummy coefficients
get_mp_pval <- function(mod) {
  co <- coef(summary(mod))
  idx <- grep("^MP", rownames(co))
  min(co[idx, "Pr(>|t|)"])
}

p_MP_A3 <- c(
  PROD = get_mp_pval(results_lm_A3_MP$PROD$model),
  LISN = get_mp_pval(results_lm_A3_MP$LISN$model),
  READ = get_mp_pval(results_lm_A3_MP$READ$model)
)

results_p_MP_A3 <- data.frame(
  Outcome = names(p_MP_A3),
  p_raw = unname(p_MP_A3),
  p_Holm = unname(p.adjust(p_MP_A3, method = "holm"))
)

cat("\n--- Multiple-testing correction: MP ---\n")
print(results_p_MP_A3)

cat("\n--- VIF ES Models ---\n")
print(vif(results_lm_A3_ES$PROD$model))
print(vif(results_lm_A3_ES$LISN$model))
print(vif(results_lm_A3_ES$READ$model))

cat("\n--- VIF MP Models ---\n")
print(vif(results_lm_A3_MP$PROD$model))
print(vif(results_lm_A3_MP$LISN$model))
print(vif(results_lm_A3_MP$READ$model))

# PENALIZED REGRESSION MODELS: RIDGE AND LASSO

fit_ridge_lasso_A3 <- function(data, task, outcome, current_predictors) {
  
  vars <- c(outcome, current_predictors)
  d <- data[complete.cases(data[, vars, drop = FALSE]), , drop = FALSE]
  d$Sex <- factor(d$Sex)
  
  X <- model.matrix(reformulate(current_predictors), data = d)[, -1, drop = FALSE]
  y <- d[[outcome]]
  
  set.seed(123)
  foldid <- sample(rep(1:10, length.out = nrow(X)))
  
  ridge <- cv.glmnet(x = X, y = y, alpha = 0, family = "gaussian", foldid = foldid, type.measure = "mse", standardize = TRUE)
  lasso <- cv.glmnet(x = X, y = y, alpha = 1, family = "gaussian", foldid = foldid, type.measure = "mse", standardize = TRUE)
  
  cat("\n==========", task, "==========\n")
  cat("N =", nrow(d), "\n")
  
  cat("\n--- RIDGE ---\n")
  cat("Lambda min:", ridge$lambda.min, "\n")
  cat("CV MSE:", min(ridge$cvm), "\n")
  
  cat("\n--- LASSO ---\n")
  cat("Lambda min:", lasso$lambda.min, "\n")
  cat("CV MSE:", min(lasso$cvm), "\n")
  
  return(list(ridge = ridge, lasso = lasso, n = nrow(d), X = X, y = y))
}

# Run penalized models for ES
results_penalized_ES <- lapply(names(outcomes), function(task) {
  fit_ridge_lasso_A3(data_project, task, outcomes[[task]], predictors_A3_ES)
})
names(results_penalized_ES) <- names(outcomes)

# Run penalized models for MP
results_penalized_MP <- lapply(names(outcomes), function(task) {
  fit_ridge_lasso_A3(data_project, task, outcomes[[task]], predictors_A3_MP)
})
names(results_penalized_MP) <- names(outcomes)

# ============================================================
# COMPARE THREE MODELS: LINEAR REGRESSION, RIDGE AND LASSO
# Repeated 10-fold cross-validation (5 repetitions)
# ============================================================

compare_predictive_A3 <- function(data, outcome, current_predictors, K = 10, repetition = 1) {
  
  vars <- c(outcome, current_predictors)
  d <- data[complete.cases(data[, vars, drop = FALSE]), , drop = FALSE]
  d$Sex <- factor(d$Sex)
  
  n <- nrow(d)
  y <- d[[outcome]]
  
  set.seed(1000 + repetition)
  folds <- sample(rep(1:K, length.out = n))
  
  pred <- matrix(NA_real_, nrow = n, ncol = 3, dimnames = list(NULL, c("Linear", "Ridge", "Lasso")))
  f <- reformulate(current_predictors, response = outcome)
  x_formula <- reformulate(current_predictors)
  
  for (k in seq_len(K)) {
    train <- d[folds != k, , drop = FALSE]
    test  <- d[folds == k, , drop = FALSE]
    
    # 1. Multiple linear regression
    fit_lm <- lm(f, data = train)
    pred[folds == k, "Linear"] <- predict(fit_lm, newdata = test)
    
    # Design matrices for Ridge and Lasso
    X_train <- model.matrix(x_formula, data = train)[, -1, drop = FALSE]
    X_test  <- model.matrix(x_formula, data = test)[, -1, drop = FALSE]
    y_train <- train[[outcome]]
    
    set.seed(2000 + repetition * 100 + k)
    inner_folds <- sample(rep(1:10, length.out = nrow(train)))
    
    # 2. Ridge regression
    fit_ridge <- cv.glmnet(X_train, y_train, alpha = 0, family = "gaussian", foldid = inner_folds, type.measure = "mse")
    pred[folds == k, "Ridge"] <- as.numeric(predict(fit_ridge, newx = X_test, s = "lambda.min"))
    
    # 3. Lasso regression
    fit_lasso <- cv.glmnet(X_train, y_train, alpha = 1, family = "gaussian", foldid = inner_folds, type.measure = "mse")
    pred[folds == k, "Lasso"] <- as.numeric(predict(fit_lasso, newx = X_test, s = "lambda.min"))
  }
  
  metrics <- do.call(rbind, lapply(colnames(pred), function(method) {
    errors <- y - pred[, method]
    data.frame(
      Outcome = outcome,
      Model = method,
      RMSE = sqrt(mean(errors^2)),
      MAE = mean(abs(errors)),
      R2_predictive = 1 - sum(errors^2) / sum((y - mean(y))^2)
    )
  }))
  metrics$Repetition <- repetition
  return(metrics)
}

# Run comparison for ES (5 repetitions)
results_repeated_ES <- do.call(rbind, lapply(1:5, function(r) {
  do.call(rbind, lapply(names(outcomes), function(task) {
    compare_predictive_A3(data_project, outcomes[[task]], predictors_A3_ES, K = 10, repetition = r)
  }))
}))

summary_repeated_ES <- results_repeated_ES %>%
  group_by(Outcome, Model) %>%
  summarise(across(c(RMSE, MAE, R2_predictive), list(mean = mean, SD = sd)), .groups = "drop")

cat("\n--- PREDICTIVE PERFORMANCE SUMMARY ANALYSIS 3: ES ---\n")
print(summary_repeated_ES)

# Run comparison for MP (5 repetitions)
results_repeated_MP <- do.call(rbind, lapply(1:5, function(r) {
  do.call(rbind, lapply(names(outcomes), function(task) {
    compare_predictive_A3(data_project, outcomes[[task]], predictors_A3_MP, K = 10, repetition = r)
  }))
}))

summary_repeated_MP <- results_repeated_MP %>%
  group_by(Outcome, Model) %>%
  summarise(across(c(RMSE, MAE, R2_predictive), list(mean = mean, SD = sd)), .groups = "drop")

cat("\n--- PREDICTIVE PERFORMANCE SUMMARY ANALYSIS 3: MP ---\n")
print(summary_repeated_MP)
----------------------------------------------------------------------------------------------

#Graphs for ES
newdata_ES <- data.frame(
  ES = seq(
    min(data_project$ES, na.rm = TRUE),
    max(data_project$ES, na.rm = TRUE),
    length.out = 100
  ),
  Age = mean(data_project$Age, na.rm = TRUE),
  Sex = factor("F", levels = levels(data_project$Sex))
)

#Models predictions
pred_prod <- predict(results_lm_A3_ES$PROD$model, newdata_ES, interval = "confidence")
pred_lisn <- predict(results_lm_A3_ES$LISN$model, newdata_ES, interval = "confidence")
pred_read <- predict(results_lm_A3_ES$READ$model, newdata_ES, interval = "confidence")

plot_ES <- rbind(
  data.frame(newdata_ES, as.data.frame(pred_prod), Task = "Production"),
  data.frame(newdata_ES, as.data.frame(pred_lisn), Task ="Listening"),
  data.frame(newdata_ES, as.data.frame(pred_read), Task = "Reading")
)

ggplot(plot_ES, aes(x = ES, y = fit, color = Task, fill = Task)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.12, color = NA) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c(
    "Production" = "#FFFF00",
    "Listening" = "#FF3399",
    "Reading" = "#FF9900"
  )) +
  scale_fill_manual(values = c(
    "Production" ="#FFFF00",
    "Listening" = "#FF3399",
    "Reading" = "#FF9900"
  )) +
  labs(
    title = "Edinburgh score and language lateralisation",
    subtitle = "Adjusted predictions for Age, Sex and Educational level (95% confidence intervals)",
    x = "Edinburgh score (ES)",
    y = "Functional language asymmetry",
    color = "Task",
    fill = "Task"
  ) +
  theme_minimal()+
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 12)
  )

# ggsave(
#   "analysis3_regressionES.png",
#   width = 10,
#   height = 7,
#   dpi = 300
# )

#Graphs for MP
newdata_MP <- expand.grid(
  MP = levels(data_project$MP),
  Age = mean(data_project$Age, na.rm = TRUE),
  Sex = factor("F", levels = levels(data_project$Sex))
)

# Calcolo le predizioni dei tre modelli
pred_prod <- predict(results_lm_A3_MP$PROD$model, newdata_MP, interval = "confidence")
pred_lisn <- predict(results_lm_A3_MP$LISN$model, newdata_MP, interval = "confidence")
pred_read <- predict(results_lm_A3_MP$READ$model, newdata_MP, interval = "confidence")

plot_MP <- rbind(
  data.frame(newdata_MP, as.data.frame(pred_prod), Task = "Production"),
  data.frame(newdata_MP, as.data.frame(pred_lisn), Task = "Listening"),
  data.frame(newdata_MP, as.data.frame(pred_read), Task = "Reading")
)

ggplot(plot_MP, aes(x = MP, y = fit, color = Task, group = Task)) +
  geom_line(
    position = position_dodge(width = 0.25),
    linewidth = 0.8
  ) +
  geom_point(
    position = position_dodge(width = 0.25),
    size = 3
  ) +
  geom_errorbar(
    aes(ymin = lwr, ymax = upr),
    position = position_dodge(width = 0.25),
    width = 0.1
  ) +
  scale_color_manual(values = c(
    "Production" = "#0000CC",
    "Listening" = "#CC3CCF",
    "Reading" = "#66CC00"
  )) +
  labs(
    title = "Manual preference and language lateralisation",
    subtitle = "Adjusted predictions for Age, Sex and Educational level",
    x = "Manual preference (MP)",
    y = "Functional language asymmetry",
    color = "Task"
  ) +
  theme_minimal()+
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 12)
  )

# ggsave(
#   "analysis3_regressionMP.png",
#   width = 10,
#   height = 7,
#   dpi = 300
# )

# GRAPHIC COMPARISON FOR ANALYSIS 3 (ES - Edinburgh Score)

# Design matrix for Ridge and Lasso (ES model)
X_new_ES <- model.matrix(
  reformulate(predictors_A3_ES),
  data = newdata_ES
)[, -1, drop = FALSE]

# Predictions of the three methods for each outcome
compare_models_ES <- function(task) {
  
  # Multiple linear regression
  pred_lm <- predict(
    results_lm_A3_ES[[task]]$model,
    newdata = newdata_ES
  )
  
  # Ridge
  pred_ridge <- as.numeric(predict(
    results_penalized_ES[[task]]$ridge,
    newx = X_new_ES,
    s = "lambda.min"
  ))
  
  # Lasso
  pred_lasso <- as.numeric(predict(
    results_penalized_ES[[task]]$lasso,
    newx = X_new_ES,
    s = "lambda.min"
  ))
  
  data.frame(
    ES = rep(newdata_ES$ES, 3),
    Prediction = c(pred_lm, pred_ridge, pred_lasso),
    Model = factor(
      rep(c("Linear regression", "Ridge", "Lasso"),
          each = nrow(newdata_ES)),
      levels = c("Linear regression", "Ridge", "Lasso")
    ),
    Task = task
  )
}

# Combine results from the three outcomes
comparison_data_ES <- do.call(
  rbind,
  lapply(names(outcomes), compare_models_ES)
)

# More readable task labels
comparison_data_ES$Task <- factor(
  comparison_data_ES$Task,
  levels = c("PROD", "LISN", "READ"),
  labels = c("Production", "Listening", "Reading")
)

# Comparative plot for ES model
ggplot(
  comparison_data_ES,
  aes(x = ES, y = Prediction, color = Model)
) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Task, scales = "free_y") +
  scale_color_manual(values = c(
    "Linear regression" = "#00ccff",
    "Ridge" = "#66ff00",
    "Lasso" = "#ff3300"
  )) +
  labs(
    title = "Comparison of regression models (Edinburgh Score)",
    subtitle = "Predictions adjusted for Age and Sex",
    x = "Edinburgh Score (ES)",
    y = "Predicted language asymmetry",
    color = "Model"
  ) +
  theme_minimal()+
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 12)
  )

# ggsave(
#   "analysis3_comparemodelsES.png",
#   width = 10,
#   height = 7,
#   dpi = 300
# )

# GRAPHIC PERFORMANCE COMPARISON: MANUAL PREFERENCE (MP) MODELS
# Point-range plot of Root Mean Squared Error (RMSE) across 10-fold CV

# Make sure Outcome has readable labels if needed (e.g., Production, Listening, Reading)
summary_repeated_MP$Outcome_Label <- factor(
  summary_repeated_MP$Outcome,
  levels = c("sentCORE_PROD_Asym", "sentCORE_LISN_Asym", "sentCORE_READ_Asym"),
  labels = c("Production", "Listening", "Reading"),
)
summary_repeated_MP$Model <- factor(
  summary_repeated_MP$Model,
  levels = c("Linear", "Ridge", "Lasso")
)

# Comparative point-range plot of predictive performance (RMSE)
ggplot(
  summary_repeated_MP,
  aes(x = Model, y = RMSE_mean, color = Model)
) +
  geom_pointrange(
    aes(ymin = RMSE_mean - RMSE_SD, ymax = RMSE_mean + RMSE_SD),
    position = position_dodge(0.5),
    size = 0.7,
    linewidth = 0.8
  ) +
  facet_wrap(~ Outcome_Label, scales = "free_y") +
  scale_color_manual(values = c(
    "Linear" = "#00ccff",
    "Ridge"  = "#66ff00",
    "Lasso"  = "#ff3300"
  )) +
  labs(
    title = "Predictive Performance Comparison (Manual Preference)",
    subtitle = "Root Mean Squared Error (RMSE) across repeated 10-fold CV",
    x = "Regression Model",
    y = "Mean RMSE (± SD)",
    color = "Model"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  )
 # 
 # ggsave(
 #   "analysis3_comparemodelsMP.png",
 #   width = 10,
 #   height = 7,
 #   dpi = 300
 # )

# ============================================================================
# 8. STATISTICAL ANALYSIS 4 - Structural Asymmetry and Language Lateralisation
# ============================================================================

# MULTIPLE LINEAR REGRESSION MODELS 

predictors_A4 <- c("GMasym", "WMasym", "CCvol", "Age", "Sex", "ES", "ITV")

outcomes <- c(
  PROD = "sentCORE_PROD_Asym",
  LISN = "sentCORE_LISN_Asym",
  READ = "sentCORE_READ_Asym"
)

analyze_lm_A4 <- function(data, task, outcome, predictors) {
  
  # Prepare complete cases and preserve original row numbers
  d <- data
  d$.row <- seq_len(nrow(d))
  vars <- c(outcome, predictors)
  d <- d[complete.cases(d[, vars, drop = FALSE]), ]
  
  # Original model
  f <- reformulate(predictors, response = outcome)
  model <- lm(f, data = d)
  
  cat("\n\n==========", task, "- A4 ORIGINAL MODEL ==========\n")
  print(summary(model))
  
  # Diagnostic plots
  par(mfrow = c(2, 2))
  plot(model, which = 1:4)
  par(mfrow = c(1, 1))
  
  # Cook's distance, correctly matched to IDs
  used_rows <- as.integer(rownames(model.frame(model)))
  used_data <- d[match(used_rows, d$.row), ]
  
  cook <- data.frame(
    ID = used_data$ID,
    row = used_data$.row,
    Cook = as.numeric(cooks.distance(model))
  )
  
  cook <- cook[order(cook$Cook, decreasing = TRUE), ]
  
  cat("\n--- Top 10 Cook's distances ---\n")
  print(head(cook, 10))
  
  # Select the three highest Cook's distances
  top3 <- head(cook, 3)
  
  cat("\n--- Participants selected for sensitivity analysis ---\n")
  print(data[top3$row, ])
  
  # Sensitivity model excluding those three observations
  d_sens <- d[!d$.row %in% top3$row, ]
  model_sens <- lm(f, data = d_sens)
  
  cat("\n--- A4 SENSITIVITY MODEL ---\n")
  print(summary(model_sens))
  
  # Save results
  list(
    model = model,
    sensitivity = model_sens,
    cook = cook,
    top3 = top3
  )
}

# Run all three tasks for A4
results_lm_A4 <- lapply(names(outcomes), function(task) {
  analyze_lm_A4(
    data = data_project,
    task = task,
    outcome = outcomes[[task]],
    predictors = predictors_A4
  )
})

names(results_lm_A4) <- names(outcomes)

#multiple testing correction has not been applied as all p-values were >0.05

cat("\n--- VIF Analysis 4 ---\n")
print(vif(results_lm_A4$PROD$model))
print(vif(results_lm_A4$LISN$model))
print(vif(results_lm_A4$READ$model))

# PENALIZED REGRESSION MODELS: RIDGE AND LASSO 
fit_ridge_lasso_A4 <- function(data, task, outcome, predictors) {
  
  vars <- c(outcome, predictors)
  d <- data[complete.cases(data[, vars, drop = FALSE]), , drop = FALSE]
  d$Sex <- factor(d$Sex)
  
  X <- model.matrix(reformulate(predictors), data = d)[, -1, drop = FALSE]
  y <- d[[outcome]]
  
  set.seed(123)
  foldid <- sample(rep(1:10, length.out = nrow(X)))
  
  ridge <- cv.glmnet(x = X, y = y, alpha = 0, family = "gaussian", foldid = foldid, type.measure = "mse", standardize = TRUE)
  lasso <- cv.glmnet(x = X, y = y, alpha = 1, family = "gaussian", foldid = foldid, type.measure = "mse", standardize = TRUE)
  
  cat("\n==========", task, "(A4) ==========\n")
  cat("N =", nrow(d), "\n")
  
  cat("\n--- RIDGE ---\n")
  cat("Lambda min:", ridge$lambda.min, "\n")
  cat("CV MSE:", min(ridge$cvm), "\n")
  
  cat("\n--- LASSO ---\n")
  cat("Lambda min:", lasso$lambda.min, "\n")
  cat("CV MSE:", min(lasso$cvm), "\n")
  
  return(list(ridge = ridge, lasso = lasso, n = nrow(d), X = X, y = y))
}

results_penalized_A4 <- lapply(names(outcomes), function(task) {
  fit_ridge_lasso_A4(
    data = data_project,
    task = task,
    outcome = outcomes[[task]],
    predictors = predictors_A4
  )
})

names(results_penalized_A4) <- names(outcomes)

# ============================================================
# COMPARE THREE MODELS: LINEAR REGRESSION, RIDGE AND LASSO 
# Repeated 10-fold cross-validation (5 repetitions)
# ============================================================
compare_predictive_A4 <- function(data, outcome, predictors, K = 10, repetition = 1) {
  
  vars <- c(outcome, predictors)
  d <- data[complete.cases(data[, vars, drop = FALSE]), , drop = FALSE]
  d$Sex <- factor(d$Sex)
  
  n <- nrow(d)
  y <- d[[outcome]]
  
  set.seed(1000 + repetition)
  folds <- sample(rep(1:K, length.out = n))
  
  pred <- matrix(NA_real_, nrow = n, ncol = 3, dimnames = list(NULL, c("Linear", "Ridge", "Lasso")))
  f <- reformulate(predictors, response = outcome)
  x_formula <- reformulate(predictors)
  
  for (k in seq_len(K)) {
    train <- d[folds != k, , drop = FALSE]
    test  <- d[folds == k, , drop = FALSE]
    
    # 1. Multiple linear regression
    fit_lm <- lm(f, data = train)
    pred[folds == k, "Linear"] <- predict(fit_lm, newdata = test)
    
    # Design matrices for Ridge and Lasso
    X_train <- model.matrix(x_formula, data = train)[, -1, drop = FALSE]
    X_test  <- model.matrix(x_formula, data = test)[, -1, drop = FALSE]
    y_train <- train[[outcome]]
    
    set.seed(2000 + repetition * 100 + k)
    inner_folds <- sample(rep(1:10, length.out = nrow(train)))
    
    # 2. Ridge regression
    fit_ridge <- cv.glmnet(X_train, y_train, alpha = 0, family = "gaussian", foldid = inner_folds, type.measure = "mse")
    pred[folds == k, "Ridge"] <- as.numeric(predict(fit_ridge, newx = X_test, s = "lambda.min"))
    
    # 3. Lasso regression
    fit_lasso <- cv.glmnet(X_train, y_train, alpha = 1, family = "gaussian", foldid = inner_folds, type.measure = "mse")
    pred[folds == k, "Lasso"] <- as.numeric(predict(fit_lasso, newx = X_test, s = "lambda.min"))
  }
  
  metrics <- do.call(rbind, lapply(colnames(pred), function(method) {
    errors <- y - pred[, method]
    data.frame(
      Outcome = outcome,
      Model = method,
      RMSE = sqrt(mean(errors^2)),
      MAE = mean(abs(errors)),
      R2_predictive = 1 - sum(errors^2) / sum((y - mean(y))^2)
    )
  }))
  metrics$Repetition <- repetition
  return(metrics)
}

# Run the comparison for all outcomes, 5 repetitions 
results_repeated_A4 <- list()

for (r in 1:5) {
  results_repeated_A4[[r]] <- do.call(
    rbind,
    lapply(names(outcomes), function(task) {
      compare_predictive_A4(
        data = data_project,
        outcome = outcomes[[task]],
        predictors = predictors_A4,
        K = 10,
        repetition = r
      )
    })
  )
}

results_repeated_A4 <- do.call(rbind, results_repeated_A4)

# Summarise mean and standard deviation across repetitions
summary_repeated_A4 <- results_repeated_A4 %>%
  group_by(Outcome, Model) %>%
  summarise(
    RMSE_mean = mean(RMSE),
    RMSE_SD = sd(RMSE),
    MAE_mean = mean(MAE),
    MAE_SD = sd(MAE),
    R2_mean = mean(R2_predictive),
    R2_SD = sd(R2_predictive),
    .groups = "drop"
  )

cat("\n--- ANALYSIS 4 (A4) PREDICTIVE PERFORMANCE SUMMARY ---\n")
print(summary_repeated_A4)

------------#second part of the analysis------------

#To investigate whether structural brain metrics vary across language dominance 
#profiles while controlling for potential confounders

table(data_project$DISSOCIATIONS_2categories)
table(data_project$Classif_MMHC)

cont_vars <- c("Age", "ITV", "ES")
cat_vars  <- c("Sex")

# 1. TEST for Classif_MMHC (3 groups)
cat("\n=== BILANCIAMENTO PER: Classif_mmhc (3 gruppi) ===\n")
for (v in cont_vars) {
  res_anova <- aov(as.formula(paste(v, "~ Classif_MMHC")), data = data_project)
  p_val <- summary(res_anova)[[1]][["Pr(>F)"]][1]
  cat(sprintf("[ANOVA] %-10s vs Classif_MMHC: p = %.4f %s\n", 
              v, p_val, if(p_val < 0.05) "--> SBILANCIATO" else "(Bilanciato)"))
}
tbl <- table(data_project$Sex, data_project$Classif_MMHC)
res_chi <- chisq.test(tbl)
cat(sprintf("[Chi2]  %-10s vs Classif_mmhc: p = %.4f %s\n", 
            "Sex", res_chi$p.value, if(res_chi$p.value < 0.05) "--> SBILANCIATO" else "(Bilanciato)"))

# 2. TEST for DISSOCIATIONS_2categories (2 groups - Congruent vs Crossed)
cat("\n=== BILANCIAMENTO PER: DISSOCIATIONS_2categories (2 gruppi) ===\n")
for (v in cont_vars) {
  res_ttest <- t.test(as.formula(paste(v, "~ DISSOCIATIONS_2categories")), data = data_project)
  p_val <- res_ttest$p.value
  cat(sprintf("[t-test] %-10s vs dissociation: p = %.4f %s\n", 
              v, p_val, if(p_val < 0.05) "--> SBILANCIATO" else "(Bilanciato)"))
}
tbl2 <- table(data_project$Sex, data_project$DISSOCIATIONS_2categories)
res_chi2 <- chisq.test(tbl2)
cat(sprintf("[Chi2]   %-10s vs dissociation: p = %.4f %s\n", 
            "sex", res_chi2$p.value, if(res_chi2$p.value < 0.05) "--> SBILANCIATO" else "(Bilanciato)"))

dv_list <- c("GMasym", "WMasym", "CCvol")
p_group <- numeric(length(dv_list))
names(p_group) <- dv_list

#First Ancova loop
cat("=======================================================\n")
cat("  1. ANCOVA & FULL DIAGNOSTICS: Classif_MMHC (3 levels)\n")
cat("=======================================================\n\n")

for (v in dv_list) {
  cat(sprintf("---------------------------------------------------\n"))
  cat(sprintf("DEPENDENT VARIABLE: %s\n", v))
  cat(sprintf("---------------------------------------------------\n"))
  
  # Main ANCOVA Model
  formula_str <- paste(v, "~ Classif_MMHC + Age + ES + Sex + ITV")
  base_model <- lm(as.formula(formula_str), data = data_project)
  
  # ANOVA Table
  tab_anova <- car::Anova(base_model, type = 2)
  print(tab_anova)
  p_group[v] <- tab_anova["Classif_MMHC", "Pr(>F)"]
  
  cat("\n--- ASSUMPTION CHECKS & DIAGNOSTICS ---\n")
  
  # Normality of Residuals (Shapiro-Wilk)
  shapiro_res <- shapiro.test(residuals(base_model))
  cat(sprintf("1. Residual Normality (Shapiro-Wilk): W = %.3f, p = %.4f %s\n",
              shapiro_res$statistic, shapiro_res$p.value,
              if(shapiro_res$p.value > 0.05) "(OK)" else "(Slight deviation)"))
  
  # Homogeneity of Slopes
  formula_slope <- paste(v, "~ Classif_MMHC * (ES + Age + ITV) + Sex")
  modello_slope <- lm(as.formula(formula_slope), data = data_project)
  
  car::Anova(modello_slope, type = 2)
  cat("\n--- Single interactions test for ", v, "---\n")
  print(car::Anova(modello_slope, type = 2))
  
  test_slope <- anova(base_model, modello_slope)
  
  p_slope <- test_slope$`Pr(>F)`[2]
  
  cat(sprintf(
    "2. Homogeneity of Slopes (all continuous covariates): p = %.4f %s\n",
    p_slope,
    if (!is.na(p_slope) && p_slope >= 0.05)
      "(No evidence of heterogeneous slopes)"
    else
      "(Warning: possible heterogeneous slopes)"
  ))
  
  #Test Levene
  model_data <- model.frame(base_model)
  levene_res <- car::leveneTest(residuals(base_model) ~ model_data$Classif_MMHC)
  cat(sprintf("4. Homogeneity of Variance (Levene): F = %.3f, p = %.4f %s\n",
              levene_res$`F value`[1], levene_res$`Pr(>F)`[1],
              if(levene_res$`Pr(>F)`[1] > 0.05) "(OK)" else "(Warning: heteroscedasticity detected)"))
  
  # Influential observations: Cook's distance
  cooks <- cooks.distance(base_model)
  
  n <- nobs(base_model)
  soglia_cook <- 4 / n
  
  casi_influenti <- which(cooks > soglia_cook)
  
  cat(sprintf(
    "3. Cook's Distance: max = %.4f; soglia 4/n = %.4f; casi segnalati = %d\n",
    max(cooks, na.rm = TRUE),
    soglia_cook,
    length(casi_influenti)
  ))
  
  if (length(casi_influenti) > 0) {
    print(data.frame(
      Riga = casi_influenti,
      Cook = cooks[casi_influenti]
    ))
  }
  # Visual Diagnostic Plots 
   par(mfrow = c(2, 2)); plot(modello); par(mfrow = c(1, 1))
  
  cat("\n\n")
}

# Holm correction across the three primary ANCOVA tests
p_holm <- p.adjust(p_group, method = "holm")

results_holm <- data.frame(
  Variabile = dv_list,
  p_orig = unname(p_group),
  p_Holm = unname(p_holm)
)
print(results_holm, row.names = FALSE)

#Second Ancova loop
cat("=======================================================\n")
cat("  2. ANCOVA & FULL DIAGNOSTICS: DISSOCIATIONS_2categories\n")
cat("=======================================================\n\n")

for (v in dv_list) {
  cat(sprintf("---------------------------------------------------\n"))
  cat(sprintf("DEPENDENT VARIABLE: %s\n", v))
  cat(sprintf("---------------------------------------------------\n"))
  
  formula_str <- paste(v, "~ DISSOCIATIONS_2categories + Age + ES + Sex + ITV")
  base_model <- lm(as.formula(formula_str), data = data_project)
  
  print(car::Anova(base_model, type = 2))
  cat("\n--- ASSUMPTION CHECKS & DIAGNOSTICS ---\n")
  
  # Shapiro-Wilk (normality of residuals)
  shapiro_res <- shapiro.test(residuals(base_model))
  cat(sprintf("1. Residual Normality (Shapiro-Wilk): W = %.3f, p = %.4f %s\n",
              shapiro_res$statistic, shapiro_res$p.value,
              if(shapiro_res$p.value > 0.05) "(OK)" else "(Slight deviation)"))
  
  # Slope Homogeneity
  formula_slope <- paste(v, "~ DISSOCIATIONS_2categories * (ES + Age + ITV) + Sex")
  modello_slope <- lm(as.formula(formula_slope), data = data_project)
  
  car::Anova(modello_slope, type = 2)
  cat("\n--- Test delle singole interazioni per", v, "---\n")
  print(car::Anova(modello_slope, type = 2))
  
  test_slope <- anova(base_model, modello_slope)
  p_slope <- test_slope$`Pr(>F)`[2]
  
  cat(sprintf(
    "2. Homogeneity of Slopes (all continuous covariates): p = %.4f %s\n",
    p_slope,
    if (!is.na(p_slope) && p_slope >= 0.05)
      "(No evidence of heterogeneous slopes)"
    else
      "(Warning: possible heterogeneous slopes)"
  ))
  
  #Test Levene (variance homogeneity)
  model_data <- model.frame(base_model)
  levene_res <- car::leveneTest(residuals(base_model) ~ model_data$DISSOCIATIONS_2categories)
  cat(sprintf("4. Homogeneity of Variance (Levene): F = %.3f, p = %.4f %s\n",
              levene_res$`F value`[1], levene_res$`Pr(>F)`[1],
              if(levene_res$`Pr(>F)`[1] > 0.05) "(OK)" else "(Warning: heteroscedasticity detected)"))

  # Influential observations: Cook's distance
  cooks <- cooks.distance(base_model)
  
  n <- nobs(base_model)
  soglia_cook <- 4 / n
  
  casi_influenti <- which(cooks > soglia_cook)
  
  cat(sprintf(
    "3. Cook's Distance: max = %.4f; soglia 4/n = %.4f; casi segnalati = %d\n",
    max(cooks, na.rm = TRUE),
    soglia_cook,
    length(casi_influenti)
  ))
  
  if (length(casi_influenti) > 0) {
    print(data.frame(
      Riga = casi_influenti,
      Cook = cooks[casi_influenti]
    ))
  }
  # Visual Diagnostic Plots 
  # par(mfrow = c(2, 2)); plot(modello); par(mfrow = c(1, 1))
}

#Dedicated analysis for the significant GMasym and WMasym models
cat("=======================================================\n")
cat("  POST-HOC ANALYSIS FOR GMasym (Classif_MMHC)\n")
cat("=======================================================\n\n")

mod_gm <- lm(GMasym ~ Classif_MMHC + Age + ES + Sex + ITV, data = data_project)

# Estimated Marginal Means (Covariate-adjusted means)
emm_gm <- emmeans(mod_gm, ~ Classif_MMHC)
print(emm_gm)
cat("\n")

# Pairwise contrasts with Tukey adjustment
pairs_gm <- pairs(emm_gm, adjust = "tukey")
print(pairs_gm)

cat("=======================================================\n")
cat("  POST-HOC ANALYSIS FOR WMasym (Classif_MMHC)\n")
cat("=======================================================\n\n")

#fit the model with explicit interaction for WMasym (as suggested by slope homogeneity tests)
mod_wm_int <- lm(WMasym ~ Classif_MMHC * Age + ES + Sex + ITV, data = data_project)

#estimate marginal means accounting for the interaction with Age
# (emmeans evaluates the groups at specific or mean levels of age)
# Valori rappresentativi dell'età
eta_media <- mean(data_project$Age, na.rm = TRUE)
eta_sd <- sd(data_project$Age, na.rm = TRUE)

emm_wm_int <- emmeans(
  mod_wm_int,
  ~ Classif_MMHC | Age,
  at = list(
    Age = c(
      eta_media - eta_sd,
      eta_media,
      eta_media + eta_sd
    )
  )
)

#comparisons between groups at each age value
pairs(emm_wm_int, adjust = "tukey")
#compare the regression slopes (trends) of Age across Classif_MMHC groups
trend_age <- emtrends(mod_wm_int, pairwise ~ Classif_MMHC, var = "Age")
print(trend_age)

---------------------------------------------------------------------------------------

#GRAPHS FOR MULTIPLE REGRESSIONS
# A. Grey matter (GMasym)
valid_gm <- na.omit(data_project$GMasym)
newdata_gm <- data.frame(
  GMasym = seq(min(valid_gm), max(valid_gm), length.out = 100),
  WMasym = mean(data_project$WMasym, na.rm = TRUE),
  CCvol = mean(data_project$CCvol, na.rm = TRUE),
  Age = mean(data_project$Age, na.rm = TRUE),
  Sex = data_project$Sex[!is.na(data_project$Sex)][1],
  ES = mean(data_project$ES, na.rm = TRUE),
  ITV= mean(data_project$ITV, na.rm = TRUE)
)

plot_data_gm <- rbind(
  data.frame(newdata_gm, as.data.frame(predict(results_lm_A4$PROD$model, newdata_gm, interval = "confidence")), Task = "Production"),
  data.frame(newdata_gm, as.data.frame(predict(results_lm_A4$LISN$model, newdata_gm, interval = "confidence")), Task = "Listening"),
  data.frame(newdata_gm, as.data.frame(predict(results_lm_A4$READ$model, newdata_gm, interval = "confidence")), Task = "Reading")
)

# B. White matter (WMasym)
valid_wm <- na.omit(data_project$WMasym)
newdata_wm <- data.frame(
  WMasym = seq(min(valid_wm), max(valid_wm), length.out = 100),
  GMasym = mean(data_project$GMasym, na.rm = TRUE),
  CCvol = mean(data_project$CCvol, na.rm = TRUE),
  Age = mean(data_project$Age, na.rm = TRUE),
  Sex = data_project$Sex[!is.na(data_project$Sex)][1],
  ES = mean(data_project$ES, na.rm = TRUE),
  ITV= mean(data_project$ITV, na.rm = TRUE)
)

plot_data_wm <- rbind(
  data.frame(newdata_wm, as.data.frame(predict(results_lm_A4$PROD$model, newdata_wm, interval = "confidence")), Task = "Production"),
  data.frame(newdata_wm, as.data.frame(predict(results_lm_A4$LISN$model, newdata_wm, interval = "confidence")), Task = "Listening"),
  data.frame(newdata_wm, as.data.frame(predict(results_lm_A4$READ$model, newdata_wm, interval = "confidence")), Task = "Reading")
)

# C. Corpus callosum (CCvol)
valid_cc <- na.omit(data_project$CCvol)
newdata_cc <- data.frame(
  CCvol = seq(min(valid_cc), max(valid_cc), length.out = 100),
  WMasym = mean(data_project$WMasym, na.rm = TRUE),
  GMasym = mean(data_project$GMasym, na.rm = TRUE),
  Age = mean(data_project$Age, na.rm = TRUE),
  Sex = data_project$Sex[!is.na(data_project$Sex)][1],
  ES = mean(data_project$ES, na.rm = TRUE),
  ITV= mean(data_project$ITV, na.rm = TRUE)
)

plot_data_cc <- rbind(
  data.frame(newdata_cc, as.data.frame(predict(results_lm_A4$PROD$model, newdata_cc, interval = "confidence")), Task = "Production"),
  data.frame(newdata_cc, as.data.frame(predict(results_lm_A4$LISN$model, newdata_cc, interval = "confidence")), Task = "Listening"),
  data.frame(newdata_cc, as.data.frame(predict(results_lm_A4$READ$model, newdata_cc, interval = "confidence")), Task = "Reading")
)

#create three individual graphs 
task_colors <- c("Production" = "#3399FF", "Listening" = "#00CC66", "Reading" = "#ff6699")

p_gm <- ggplot(plot_data_gm, aes(x = GMasym, y = fit, color = Task, fill = Task)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = task_colors) +
  scale_fill_manual(values = task_colors) +
  labs(x = "GM Asymmetry (GMasym)", y = "Predicted Lateralisation") +
  theme_minimal(base_size = 11)

p_wm <- ggplot(plot_data_wm, aes(x = WMasym, y = fit, color = Task, fill = Task)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = task_colors) +
  scale_fill_manual(values = task_colors) +
  labs(x = "WM Asymmetry (WMasym)", y = "") +
  theme_minimal(base_size = 11)

p_cc <- ggplot(plot_data_cc, aes(x = CCvol, y = fit, color = Task, fill = Task)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = task_colors) +
  scale_fill_manual(values = task_colors) +
  labs(x = "Corpus Callosum Vol (CCvol)", y = "") +
  theme_minimal(base_size = 11)

#merge graphs with patchwork
figure_multivariate <- (p_gm | p_wm | p_cc) +
  plot_layout(guides = "collect") + 
  plot_annotation(
    title = "Macro-structural Predictors and Task-based Language Lateralisation",
    subtitle = "Adjusted predictions for Age, Sex, ES and InTracranial volume",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  ) &
  theme(legend.position = "right")

print(figure_multivariate)

  ggsave(
    "analysis4_regression.png",
    width = 12,
    height = 7,
    dpi = 300
  )

# GRAPHIC MODEL COMPARISON FOR ANALYSIS 4

# Funzione generica per creare il dataframe di confronto per un determinato predittore
generate_data_for_predictor <- function(predictor_name) {
  
  compare_models_A4 <- function(task) {
    
    # 1. Creiamo la sequenza per il predittore scelto (100 punti)
    x_seq <- seq(
      min(data_project[[predictor_name]], na.rm = TRUE),
      max(data_project[[predictor_name]], na.rm = TRUE),
      length.out = 100
    )
    
    # 2. Creiamo newdata_A4 assegnando la sequenza al predittore e replicando le medie per gli altri
    newdata_A4 <- data.frame(
      GMasym = rep(mean(data_project$GMasym, na.rm = TRUE), 100),
      WMasym = rep(mean(data_project$WMasym, na.rm = TRUE), 100),
      CCvol  = rep(mean(data_project$CCvol, na.rm = TRUE), 100),
      Age    = rep(mean(data_project$Age, na.rm = TRUE), 100),
      Sex    = factor(rep("F", 100), levels = levels(data_project$Sex)),
      ES     = rep(mean(data_project$ES, na.rm = TRUE), 100),
      ITV    = rep(mean(data_project$ITV, na.rm = TRUE), 100)
    )
    
    # Sostituisce la colonna del predittore corrente con la sequenza esatta
    newdata_A4[[predictor_name]] <- x_seq
    
    # 3. Matrice di disegno per i modelli penalizzati
    X_new_A4 <- model.matrix(
      reformulate(predictors_A4),
      data = newdata_A4
    )[, -1, drop = FALSE]
    
    # 4. Predizioni dei tre modelli
    pred_lm <- predict(
      results_lm_A4[[task]]$model,
      newdata = newdata_A4
    )
    
    pred_ridge <- as.numeric(predict(
      results_penalized_A4[[task]]$ridge,
      newx = X_new_A4,
      s = "lambda.min"
    ))
    
    pred_lasso <- as.numeric(predict(
      results_penalized_A4[[task]]$lasso,
      newx = X_new_A4,
      s = "lambda.min"
    ))
    
    # 5. Restituisce il dataframe formattato
    data.frame(
      X_var = newdata_A4[[predictor_name]],
      Prediction = c(pred_lm, pred_ridge, pred_lasso),
      Model = factor(
        rep(c("Linear regression", "Ridge", "Lasso"),
            each = nrow(newdata_A4)),
        levels = c("Linear regression", "Ridge", "Lasso")
      ),
      Task = task,
      PredictorName = predictor_name
    )
  }
  
  comparison_data_A4 <- do.call(
    rbind,
    lapply(names(outcomes), compare_models_A4)
  )
  
  comparison_data_A4$Task <- factor(
    comparison_data_A4$Task,
    levels = c("PROD", "LISN", "READ"),
    labels = c("Production", "Listening", "Reading")
  )
  
  return(comparison_data_A4)
}

# Generiamo i dataset per i tre predittori
df_gm <- generate_data_for_predictor("GMasym")
df_wm <- generate_data_for_predictor("WMasym")
df_cc <- generate_data_for_predictor("CCvol")

# Funzione per creare il ggplot per il singolo predittore
make_plot_for_predictor <- function(df, x_label) {
  ggplot(
    df,
    aes(x = X_var, y = Prediction, color = Model)
  ) +
    geom_line(linewidth = 1) +
    facet_wrap(~ Task, scales = "free_y") +
    scale_color_manual(values = c(
      "Linear regression" = "#00ccff",
      "Ridge" = "#66ff00",
      "Lasso" = "#ff3300"
    )) +
    labs(
      title = paste("Predictor:", x_label),
      subtitle = "Predictions adjusted for other covariates",
      x = x_label,
      y = "Predicted language asymmetry",
      color = "Model"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11)
    )
}

p_gm <- make_plot_for_predictor(df_gm, "GMasym")
p_wm <- make_plot_for_predictor(df_wm, "WMasym")
p_cc <- make_plot_for_predictor(df_cc, "CCvol")

# Unione finale con patchwork in un'unica schermata
final_grid <- (p_gm / p_wm / p_cc) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Comparison of regression models across all anatomical predictors",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      legend.position = "bottom"
    )
  )

print(final_grid)

 # ggsave(
 #   "analysis4_comparemodels.png",
 #   width = 10,
 #   height = 8,
 #   dpi = 300
 #  )

#GRAPHS FOR CLASSIFICATION ANALYSIS WITH ANCOVA
#GMasym
df_gm <- as.data.frame(confint(emm_gm))

p_gm <- ggplot(
  df_gm,
  aes(x = Classif_MMHC, y = emmean,
      color = Classif_MMHC)
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.12, linewidth = 0.7
  ) +
  scale_color_manual(values = c(
    "ATYP" = "#ff9900", "TYP_MILD" = "#33cccc", "TYP_STRONG" = "#6666ff"
  )) +
  labs(
    title = "A. Grey matter asymmetry",
    subtitle = sprintf(
      "ANCOVA: p = %.3f; Holm-adjusted p = %.3f",
      p_group["GMasym"], p_holm["GMasym"]
    ),
    x = "Language lateralisation profile",
    y = "Adjusted estimated mean (GMasym)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

#WMasym
# Griglia di età nell'intervallo osservato
eta_seq <- seq(
  min(data_project$Age, na.rm = TRUE),
  max(data_project$Age, na.rm = TRUE),
  length.out = 100
)

# Valori di riferimento delle altre covariate
newdata_wm <- expand.grid(
  Age = eta_seq,
  Classif_MMHC = levels(factor(data_project$Classif_MMHC))
)

newdata_wm$ES <- mean(data_project$ES, na.rm = TRUE)
newdata_wm$ITV <- mean(data_project$ITV, na.rm = TRUE)

newdata_wm$Sex <- factor(
  levels(factor(data_project$Sex))[1],
  levels = levels(factor(data_project$Sex))
)

#Predictions and confidence interval of 95%
pred_wm <- predict(
  mod_wm_int,
  newdata = newdata_wm,
  interval = "confidence"
)

df_wm <- cbind(newdata_wm, as.data.frame(pred_wm))

p_wm <- ggplot(
  df_wm,
  aes(
    x = Age, y = fit,
    color = Classif_MMHC,
    fill = Classif_MMHC
  )
) +
  geom_ribbon(
    aes(ymin = lwr, ymax = upr),
    alpha = 0.12, color = NA
  ) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = c(
    "ATYP" = "#ff9900", "TYP_MILD" = "#33cccc", "TYP_STRONG" = "#6666ff"
  )) +
  scale_fill_manual(values = c(
    "ATYP" = "#ff9900", "TYP_MILD" = "#33cccc", "TYP_STRONG" = "#6666ff"
  )) +
  labs(
    title = "B. White matter asymmetry",
    subtitle = "Group × Age interaction",
    x = "Age",
    y = "Predicted WMasym",
    color = "Language profile",
    fill = "Language profile"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

#merge into the same 
figure_GM_WM <- (p_gm | p_wm) +
  plot_layout(widths = c(1, 1.5)) +
  plot_annotation(
    title = "Structural asymmetry across language lateralisation profiles",
    subtitle = "Adjusted estimates with 95% confidence intervals",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14)
    )
  )

print(figure_GM_WM)

ggsave(
  "analysis4_ANCOVA.png",
  width = 11,
  height = 8,
  dpi = 300
 )
