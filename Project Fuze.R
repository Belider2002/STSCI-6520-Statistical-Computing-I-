library(IRanges)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(randomForest)
library(pROC)
library(tidyr)
library(ggplot2)
setwd("/Users/baolide/Desktop/MPS Acct/STSCI 6520 Statistical Computing I (2024FA)/Project")

# Keep only gene data on chr2
expression <- import("wgEncodeRikenCageGm12878NucleusPapPlusSignalRep1.bigWig", format = "bigWig")
expression<- expression[seqnames(expression) == "chr2"]
file_path1 <- "/Users/baolide/Desktop/MPS Acct/STSCI 6520 Statistical Computing I (2024FA)/Project/gencode.v7.annotation.gtf"
gtf_data <- import(file_path1, format = "gtf")
gtf_data <- gtf_data[seqnames(gtf_data) == "chr2"]
genes <- gtf_data[gtf_data$type == "gene"]
gene_gr <- GRanges(
  seqnames = seqnames(genes),
  ranges   = ranges(genes),
  gene_id  = genes$gene_id
)
gene_gr <- unique(gene_gr)

# Using ranges of gene on annotations to matches scores of expression and modification
matched_expression <- mergeByOverlaps(gene_gr, expression)
score_expression <-  aggregate(cbind(score, width(matched_expression$expression)) ~ gene_id, matched_expression, sum, na.rm = TRUE)
score_expression$expression_level <- score_expression$score/score_expression$V2
save(score_expression, file = "score_expression.RData")
load("score_expression.Rdata")
score_expression <- score_expression %>% select(-score, -V2)

H3k4me3 <- import("wgEncodeUwHistoneHl60H3k4me3StdAln_2Reps.norm5.rawsignal.bw", format = "bigWig")
H3k4me3<- H3k4me3[seqnames(H3k4me3) == "chr2"]
matched_H3k4me3 <- mergeByOverlaps(gene_gr, H3k4me3)
matched_H3k4me3 <- aggregate(cbind(score, width(matched_H3k4me3$H3k4me3)) ~ gene_id, matched_H3k4me3, sum, na.rm = TRUE)
matched_H3k4me3$mean.H3k4me3 <- matched_H3k4me3$score/matched_H3k4me3$V2
matched_H3k4me3 <- matched_H3k4me3 %>% select(-score, -V2)


H2az <- import("wgEncodeBroadHistoneHepg2H2azStdAln_2Reps.norm5.rawsignal.bw", format = "bigWig")
H2az<- H2az[seqnames(H2az) == "chr2"]
matched_H2az <- mergeByOverlaps(gene_gr, H2az)
matched_H2az <- aggregate(cbind(score, width(matched_H2az$H2az)) ~ gene_id, matched_H2az, sum, na.rm = TRUE)
matched_H2az$mean.H2az <- matched_H2az$score/matched_H2az$V2
matched_H2az <- matched_H2az %>% select(-score, -V2)

H3k4me2 <- import("wgEncodeBroadHistoneGm12878H3k4me2StdAln_2Reps.norm5.rawsignal.bw", format = "bigWig")
H3k4me2 <- H3k4me2[seqnames(H3k4me2) == "chr2"]
matched_H3k4me2 <- mergeByOverlaps(gene_gr, H3k4me2)
matched_H3k4me2 <- aggregate(cbind(score, width(matched_H3k4me2$H3k4me2)) ~ gene_id, matched_H3k4me2, sum, na.rm = TRUE)
matched_H3k4me2$mean.H3k4me2 <- matched_H3k4me2$score/matched_H3k4me2$V2
matched_H3k4me2 <- matched_H3k4me2 %>% select(-score, -V2)

Dnase <- import("wgEncodeOpenChromDnaseIpsAln_3Reps.norm5.rawsignal.bw", format = "bigWig")
Dnase <- Dnase[seqnames(Dnase) == "chr2"]
matched_Dnase<- mergeByOverlaps(gene_gr, Dnase)
matched_Dnase <- aggregate(cbind(score, width(matched_Dnase$Dnase)) ~ gene_id, matched_Dnase, sum, na.rm = TRUE)
matched_Dnase$mean.Dnase <- matched_Dnase$score/matched_Dnase$V2
matched_Dnase <- matched_Dnase %>% select(-score, -V2)

H3k9ac <- import("wgEncodeBroadHistoneGm12878H3k9acStdAln_2Reps.norm5.rawsignal.bw", format = "bigWig")
H3k9ac <- H3k9ac[seqnames(H3k9ac) == "chr2"]
matched_H3k9ac <- mergeByOverlaps(gene_gr, H3k9ac)
matched_H3k9ac <-  aggregate(cbind(score, width(matched_H3k9ac$H3k9ac)) ~ gene_id, matched_H3k9ac, sum, na.rm = TRUE)
matched_H3k9ac$mean.H3k9ac <- matched_H3k9ac$score/matched_H3k9ac$V2
matched_H3k9ac <- matched_H3k9ac %>% select(-score, -V2)


H4k20me1 <- import("wgEncodeBroadHistoneH1hescH4k20me1StdAln_2Reps.norm5.rawsignal.bw", format = "bigWig")
H4k20me1 <- H4k20me1[seqnames(H4k20me1) == "chr2"]
matched_H4k20me1 <- mergeByOverlaps(gene_gr, H4k20me1)
matched_H4k20me1 <-  aggregate(cbind(score, width(matched_H4k20me1$H4k20me1)) ~ gene_id, matched_H4k20me1, sum, na.rm = TRUE)
matched_H4k20me1$mean.H4k20me1<- matched_H4k20me1$score/matched_H4k20me1$V2
matched_H4k20me1 <- matched_H4k20me1 %>% select(-score, -V2)
 
TfbsK562EjundControl <- import("wgEncodeUchicagoTfbsK562EjundControlAln_3Reps.norm5.rawsignal.bw", format = "bigWig")
TfbsK562EjundControl <- TfbsK562EjundControl[seqnames(TfbsK562EjundControl) == "chr2"]
matched_TfbsK562EjundControl <- mergeByOverlaps(gene_gr, TfbsK562EjundControl)
matched_TfbsK562EjundControl <-  aggregate(cbind(score, width(matched_TfbsK562EjundControl$TfbsK562EjundControl)) ~ gene_id, matched_TfbsK562EjundControl, sum, na.rm = TRUE)
matched_TfbsK562EjundControl$mean.TfbsK562EjundControl<- matched_TfbsK562EjundControl$score/matched_TfbsK562EjundControl$V2
matched_TfbsK562EjundControl <- matched_TfbsK562EjundControl %>% select(-score, -V2)

Helas3Control <- import("wgEncodeBroadHistoneHelas3ControlStdAln_2Reps.norm5.rawsignal.bw", format = "bigWig")
Helas3Control <- Helas3Control[seqnames(Helas3Control) == "chr2"]
matched_Helas3Control <- mergeByOverlaps(gene_gr, Helas3Control)
matched_Helas3Control <-  aggregate(cbind(score, width(matched_Helas3Control$Helas3Control)) ~ gene_id, matched_Helas3Control, sum, na.rm = TRUE)
matched_Helas3Control$mean.Helas3Control<- matched_Helas3Control$score/matched_Helas3Control$V2
matched_Helas3Control <- matched_Helas3Control %>% select(-score, -V2)

H3k9me3 <- import("wgEncodeBroadHistoneK562H3k9me3StdAln_2Reps.norm5.rawsignal.bw", format = "bigWig")
H3k9me3 <- H3k9me3[seqnames(H3k9me3) == "chr2"]
matched_H3k9me3 <- mergeByOverlaps(gene_gr, H3k9me3)
matched_H3k9me3 <-  aggregate(cbind(score, width(matched_H3k9me3$H3k9me3)) ~ gene_id, matched_H3k9me3, sum, na.rm = TRUE)
matched_H3k9me3$mean.H3k9me3<- matched_H3k9me3$score/matched_H3k9me3$V2
matched_H3k9me3 <- matched_H3k9me3 %>% select(-score, -V2)

H4k20me_2 <- import("wgEncodeBroadHistoneK562H4k20me1StdAln_2Reps.norm5.rawsignal.bw", format = "bigWig")
H4k20me_2 <- H4k20me_2[seqnames(H4k20me_2) == "chr2"]
matched_H4k20me_2 <- mergeByOverlaps(gene_gr, H4k20me_2)
matched_H4k20me_2 <-  aggregate(cbind(score, width(matched_H4k20me_2$H4k20me_2)) ~ gene_id, matched_H4k20me_2, sum, na.rm = TRUE)
matched_H4k20me_2$mean.H4k20me_2<- matched_H4k20me_2$score/matched_H4k20me_2$V2
matched_H4k20me_2 <- matched_H4k20me_2 %>% select(-score, -V2)
 
df_list <- list(score_expression, matched_H3k4me2, matched_H3k4me3, matched_H3k9ac, matched_H2az, matched_Dnase, matched_H4k20me1, matched_TfbsK562EjundControl, matched_Helas3Control, matched_H3k9me3, matched_H4k20me_2)
)


merged_df <- score_expression
for (i in 2:length(df_list)) {
  merged_df <- merge(merged_df, df_list[[i]], by = "gene_id")
}

# Save only gene_id and scores
all_variables<-merged_df[,c("gene_id", "expression_level", "mean.H3k4me2", "mean.H3k4me3", "mean.H2az", "mean.Dnase", "mean.H3k9ac", 
                            "mean.H4k20me1", "mean.TfbsK562EjundControl", "mean.Helas3Control", "mean.H3k9me3", "mean.H4k20me_2")]

# add pseudocount of 0.1 on X variabels 
pseudocount <- all_variables %>%
  mutate(across(3:12, ~ . + 0.1))
# log2 transform Y varible 
pseudocount_log2 <- pseudocount %>%
  mutate(across(2:12, log2))



#-------------------------------------------------------------------------------

#Classification moedel

# Turn expression level to binary variable
pseudocount_log2$on_off <- ifelse(pseudocount_log2$expression_level >= 0.1, 1, 0 )
pseudocount_log2 <- pseudocount_log2 %>% mutate(on_off = factor(on_off, levels = c(0, 1), labels = c("off", "on")))


# Try one time Randomforest model       
n <- nrow(pseudocount_log2)
train_indicies <- sample(1:n, size = floor(0.4*n))
train_data <- pseudocount_log2 %>%
  slice(train_indicies)
test_data <- pseudocount_log2 %>% 
  slice(-train_indicies)

m_forest <- randomForest(on_off ~. -gene_id - expression_level, data = train_data,
                         importance = TRUE)

importance(m_forest)

#Plot ROC curve

pred_onoff <- predict(m_forest, test_data, type = "prob")[, 2]
actual <- test_data$on_off

roc_obj <- roc(actual, pred_onoff)
plot(roc_obj)
auc_value <- auc(roc_obj)
auc_value

#Cross validation
mtry <- ceiling(sqrt(ncol(pseudocount_log2)))
repeat_cv_random <- trainControl(method='repeatedcv', number=10, repeats=5)
forest <- train(on_off ~. -gene_id - expression_level, data=train_data, method = 'rf',
                tuneLength = 9,
                trControl = repeat_cv_random, mtric = 'Accuracy', )
forest
plot(forest)
# mtry=4 is recommended






# For reproducibility
set.seed(123)

# Number of iterations
n_iter <- 30

# Store results
auc_values <- numeric(n_iter)
var_importance_list <- list()

for (i in 1:n_iter) {
  # Random train/test split
  n <- nrow(pseudocount_log2)
  train_indicies <- sample(1:n, size = floor(0.4 * n))
  train_data <- pseudocount_log2 %>% slice(train_indicies)
  test_data <- pseudocount_log2 %>% slice(-train_indicies)
  
# Fit Random Forest model
# Adjust the formula and data as needed for your actual predictor variables
m_forest <- randomForest(on_off ~ . - gene_id - expression_level, 
                           data = train_data,
                           mtry = 4,
                           importance = TRUE)
  
# Predictions on test set
pred_onoff <- predict(m_forest, test_data, type = "prob")[, 2]
actual <- test_data$on_off
  
# Compute ROC and AUC
roc_obj <- roc(actual, pred_onoff)
auc_value <- auc(roc_obj)
auc_values[i] <- auc_value
  
# Extract variable importance
var_imp <- importance(m_forest)
  
# Store importance in a list
var_importance_list[[i]] <- var_imp
}

# Calculate the mean AUC across all iterations
mean_auc <- mean(auc_values)
mean_auc
# Average variable importance across all iterations
# First, combine all importance matrices into a single array
importance_array <- simplify2array(var_importance_list)
# importance_array will be a 3D array if multiple measures of importance are included (e.g. MeanDecreaseGini and MeanDecreaseAccuracy)
# Typically, importance(m_forest) returns a matrix with predictors in rows and importance metrics in columns.
# We can compute means across the third dimension (the iterations) using apply.
mean_importance <- apply(importance_array, c(1, 2), mean)

# mean_auc now contains the average AUC
print(paste("Mean AUC over", n_iter, "runs:", mean_auc))

# mean_importance now contains the averaged importance values for each predictor
print(mean_importance)

# Assuming you have run your for-loop code and have:
# var_importance_list: a list of importance matrices from each iteration
# Each matrix typically has variables as rows and various importance metrics as columns,
# including "MeanDecreaseGini".


# Convert the list of matrices to a single long-format data frame
df_list <- lapply(seq_along(var_importance_list), function(i) {
imp_mat <- var_importance_list[[i]]
  
# Convert the matrix to a data frame
df <- as.data.frame(imp_mat)
  
# Add a column for the variable names from row names
df$Variable <- rownames(imp_mat)
  
# Add an iteration column
df$Iteration <- i
  
df
})

var_imp_long <- do.call(rbind, df_list)

# Now var_imp_long is a data frame with columns for each importance metric, plus Variable and Iteration.
# Typically, MeanDecreaseGini is one of the columns returned by importance().
# Check the column names of var_imp_long to confirm the exact name.

# Example column names often are:
# "MeanDecreaseAccuracy" and "MeanDecreaseGini"
# We'll assume "MeanDecreaseGini" is the correct column name.

# Create a box plot of MeanDecreaseGini for each variable
ggplot(var_imp_long, aes(x = Variable, y = MeanDecreaseGini)) +
  geom_boxplot() +
  coord_flip() +  # Flip coordinates to make variables more readable on the y-axis
  labs(title = "Mean Decrease Gini Importance Distribution",
       x = "Variables",
       y = "Mean Decrease Gini") +
  theme_minimal()
  










#-------------------------------------------------------------------------------
#Linear Regression
library(relaimpo)  # for relative importance

# Set seed for reproducibility
set.seed(123)

# Number of iterations
n_iter <- 30

# Store results
cor_values <- numeric(n_iter)
rmse_values <- numeric(n_iter)
rel_importance_list <- list()

for (i in 1:n_iter) {
# Random train/test split
n <- nrow(pseudocount_log2)
train_indices <- sample(seq_len(n), size = floor(0.4 * n))
train_data <- pseudocount_log2[train_indices, ]
test_data <- pseudocount_log2[-train_indices, ]
  
# Fit linear regression model
# Adjust the formula to exclude non-numeric or unwanted columns
# For example, if you want to predict expression_level from all other numeric columns 
# except gene_id and on_off:
predictor_cols <- setdiff(names(pseudocount_log2), c("gene_id", "on_off", "expression_level"))
formula_str <- paste("expression_level ~", paste(predictor_cols, collapse = " + "))
  
lm_model <- lm("expression_level ~ mean.H3k4me2 + mean.H3k4me3 + mean.H2az + mean.Dnase + mean.H3k9ac + mean.H4k20me1 + mean.TfbsK562EjundControl + mean.Helas3Control + mean.H3k9me3 + mean.H4k20me_2"
                 , data = train_data)
  
# Predict on test set
predictions <- predict(lm_model, newdata = test_data)
  
# Calculate Pearson correlation
actual <- test_data$expression_level
cor_values[i] <- cor(predictions, actual, use = "complete.obs")
  
# Calculate RMSE
rmse_values[i] <- sqrt(mean((predictions - actual)^2, na.rm = TRUE))
  
# Compute relative importance
# lmg metric is a common choice for linear models
rel_imp <- calc.relimp(lm_model, type = "lmg", rela = TRUE)
# rel_imp$lmg gives a named vector of relative importances (sum to 1)
  
rel_importance_list[[i]] <- rel_imp$lmg
}

# Compute average correlation and RMSE
mean_cor <- mean(cor_values)
mean_rmse <- mean(rmse_values)

cat("Average Pearson Correlation:", mean_cor, "\n")
cat("Average RMSE:", mean_rmse, "\n")

# Combine relative importance across iterations
# First ensure all vectors have the same ordering of predictors
predictors <- names(rel_importance_list[[1]])
rel_mat <- do.call(rbind, rel_importance_list)  # each row is one iteration
colnames(rel_mat) <- predictors

mean_rel_importance <- colMeans(rel_mat)

# mean_rel_importance now holds the average relative contribution to R² of each predictor

# Plot the relative importance(Bar Chart)
# Convert to a data frame for ggplot
rel_imp_df <- data.frame(Predictor = names(mean_rel_importance),
                         MeanRelativeImportance = mean_rel_importance)

ggplot(rel_imp_df, aes(x = reorder(Predictor, MeanRelativeImportance), 
                       y = MeanRelativeImportance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Average Relative Contribution to R²",
       x = "Predictors",
       y = "Mean Relative Importance (R²)") +
  theme_minimal()

# Box Plot
# Combine all relative importance vectors into a data frame
var_imp_df <- do.call(rbind, rel_importance_list)  
var_imp_df <- as.data.frame(var_imp_df)

# Add iteration column
var_imp_df$Iteration <- 1:nrow(var_imp_df)

# Convert to long format
var_imp_long <- var_imp_df %>%
  pivot_longer(cols = -Iteration, names_to = "Variable", values_to = "Importance")

# Create boxplot of variable importance
ggplot(var_imp_long, aes(x = Importance, y = Variable)) +
  geom_boxplot() +
  labs(title = "Variable Importance in Regression",
       x = "Relative Contribution to R²",
       y = "Variables") +
  theme_minimal() 
 

