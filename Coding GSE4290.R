# =========================
# 1. INSTALL & LOAD PACKAGE
# =========================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) BiocManager::install("AnnotationDbi")
if (!requireNamespace("hgu133plus2.db", quietly = TRUE)) BiocManager::install("hgu133plus2.db")

if (!requireNamespace("randomForest", quietly = TRUE)) install.packages("randomForest")
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
if (!requireNamespace("caret", quietly = TRUE)) install.packages("caret")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(GEOquery)
library(limma)
library(randomForest)
library(pROC)
library(caret)
library(ggplot2)
library(AnnotationDbi)
library(hgu133plus2.db)

# =========================
# 2. AMBIL DATA GEO
# =========================

gset <- getGEO("GSE4290", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]

dim(exprs(gset))
head(pData(gset))
unique(pData(gset)$`Histopathological diagnostic:ch1`)

# =========================
# 3. BUAT LABEL KELAS
# =========================

diagnosis <- pData(gset)$`Histopathological diagnostic:ch1`

group <- ifelse(
  grepl("glioblastoma", diagnosis, ignore.case = TRUE),
  "GBM",
  ifelse(
    grepl("non-tumor", diagnosis, ignore.case = TRUE),
    "Normal",
    NA
  )
)

table(group, useNA = "ifany")

# =========================
# 4. FILTER SAMPEL
# =========================

keep <- !is.na(group)

gset_filtered <- gset[, keep]
group_filtered <- group[keep]

table(group_filtered)

# =========================
# 5. AMBIL EXPRESSION MATRIX
# =========================

ex <- exprs(gset_filtered)
dim(ex)

# =========================
# 6. CEK LOG2 TRANSFORM
# =========================

qx <- quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE)
qx

needs_log <- (qx[5] > 100) || (qx[6] - qx[1] > 50)

if (needs_log) {
  ex_log2 <- log2(ex + 1)
  cat("Data di-log2 transform\n")
} else {
  ex_log2 <- ex
  cat("Data sudah tampak dalam skala log\n")
}

# =========================
# 7. DESIGN MATRIX
# =========================

group_factor <- factor(group_filtered)
group_factor <- relevel(group_factor, ref = "Normal")

design <- model.matrix(~ group_factor)
head(design)

# =========================
# 8. LIMMA
# =========================

fit <- lmFit(ex_log2, design)
fit2 <- eBayes(fit)

results <- topTable(
  fit2,
  coef = "group_factorGBM",
  number = Inf,
  sort.by = "P"
)

head(results)
dim(results)

# =========================
# 9. FILTER DEG
# =========================

deg <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]
dim(deg)

# =========================
# 10. AMBIL TOP 100
# =========================

deg_sorted <- deg[order(deg$adj.P.Val), ]
top_n <- min(100, nrow(deg_sorted))
top100 <- deg_sorted[1:top_n, ]

dim(top100)

# =========================
# 11. SIAPKAN DATA ML
# =========================

ml_data <- ex_log2[rownames(top100), ]
ml_data_t <- t(ml_data)

ml_df <- as.data.frame(ml_data_t)
ml_df$Group <- factor(group_filtered)

colnames(ml_df) <- make.names(colnames(ml_df))

dim(ml_df)
table(ml_df$Group)

# =========================
# 12. TRAIN / TEST SPLIT
# =========================

set.seed(123)

train_index <- sample(seq_len(nrow(ml_df)), size = 0.8 * nrow(ml_df))
train_data <- ml_df[train_index, ]
test_data  <- ml_df[-train_index, ]

dim(train_data)
dim(test_data)

# =========================
# 13. LOGISTIC REGRESSION
# =========================

model_glm <- glm(Group ~ ., data = train_data, family = binomial)

glm_prob <- predict(model_glm, newdata = test_data, type = "response")
glm_pred <- ifelse(glm_prob > 0.5, "GBM", "Normal")
glm_pred <- factor(glm_pred, levels = levels(test_data$Group))

table(Predicted = glm_pred, Actual = test_data$Group)

glm_acc <- mean(glm_pred == test_data$Group)
glm_acc

# =========================
# 14. RANDOM FOREST
# =========================

set.seed(123)

rf_model <- randomForest(
  Group ~ .,
  data = train_data,
  ntree = 500,
  importance = TRUE
)

rf_pred <- predict(rf_model, newdata = test_data)

table(Predicted = rf_pred, Actual = test_data$Group)

rf_acc <- mean(rf_pred == test_data$Group)
rf_acc

# =========================
# 15. FEATURE IMPORTANCE
# =========================

importance_values <- importance(rf_model)

if ("MeanDecreaseGini" %in% colnames(importance_values)) {
  importance_metric <- importance_values[, "MeanDecreaseGini"]
} else {
  importance_metric <- importance_values[, 1]
}

importance_df <- data.frame(
  Gene = rownames(importance_values),
  Importance = importance_metric
)

importance_df <- importance_df[order(-importance_df$Importance), ]

head(importance_df, 10)

varImpPlot(rf_model, n.var = 20)

# =========================
# 16. ROC & AUC
# =========================

rf_prob <- predict(rf_model, newdata = test_data, type = "prob")[, "GBM"]

roc_obj <- roc(response = test_data$Group, predictor = rf_prob)

plot(roc_obj)
auc_value <- auc(roc_obj)
auc_value

# =========================
# 17. CROSS-VALIDATION
# =========================

set.seed(123)

control <- trainControl(
  method = "cv",
  number = 5
)

rf_cv <- train(
  Group ~ .,
  data = ml_df,
  method = "rf",
  trControl = control,
  ntree = 500
)

rf_cv

# =========================
# 18. VOLCANO PLOT
# =========================

results$Significant <- ifelse(
  results$adj.P.Val < 0.05 & abs(results$logFC) > 1,
  "Yes",
  "No"
)

ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
  geom_point(alpha = 0.6) +
  theme_minimal()

# =========================
# 19. SIMPAN HASIL
# =========================

write.csv(deg, "DEG_results.csv", row.names = TRUE)
write.csv(top100, "Top100_genes.csv", row.names = TRUE)
write.csv(importance_df, "Feature_importance.csv", row.names = FALSE)

# =========================
# 20. MAPPING PROBE -> GENE SYMBOL
# =========================

annotation(gset_filtered)

importance_df$ProbeID <- gsub("^X", "", importance_df$Gene)

top10_probes <- head(importance_df$ProbeID, 10)

gene_symbols <- mapIds(
  hgu133plus2.db,
  keys = top10_probes,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

gene_symbols

importance_df$GeneSymbol <- mapIds(
  hgu133plus2.db,
  keys = importance_df$ProbeID,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

head(importance_df, 10)