# Differential Gene Expression (DEG) + Machine Learning Classification (GBM vs Normal) — GSE4290  
## limma DEG Analysis & ML-based Classification (LogReg vs Random Forest)

## 📌 Overview
This project analyzes gene expression differences between glioblastoma (GBM) and non-tumor (Normal) brain tissues using the public GEO dataset GSE4290. Differentially Expressed Genes (DEGs) are identified with limma, then the Top 100 DEGs are used as features to test whether a machine learning model can classify samples as GBM or Normal using Logistic Regression (baseline) and Random Forest (main model).

## 🧬 Dataset
- Accession: GSE4290  
- Platform: GPL570 (Affymetrix Human Genome U133 Plus 2.0 Array)  
- Initial samples: 180  
- After filtering: 100 samples (77 GBM, 23 Normal)  
- Comparison: GBM vs Normal (non-tumor)

Note: The expression data are fetched directly from GEO using `GEOquery::getGEO()`, so storing the raw dataset in this repository is optional.

## 🔬 Analysis Workflow
1. Retrieve data from GEO (GSE4290)
2. Assign class labels from histopathological metadata (GBM vs Normal)
3. Filter samples (keep only GBM and non-tumor)
4. Apply log2 transformation (to stabilize expression distributions)
5. Perform differential expression analysis with limma (lmFit + eBayes)
6. Select significant DEGs using thresholds:
   - adjusted p-value < 0.05
   - |logFC| > 1
7. Feature selection: Top 100 genes ranked by adjusted p-value
8. Build ML dataset:
   - Features: expression levels of Top 100 genes
   - Label: GBM / Normal
9. Split data into train/test (80/20)
10. Train models:
    - Logistic Regression (baseline)
    - Random Forest (main)
11. Evaluate performance:
    - Confusion matrix & accuracy
    - ROC curve & AUC
    - 5-fold cross-validation
12. Interpret results:
    - Random Forest feature importance
    - Map probe IDs → gene symbols

## 📊 Key Results
- Significant DEGs (adj p < 0.05, |logFC| > 1): 7,718
- Logistic Regression:
  - Test accuracy ~0.65 (baseline; can be unstable with high-dimensional correlated features)
- Random Forest:
  - Test accuracy = 1.00 on one split
  - AUC = 1.00 on one split
  - 5-fold cross-validation accuracy ~0.95–0.96 (more realistic estimate)

## 🧠 Interpretation Highlights
- Thousands of DEGs indicate a strong transcriptomic shift between GBM and non-tumor tissues.
- Random Forest captures a robust separation using the Top 100 DEG feature set.
- Feature importance highlights candidate genes that contribute most to classification and may serve as early biomarker candidates (requires external validation).

## 🧬 One-line Pipeline
GSE4290 -> Labeling -> Filtering -> Log2 -> Limma -> DEG -> Top100 -> DatasetML -> Split -> LogReg -> RandomForest -> ROC/AUC -> CV -> Importance -> GeneMapping -> Export

## 📂 Repository Structure
- `scripts/` : main analysis script  
- `results/` : output tables (DEG, Top100, feature importance)  
- `figures/` : volcano plot and feature importance plot  
- `report/`  : short report (PDF)

## ⚙️ Tools & Packages
- R / RStudio
- GEOquery
- limma
- ggplot2
- randomForest
- pROC
- caret
- hgu133plus2.db / AnnotationDbi

## ▶️ How to Run
1. Open R/RStudio
2. Run: `scripts/01_gse4290_deg_ml.R`
3. Outputs will be saved to `results/` and `figures/`

## 🖼️ Visualizations
![Volcano Plot](figures/volcano_plot.png)  
![Random Forest Feature Importance](figures/rf_feature_importance.png)

## 🎥 Video Presentation
Google Drive link: <https://drive.google.com/file/d/1qBxmjMEWlX65crPjVVKVZZbpyLUJBTfI/view?usp=sharing>

## 📚 References
1. Sun L, Hui AM, Su Q, et al. Neuronal and glioma-derived stem cell factor induces angiogenesis within the brain. *Cancer Cell*. 2006;9(4):287-300. doi:10.1016/j.ccr.2006.03.003.  
2. The Cancer Genome Atlas Research Network. Comprehensive genomic characterization defines human glioblastoma genes and core pathways. *Nature*. 2008;455(7216):1061-1068. doi:10.1038/nature07385.  
3. Verhaak RGW, Hoadley KA, Purdom E, et al. Integrated genomic analysis identifies clinically relevant subtypes of glioblastoma characterized by abnormalities in PDGFRA, IDH1, EGFR, and NF1. *Cancer Cell*. 2010;17(1):98-110. doi:10.1016/j.ccr.2009.12.020.
