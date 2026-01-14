# Interpretable Machine Learning for Characterizing Histone Modifications

**Author:** Dylan Setiawan  
**Project Type:** Senior Thesis (Computational Biology / Machine Learning)

---

## Overview

This repository contains the full codebase, analysis, and results for my senior thesis project focused on **interpretable machine learning for characterizing the regulatory effects of histone modifications on gene expression**.

The goal of this project was to build a **reproducible, end-to-end bioinformatics and machine learning pipeline** that:

- Integrates **histone modification sequencing data** with RNA-seq gene expression  
- Trains and compares multiple ML models (tree-based, kernel-based, and neural networks)  
- Uses **SHAP (SHapley Additive exPlanations)** to interpret feature importance  
- Evaluates whether learned feature effects align with **known biological activating vs. repressive histone functions**  
- Generalizes to the prediction of **previously uncharacterized histone modifications**

The final pipeline achieves **~90% accuracy** in correctly identifying the functional direction (activating vs. repressive) of histone marks, demonstrating that interpretable ML can recover biologically meaningful regulatory signals from high-dimensional sequencing data.

---

## Key Contributions

- Designed a **modular ML + bioinformatics pipeline** from raw sequencing matrices to interpretable biological conclusions  
- Implemented **robust preprocessing** (log transforms, outlier detection, quantile filtering)  
- Compared multiple models:
  - XGBoost (GPU-accelerated)
  - LightGBM
  - Support Vector Regression
  - k-Nearest Neighbors (GPU)
  - Multilayer Perceptron (PyTorch)
- Developed a **custom SHAP analysis framework** to:
  - Quantify feature effects at high values  
  - Compare model explanations against known histone biology  
- Performed **feature validation and leave-one-out analysis**  
- Organized and cleaned a **bioinformatics shell-based sequencing pipeline**

---


## Pipeline Summary

1. **Data Integration**
   - Merge histone modification signal matrices with RNA-seq expression (FPKM)
   - Align features by gene ID

2. **Exploratory Data Analysis**
   - Distribution analysis
   - Log transformations
   - Outlier detection (Isolation Forest + quantile filtering)

3. **Model Training**
   - Train multiple ML models under consistent preprocessing
   - GPU acceleration where applicable

4. **Model Selection & Optimization**
   - Hyperparameter tuning with Optuna
   - Performance comparison across models

5. **Interpretability & Biological Validation**
   - SHAP value computation
   - High-value SHAP analysis
   - Agreement testing against known activating/repressive histone marks

---

## Technologies Used

- **Python** (pandas, NumPy, scikit-learn)
- **GPU-accelerated ML** (XGBoost, cuML, LightGBM)
- **Deep Learning** (PyTorch)
- **Explainable AI** (SHAP)
- **Hyperparameter Optimization** (Optuna)
- **Bioinformatics Pipelines** (Bash / shell scripting)
- **Visualization** (matplotlib, seaborn)

---

## Results Highlights

- Gradient-boosted tree models (XGBoost) provided the best balance of:
  - Predictive performance
  - Computational efficiency
  - Interpretability
- Learned SHAP feature effects aligned with known biology for the majority of histone marks
- Demonstrated feasibility of **predicting regulatory direction for unknown histone modifications**

---

## Thesis Report

The full written thesis, including background, methodology, results, and discussion, is included in this repository:
```
report/SetiawanDylanSpring2025HonorsThesis.pdf
```


