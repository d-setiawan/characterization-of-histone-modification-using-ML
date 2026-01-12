# %%
#Import Libraries and Data

import lightgbm as lgb
import pandas as pd
import numpy as np
df = pd.read_csv("cleaned_transformed_histone_dataset_categorical.tsv", sep="\t")

from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(
    df.drop(columns=["Gene Expression (FPKM)_log"]),
    df["Gene Expression (FPKM)_log"],
    test_size=0.2,
    random_state=42
)

# %%
#Set LightGBM parameters and train model
params = {
    'objective': 'regression',         # Regression task
    'metric': 'rmse',                  # Root Mean Squared Error
    'boosting_type': 'gbdt',           # Gradient Boosting Decision Tree          
    'num_leaves': 31,                  # Max leaves in a tree
    'learning_rate': 0.05,             # Learning rate
    'feature_fraction': 0.8,           # Fraction of features to use in each iteration
    'verbose': -1                      # Suppress output
}

train_data = lgb.Dataset(X_train, label=y_train)
test_data = lgb.Dataset(X_test, label=y_test, reference=train_data)

model = lgb.train(
    params,
    train_data,
    num_boost_round=500,
    valid_sets=[train_data, test_data],
    valid_names=['train', 'test']
)


# %%
from sklearn.metrics import root_mean_squared_error
# Make predictions
y_pred = model.predict(X_test)

# Calculate RMSE
rmse = root_mean_squared_error(y_test, y_pred)
print(f"Test RMSE: {rmse:.4f}")

# %%
#Evaluate Feature Importance through built in model functions
import matplotlib.pyplot as plt

# Plot feature importance
lgb.plot_importance(model, max_num_features=20, importance_type='gain', figsize=(10, 6))
plt.title("Feature Importance")
plt.show()

# %%
#Analyze Shap Values and results

import shap

# Calculate SHAP values
explainer = shap.Explainer(model)
shap_values = explainer.shap_values(X_test)

shap.summary_plot(shap_values, X_test)

import shap_analysis

function_map = {
    "Features": ['H3K4me1', 'H3K9me2', 'H3K4me3', 'H3K36me3',
       'H4K5Ac', 'H2A.Z.11', 'H3K27me3_log',
       'H3K9Ac_log', 'cpg_percentage_log', 'chh_percentage_log',
       'H4K20me1_log', 'chg_percentage_log', 'H2A.W.7_log', 'H3K9me1_log',
       'H2A.W.6_log', 'H3Ac_log', 'H3K9K14Ac_log'],
    "Known Function": ['Activating', 'Repressive', 'Activating', 'Activating', 'Activating', 'Activating',
                       'Repressive', 'Activating', 'Repressive','Repressive', 'Repressive', 'Repressive', 'Repressive',
                       'Repressive', 'Repressive', 'Activating', 'Activating']}

analyzer = shap_analysis.SHAPAnalyzer(X_test, shap_values, function_map)

analyzer.calculate_high_value_shap_means()

results_df, total_score, result_summary = analyzer.get_results()

print(results_df)
print(f"Total Score: {total_score}")
print(f"Result Summary: {result_summary}")


