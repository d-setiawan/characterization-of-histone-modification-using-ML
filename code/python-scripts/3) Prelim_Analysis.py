# %%
import numpy as np
import pandas as pd
import xgboost as xgb
import cuml

# %%
# Load the cleaned and transformed dataset
df = pd.read_csv("cleaned_transformed_histone_dataset_categorical.tsv", sep="\t")
df = df.dropna(axis=0, how='any')


# %%
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score

X = df.drop(columns=["Gene Expression (FPKM)_log"])
y = df["Gene Expression (FPKM)_log"]

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

# %%
#Train XGBoost model

# Convert to DMatrix, XGBoost's native GPU supported Data Format
dtrain = xgb.DMatrix(X_train, label=y_train)
dtest = xgb.DMatrix(X_test, label=y_test)

#Set basic parameters for XGBoost
params = {
    'booster': 'gbtree',
    'device': 'cuda',
    'objective': 'reg:squarederror',
    'eval_metric': 'rmse',
    'tree_method': 'gpu_hist',
    'verbosity': 1,
}

# Provide watchlist for validation tracking and train model
watchlist = [(dtrain, 'train'), (dtest, 'test')]
bst = xgb.train(params, dtrain, num_boost_round=100, evals= watchlist, early_stopping_rounds=50)

y_pred = bst.predict(dtest)

#Calculate r2_score as metric
r2 = r2_score(y_test, y_pred)
print(r2)

# %%
import shap
#Use Shap to explain model throgh feature importance and effects

explainer = shap.TreeExplainer(bst)
shap_values = explainer.shap_values(X_test)
shap.summary_plot(shap_values, X_test)

# %%
# Define the method to calculate mean SHAP values for high feature values
def calculate_high_value_shap_means(X_test, shap_values, quantile=0.75):
    """
    Calculate the mean SHAP values for high feature values (top quantile) across all features.

    Parameters:
    - X_test (pd.DataFrame): The feature matrix used for predictions.
    - shap_values (np.ndarray): The SHAP values corresponding to the features in X_test.
    - quantile (float): The quantile to define "high" values (default is 0.75, i.e., top 25%).

    Returns:
    - pd.DataFrame: A DataFrame of mean SHAP values for high feature values, sorted by impact.
    """
    # Define the threshold for "high" values
    high_value_thresholds = X_test.quantile(quantile)

    # Ensure indices of SHAP values and X_test match
    shap_values_df = pd.DataFrame(shap_values, columns=X_test.columns)
    shap_values_df = shap_values_df.set_index(X_test.index)

    # Initialize a dictionary to store mean SHAP values for high feature values
    high_value_shap_means = {}

    # Loop through all features
    for feature in X_test.columns:
        # Filter rows where the feature value is in the top quantile (high values)
        high_group = shap_values_df.loc[X_test[feature] > high_value_thresholds[feature], feature]
        
        # Calculate the mean SHAP value for high feature values
        high_value_shap_means[feature] = high_group.mean()

    # Convert results to a DataFrame for easier interpretation
    high_value_shap_means_df = pd.DataFrame.from_dict(
        high_value_shap_means, orient='index', columns=['Mean SHAP Value (High)']
    )
    high_value_shap_means_df = high_value_shap_means_df.sort_values(by='Mean SHAP Value (High)', ascending=False)

    return high_value_shap_means_df

# Example usage
high_value_shap_means_df = calculate_high_value_shap_means(X_test, shap_values)
print("XGBoost Prelim Analysis Results")
print(high_value_shap_means_df)



