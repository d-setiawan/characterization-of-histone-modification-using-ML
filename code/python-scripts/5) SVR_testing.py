# %%
#Import libraries and data
import cudf as cd
import cupy as cp

df = cd.read_csv("cleaned_transformed_histone_dataset_categorical.tsv", sep="\t")

from cuml.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(df.drop(columns = ['Gene Expression (FPKM)_log']), df['Gene Expression (FPKM)_log'], test_size=0.2, random_state=42)

# %%
#Scale data since SVR is sensitive to scale
from cuml.preprocessing import StandardScaler

scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

X_train_scaled.columns = X_train.columns
X_test_scaled.columns = X_test.columns

# %%
#Train and evaluate model
from cuml.svm import SVR, LinearSVR

model = SVR(kernel = 'linear')
model.fit(X_train_scaled, y_train)

y_pred = model.predict(X_test_scaled)

from cuml.metrics import r2_score, mean_squared_error

r2 = r2_score(y_test, y_pred)
rmse = mean_squared_error(y_test, y_pred)

print(f"R-squared: {r2}")
print(f"Root Mean Squared Error: {rmse}")

# %%
#Calculate SHAP values
import shap
import numpy as np

#Sampling of data is needed because shap value calculation is computationally expensive
background_indices = np.random.choice(X_train_scaled.shape[0], size=100, replace=False)
background_sample = X_train_scaled.iloc[background_indices]
background_sample_pd = background_sample.to_pandas()

sample_indices = np.random.choice(X_test_scaled.shape[0], size=100, replace=False)
sample = X_test_scaled.iloc[sample_indices]
sample_pd = sample.to_pandas()

#Calculate Shap Values and Plot
explainer = shap.Explainer(model.predict, background_sample_pd)
shap_values = explainer.shap_values(sample_pd)

shap.summary_plot(shap_values, sample_pd)

# %%
#Analyze Shap Results
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

analyzer = shap_analysis.SHAPAnalyzer(sample.to_pandas(), shap_values, function_map)

analyzer.calculate_high_value_shap_means()

results_df, total_score, result_summary = analyzer.get_results()

print(results_df)
print(f"Total Score: {total_score}")
print(f"Result Summary: {result_summary}")
