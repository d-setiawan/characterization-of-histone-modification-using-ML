# %%

#Import libraries and dataset
import cudf as cd
import cupy as cp
from cuml.model_selection import train_test_split
from cuml.neighbors import KNeighborsRegressor

df = cd.read_csv("cleaned_transformed_histone_dataset_categorical.tsv", sep="\t")

X = df.drop(columns=["Gene Expression (FPKM)_log", "H4K5Ac", "H2A.Z.11"])
y = df["Gene Expression (FPKM)_log"]

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

# %%
#Scale data due to nature of KNN
from cuml.preprocessing import StandardScaler

scaler = StandardScaler()

X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)
X_train_scaled.columns = X_train.columns
X_test_scaled.columns = X_test.columns

# %%
knn_model = KNeighborsRegressor(
    n_neighbors=8,
    # weights="distance",
    algorithm="ivfflat",
    # p=1,
    # batch_size=512
)
knn_model.fit(X_train_scaled, y_train)

y_pred = knn_model.predict(X_test_scaled)

from cuml.metrics import r2_score, mean_squared_error

r2 = r2_score(y_test, y_pred)
mse = mean_squared_error(y_test, y_pred)
print(f"R-squared: {r2}")
print(f"Root Mean Squared Error: {mse}")

# %%
from shap_analysis import SHAPAnalyzer, function_map
import shap
import numpy as np

#Train model and get shap values
def predict_knn(data):
    data_cudf = cd.DataFrame(data)
    return knn_model.predict(data_cudf).to_numpy()

X_train_pd = X_train_scaled.to_pandas()
X_test_pd = X_test_scaled.to_pandas()

background_indices = np.random.choice(X_train_pd.shape[0], size=100, replace=False)
background = X_train_pd.iloc[background_indices]

sample_indices = np.random.choice(X_test_pd.shape[0], size=100, replace=False)
sample = X_test_pd.iloc[sample_indices]


explainer = shap.KernelExplainer(predict_knn, background)
shap_values = explainer.shap_values(sample)


# %%
#Analyze Shap Values
analyzer = SHAPAnalyzer(sample, shap_values, function_map)

analyzer.calculate_high_value_shap_means()

results_df, total_score, result_summary = analyzer.get_results()

print(f"Total Score: {total_score}")
print(f"Results Summary: {result_summary}")


