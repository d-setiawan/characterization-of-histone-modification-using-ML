# %%
# Import Libraries and Data
import numpy as np
import pandas as pd
import xgboost as xgb
import optuna
from shap_analysis import SHAPAnalyzer, function_map
from sklearn.model_selection import train_test_split
import shap

df = pd.read_csv("cleaned_transformed_histone_dataset_categorical.tsv", sep="\t")

X = df.drop(columns = ['Gene Expression (FPKM)_log', 'H3K9me2'])
y = df['Gene Expression (FPKM)_log']

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# %%
# Define and run Optuna trial

# Optuna trial defined
def objective(trial):
    # Define the hyperparameters to tune
    params = {
        'booster': 'gbtree',
        'objective': 'reg:squarederror',
        'eval_metric': 'rmse',
        'tree_method': 'hist',
        'device': 'cuda',
        'learning_rate': trial.suggest_float("learning_rate", 0.001, 0.05),
        'max_depth': trial.suggest_int("max_depth", 10, 30),
        'subsample': trial.suggest_float("subsample", 0.4, 1.0),
        'gamma': trial.suggest_float("gamma", 0.0, 1.0),
        'colsample_bytree': trial.suggest_float("colsample_bytree", 0.5, 1.0),
        'lambda': trial.suggest_float("lambda", 1e-3, 15),
        'alpha': trial.suggest_float("alpha", 1e-3, 15),
    }

    # Prepare the DMatrix
    dtrain = xgb.DMatrix(X_train, label=y_train)
    dvalid = xgb.DMatrix(X_test, label=y_test)

    # Train the model
    model = xgb.train(params, dtrain, num_boost_round=100, 
                      evals=[(dvalid, 'validation')],
                      early_stopping_rounds=10, 
                      verbose_eval=False)

    # Predict on validation data
    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X_test)

    # Analyze the SHAP values
    analyzer = SHAPAnalyzer(X_test, shap_values, function_map)

    analyzer.calculate_high_value_shap_means()

    results_df, total_score, result_summary = analyzer.get_results()

    total_mismatches = result_summary.get("Mismatch", 0)
    total_matches = result_summary.get("Match", 0)
    accuracy = total_matches / (total_matches + total_mismatches) if (total_matches + total_mismatches) > 0 else 0

    # Define the objective score as a combination of accuracy and score of high value SHAP means
    alpha = 1.0  # Weight for accuracy
    beta = 1.0   # Weight for mismatch count
    objective_score = alpha * accuracy + beta * total_score
    
    trial.set_user_attr("result_summary", result_summary)
    trial.set_user_attr("results_df", results_df)

    print(f"\n--- Trial {trial.number} ---")
    print(f"Params: {params}")
    print(f"Accuracy: {accuracy:.4f}")
    print(f"Total Matches: {total_matches}, Total Mismatches: {total_mismatches}")
    print(f"Objective Score: {objective_score:.4f}")

    return objective_score

# Create a study object and optimize
study = optuna.create_study(direction='maximize')
study.optimize(objective, n_trials=100)

# Output the best parameters
print("Best parameters:", study.best_params)
print("Best RMSE:", study.best_value)

# %%
#Analyze the output of the best trial
best_trial = study.best_trial

results_df = best_trial.user_attrs["results_df"]
results_df
