# %%
# Import Libraries
import xgboost as xgb
import shap
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from shap_analysis import SHAPAnalyzer, function_map
import optuna

# %%
# Load Data
df = pd.read_csv("cleaned_transformed_histone_dataset_categorical.tsv", sep="\t")

X = df.drop(columns = ['Gene Expression (FPKM)_log'])
y = df['Gene Expression (FPKM)_log']

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)


# %%
#Make a function to evaluate the model based on the optuna study
def evaluate_model(X_train, X_test, y_train, y_test, params=None, optuna_study=None):
    """
    Evaluate the model either using the provided parameters or the Optuna study's best trial.
    
    Parameters:
    - X_train (pd.DataFrame): Training feature data.
    - X_test (pd.DataFrame): Test feature data.
    - y_train (pd.Series): Training target data.
    - y_test (pd.Series): Test target data.
    - params (dict): Parameters for the XGBoost model. If None, use `optuna_study`.
    - optuna_study (optuna.Study): Optuna study object containing the best trial.
    
    Returns:
    - model: Trained XGBoost model.
    - total_score (float): The evaluation score based on SHAP analysis.
    - result_summary (dict): Summary of matches and mismatches.
    - results_df (pd.DataFrame): Detailed SHAP results for features.
    """
    # Determine parameters: Use the best trial from Optuna if no parameters are provided
    if optuna_study is not None and params is None:
        params = optuna_study.best_trial.params
        params.update({
            'booster': 'gbtree',
            'objective': 'reg:squarederror',
            'eval_metric': 'rmse',
            'tree_method': 'hist',
            'device': 'cuda'
        })

    # Prepare the data
    dtrain = xgb.DMatrix(X_train, label=y_train)
    dtest = xgb.DMatrix(X_test, label=y_test)

    # Train the model
    model = xgb.train(
        params, 
        dtrain, 
        num_boost_round=100, 
        evals=[(dtest, 'validation')],
        early_stopping_rounds=10, 
        verbose_eval=False
    )

    # SHAP analysis
    background_indices = np.random.choice(X_train.shape[0], size=200, replace=False)
    background_sample = X_train.iloc[background_indices]

    sample_indices = np.random.choice(X_test.shape[0], size=200, replace=False)
    sample = X_test.iloc[sample_indices]

    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X_test)

    # Use the SHAPAnalyzer to calculate metrics
    analyzer = SHAPAnalyzer(X_test, shap_values, function_map)
    analyzer.calculate_high_value_shap_means()

    results_df, total_score, result_summary = analyzer.get_results()

    print(f"Results Summary: {result_summary}")

    return model, total_score, result_summary, results_df

# %%
#Make a function to perform feature selection with Optuna involved
def feature_selection_with_optuna(X_train, X_test, y_train, y_test, function_map):
    """
    Perform feature selection using Optuna to maximize the accuracy of the model.

    Parameters:
    - X_train (pd.DataFrame): Training feature data.
    - X_test (pd.DataFrame): Test feature data.
    - y_train (pd.Series): Training target data.
    - y_test (pd.Series): Test target data.
    - function_map (dict): Map of feature names to their corresponding functions for SHAP analysis.

    Returns:
    - selected_features (list): List of selected features.
    - study (optuna.Study): Optuna study object containing the best trial.
    """
    remaining_features = X_train.columns.tolist()
    selected_features = []

    # Define the Objective (Make Optuna study object)
    def optuna_objective(trial, X_train, X_test, y_train, y_test):
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

        dtrain = xgb.DMatrix(X_train, label=y_train)
        dvalid = xgb.DMatrix(X_test, label=y_test)

        model = xgb.train(params, dtrain, num_boost_round=100,
                          evals=[(dvalid, 'validation')],
                          early_stopping_rounds=10,
                          verbose_eval=False)

        explainer = shap.TreeExplainer(model)
        shap_values = explainer.shap_values(X_test)

        analyzer = SHAPAnalyzer(X_test, shap_values, function_map)
        analyzer.calculate_high_value_shap_means()

        results_df, total_score, result_summary = analyzer.get_results()
        total_mismatches = result_summary.get("Mismatch", 0)
        total_matches = result_summary.get("Match", 0)
        accuracy = total_matches / (total_matches + total_mismatches) if (total_matches + total_mismatches) > 0 else 0

        # Define the objective score
        alpha = 1.0  # Weight for accuracy
        beta = 1.0   # Weight for mismatch count
        objective_score = alpha * accuracy
        
        trial.set_user_attr("result_summary", result_summary)
        trial.set_user_attr("results_df", results_df)
        
        return objective_score

    # Perform feature selection after using the best optuna trial
    while remaining_features:
        print(f"Remaining Features: {remaining_features}")
        
        # Run Optuna for the current subset of features
        study = optuna.create_study(direction='maximize')
        study.optimize(lambda trial: optuna_objective(trial, 
                                                      X_train[remaining_features], 
                                                      X_test[remaining_features], 
                                                      y_train, y_test), 
                       n_trials=40)
        
        best_params = study.best_params
        print(f"Best Parameters: {best_params}")

        # Evaluate model with the current subset
        reduced_X_train = X_train[remaining_features]
        reduced_X_test = X_test[remaining_features]

        model, total_score, result_summary, results_df = evaluate_model(
            reduced_X_train, reduced_X_test, y_train, y_test, optuna_study=study
)

        # Identify mismatched features
        mismatch_features = results_df[results_df["Result"] == "Mismatch"]

        #If features all match, break the loop
        if mismatch_features.empty:
            selected_features = remaining_features
            break

        # Sort mismatching features and remove the highest impact wrong feature
        mismatch_to_remove = mismatch_features.sort_values(
            by="Mean SHAP Value (High)", ascending=False
        ).iloc[0]["Features"]

        print(f"Removing mismatch feature: {mismatch_to_remove}")
        remaining_features.remove(mismatch_to_remove)

    # Return the selected features and the best optuna study
    return selected_features, study

# %%
# Perform Leave one out cross validation with features using feature selection and optuna
optuna.logging.set_verbosity(optuna.logging.WARNING)

total_matches = 0
total_mismatches = 0

features = function_map["Features"]
known_functions = function_map["Known Function"]

feature_scores = {}

# Loops through all features to perform LOO
for i, feature_to_leave_out in enumerate(features):
    print(f"\n--- Performing LOO for feature {feature_to_leave_out} ---\n")

    # Create reduced X_train and X_test for optuna optimization
    reduced_X_train = X_train.drop(columns = [feature_to_leave_out])
    reduced_X_test = X_test.drop(columns = [feature_to_leave_out])

    # Perform feature selection with optuna
    selected_features, study = feature_selection_with_optuna(reduced_X_train, reduced_X_test, y_train, y_test, function_map)

    # Create final X_train and X_test by adding back the left out feature
    final_X_train = X_train[selected_features + [feature_to_leave_out]]
    final_X_test = X_test[selected_features + [feature_to_leave_out]]

    # Evaluate model with the left out feature added back
    print(f"\n--- Evaluating model with feature {feature_to_leave_out} added back ---\n")

    model, total_score, result_summary, results_df = evaluate_model(
        final_X_train, final_X_test, y_train, y_test, optuna_study=study
        )

    # Get the known function for the left out feature
    known_function = known_functions[i]

    # Get the result for the left out feature
    feature_result = results_df[results_df["Features"] == feature_to_leave_out]

    # If the feature result is not empty, get the match result and update the counters
    if not feature_result.empty:
        match_result = feature_result["Result"].values[0]
        print(f"Feature: {feature_to_leave_out}, Known Function: {known_function}, Result: {match_result}")

        if match_result == "Match":
            total_matches += 1
        else:
            total_mismatches += 1

        # Store the result in the feature scores dictionary
        feature_scores[feature_to_leave_out] = {
            "Total Score": total_score,
            "Known Function": known_function,
            "Result": match_result
        }

# Calculate the total number of features and accuracy
total_features = total_matches + total_mismatches
accuracy = total_matches / total_features if total_features > 0 else 0

# Print the feature analysis results
print("\n--- Feature Analysis Results ---")

for feature, score in feature_scores.items():
    print(f"Feature: {feature}, Total Score: {score['Total Score']}, Known Function: {score['Known Function']}, Result: {score['Result']}")

# Print the overall accuracy
print("\n--- Overall Accuracy ---")
print(f"Total Matches: {total_matches}")
print(f"Total Mismatches: {total_mismatches}")
print(f"Accuracy: {accuracy:.2f}")

# %%
#Perform Leave one out cross validation with optuna but without feature selection
optuna.logging.set_verbosity(optuna.logging.INFO)

total_matches = 0
total_mismatches = 0

features = function_map["Features"]
known_functions = function_map["Known Function"]

feature_scores = {}

for i, feature_to_leave_out in enumerate(features):
    print(f"\n--- Performing LOO for feature {feature_to_leave_out} ---\n")

    #Features for optuna optimization
    reduced_X_train = X_train.drop(columns = [feature_to_leave_out])
    reduced_X_test = X_test.drop(columns = [feature_to_leave_out])

    #Optuna Optimization
    def optuna_objective(trial, reduced_X_train, reduced_X_test, y_train, y_test):
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

        dtrain = xgb.DMatrix(reduced_X_train, label=y_train)
        dvalid = xgb.DMatrix(reduced_X_test, label=y_test)

        model = xgb.train(params, dtrain, num_boost_round=100,
                          evals=[(dvalid, 'validation')],
                          early_stopping_rounds=10,
                          verbose_eval=False)

        explainer = shap.TreeExplainer(model)
        shap_values = explainer.shap_values(reduced_X_test)

        analyzer = SHAPAnalyzer(reduced_X_test, shap_values, function_map)
        analyzer.calculate_high_value_shap_means()

        results_df, total_score, result_summary = analyzer.get_results()
        total_mismatches = result_summary.get("Mismatch", 0)
        total_matches = result_summary.get("Match", 0)
        accuracy = total_matches / (total_matches + total_mismatches) if (total_matches + total_mismatches) > 0 else 0

        # Define the objective score
        alpha = 1.0  # Weight for accuracy
        beta = 1.0   # Weight for mismatch count
        objective_score = alpha * accuracy
        
        trial.set_user_attr("result_summary", result_summary)
        trial.set_user_attr("results_df", results_df)
        
        return objective_score
    
    # Create a study object and optimize
    study = optuna.create_study(direction='maximize')
    study.optimize(lambda trial: optuna_objective(trial, 
                                                    reduced_X_train, 
                                                    reduced_X_test, 
                                                    y_train, y_test), 
                    n_trials=40)

    # Evaluate the model on the feature that was left out
    print(f"\n--- Evaluating model with feature {feature_to_leave_out} added back ---\n")

    model, total_score, result_summary, results_df = evaluate_model(
        X_train, X_test, y_train, y_test, optuna_study=study
        )

    # Get the known function for the left out feature
    known_function = known_functions[i]
    feature_result = results_df[results_df["Features"] == feature_to_leave_out]

    # If the feature result is not empty, get the match result (whether feature was correctly classified or not) and update the counters
    if not feature_result.empty:
        match_result = feature_result["Result"].values[0]
        print(f"Feature: {feature_to_leave_out}, Known Function: {known_function}, Result: {match_result}")

        if match_result == "Match":
            total_matches += 1
        else:
            total_mismatches += 1

        feature_scores[feature_to_leave_out] = {
            "Total Score": total_score,
            "Known Function": known_function,
            "Result": match_result
        }

# Print metrics
total_features = total_matches + total_mismatches
accuracy = total_matches / total_features if total_features > 0 else 0

print("\n--- Feature Analysis Results ---")

for feature, score in feature_scores.items():
    print(f"Feature: {feature}, Total Score: {score['Total Score']}, Known Function: {score['Known Function']}, Result: {score['Result']}")

print("\n--- Overall Accuracy ---")
print(f"Total Matches: {total_matches}")
print(f"Total Mismatches: {total_mismatches}")
print(f"Accuracy: {accuracy:.2f}")


