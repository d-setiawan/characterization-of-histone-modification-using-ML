# Define the method to calculate mean SHAP values for high feature values
import pandas as pd
import numpy as np

# Mapping of features to known biological functions
#  - 'Activating' marks are associated with active gene expression
#  - 'Repressive' marks are associated with repressed gene expression
function_map = {
    "Features": ['H3K4me1', 'H3K9me2', 'H3K4me3', 'H3K36me3',
       'H4K5Ac', 'H3K27me3_log',
       'H3K9Ac_log', 'cpg_percentage_log', 'chh_percentage_log',
       'H4K20me1_log', 'chg_percentage_log', 'H2A.W.7_log', 'H3K9me1_log',
       'H2A.W.6_log', 'H3Ac_log', 'H3K9K14Ac_log'],
    "Known Function": ['Activating', 'Repressive', 'Activating', 'Activating', 'Activating',
                       'Repressive', 'Activating', 'Repressive','Repressive', 'Repressive', 'Repressive', 'Repressive',
                       'Repressive', 'Repressive', 'Activating', 'Activating']}

#Python class object to analyze SHAP values and create summary statistics
class SHAPAnalyzer:
    def __init__(self, X_test, shap_values, function_map, quantile=0.75):
        """
        Initialize the SHAPAnalyzer with the dataset, SHAP values, and configuration.

        Parameters:
        - X_test (pd.DataFrame): The feature matrix used for predictions.
        - shap_values (np.ndarray): The SHAP values corresponding to the features in X_test.
        - function_map (dict): A dictionary mapping feature names to their function (Activating/Repressive).
        - quantile (float): The quantile to define "high" values (default is 0.75, i.e., top 25%).
        """
        self.X_test = X_test
        self.shap_values = shap_values
        self.function_map = function_map
        self.quantile = quantile
        self.results_df = None
        self.total_score = None
        self.result_summary = None


    def calculate_high_value_shap_means(self):
        # Validate SHAP values and features
        """
        Calculate mean SHAP values for high feature values (top quantile) across all features.

        Parameters:
        - X_test (pd.DataFrame): The feature matrix used for predictions.
        - shap_values (np.ndarray): The SHAP values corresponding to the features in X_test.
        - quantile (float): The quantile to define "high" values (default is 0.75, i.e., top 25%).

        Returns:
        - pd.DataFrame: A DataFrame of mean SHAP values for high feature values, sorted by impact.
        """
        if self.shap_values.shape[1] != self.X_test.shape[1]:
            raise ValueError("Mismatch between SHAP values and input features.")
        # if set(self.X_test.columns) != set(self.function_map['Features']):
        #     raise ValueError("Mismatch between X_test features and function_map.")

        # Define the threshold for "high" values
        high_value_thresholds = self.X_test.quantile(self.quantile)

        # Create DataFrame for SHAP values
        try:
            shap_values_df = pd.DataFrame(self.shap_values, columns=self.X_test.columns)
            shap_values_df = shap_values_df.set_index(self.X_test.index)
        except Exception as e:
            raise ValueError("Error creating SHAP values DataFrame.") from e

        # Calculate mean SHAP values for high feature values
        high_value_shap_means = {}
        for feature in self.X_test.columns:
            high_group = shap_values_df.loc[self.X_test[feature] > high_value_thresholds[feature], feature]
            high_value_shap_means[feature] = high_group.mean()

        # Convert results to DataFrame
        results_df = pd.DataFrame.from_dict(
            high_value_shap_means, orient='index', columns=['Mean SHAP Value (High)']
        ).reset_index().rename(columns={'index': 'Features'})

        # Add Known Function column
        results_df = pd.merge(
            results_df, pd.DataFrame(self.function_map), how='inner', on='Features'
        )

        # Map functions and calculate error
        results_df['Function Value'] = results_df['Known Function'].map({'Activating': 1, 'Repressive': -1})
        results_df['Accuracy Score'] = results_df['Mean SHAP Value (High)'] * results_df['Function Value']

        # Compare SHAP results
        results_df['Result'] = results_df.apply(self.compare_shap_function, axis=1)
        self.results_df = results_df
        self.total_score = results_df['Accuracy Score'].sum()
        self.result_summary = results_df['Result'].value_counts().to_dict()

    @staticmethod
    def compare_shap_function(row):
        """Determine if the SHAP value aligns with the known function."""
        if (row['Mean SHAP Value (High)'] > 0 and row['Known Function'] == "Activating") or \
           (row['Mean SHAP Value (High)'] < 0 and row['Known Function'] == "Repressive"):
            return "Match"
        else:
            return "Mismatch"

    def get_results(self):
        """
        Get the detailed results DataFrame, total score, and match/mismatch summary.

        Returns:
        - pd.DataFrame: The detailed results DataFrame.
        - float: The total score (sum of all error values).
        - dict: The summary of matches and mismatches.
        """
        if self.results_df is None or self.total_score is None or self.result_summary is None:
            raise ValueError("You must call calculate_high_value_shap_means() before getting results.")
        return self.results_df, self.total_score, self.result_summary


