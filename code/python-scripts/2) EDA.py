# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# %%
#Load in the dataset from previous script and drop index columns and NA's
df = pd.read_csv('Histone_Modification_Characterization.tsv', sep='\t')
df_cleaned = df.drop(columns=['ID'])
df_cleaned = df_cleaned.dropna(axis=0, how='any')
df_cleaned.describe()

# %%
#Add feature column expression category to indicate whether gene is expressed or not
df_cleaned['expression_category'] = df_cleaned['Gene Expression (FPKM)'].apply(lambda x: 0 if x == 0 else 1)


# %%
#Plot histograms of all features below the .99 quantile
#This is done to ignore outliers

#Get the names of the features
features = df_cleaned.drop(columns=['expression_category']).columns
filtered_df = df_cleaned.copy()

#Create subplots for histograms
fig, axes = plt.subplots(nrows=8, ncols=3, figsize=(20, 30))
axes = axes.flatten()

#Loop through each feature
for i, feature in enumerate(features):
    #Get the upper limit for the feature (i.e. the value below the 99th percentile)
    upper_limit = df_cleaned[feature].quantile(0.99)
    
    #Plot the histogram
    sns.histplot(x=df_cleaned[feature][df[feature] < upper_limit], bins=50, kde = True, ax=axes[i])
    
    #Add labels and title
    axes[i].set_title(feature)
    axes[i].set_xlabel('Value')
    axes[i].set_ylabel('Frequency')

#Adjust layout and show the plot
plt.tight_layout(pad=3.0) 
plt.show()


# %%
#Determine which features should be log transformed based on above histograms
columns_to_log_transform = ['H3K27me3', 'H3K9Ac', 'cpg_percentage', 'H3K9K14Ac', 
                            'H3Ac', 'chh_percentage', 'H4K20me1', 'chg_percentage',
                            'H2A.W.7', 'H3K9me1', 'H2A.W.6', 'H2A.Z.9', 'Gene Expression (FPKM)']

#Log transform the features
df_transformed = df_cleaned.copy()
for feature in columns_to_log_transform:
    df_transformed[f'{feature}_log'] = np.log(df_transformed[feature]+0.0001)
    df_transformed.drop(columns=[feature], inplace=True)

#Plot the features again like above but now also ignoring the bottom 1%
features = df_transformed.columns

fig, axes = plt.subplots(nrows=8, ncols=3, figsize=(20, 30))
axes = axes.flatten()

for i, feature in enumerate(features):
    upper_limit = df_transformed[feature].quantile(0.99)
    lower_limit = df_transformed[feature].quantile(0.01)
    sns.histplot(x=df_transformed[feature][(df_transformed[feature] < upper_limit) & (df_transformed[feature] > lower_limit)], bins=50, kde = True, ax=axes[i])
    axes[i].set_title(feature)
    axes[i].set_xlabel('Value')
    axes[i].set_ylabel('Frequency')

plt.tight_layout(pad=3.0) 
plt.show()

# %%
#Implement Isolation Forest, an outlier detection algorithm that is based on tree clustering

from sklearn.ensemble import IsolationForest

#Fit data and predict outliers
iso_forest = IsolationForest(contamination=0.05, random_state=42)
iso_forest.fit(df_transformed)

df_transformed['is outlier'] = iso_forest.predict(df_transformed)


#Plot outliers on histograms
fig, axes = plt.subplots(nrows=8, ncols=3, figsize=(20, 30))
axes = axes.flatten()

for i, feature in enumerate(df_transformed.columns[:-1]):
    sns.histplot(data = df_transformed, x = feature, 
                 hue='is outlier', bins = 50, kde = True, palette = {1: 'blue', -1: 'red'}, 
                 ax=axes[i])
    axes[i].set_title(f'Distribution of {feature}') 
    axes[i].set_xlabel(feature)
    axes[i].set_ylabel('Frequency')


plt.tight_layout()
plt.show()

# %%
#Drop outliers column from Isolation forest
df_transformed_cleaned = df_transformed[df_transformed['is outlier'] == 1].drop(columns=['is outlier'])
# %%
#Replot Histogram after outlier removal
features = df_transformed_cleaned.drop(columns=['expression_category']).columns

fig, axes = plt.subplots(nrows=7, ncols=3, figsize=(20, 30))
axes = axes.flatten()

for i, feature in enumerate(features):
    upper_limit = df_transformed_cleaned[feature].quantile(0.999)
    lower_limit = df_transformed_cleaned[feature].quantile(0.001)
    sns.histplot(x=df_transformed_cleaned[feature][(df_transformed[feature] < upper_limit) & (df_transformed[feature] > lower_limit)], bins=50, kde = True, ax=axes[i])
    axes[i].set_title(feature)
    axes[i].set_xlabel('Value')
    axes[i].set_ylabel('Frequency')

plt.tight_layout(pad=3.0) 
plt.show()

# %%
#Apply quantile mask to further remove outliers based on visual inspection
lower_quantile = 0.001
upper_quantile = 0.999

# Initialize a mask with True values for all indices
mask = pd.Series(True, index=df_transformed_cleaned.index)

# Iterate over each feature to apply quantile limits
for feature in df_transformed_cleaned.columns:
    # Calculate lower and upper quantile limits
    lower_limit = df_transformed_cleaned[feature].quantile(lower_quantile)
    upper_limit = df_transformed_cleaned[feature].quantile(upper_quantile)
    # Update the mask to retain values within the quantile limits
    mask &= (df_transformed_cleaned[feature] <= upper_limit) & (df_transformed_cleaned[feature] >= lower_limit)

# Apply the mask to filter the DataFrame based on quantile limits
df_transformed_cleaned_quantile = df_transformed_cleaned[mask]

# %%
# Save the cleaned and transformed dataset to a TSV file
df_transformed_cleaned_quantile.to_csv('cleaned_transformed_histone_dataset_categorical.tsv', sep='\t', index = False)
