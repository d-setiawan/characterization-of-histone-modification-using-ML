# %%
#Import all necessary libraries
import os
import pandas as pd
import numpy as np

print(os.getcwd())
# %%
#Set the folder path
folder_path = "/home/dgsetiawan/MachineLearning/ZhongLab/HistoneModificationInterpretation/Data/Matrices/"

#Create a list to store the dataframes
dataframes = []

#Loop through all files in folder
for file_name in os.listdir(folder_path):
    #Check if the file ends with ".mat"
    if file_name.endswith(".mat"):
        #Get the full file path
        file_path = os.path.join(folder_path, file_name)
        print(f"Processing file:{file_name}")

        #Load data
        df = pd.read_csv(file_path, sep="\t", skiprows=1, header=None)

        #Select only the ID and value columns
        filtered_df = df[[3,6]]

        #Extract ID from column
        filtered_df[3] = filtered_df[3].str.extract(r'ID=([^;]+)')

        #Rename the columns
        filtered_df.columns = ['ID', os.path.splitext(os.path.basename(file_name))[0]]

        #Add the dataframe to the list
        dataframes.append(filtered_df)

#Merge all the dataframes using ID's as key, only merging when gene ID's exist
merged_df = dataframes[0]
for df in dataframes[1:]:
    merged_df = pd.merge(merged_df, df, on='ID', how='outer')

# %%
#Load the target data (RNA seq measured in FPKM) and rename columns
target_df = pd.read_csv("/home/dgsetiawan/MachineLearning/ZhongLab/HistoneModificationInterpretation/Data/Matrices/SRR7405039.mat.genes.results", sep="\t")
target_df = target_df[['gene_id', 'FPKM']]
target_df.columns = ['ID', 'Gene Expression (FPKM)']

# %%
#Merge the histone modification data with the target data (RNA seq measured in FPKM)
#Inner join to only include genes with measured expression
merged_df = pd.merge(merged_df, target_df, on='ID', how='inner')

#Save the merged dataframe to disk as a tab-separated file
merged_df.to_csv("Histone_Modification_Characterization.tsv", index=False, sep="\t")
