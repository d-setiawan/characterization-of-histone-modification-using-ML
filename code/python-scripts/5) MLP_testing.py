# %%
#Import Libraries
import torch
import torch.nn as nn
import torch.optim as optim

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import root_mean_squared_error

import numpy as np
import pandas as pd

# Check if a GPU is available
if torch.cuda.is_available():
    print("GPU is available!")
    print(f"CUDA version: {torch.version.cuda}")
    print(f"GPU Name: {torch.cuda.get_device_name(0)}")
else:
    print("GPU is not available. Running on CPU.")


# %%
#Set device to gpu
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
device

# %%
#Load data, Scale data becauase nn's are sensitive to scale
df = pd.read_csv('cleaned_transformed_histone_dataset_categorical.tsv', sep='\t')

X_train, X_test, y_train, y_test = train_test_split(df.drop(columns=['Gene Expression (FPKM)_log']), df['Gene Expression (FPKM)_log'], test_size=0.2, random_state=42)

scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

X_train_tensor = torch.tensor(X_train_scaled, dtype=torch.float32)
y_train_tensor = torch.tensor(y_train.values, dtype=torch.float32)
X_test_tensor = torch.tensor(X_test_scaled, dtype=torch.float32)
y_test_tensor = torch.tensor(y_test.values, dtype=torch.float32)

#Move to GPU
X_train_tensor = X_train_tensor.to(device)
y_train_tensor = y_train_tensor.to(device)
X_test_tensor = X_test_tensor.to(device)
y_test_tensor = y_test_tensor.to(device)


# %%
#Make the model, after testing, a wide and deep model worked best
class MLPRegressor(nn.Module):
    def __init__(self, input_dim):
        super(MLPRegressor, self).__init__()
        self.network = nn.Sequential(
            nn.Linear(input_dim, 256),  # Wider input layer
            nn.Tanh(),
            nn.Linear(256, 256),  # Wider hidden layers
            nn.Tanh(),
            nn.Linear(256, 128),  # Transition to slightly narrower layers
            nn.Tanh(),
            nn.Linear(128, 128),
            nn.Tanh(),
            nn.Linear(128, 64),
            nn.Tanh(),
            nn.Linear(64, 64),
            nn.Tanh(),
            nn.Linear(64, 32),  # Gradual narrowing
            nn.Tanh(),
            nn.Linear(32, 32),
            nn.Tanh(),
            nn.Linear(32, 1)  # Output layer
        )

    def forward(self, x):
        return self.network(x)

# %%
#Move the model to GPU and set the loss function and optimizer
input_dim = X_train.shape[1]
model = MLPRegressor(input_dim)

#Move to GPU
model = model.to(device)

criterion = nn.MSELoss()

optimizer = optim.Adam(model.parameters(), lr=0.0005)

print(next(model.parameters()).device)  # Output: cuda:0 if on GPU
print(X_train_tensor.device)  # Output: cuda:0 if on GPU


# %%
#Train the model

epochs = 200
for epoch in range(epochs):
    model.train()  # Set model to training mode
    optimizer.zero_grad()  # Clear gradients
    
    # Forward pass for training
    y_pred_train = model(X_train_tensor)
    train_loss = criterion(y_pred_train, y_train_tensor)  # Compute training loss

    # Backward pass
    train_loss.backward()  # Compute gradients
    optimizer.step()  # Update weights

    # Evaluate on test data
    model.eval()  # Set model to evaluation mode
    with torch.no_grad():  # Disable gradient computation for test data
        y_pred_test = model(X_test_tensor)
        test_loss = criterion(y_pred_test, y_test_tensor)  # Compute test loss

    # Print losses every epoch
    if epoch % 10 == 0:
        print(f"Epoch {epoch + 1}/{epochs}, Training Loss: {train_loss.item():.4f}, Test Loss: {test_loss.item():.4f}")

# %%
#Evaluate model

# Switch to evaluation mode
model.eval()

# Make predictions
with torch.no_grad():  # Disable gradient computation for evaluation
    y_pred = model(X_test_tensor)

y_pred_numpy = y_pred.cpu().numpy()

# Compute RMSE
rmse = root_mean_squared_error(y_test, y_pred_numpy)
print(f"Test RMSE: {rmse:.4f}")


# %%
#Calculate SHAP values and make plots
import shap

#Sampling of data is needed because shap value calculation is computationally expensive
X_train_df = pd.DataFrame(X_train_tensor.cpu().numpy(), columns=X_train.columns)
X_test_df = pd.DataFrame(X_test_tensor.cpu().numpy(), columns=X_train.columns)

background_indices = np.random.choice(X_train.shape[0], size=100, replace=False)
background_sample_tensor = X_train_tensor[background_indices]
background_sample = X_train_df.iloc[background_indices]

sample_indices = np.random.choice(X_test.shape[0], size=100, replace=False)
sample_tensor = X_test_tensor[sample_indices]
sample = X_test_df.iloc[sample_indices]

#Calculate Shap Values
explainer = shap.DeepExplainer(model, background_sample_tensor)
shap_values = explainer.shap_values(sample_tensor)

shap.summary_plot(shap_values, sample)

# %%
#Analyze the Shap values
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

analyzer = shap_analysis.SHAPAnalyzer(background_sample, shap_values, function_map)

analyzer.calculate_high_value_shap_means()

results_df, total_score, result_summary = analyzer.get_results()

print(results_df)
print(f"Total Score: {total_score}")
print(f"Result Summary: {result_summary}")
