# Imports
import pandas as pd
import os

from sklearn.preprocessing import RobustScaler,MinMaxScaler,StandardScaler
from sklearn.model_selection import train_test_split, StratifiedKFold, GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.utils import resample
from sklearn.metrics import confusion_matrix,classification_report,accuracy_score,roc_auc_score

import joblib 

# Load matrix
# Get the current directory of the script
current_directory = os.path.dirname(__file__)

# This is you execute it in Jupyter Notebook:
# current_directory = os.getcwd() 

# Navigate up one directory to mainfolder
mainfolder_directory = os.path.dirname(current_directory)

# Construct the path to the full_matrix.csv file in the data folder
full_matrix_path = os.path.join(mainfolder_directory, 'data', 'full_matrix.csv')

# Read the CSV file into a DataFrame
full_matrix = pd.read_csv(full_matrix_path)

# Select all features for the prediction
X = full_matrix.iloc[:, 1:]

# Select the first column (binding_site), which will be compared with the result of the prediction
Y = full_matrix.iloc[:, 0]

# Separate majority and minority classes
X_majority = X[Y==0]
X_minority = X[Y==1]

# Upsample minority class
X_minority_upsampled = resample(X_minority,
                                replace=True, # Sample with replacement
                                n_samples=len(X_majority), # Match majority class
                                random_state=42) # Reproducible results

# Combine majority class with upsampled minority class
X_upsampled = pd.concat([X_majority, X_minority_upsampled])

# Create new target variable
Y_upsampled = pd.Series([0]*len(X_majority) + [1]*len(X_minority_upsampled))

# Step 1: Split the upsampled dataset into training set and test set
X_train, X_test, Y_train, Y_test = train_test_split(X_upsampled, Y_upsampled, test_size=0.3, random_state=42)

# Step 2: Define the parameter grid
param_grid = {
    'n_estimators': [50, 100, 200],
    'max_depth': [None, 10, 20],
    'min_samples_split': [2, 5, 10],
    'min_samples_leaf': [1, 2, 4],
    'class_weight': [{0: 1, 1: 1000}]
}

# Step 3: Create a RandomForestClassifier
rfc = RandomForestClassifier()

# Step 4: Create a GridSearchCV object
grid_search = GridSearchCV(rfc, param_grid, cv=5, scoring='accuracy', n_jobs=-1, verbose=2)

# Step 5: Fit the GridSearchCV object on the training data
grid_search.fit(X_train, Y_train)

# Step 6: Extract and print the best parameters
best_params = grid_search.best_params_
print("Best Hyperparameters:", best_params)

# Step 7: Evaluate the best model on the test set
best_model = grid_search.best_estimator_
Y_pred = best_model.predict(X_test)

# Print classification report and accuracy
print("Classification Report:")
print(classification_report(Y_test, Y_pred))
print("Accuracy:", accuracy_score(Y_test, Y_pred))

# Save model as a file (uncomment to reproduce it)
# best_model_filename = 'BSmodel.pkl'
# joblib.dump(best_model, best_model_filename)