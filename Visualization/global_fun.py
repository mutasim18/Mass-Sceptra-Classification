import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report, accuracy_score, confusion_matrix, roc_auc_score
import xgboost as xgb
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import seaborn as sns
import os
import joblib
import warnings
warnings.filterwarnings('ignore')

# Function to read mass spectra data from CSV
def load_spectrum(filepath):
    """Reads CSV file containing mass spectrum data
    Returns: m/z values and intensity values"""
    try:
        spectrum = pd.read_csv(filepath, skiprows=1)  # Skip header row
        mz = spectrum['X(Thompsons)'].values         # Mass/charge values
        intensities = spectrum['Y(Counts)'].values    # Intensity values
        return mz, intensities
    except Exception as e:
        if 'skiprows' not in str(e):
            print(f"Error loading {filepath}: {str(e)}")
        return None, None

# Function to preprocess and align two spectra
def align_spectra(mz1, intensities1, mz2, intensities2, n_points=1000):
    """Aligns two spectra for comparison by:
    1. Creating common m/z axis
    2. Interpolating intensities
    3. Normalizing intensities to [0,1] range"""
    try:
        # Find overlapping m/z range
        min_mz = max(np.min(mz1), np.min(mz2))
        max_mz = min(np.max(mz1), np.max(mz2))
        
        # Create uniform m/z axis with 1000 points
        common_mz = np.linspace(min_mz, max_mz, n_points)
        
        # Interpolate intensities to match common m/z axis
        f1 = interp1d(mz1, intensities1, kind='linear', bounds_error=False, fill_value=0)
        f2 = interp1d(mz2, intensities2, kind='linear', bounds_error=False, fill_value=0)
        
        int1_aligned = f1(common_mz)
        int2_aligned = f2(common_mz)
        
        # Normalize intensities to [0,1] range
        max_int1 = np.max(int1_aligned)
        max_int2 = np.max(int2_aligned)
        
        if max_int1 > 0:
            int1_aligned = int1_aligned / max_int1
        if max_int2 > 0:
            int2_aligned = int2_aligned / max_int2
        
        return common_mz, int1_aligned, int2_aligned
        
    except Exception as e:
        print(f"Error aligning spectra: {str(e)}")
        return None, None, None

# Function to calculate cosine similarity only
def calculate_similarity_features(mz, int1, int2):
    """Calculates cosine similarity between spectra vectors"""
    try:
        # Calculate cosine similarity (angle between vectors)
        cosine_sim = np.dot(int1, int2) / (np.linalg.norm(int1) * np.linalg.norm(int2))
        
        features = {
            'cosine_similarity': cosine_sim,
        }
        
        # Check for invalid values
        if not all(np.isfinite(v) for v in features.values()):
            return None
            
        return features
        
    except Exception as e:
        print(f"Error calculating features: {str(e)}")
        return None

# Function to check if compound files exist
def check_compound_existence(base_path, compound_num):
    """Checks if at least one spectrum file exists for a compound"""
    for i in range(1, 11):
        if os.path.exists(f"{base_path}/{compound_num:02d}-{i:02d}.csv"):
            return True
    return False

# Main function to generate comparison dataset
def generate_balanced_comparisons(base_path, debug=False):
    """Creates balanced dataset of spectrum comparisons:
    - Positive class (y=1): Same compound comparisons
    - Negative class (y=0): Different compound comparisons"""
    X = []  # Features
    y = []  # Labels
    comparisons = []  # File pairs
    same_count = 0    # Count of same-compound pairs
    diff_count = 0    # Count of different-compound pairs
    
    # Find all existing compounds
    existing_compounds = []
    for compound_i in range(1, 29):  # CM1 has 28 compounds
        if check_compound_existence(base_path, compound_i):
            existing_compounds.append(compound_i)
    
    if debug:
        print(f"\nFound {len(existing_compounds)} compounds")
        
    # Generate same-compound comparisons
    for compound_i in existing_compounds:
        if debug:
            print(f"Processing same-compound comparisons for compound {compound_i}")
            
        # Compare first 5 spectra of each compound
        for i in range(1, 5):
            for j in range(i + 1, 6):
                file1 = f"{base_path}/{compound_i:02d}-{i:02d}.csv"
                file2 = f"{base_path}/{compound_i:02d}-{j:02d}.csv"
                
                if not os.path.exists(file1) or not os.path.exists(file2):
                    continue
                    
                features = calculate_similarity_features(
                    *align_spectra(*load_spectrum(file1), *load_spectrum(file2))
                )
                
                if features is not None:
                    X.append(list(features.values()))
                    y.append(1)  # Same compound
                    comparisons.append((file1, file2))
                    same_count += 1
    
    # Generate different-compound comparisons
    for idx, compound_i in enumerate(existing_compounds[:-1]):
        if debug:
            print(f"Processing different-compound comparisons for compound {compound_i}")
            
        # Compare with next 2 compounds for balance
        for compound_j in existing_compounds[idx+1:min(idx+3, len(existing_compounds))]:
            for i in range(1, 6):
                for j in range(1, 6):
                    file1 = f"{base_path}/{compound_i:02d}-{i:02d}.csv"
                    file2 = f"{base_path}/{compound_j:02d}-{j:02d}.csv"
                    
                    if not os.path.exists(file1) or not os.path.exists(file2):
                        continue
                    
                    features = calculate_similarity_features(
                        *align_spectra(*load_spectrum(file1), *load_spectrum(file2))
                    )
                    
                    if features is not None:
                        X.append(list(features.values()))
                        y.append(0)  # Different compounds
                        comparisons.append((file1, file2))
                        diff_count += 1
    
    # Print dataset statistics
    print("\nComparison Statistics:")
    print(f"Same compound comparisons: {same_count}")
    print(f"Different compound comparisons: {diff_count}")
    print(f"Class balance ratio (same/diff): {same_count/diff_count:.2f}")
    
    if len(X) == 0:
        print("No valid comparisons generated!")
        return None, None, None
        
    return np.array(X), np.array(y), comparisons

# Function to visualize confusion matrix
def plot_confusion_matrix(y_true, y_pred):
    """Creates and plots confusion matrix for model evaluation"""
    cm = confusion_matrix(y_true, y_pred)
    plt.figure(figsize=(8, 6))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues',
                xticklabels=['Different', 'Same'],
                yticklabels=['Different', 'Same'])
    plt.title('Confusion Matrix')
    plt.ylabel('True Label')
    plt.xlabel('Predicted Label')
    plt.tight_layout()
    plt.show()
    return cm

# Dataset preparation
def prepare_dataset(base_path, test_size=0.4, random_state=42):
    """Generates dataset, splits it into training and testing sets, and scales the features."""
    X, y, comparisons = generate_balanced_comparisons(base_path, debug=True)
    
    if X is None:
        print("No data generated. Exiting...")
        return None, None, None, None
    
    # Split into training and test sets
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=random_state, stratify=y
    )
    
    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    return X_train_scaled, X_test_scaled, y_train, y_test, scaler

# Initialize XGBoost model
def initialize_xgboost_model(): 
    """Initializes the XGBoost model with optimized parameters."""
    return xgb.XGBClassifier(
        max_depth=4,
        min_child_weight=2,
        gamma=0.1,
        subsample=0.9,
        colsample_bytree=0.9,
        learning_rate=0.005,
        n_estimators=300,
        scale_pos_weight=4,
        reg_lambda=0.8,
        reg_alpha=0.1,
        objective='binary:logistic',
        random_state=42,
        use_label_encoder=False
    )

# Cross-validation
def perform_cross_validation(model, X_train, y_train):
    """Performs cross-validation and prints the results."""
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    cv_scores = cross_val_score(model, X_train, y_train, cv=cv, scoring='accuracy')
    print("\nCross-validation scores:", cv_scores)
    print(f"Mean CV score: {cv_scores.mean():.4f} (+/- {cv_scores.std() * 2:.4f})")
    return cv_scores

# Model train and evaluation
def train_and_evaluate(model, X_train, y_train, X_test, y_test):
    """Trains the model and evaluates its performance on the test set."""
    print("\nTraining model...")
    model.fit(X_train, y_train)
    
    y_pred = model.predict(X_test)
    y_prob = model.predict_proba(X_test)[:, 1]
    
    print("\nTest Set Performance:")
    print(f"Accuracy: {accuracy_score(y_test, y_pred):.4f}")
    print(f"ROC AUC: {roc_auc_score(y_test, y_prob):.4f}")
    print("\nClassification Report:")
    print(classification_report(y_test, y_pred))
    
    return y_pred, y_prob

# feature importance
def plot_feature_importance(model, feature_names):
    """Plots the importance of features."""
    importance = model.feature_importances_
    
    # Create a DataFrame for feature importance
    importance_df = pd.DataFrame({
        'Feature': feature_names,
        'Importance': importance
    }).sort_values('Importance', ascending=False)

    # Plot the feature importance
    plt.figure(figsize=(10, 6))
    sns.barplot(x='Importance', y='Feature', data=importance_df, palette='viridis')
    plt.title('Feature Importance')
    plt.xlabel('Importance Score')
    plt.ylabel('Features')
    plt.tight_layout()
    plt.show()
    
    # Print feature importance values
    print("\nFeature Importance Values:")
    print(importance_df.to_string(index=False))  # Print DataFrame without the index

# correlation
def plot_correlation(y_test, y_prob):
    """Plots the correlation between predicted probabilities and actual labels."""
    correlation = np.corrcoef(y_test, y_prob)[0, 1]
    print(f"\nCorrelation between predicted probabilities and actual labels: {correlation:.4f}")

    plt.figure(figsize=(8, 6))
    plt.scatter(y_test, y_prob, alpha=0.5, edgecolor='k')
    plt.title(f"Correlation: {correlation:.4f}", fontsize=16)
    plt.xlabel("True Labels", fontsize=14)
    plt.ylabel("Predicted Probabilities", fontsize=14)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()