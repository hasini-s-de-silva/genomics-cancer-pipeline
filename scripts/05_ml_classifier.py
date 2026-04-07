import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.preprocessing import LabelEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    roc_auc_score,
    roc_curve,
    confusion_matrix,
    ConfusionMatrixDisplay
)

# -----------------------------
# Load data
# -----------------------------
expr = pd.read_csv("data/processed/expression_filtered.csv")
metadata = pd.read_csv("data/processed/metadata_matched.csv")

print("Expression shape:", expr.shape)
print("Metadata shape:", metadata.shape)

# -----------------------------
# Prepare expression matrix
# -----------------------------
expr = expr.set_index("gene")
X = expr.T.copy()

# Clean sample IDs
X.index = (
    X.index.astype(str)
    .str.replace(".", "-", regex=False)
    .str.replace(r"^X", "", regex=True)
    .str[:15]
)

metadata["sample_id"] = (
    metadata["sample_id"].astype(str)
    .str.replace(".", "-", regex=False)
    .str.replace(r"^X", "", regex=True)
    .str[:15]
)

# Match metadata to expression samples
metadata = metadata.set_index("sample_id").loc[X.index].reset_index()

# Labels
y = metadata["sample_type"]

print("X shape before split:", X.shape)
print("y distribution:")
print(y.value_counts())

# -----------------------------
# Encode labels
# -----------------------------
le = LabelEncoder()
y_encoded = le.fit_transform(y)

print("Label mapping:", dict(zip(le.classes_, le.transform(le.classes_))))

# -----------------------------
# Train/test split FIRST
# -----------------------------
X_train, X_test, y_train, y_test = train_test_split(
    X,
    y_encoded,
    test_size=0.2,
    random_state=42,
    stratify=y_encoded
)

print("Train shape before feature selection:", X_train.shape)
print("Test shape before feature selection:", X_test.shape)

# -----------------------------
# Feature selection on TRAINING data only
# -----------------------------
y_train_series = pd.Series(y_train, index=X_train.index, name="label")
train_means = X_train.groupby(y_train_series).mean()

if train_means.shape[0] != 2:
    raise ValueError("Training data does not contain both classes.")

mean_diff = (train_means.iloc[1] - train_means.iloc[0]).abs()

# Select top 100 genes
top_genes = mean_diff.sort_values(ascending=False).head(100).index.tolist()

print("Number of top genes selected from training data:", len(top_genes))
print("Top 10 selected genes:", top_genes[:10])

# Subset train/test using selected genes
X_train = X_train[top_genes]
X_test = X_test[top_genes]

print("Train shape after feature selection:", X_train.shape)
print("Test shape after feature selection:", X_test.shape)

# -----------------------------
# Cross-validation setup
# -----------------------------
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# -----------------------------
# Logistic Regression
# -----------------------------
log_model = LogisticRegression(max_iter=5000, random_state=42)

# Cross-validation scores
log_cv_acc = cross_val_score(log_model, X_train, y_train, cv=cv, scoring="accuracy")
log_cv_auc = cross_val_score(log_model, X_train, y_train, cv=cv, scoring="roc_auc")

# Fit on full training set
log_model.fit(X_train, y_train)

y_pred_log = log_model.predict(X_test)
y_prob_log = log_model.predict_proba(X_test)[:, 1]

log_acc = accuracy_score(y_test, y_pred_log)
log_auc = roc_auc_score(y_test, y_prob_log)

print("\n=== Logistic Regression ===")
print("Test Accuracy:", round(log_acc, 4))
print("Test ROC-AUC:", round(log_auc, 4))
print("CV Accuracy Mean:", round(log_cv_acc.mean(), 4))
print("CV Accuracy SD:", round(log_cv_acc.std(), 4))
print("CV ROC-AUC Mean:", round(log_cv_auc.mean(), 4))
print("CV ROC-AUC SD:", round(log_cv_auc.std(), 4))
print(classification_report(y_test, y_pred_log, target_names=le.classes_))

# -----------------------------
# Random Forest
# -----------------------------
rf_model = RandomForestClassifier(
    n_estimators=300,
    random_state=42,
    class_weight="balanced"
)

# Cross-validation scores
rf_cv_acc = cross_val_score(rf_model, X_train, y_train, cv=cv, scoring="accuracy")
rf_cv_auc = cross_val_score(rf_model, X_train, y_train, cv=cv, scoring="roc_auc")

# Fit on full training set
rf_model.fit(X_train, y_train)

y_pred_rf = rf_model.predict(X_test)
y_prob_rf = rf_model.predict_proba(X_test)[:, 1]

rf_acc = accuracy_score(y_test, y_pred_rf)
rf_auc = roc_auc_score(y_test, y_prob_rf)

print("\n=== Random Forest ===")
print("Test Accuracy:", round(rf_acc, 4))
print("Test ROC-AUC:", round(rf_auc, 4))
print("CV Accuracy Mean:", round(rf_cv_acc.mean(), 4))
print("CV Accuracy SD:", round(rf_cv_acc.std(), 4))
print("CV ROC-AUC Mean:", round(rf_cv_auc.mean(), 4))
print("CV ROC-AUC SD:", round(rf_cv_auc.std(), 4))
print(classification_report(y_test, y_pred_rf, target_names=le.classes_))

# -----------------------------
# Save performance summary
# -----------------------------
performance_df = pd.DataFrame({
    "Model": ["Logistic Regression", "Random Forest"],
    "Test_Accuracy": [log_acc, rf_acc],
    "Test_ROC_AUC": [log_auc, rf_auc],
    "CV_Accuracy_Mean": [log_cv_acc.mean(), rf_cv_acc.mean()],
    "CV_Accuracy_SD": [log_cv_acc.std(), rf_cv_acc.std()],
    "CV_ROC_AUC_Mean": [log_cv_auc.mean(), rf_cv_auc.mean()],
    "CV_ROC_AUC_SD": [log_cv_auc.std(), rf_cv_auc.std()]
})

performance_df.to_csv("results/tables/ml_model_performance.csv", index=False)

# -----------------------------
# ROC curve
# -----------------------------
fpr_log, tpr_log, _ = roc_curve(y_test, y_prob_log)
fpr_rf, tpr_rf, _ = roc_curve(y_test, y_prob_rf)

plt.figure(figsize=(7, 5))
plt.plot(fpr_log, tpr_log, label=f"Logistic Regression (AUC = {log_auc:.3f})")
plt.plot(fpr_rf, tpr_rf, label=f"Random Forest (AUC = {rf_auc:.3f})")
plt.plot([0, 1], [0, 1], linestyle="--")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve for Tumor vs Normal Classification")
plt.legend()
plt.tight_layout()
plt.savefig("results/figures/roc_curve.png", dpi=300)
plt.close()

# -----------------------------
# Confusion matrix (Random Forest)
# -----------------------------
cm = confusion_matrix(y_test, y_pred_rf)

disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=le.classes_)
fig, ax = plt.subplots(figsize=(5, 5))
disp.plot(ax=ax, colorbar=False)
plt.title("Confusion Matrix, Random Forest")
plt.tight_layout()
plt.savefig("results/figures/confusion_matrix_rf.png", dpi=300)
plt.close()

# -----------------------------
# Feature importance (Random Forest)
# -----------------------------
feature_importance = pd.DataFrame({
    "gene": X_train.columns,
    "importance": rf_model.feature_importances_
}).sort_values("importance", ascending=False)

feature_importance.to_csv("results/tables/rf_feature_importance.csv", index=False)

top_features = feature_importance.head(20)

plt.figure(figsize=(8, 6))
plt.barh(top_features["gene"][::-1], top_features["importance"][::-1])
plt.xlabel("Importance")
plt.ylabel("Gene")
plt.title("Top 20 Important Genes (Random Forest)")
plt.tight_layout()
plt.savefig("results/figures/rf_top20_genes.png", dpi=300)
plt.close()

print("\nSaved outputs:")
print("- results/tables/ml_model_performance.csv")
print("- results/tables/rf_feature_importance.csv")
print("- results/figures/roc_curve.png")
print("- results/figures/confusion_matrix_rf.png")
print("- results/figures/rf_top20_genes.png")