#!/Users/bm0211/miniconda3/bin/python

import pandas as pd
import numpy as np
import math

# Load the data
df = pd.read_csv("RAW_Data_for_NMR_ANALYSIS.csv", low_memory=False)

# Convert key variables to numeric (coerce errors to NaN)
numeric_cols = [
    'Age',
    'Total.cholesterol.mg.dL.',
    'HDL.cholesterol.mg.dL',
    'Blood.pressure.mm.Hg.systolic.'
]
for col in numeric_cols:
    df[col] = pd.to_numeric(df[col], errors='coerce')

# Clean relevant string variables
df['Sex'] = df['Sex'].str.strip().str.capitalize()
df['Tobacco.yes.no'] = df['Tobacco.yes.no'].astype(str).str.strip().str.lower()
df['Type.I.diabetes'] = df['Type.I.diabetes'].astype(str).str.strip().str.lower()
df['Type.II.diabetes'] = df['Type.II.diabetes'].astype(str).str.strip().str.lower()
df['medication_yes_no'] = df['medication_yes_no'].astype(str).str.strip().str.lower()

# Drop rows with missing required values
df = df.dropna(subset=numeric_cols + ['Sex', 'Tobacco.yes.no', 'medication_yes_no'])

# Define the Framingham risk function
def compute_risk(row):
    # Coefficients by sex
    if row['Sex'] == 'Male':
        coeffs = {
            'ln_age': 3.06117,
            'ln_total_chol': 1.12370,
            'ln_hdl': -0.93263,
            'ln_sbp_treated': 1.99881,
            'ln_sbp_untreated': 1.93303,
            'smoker': 0.65451,
            'diabetic': 0.57367,
            'mean': 23.9802,
            'baseline_survival': 0.88936
        }
    else:
        coeffs = {
            'ln_age': 2.32888,
            'ln_total_chol': 1.20904,
            'ln_hdl': -0.70833,
            'ln_sbp_treated': 2.82263,
            'ln_sbp_untreated': 2.76157,
            'smoker': 0.52873,
            'diabetic': 0.69154,
            'mean': 26.1931,
            'baseline_survival': 0.95012
        }

    # Protect against log(0) or negative input
    age = max(row['Age'], 1)
    total_chol = max(row['Total.cholesterol.mg.dL.'], 1)
    hdl = max(row['HDL.cholesterol.mg.dL'], 1)
    sbp = max(row['Blood.pressure.mm.Hg.systolic.'], 1)

    # Log-transform
    ln_age = np.log(age)
    ln_total_chol = np.log(total_chol)
    ln_hdl = np.log(hdl)
    ln_sbp = np.log(sbp)

    # Choose SBP coefficient
    ln_sbp_coef = coeffs['ln_sbp_treated'] if row['medication_yes_no'] == 'hypertension' else coeffs['ln_sbp_untreated']

    # Smoking
    smoker = 1 if row['Tobacco.yes.no'] == 'yes' else 0

    # Diabetes
    diabetic = 1 if row['Type.I.diabetes'] == 'yes' or row['Type.II.diabetes'] == 'yes' else 0

    # Framingham risk score
    sum_beta_x = (
        coeffs['ln_age'] * ln_age +
        coeffs['ln_total_chol'] * ln_total_chol +
        coeffs['ln_hdl'] * ln_hdl +
        ln_sbp_coef * ln_sbp +
        coeffs['smoker'] * smoker +
        coeffs['diabetic'] * diabetic
    )

    risk_score = 1 - (coeffs['baseline_survival'] ** np.exp(sum_beta_x - coeffs['mean']))
    return round(risk_score * 100, 2)

# Apply model
df['Framingham_10yr_CVD_Risk'] = df.apply(compute_risk, axis=1)

# Preview output
print(df[['Sex', 'Age', 'Framingham_10yr_CVD_Risk']].head())

# Optional: Save to CSV
df.to_csv("RAW_Data_for_NMR_ANALYSIS_python_ran.csv", index=False)

