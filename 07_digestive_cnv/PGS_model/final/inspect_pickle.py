"""Analyze digestive-enzyme copy numbers and haplotypes."""

import os
import pickle
import pandas as pd

def safe_load_pickle(path):
    try:
        with open(path, "rb") as f:
            return pickle.load(f)
    except Exception as e:
        print(f"❌ Failed to load {path}: {e}")
        return None

correlations = []

for i in range(1, 103):
    filename = f"iteration_{i}_validation.pkl"
    filepath = os.path.join(os.getcwd(), filename)
    data = safe_load_pickle(filepath)
    if data is not None:
        correlations.append({
            "iteration": i,
            "train_corr": data["train"]["correlation"],
            "test_corr": data["test"]["correlation"],
            "malay_corr": data["malay"]["correlation"],
            "bajau_corr": data["bajau"]["correlation"],
        })

df = pd.DataFrame(correlations)
df.set_index("iteration", inplace=True)

print(df.head())
df.to_csv("../correlations_summary.csv")  #  Save results
