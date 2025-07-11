import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# === Parameters ===
input_file = "results/publication_ready/SRR32229173.top_functional.clinvar_merged.tsv"
chip_genes = {"DNMT3A", "TET2", "ASXL1", "JAK2", "TP53", "SF3B1", "SRSF2", "U2AF1", "CBL", "IDH2", "IDH1"}

# === Load merged table ===
df = pd.read_csv(input_file, sep="\t")

# --- Ensure AF column is numeric and drop missing ---
if "AF" not in df.columns:
    df["AF"] = np.nan  # Placeholder if missing

df["AF"] = pd.to_numeric(df["AF"], errors="coerce")
df = df.dropna(subset=["AF"])

# Binning
bins = [0, 0.05, 0.10, 0.20, 1.01]
labels = ["<5%", "5–10%", "10–20%", ">20%"]
df["AF_bin"] = pd.cut(df["AF"], bins=bins, labels=labels, right=False)

# Filter for CHIP gene hits
df_chip = df[df["Effect"].notnull() & df["CHROM"].notnull() & df["ALT"].notnull()]
df_chip["Gene"] = df_chip["Effect"]  # Placeholder if no gene info
df_chip = df_chip[df_chip["Gene"].isin(chip_genes)]

# Count by bin
summary = df_chip.groupby(["Gene", "AF_bin"]).size().reset_index(name="Count")
summary.to_csv("results/publication_ready/vaf_clone_bins.tsv", sep="\t", index=False)

# Plot
plt.figure(figsize=(10, 5))
sns.barplot(data=summary, x="AF_bin", y="Count", hue="Gene", palette="tab10")
plt.title("CHIP Gene Variant Burden by Clonal VAF Bin")
plt.xlabel("Variant Allele Frequency (VAF)")
plt.ylabel("Variant Count")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig("results/publication_ready/vaf_clone_bins.png")
plt.close()

print("✅ VAF clone bin summary and plot saved to: publication_ready/")
