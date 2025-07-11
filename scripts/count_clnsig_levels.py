import re
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

vcf_path = "results/publication_ready/SRR32229173.mutect2.clinvar.vcf"
counts = Counter()

with open(vcf_path) as f:
    for line in f:
        if line.startswith("#"): continue
        match = re.search(r'CLNSIG=([^;]+)', line)
        if match:
            vals = match.group(1).split('|')[0].split(',')
            for val in vals:
                counts[val.strip()] += 1

# Save TSV
df = pd.DataFrame(counts.items(), columns=["CLNSIG", "Count"]).sort_values("Count", ascending=False)
df.to_csv("results/publication_ready/SRR32229173.clinvar.clnsig_counts.tsv", sep="\t", index=False)

# Barplot
plt.figure(figsize=(8,5))
sns.barplot(data=df, x="CLNSIG", y="Count", palette="Set2")
plt.title("ClinVar CLNSIG Variant Counts")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig("results/publication_ready/SRR32229173.clinvar.clnsig_counts.png")
print("✅ CLNSIG counts saved + plotted.")
