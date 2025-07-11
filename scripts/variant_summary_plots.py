import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib_venn import venn2
from collections import Counter

# === Load AF TSVs ===
mutect_df = pd.read_csv("results/publication_ready/SRR32229173.mutect2.tsv", sep='\t', header=None,
                        names=["CHROM", "POS", "REF", "ALT", "AF"])
lofreq_df = pd.read_csv("results/publication_ready/SRR32229173.lofreq.tsv", sep='\t', header=None,
                        names=["CHROM", "POS", "REF", "ALT", "AF"])

# Convert AF to numeric
mutect_df["AF"] = pd.to_numeric(mutect_df["AF"], errors='coerce')
lofreq_df["AF"] = pd.to_numeric(lofreq_df["AF"], errors='coerce')

# === AF Histogram ===
plt.figure(figsize=(10, 5))
sns.histplot(mutect_df["AF"], bins=50, kde=True, label="Mutect2", color="skyblue", stat="density")
sns.histplot(lofreq_df["AF"], bins=50, kde=True, label="LoFreq", color="orange", stat="density")
plt.xlabel("Allele Frequency (AF)")
plt.ylabel("Density")
plt.title("AF Distribution: Mutect2 vs LoFreq")
plt.legend()
plt.tight_layout()
plt.savefig("results/publication_ready/af_distribution.png")
plt.close()

# === Venn Diagram ===
def load_positions(filepath):
    with open(filepath) as f:
        return set(line.strip() for line in f if line.strip())

mutect_pos = load_positions("results/publication_ready/mutect_pos.txt")
lofreq_pos = load_positions("results/publication_ready/lofreq_pos.txt")

plt.figure(figsize=(5, 5))
venn2([mutect_pos, lofreq_pos], set_labels=("Mutect2", "LoFreq"))
plt.title("Variant Position Overlap")
plt.tight_layout()
plt.savefig("results/publication_ready/variant_overlap_venn.png")
plt.close()

# === CHIP Gene Hit Count ===
def extract_genes_from_vcf(vcf_path):
    genes = []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            if "ANN=" in line:
                try:
                    ann_field = line.strip().split("ANN=")[1].split(";")[0]
                    gene = ann_field.split("|")[3]
                    if gene:
                        genes.append(gene)
                except IndexError:
                    continue
    return genes

chip_mutect_genes = extract_genes_from_vcf("results/publication_ready/SRR32229173.mutect2.chip_variants.vcf")
chip_lofreq_genes = extract_genes_from_vcf("results/publication_ready/SRR32229173.lofreq.chip_variants.vcf")

# Combine and count
chip_counts = Counter(chip_mutect_genes + chip_lofreq_genes)
df_chip = pd.DataFrame(chip_counts.items(), columns=["Gene", "Count"]).sort_values("Count", ascending=False)

plt.figure(figsize=(10, 5))
sns.barplot(data=df_chip, x="Gene", y="Count", palette="mako")
plt.xticks(rotation=45, ha="right")
plt.title("CHIP Gene Hits (Mutect2 + LoFreq)")
plt.tight_layout()
plt.savefig("results/publication_ready/chip_gene_hits.png")
plt.close()

print("✅ Plots saved to results/publication_ready/")
