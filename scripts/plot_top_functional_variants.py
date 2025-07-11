import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

VCF_DIR = "results/publication_ready"
FILES = [
    "SRR32229173.mutect2.top_functional.vcf",
    "SRR32229173.lofreq.top_functional.vcf"
]

def parse_vcf_effects(vcf_file):
    effects = []
    genes = []
    with open(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            if "ANN=" in line:
                try:
                    ann_field = line.split("ANN=")[1].split(";")[0]
                    annotations = ann_field.split(",")
                    for ann in annotations:
                        fields = ann.split("|")
                        if len(fields) > 4:
                            effects.append(fields[1])
                            genes.append(fields[3])
                except IndexError:
                    continue
    return effects, genes

plot_dir = os.path.join(VCF_DIR, "plots")
os.makedirs(plot_dir, exist_ok=True)

for file in FILES:
    tag = "mutect2" if "mutect2" in file else "lofreq"
    print(f"🔍 Parsing {file} ...")

    effects, genes = parse_vcf_effects(os.path.join(VCF_DIR, file))

    # --- Plot 1: Effect type barplot ---
    eff_df = pd.DataFrame(Counter(effects).items(), columns=["Effect", "Count"])
    eff_df = eff_df.sort_values("Count", ascending=False)

    plt.figure(figsize=(10, 5))
    sns.barplot(data=eff_df, x="Effect", y="Count", palette="rocket")
    plt.title(f"Top Functional Variant Types ({tag})")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(f"{plot_dir}/{tag}_effect_barplot.png")
    plt.close()

    # --- Plot 2: CHIP gene counts ---
    CHIP_GENES = {"DNMT3A","TET2","ASXL1","JAK2","TP53","SF3B1","SRSF2","U2AF1","CBL","IDH2","IDH1"}
    chip_hits = [g for g in genes if g in CHIP_GENES]
    chip_df = pd.DataFrame(Counter(chip_hits).items(), columns=["Gene", "Count"]).sort_values("Count", ascending=False)

    if not chip_df.empty:
        plt.figure(figsize=(8, 4))
        sns.barplot(data=chip_df, x="Gene", y="Count", palette="mako")
        plt.title(f"CHIP Gene Hits ({tag})")
        plt.tight_layout()
        plt.savefig(f"{plot_dir}/{tag}_chip_gene_hits.png")
        plt.close()

print("✅ Graphs saved to:", plot_dir)
