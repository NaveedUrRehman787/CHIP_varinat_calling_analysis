import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
import os

fig_dir = "results/publication_ready"
fig, axs = plt.subplots(3, 2, figsize=(13, 12))
gs = GridSpec(3, 2, figure=fig)

# A: AF Distribution
img = mpimg.imread(f"{fig_dir}/af_distribution.png")
axs[0, 0].imshow(img)
axs[0, 0].axis("off")
axs[0, 0].set_title("A. Allele Frequency Distribution")

# B: Venn
img = mpimg.imread(f"{fig_dir}/variant_overlap_venn.png")
axs[0, 1].imshow(img)
axs[0, 1].axis("off")
axs[0, 1].set_title("B. Overlapping Variant Positions")

# C: Functional Effects
img = mpimg.imread(f"{fig_dir}/plots/mutect2_effect_barplot.png")
axs[1, 0].imshow(img)
axs[1, 0].axis("off")
axs[1, 0].set_title("C. Functional Effects (Mutect2)")

# D: CHIP gene hits
img = mpimg.imread(f"{fig_dir}/plots/mutect2_chip_gene_hits.png")
axs[1, 1].imshow(img)
axs[1, 1].axis("off")
axs[1, 1].set_title("D. CHIP Gene Hits")

# E: ClinVar Significance
img = mpimg.imread(f"{fig_dir}/SRR32229173.clinvar.clnsig_counts.png")
axs[2, 0].imshow(img)
axs[2, 0].axis("off")
axs[2, 0].set_title("E. ClinVar Clinical Significance (CLNSIG)")

# Empty panel
axs[2, 1].axis("off")

plt.tight_layout()
plt.savefig(f"{fig_dir}/figures_summary.pdf")
print("✅ Publication figure saved: figures_summary.pdf")
