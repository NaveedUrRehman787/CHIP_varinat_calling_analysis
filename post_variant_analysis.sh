#!/bin/bash

# ============================
# Post-Variant Analysis Script
# ============================

set -e
mkdir -p results/publication_ready
SAMPLE="SRR32229173"
REF="ref/hg38.fa"
OUTDIR="results/publication_ready"
INVCF_MUTECT="results/${SAMPLE}.mutect2.filtered.vcf"
VCF_MUTECT="$OUTDIR/${SAMPLE}.mutect2.filtered.highconf.vcf"
VCF_LOFREQ="results/${SAMPLE}.lofreq.vcf"

# Load conda shell support
eval "$(conda shell.bash hook)"
JAR_PATH=$(find $CONDA_PREFIX -name "snpEff.jar" | head -n 1)

# === 0. Pre-filter Mutect2 VCF ===
echo "🔧 Filtering Mutect2 VCF for AF > 0.01 and DP > 10 ..."
conda activate bcftools_env >/dev/null
bcftools filter -i 'FORMAT/DP>10 && FORMAT/AF>0.01' "$INVCF_MUTECT" > "$VCF_MUTECT"

# === 1. Annotate using SnpEff ===
echo "🧬 Annotating VCFs using SnpEff with 10G memory..."
conda activate variant_calling >/dev/null
java -Xmx10g -jar "$JAR_PATH" -noStats -noLog -v hg38 "$VCF_MUTECT" > "$OUTDIR/${SAMPLE}.mutect2.ann.vcf"
java -Xmx4g  -jar "$JAR_PATH" -noStats -noLog -v hg38 "$VCF_LOFREQ" > "$OUTDIR/${SAMPLE}.lofreq.ann.vcf"

# === 2. Convert VCFs to TSV ===
echo "📁 Extracting tables from VCFs..."
conda activate bcftools_env >/dev/null
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FORMAT/AF\n' "$VCF_MUTECT" > "$OUTDIR/${SAMPLE}.mutect2.tsv"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' "$VCF_LOFREQ" > "$OUTDIR/${SAMPLE}.lofreq.tsv"

# === 3. Filter CHIP-related genes ===
echo "🧬 Filtering CHIP genes..."
CHIP_GENES="DNMT3A|TET2|ASXL1|JAK2|TP53|SF3B1|SRSF2|U2AF1|CBL|IDH2|IDH1"
grep -Ei "$CHIP_GENES" "$OUTDIR/${SAMPLE}.mutect2.ann.vcf" > "$OUTDIR/${SAMPLE}.mutect2.chip_variants.vcf" || true
grep -Ei "$CHIP_GENES" "$OUTDIR/${SAMPLE}.lofreq.ann.vcf" > "$OUTDIR/${SAMPLE}.lofreq.chip_variants.vcf" || true

# === 4. Compare variant positions ===
echo "🔍 Comparing LoFreq vs Mutect2 variants..."
grep -v "^#" "$VCF_MUTECT" | cut -f1,2 | sort > "$OUTDIR/mutect_pos.txt"
grep -v "^#" "$VCF_LOFREQ" | cut -f1,2 | sort > "$OUTDIR/lofreq_pos.txt"

comm -12 "$OUTDIR/mutect_pos.txt" "$OUTDIR/lofreq_pos.txt" > "$OUTDIR/overlap_positions.txt"
comm -23 "$OUTDIR/mutect_pos.txt" "$OUTDIR/lofreq_pos.txt" > "$OUTDIR/mutect_unique.txt"
comm -13 "$OUTDIR/mutect_pos.txt" "$OUTDIR/lofreq_pos.txt" > "$OUTDIR/lofreq_unique.txt"

# === 5. Summary ===
echo "✅ Done. Publication-ready outputs:"
ls -lh "$OUTDIR"
