
#!/bin/bash

# =============================
# ✅ Reproducible Variant Calling Pipeline
# From FASTQ to Annotated VCFs
# =============================

# === 0. Setup and Conda Environments ===
conda activate variant_calling

# === 1. BWA Index ===
bwa index ref/hg38.fa

# === 2. Align Paired-End Reads ===
bwa mem -t 8 ref/hg38.fa data/SRR32229173_1.fastq data/SRR32229173_2.fastq > results/SRR32229173.sam

# === 3. Convert SAM to BAM, Sort ===
conda activate samtools_env
samtools view -@ 8 -bS results/SRR32229173.sam > results/SRR32229173.bam
samtools sort -@ 8 -o results/SRR32229173.sorted.bam results/SRR32229173.bam
samtools index results/SRR32229173.sorted.bam
conda activate variant_calling

# === 4. Mark Duplicates with Picard ===
picard MarkDuplicates   I=results/SRR32229173.sorted.bam   O=results/SRR32229173.dedup.bam   M=results/SRR32229173.dup_metrics.txt   CREATE_INDEX=true

# === 5. Add Read Group (optional, for Mutect2) ===
gatk AddOrReplaceReadGroups   -I results/SRR32229173.dedup.bam   -O results/SRR32229173.rg.bam   -ID ID -LB lib -PL illumina -PU unit -SM SRR32229173   --CREATE_INDEX true

# === 6. LoFreq Indel Quality + Call ===
conda activate samtools_env
lofreq indelqual --dindel -f ref/hg38.fa -o results/SRR32229173.indelqual.bam results/SRR32229173.dedup.bam
samtools index results/SRR32229173.indelqual.bam
lofreq call-parallel --pp-threads 8 -f ref/hg38.fa -o results/SRR32229173.lofreq.vcf results/SRR32229173.indelqual.bam
conda activate variant_calling

# === 7. Mutect2 Variant Calling + Filtering ===
gatk Mutect2   -R ref/hg38.fa   -I results/SRR32229173.rg.bam   -L ref/covered_regions.bed   -O results/SRR32229173.mutect2.unfiltered.vcf

gatk FilterMutectCalls   -R ref/hg38.fa   -V results/SRR32229173.mutect2.unfiltered.vcf   -O results/SRR32229173.mutect2.filtered.vcf

# === 8. Optional High Confidence Filtering ===
cp results/SRR32229173.mutect2.filtered.vcf results/publication_ready/SRR32229173.mutect2.filtered.highconf.vcf

# === 9. Variant Annotation with SnpEff ===
java -Xmx8g -jar $(which snpEff) -v hg38 results/publication_ready/SRR32229173.mutect2.filtered.highconf.vcf > results/publication_ready/SRR32229173.mutect2.ann.vcf
java -Xmx8g -jar $(which snpEff) -v hg38 results/SRR32229173.lofreq.vcf > results/publication_ready/SRR32229173.lofreq.ann.vcf

# === 10. Generate Summary Report ===
bash summary_report.sh

# === 11. Post-Variant Analysis ===
bash post_variant_analysis.sh

# === 12. Plot Functional Variant Types ===
python plot_top_functional_variants.py

# === 13. ClinVar Annotation ===
bgzip -c results/publication_ready/SRR32229173.mutect2.filtered.highconf.vcf > results/publication_ready/SRR32229173.mutect2.filtered.highconf.vcf.gz
tabix -p vcf results/publication_ready/SRR32229173.mutect2.filtered.highconf.vcf.gz

bcftools annotate \
  -a ref/databases/clinvar.vcf.gz \
  -c CHROM,POS,REF,ALT,INFO/CLNSIG,INFO/CLNREVSTAT \
  -h <(echo "##INFO=<ID=CLNSIG,Number=.,Type=String,Description='Clinical significance from ClinVar'>") \
  -h <(echo "##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description='Review status from ClinVar'>") \
  -o results/publication_ready/SRR32229173.mutect2.clinvar.vcf \
  -O v results/publication_ready/SRR32229173.mutect2.filtered.highconf.vcf.gz

# === 14. Extract Top Functional Variants ===
grep -Ev "^#" results/publication_ready/SRR32229173.mutect2.ann.vcf | grep -E "missense_variant|stop_gained|frameshift_variant|start_lost|splice_acceptor|splice_donor" > results/publication_ready/SRR32229173.mutect2.top_functional.vcf

# === 15. Analyze ClinVar Significance ===
python count_clnsig_levels.py

# === 16. Merge Functional + ClinVar ===
python merge_topfx_clinvar.py

# === 17. Clonal Burden Analysis ===
python analyze_clonal_structure.py

# === 18. Publication-Ready Figure Layout ===
python generate_figure_layout.py

# === 19. Save Conda Environments ===
conda env export -n variant_calling > environment_variant_calling.yaml
conda env export -n samtools_env > environment_samtools.yaml

# === 20. Package Output ===
tar -czf chip_variants_analysis_results.tar.gz results/publication_ready *.py *.sh *.yaml
