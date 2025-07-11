#!/bin/bash

# ==============================
# Variant Calling Pipeline (Somatic Panel Data)
# Using: GATK Mutect2 + LoFreq
# TMUX-ready script
# ==============================

# === CONFIG ===
SAMPLE="SRR32229173"
REF="ref/hg38.fa"
FASTQ1="data/${SAMPLE}_1.fastq"
FASTQ2="data/${SAMPLE}_2.fastq"
THREADS=8
OUTDIR="results"

mkdir -p $OUTDIR logs

# === 1. Alignment ===
echo "🧬 [1] Aligning reads with BWA..."
bwa mem -t $THREADS $REF $FASTQ1 $FASTQ2 2> logs/bwa.err | \
    samtools view -bS - > $OUTDIR/${SAMPLE}.bam

# === 2. Sort BAM ===
echo "📂 [2] Sorting BAM..."
samtools sort -@ $THREADS -o $OUTDIR/${SAMPLE}.sorted.bam $OUTDIR/${SAMPLE}.bam
samtools index $OUTDIR/${SAMPLE}.sorted.bam

# === 3. Mark Duplicates ===
echo "🔁 [3] Marking Duplicates..."
picard MarkDuplicates \
    I=$OUTDIR/${SAMPLE}.sorted.bam \
    O=$OUTDIR/${SAMPLE}.dedup.bam \
    M=$OUTDIR/${SAMPLE}.dup_metrics.txt \
    CREATE_INDEX=true \
    2> logs/picard.err

# === 4. GATK Mutect2 Variant Calling ===
echo "🧪 [4] GATK Mutect2 Variant Calling..."
gatk Mutect2 \
    -R $REF \
    -I $OUTDIR/${SAMPLE}.dedup.bam \
    -O $OUTDIR/${SAMPLE}.mutect2.unfiltered.vcf \
    --max-mnp-distance 0 \
    2> logs/mutect2.err

gatk FilterMutectCalls \
    -V $OUTDIR/${SAMPLE}.mutect2.unfiltered.vcf \
    -O $OUTDIR/${SAMPLE}.mutect2.filtered.vcf \
    2> logs/filter.err

# === 5. LoFreq Variant Calling ===
echo "⚡ [5] LoFreq Variant Calling..."
# Indel quality addition
lofreq indelqual \
    --dindel \
    -f $REF \
    -o $OUTDIR/${SAMPLE}.indelqual.bam \
    $OUTDIR/${SAMPLE}.dedup.bam \
    2> logs/lofreq_indel.err

samtools index $OUTDIR/${SAMPLE}.indelqual.bam

# Variant calling
lofreq call-parallel --pp-threads $THREADS \
    -f $REF \
    -o $OUTDIR/${SAMPLE}.lofreq.vcf \
    $OUTDIR/${SAMPLE}.indelqual.bam \
    2> logs/lofreq_call.err

echo "✅ Pipeline completed. Results in $OUTDIR"
