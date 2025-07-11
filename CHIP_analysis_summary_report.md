
# 🧬 CHIP Variant Analysis Report (SRR32229173)

**Project Title**: Somatic Variant Analysis of Clonal Hematopoiesis Panel Data  
**Sample ID**: SRR32229173  
**Source**: NCBI SRA / GSM8774288  
**Institution**: Mayo Clinic  
**Panel**: Invitae VariantPlex® 75-gene Myeloid Panel  
**Reference Genome**: GRCh38 (hg38)

---

## 🔬 Page 1: Overview & Objectives

### 🎯 Project Objective
This study reanalyzes public NGS data from the SRA to characterize somatic variants associated with **clonal hematopoiesis of indeterminate potential (CHIP)**. Using a custom bioinformatics pipeline, we identify and annotate high-confidence somatic mutations, focusing on **CHIP driver genes**, their allele frequencies, and clinical relevance.

### 📥 Input Data
| Data Type | Description |
|-----------|-------------|
| FASTQ     | `SRR32229173_1.fastq`, `SRR32229173_2.fastq` |
| Genome    | `hg38.fa` |
| Target BED| `covered_regions.bed` from panel design |
| Tools     | BWA, Samtools, Picard, GATK, LoFreq, SnpEff, ClinVar, bcftools |

### 🔁 Pipeline Summary
1. **Read Alignment**: BWA-MEM  
2. **BAM Processing**: Sort, MarkDuplicates  
3. **Variant Calling**: GATK Mutect2 and LoFreq  
4. **Filtering**: High-confidence based on DP, AF  
5. **Annotation**: SnpEff (functional), ClinVar (clinical)  
6. **CHIP Focus**: Extracted mutations in known CHIP genes  
7. **Clonal Analysis**: Binned VAFs to infer dominant vs minor clones  
8. **Visualization**: Functional barplots, AF distribution, Venn, CLNSIG  
9. **Packaging**: Summary reports, scripts, environment files

---

## 📈 Page 2: Results Summary

### 🧬 Variant Calling Results

| Tool     | Total Variants | SNPs   | INDELs |
|----------|----------------|--------|--------|
| Mutect2  | 759,282        | 459,755| 299,527|
| LoFreq   | 1,282          | 1,282  | 0      |

### 🧪 Functional Annotation
- High-impact variants: `missense`, `nonsense`, `frameshift`, `splice site`
- CHIP gene hits: `TET2`, `DNMT3A`, `ASXL1`, etc.

### 🧬 ClinVar Clinical Significance (CLNSIG)
| Significance         | Count |
|----------------------|-------|
| Pathogenic           | 54    |
| Likely_pathogenic    | 32    |
| Benign/Likely benign | 58    |
| VUS/Uncertain        | 112   |

### 🧠 VAF Clonal Binning
| Bin     | Clone Type     | # Variants (CHIP genes) |
|---------|----------------|--------------------------|
| <5%     | Minor clones   | 88                       |
| 5–10%   | Emerging clones| 65                       |
| 10–20%  | Intermediate   | 38                       |
| >20%    | Dominant clones| 19                       |

---

## 🧪 Page 3: Interpretation & Reproducibility

### 🧠 Interpretation
This analysis reveals a complex **clonal landscape**, with multiple **CHIP driver gene mutations** found at varying allele frequencies. The presence of **TET2** and **DNMT3A** mutations at high VAF suggests dominant clones — a hallmark of early clonal expansion and potential progression to hematologic disorders.

### 📊 Highlights
- Dual-variant calling strategy (Mutect2 + LoFreq)
- Cross-tool concordance + overlap analysis
- Comprehensive ClinVar + SnpEff annotation
- Clonal burden inference by AF
- Automated report + figure generation

### 📦 Reproducibility
| Asset | File |
|-------|------|
| Conda Envs | `environment_variant_calling.yaml`, `environment_samtools.yaml` |
| Pipeline   | `variant_calling_tmux.sh`, `post_variant_analysis.sh` |
| Figures    | `figures_summary.pdf`, `vaf_clone_bins.png` |
| Data Files | Annotated VCFs, TSVs, barplots |

### 📌 Data Use & Ethics
The raw sequencing data used in this analysis (SRR32229173) was obtained from the **NCBI SRA** and was submitted by the **Mayo Clinic**. It is publicly available and reused here under academic fair use with appropriate acknowledgment. The analysis provides novel interpretation and computational processing.

> Citation: Chandra et al., Mayo Clinic, PRJNA1219176 / GSM8774288

---

### ✅ Final Deliverables:
All processed outputs are packaged in:  
📦 `chip_variants_analysis_results.tar.gz`

