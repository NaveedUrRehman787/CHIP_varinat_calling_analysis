import pandas as pd

# Load annotated VCF with ClinVar
clinvar_vcf = "results/publication_ready/SRR32229173.mutect2.clinvar.vcf"
topfx_vcf = "results/publication_ready/SRR32229173.mutect2.top_functional.vcf"

def extract_core_info(vcf):
    records = []
    with open(vcf) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split('\t')
            chrom, pos, ref, alt, info = parts[0], parts[1], parts[3], parts[4], parts[7]
            clnsig = "-"
            if "CLNSIG=" in info:
                clnsig = info.split("CLNSIG=")[1].split(";")[0]
            records.append((chrom, pos, ref, alt, clnsig))
    df = pd.DataFrame(records, columns=["CHROM", "POS", "REF", "ALT", "CLNSIG"])
    return df

def extract_topfx(vcf):
    records = []
    with open(vcf) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split('\t')
            chrom, pos, ref, alt, info = parts[0], parts[1], parts[3], parts[4], parts[7]
            effect = "-"
            if "ANN=" in info:
                ann = info.split("ANN=")[1].split(";")[0]
                effect = ann.split("|")[1]
            records.append((chrom, pos, ref, alt, effect))
    df = pd.DataFrame(records, columns=["CHROM", "POS", "REF", "ALT", "Effect"])
    return df

df_clinvar = extract_core_info(clinvar_vcf)
df_topfx = extract_topfx(topfx_vcf)

# Merge
merged = pd.merge(df_topfx, df_clinvar, on=["CHROM", "POS", "REF", "ALT"], how="inner")
merged.to_csv("results/publication_ready/SRR32229173.top_functional.clinvar_merged.tsv", sep="\t", index=False)

print("✅ Merged file saved as: SRR32229173.top_functional.clinvar_merged.tsv")
