# Data download from GSE159929 and reanalysis with CellRanger 6.1.2 with code 'SCP5P-PE'

# Call genetic varints for each library  and then merge all vcf together
bcftools mpileup -Ou -f /home/projects/ku_00009/data/TranscriptScar/pipeline/refdata-gex-GRCh38-2020-A/fasta/genome.fa possorted_genome_bam.bam | bcftools call -vmO z -o in1.vcf.gz
bcftools index in1.vcf.gz

## Merge all variants into one VCF file
bcftools merge --force-samples -O z -o var.vcf.gz in1.vcf.gz ...
bcftools norm -f  /home/projects/ku_00009/data/TranscriptScar/pipeline/refdata-gex-GRCh38-2020-A/fasta/genome.fa var.vcf.gz -O z -o norm.vcf.gz

# Annotate features
PISA -gtf /home/projects/ku_00009/data/TranscriptScar/pipeline/refdata-gex-GRCh38-2020-A/genes/genes.gtf -vcf norm.vcf.gz -ref-alt possorted_genome_bam.bam -o anno.bam

# Count features
mkdir -p exp
mkdir -p var
PISA count -tag CB -umi UB -anno-tag GN -outdir exp anno.bam
PISA count -tag CB -umi UB -anno-tag VF -outdir var anno.bam

