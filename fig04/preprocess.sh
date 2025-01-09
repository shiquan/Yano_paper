# scRNA-seq BAM files download from GEO websites with accession: SRR11832853 (P8_cSCC_scRNA_1), SRR11832854 (P8_normal_scRNA_1), SRR11832855 (P8_cSCC_scRNA_2), SRR11832856 (P8_normal_scRNA_1)

# WES BAM files download from GEO with accession: SRR11870482 (P8_normal_WES), SRR11870481 (P8_cSCC_WES)
# Index bam files
# Here, we use a single file for demonstration purposes.
sambamba index SRR11832853.bam

# call genetic varints with bcftools
bcftools mpileup -Ou -f /home/projects/ku_00009/data/TranscriptScar/pipeline/ensemble_genome/Homo_sapiens.GRCh37.dna.toplevel.fa SRR11832853.bam | bcftools call -vmO z -o in1.vcf.gz
bcftools index in1.vcf.gz

# For WES data, we start from fastq files, follow protocol https://www.htslib.org/workflow/wgs-call.html
bwa index /home/projects/ku_00009/data/TranscriptScar/pipeline/ensemble_genome/Homo_sapiens.GRCh37.dna.toplevel.fa
bwa mem -R -R '@RG\tID:P8_cSCC_WES\tSM:P8_cSCC_WES' /home/projects/ku_00009/data/TranscriptScar/pipeline/ensemble_genome/Homo_sapiens.GRCh37.dna.toplevel.fa SRR11870481_1.fastq SRR11870481_2.fastq > aln1.sam
bwa mem -R -R '@RG\tID:P8_normal_WES\tSM:P8_normal_WES' /home/projects/ku_00009/data/TranscriptScar/pipeline/ensemble_genome/Homo_sapiens.GRCh37.dna.toplevel.fa SRR11870482_1.fastq SRR11870482_2.fastq > aln2.sam

## Note: The following step process for each sample seperately
#
#
samtools fixmate -O bam aln.sam fixmated.bam
samtools sort -O bam -o sorted.bam fixmated.bam
samtools rmdup sorted.bam rmdup.bam
samtools index rmdup.bam
bcftools mpileup -Ou -f /home/projects/ku_00009/data/TranscriptScar/pipeline/ensemble_genome/Homo_sapiens.GRCh37.dna.toplevel.fa rmdup.bam | bcftools call -vmO z -o WES_in1.vcf.gz
bcftools index WES_in1.vcf.gz

## Merge all variants into one VCF file
bcftools merge --force-samples -O z -o var.vcf.gz in1.vcf.gz ...
bcftools norm -f /home/projects/ku_00009/data/TranscriptScar/pipeline/ensemble_genome/Homo_sapiens.GRCh37.dna.toplevel.fa var.vcf.gz -O z -o norm.vcf.gz


# Annotate features
## Annotion file download from ftp://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz
PISA -gtf Homo_sapiens.GRCh37.87.chr.gtf.gz -vcf norm.vcf.gz -ref-alt SRR11832853.bam -o anno.bam

# Count features
mkdir -p exp
mkdir -p var
PISA count -tag CB -umi UB -anno-tag GN -outdir exp anno.bam
PISA count -tag CB -umi UB -anno-tag VF -outdir var anno.bam


# Megre WES bam files, use to distinguish normal and tumor only genetic variants
samtools merge merge.bam cSCC.bam normal.bam
PISA -vcf norm.vcf.gz -ref-alt -is merge.bam -o WES_anno.bam
mkdir -p WES_var
PISA count -tag RG -anno-tag VF -outdir WES_var WES_anno.bam
