# Download raw data

## Since the original study provided the BAM files, we only need to download them from the server.
wget -c -i 10Xv2_files.txt

## unpack the BAM files one by one
tar xvf xxx.tar

## And index these BAM files for downstream visulization, because Yano::TrackPlot requires the input BAMs are sorted and indexed
## Here we have 10 BAM files in total, we need index and annotate them one by one, I only give an example for one file.
sambamba index in.bam

# Annotation
## PISA anno enable the strand sensitive annotation in default.
## -exon is used because we want annotate the exon and junction information rather than gene and transcript only
## -psi is used to annotate exon exluded reads for fig03, benchmark part
PISA anno -gtf /home/projects/ku_00009/data/TranscriptScar/pipeline/refdata-gex-mm10-2020-A/genes/genes.gtf -psi -exon -o anno.bam merged.bam

# Feature counting
mkdir -p exp
mkdir -p exon
mkdir -p junction
mkdir -p exclude

## generate gene counts
PISA count -tag CB -umi UB -anno-tag GN -outdir exp anno.bam
## generate exon counts
PISA count -tag CB -umi UB -anno-tag EX -outdir exon anno.bam
## generate junction counts
PISA count -tag CB -umi UB -anno-tag JC -outdir junction anno.bam
## generate exon exlude counts
PISA count -tag CB -umi UB -anno-tag ER -outdir exclude anno.bam
