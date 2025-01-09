# Download Raw FATSQs

## This step downloads all FASTQ files from the Nemoarchive server.
## Since the FASTQ files are split by cell (each cell has two files: Read 1 and Read 2) and packaged into tar files,
## we need to unpack them and perform alignment for each cell individually.

wget -i Smartseq2_files.txt

## FOr each cell
tar xvf xxx.tar

# Alignment

## We use STAR to align all reads to the mm10 reference genome.
## Below is an example command for processing a single cell.
/home/projects/ku_00009/data/TranscriptScar/pipeline/cellranger-6.1.2/lib/bin/STAR --genomeDir /home/projects/ku_00009/data/TranscriptScar/pipeline/refdata-gex-mm10-2020-A/star/ --readFilesIn SM-GE935_S384_E1-50_R1.fastq.gz SM-GE935_S384_E1-50_R2.fastq.gz --readFilesCommand zcat --outSAMattrRGline ID:SM-GE935_S384 --outSAMtype BAM SortedByCoordinate

# Merge all BAMs into one
## After alignment, the BAM files are sorted by coordinates.
## Then, we merge all BAM files into a single file for easier annotation and downstream visualization.

## we need first index the BAM files one by one
sambamba index in.bam

## In practice, the merging process is divided into multiple steps.
## Approximately 100 cells are merged into a temporary file, and then these temporary files are merged together.

sambamba merge merged.bam in1.bam in2.bam ..

# Annotation

## -is is used for ignore the strand information during annotation, because Smartseq2 is not strand sensitive
## -exon is used because we want annotate the exon and junction information rather than gene and transcript only
## -psi is used to annotate exon exluded reads for fig03, benchmark part
PISA anno -gtf /home/projects/ku_00009/data/TranscriptScar/pipeline/refdata-gex-mm10-2020-A/genes/genes.gtf -is -psi -exon -o anno.bam merged.bam

# Feature Counting
mkdir -p exp
mkdir -p exon
mkdir -p junction
mkdir -p exclude

## generate gene counts
PISA count -tag RG -outdir exp -anno-tag GN anno.bam

## generate exon counts
PISA count -tag RG -outdir exon -anno-tag EX anno.bam

## generate junction counts
PISA count -tag RG -outdir junction -anno-tag JC anno.bam

## generate exlude reads counts
PISA count -tag RG -outdir exclude -anno-tag ER anno.bam
