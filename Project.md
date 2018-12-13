# WGCNA ANALYSIS
This markdown file will provide a brief description on Weighted Correlation Network Analysis on high-dimensional data.
## STEP 1: Downloading the FASTQ files and Reference Genome.
> To download the .fastq or .fa files from the server, one could use filezilla from the ionode. 
<pre><code>qrsh -q ionode
mkdir ./fastq_files
filezilla
gunzip *.gz
</code></pre>
> The genome could be accessed from the ftp website of [ensembl](http://metazoa.ensembl.org/info/website/ftp/index.html) or [gencode](https://www.gencodegenes.org/) website or [ucsc](http://hgdownload.cse.ucsc.edu/downloads.html) genome website. Once downloaded you would want to concatenate the reference genome such that the chromosomes are in the correct order starting from chromosome 1.
<pre><code>mkdir ./star_genome
wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/*
gunzip *.gz
echo "$(ls chr*.fa | sort -V | grep -vP 'chr[^X|Y|\d]'; ls chr*.fa | sort -V | grep -vP 'chr[\d|X|Y]')" | xargs cat > refgenome.fa
</code></pre>
## STEP2: Indexing the genome for STAR
> Indexing is the first step for STAR alignment.
<pre><code>#! /bin/bash/
#$ -N Star_genome_generation
#$ -pe openmp 30
#$ -q bio,epyc,pub8i

module load STAR/2.5.2a
STAR --runMode genomeGenerate --genomeDir /bio/neelakss/star_genome --genomeFastaFiles /bio/neelakss/star_genome/refgenome.fa --runThreadN 24
</code></pre>
### STEP 3: Alignment of fastq files to the refrence genome
> Paired end alignemnt is shown in the code below, however one could refer to the manual. Changing the fastq files in the given code for each sample will give you the alignment for all.
<pre><code>#! /bin/bash/
#$ -N Star_alignment_A1
#$ -pe openmp 10
#$ -q bio,epyc,pub8i

module load STAR/2.5.2a
STAR --genomeDir /bio/neelakss/star-genome  --runThreadN 24 --readFilesIn /bio/neelakss/fastq_files/A1_S1_R1.fastq /bio/neelakss/fastq_files/A1_S1_R2.fastq --outFileNamePrefix A1 --outSAMtype BAM Unsorted SortedByCoordinate
</code></pre>
### STEP 4: Picard Matrix
> As the data is RNA-Seq we need to perform a linear regession on the TPM matrix using the PC values calculated from the alignmnent metrics run on the bam files.
<pre><code>#! /bin/bash/
#$ -N BuildBai_A1
#$ -pe openmp 6
#$ -q epyc,bio,pub8i

module load picard-tools/1.130

java -jar /data/apps/picard-tools/1.130/picard.jar BuildBamIndex I=/bio/neelakss/picard_tools/buildbai/A1Aligned.sortedByCoord.out.bam
</code></pre>
<pre><code>#! /bin/bash/
#$ -N PCA_Req_Data_A1
#$ -pe openmp 6
#$ -q epyc,bio,pub8i

module load picard-tools/1.130

java -jar /data/apps/picard-tools/1.130/picard.jar CollectAlignmentSummaryMetrics\
 R=refgenome.fa\
 I=/pub/neelakss/picard_tools/buildbai/A1Aligned.sortedByCoord.out.bam\
 O=output.Alignment_Summary_A1.txt

java -jar /data/apps/picard-tools/1.130/picard.jar CollectRnaSeqMetrics \
      I=/pub/neelakss/picard_tools/buildbai/A1Aligned.sortedByCoord.out.bam \
      O=output_A1.RNA_Metrics \
      REF_FLAT=refFlat.txt \
      STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
      RIBOSOMAL_INTERVALS=ribosomal.interval_list

java -jar /data/apps/picard-tools/1.130/picard.jar CollectGcBiasMetrics \
      I=/pub/neelakss/picard_tools/buildbai/A1Aligned.sortedByCoord.out.bam \
      O=gc_bias_metrics_A1.txt \
      CHART=gc_bias_metrics_A1.pdf \
      S=summary_metrics_A1.txt \
      R=refgenome.fa

java -jar /data/apps/picard-tools/1.130/picard.jar MarkDuplicates \
      I=/pub/neelakss/picard_tools/buildbai/A1Aligned.sortedByCoord.out.bam \
      O=marked_duplicates_A1.bam \
      M=marked_dup_metrics_A1.txt

java -jar /data/apps/picard-tools/1.130/picard.jar AddOrReplaceReadGroups \
      I=/pub/neelakss/picard_tools/buildbai/A1Aligned.sortedByCoord.out.bam \
      O=A1Aligned.sortedByCoord.out.read.bam \
      RGID=HMCVCBGX5.11101.1 \
      RGLB=lib1 \
      RGPL=illumina \
      RGPU=TAAGGC+TCTTAC \
      RGSM=Wildtype1

java -jar /data/apps/picard-tools/1.130/picard.jar CollectGcBiasMetrics \
      I=A1Aligned.sortedByCoord.out.read.bam\
      O=gc_bias_metrics_A1.txt \
      CHART=gc_bias_metrics_A1.pdf \
      S=summary_metrics_A1.txt \
      R=refgenome.fa
</pre></code>
<pre><code>#! /bin/bash/
#$ -N reading_metrics
#$ -pe openmp 6
#$ -q epyc,bio,pub8i


awk ' FNR== 7 ' summary_metrics_A1.txt > GC_Header.txt

awk ' FNR== 8 ' summary_metrics_A1.txt summary_metrics_A2.txt \
summary_metrics_A3.txt summary_metrics_A4.txt \
summary_metrics_A5.txt summary_metrics_A6.txt \
summary_metrics_B1.txt summary_metrics_B2.txt \
summary_metrics_B3.txt summary_metrics_B4.txt \
summary_metrics_B5.txt summary_metrics_B6.txt \
summary_metrics_C1.txt summary_metrics_C2.txt \
summary_metrics_C3.txt summary_metrics_C4.txt \
summary_metrics_C5.txt summary_metrics_C6.txt \
summary_metrics_D1.txt summary_metrics_D2.txt \
summary_metrics_D3.txt summary_metrics_D4.txt \
summary_metrics_D5.txt summary_metrics_D6.txt > GC_Metrics.txt

awk ' FNR== 7 ' output_A1.RNA_Metrics > RNA_Header.txt

awk ' FNR== 8 ' output_A1.RNA_Metrics output_A2.RNA_Metrics \
output_A3.RNA_Metrics output_A4.RNA_Metrics \
output_A5.RNA_Metrics output_A6.RNA_Metrics \
output_B1.RNA_Metrics output_B2.RNA_Metrics \
output_B3.RNA_Metrics output_B4.RNA_Metrics \
output_B5.RNA_Metrics output_B6.RNA_Metrics \
output_C1.RNA_Metrics output_C2.RNA_Metrics \
output_C3.RNA_Metrics output_C4.RNA_Metrics \
output_C5.RNA_Metrics output_C6.RNA_Metrics \
output_D1.RNA_Metrics output_D2.RNA_Metrics \
output_D3.RNA_Metrics output_D4.RNA_Metrics \
output_D5.RNA_Metrics output_D6.RNA_Metrics > RNA_Metrics.txt

awk ' FNR== 7 ' marked_dup_metrics_A1.txt > MDUP_Header.txt

awk ' FNR== 8 ' marked_dup_metrics_A1.txt marked_dup_metrics_A2.txt \
marked_dup_metrics_A3.txt marked_dup_metrics_A4.txt \
marked_dup_metrics_A5.txt marked_dup_metrics_A6.txt \
marked_dup_metrics_B1.txt marked_dup_metrics_B2.txt \
marked_dup_metrics_B3.txt marked_dup_metrics_B4.txt \
marked_dup_metrics_B5.txt marked_dup_metrics_B6.txt \
marked_dup_metrics_C1.txt marked_dup_metrics_C2.txt \
marked_dup_metrics_C3.txt marked_dup_metrics_C4.txt \
marked_dup_metrics_C5.txt marked_dup_metrics_C6.txt \
marked_dup_metrics_D1.txt marked_dup_metrics_D2.txt \
marked_dup_metrics_D3.txt marked_dup_metrics_D4.txt \
marked_dup_metrics_D5.txt marked_dup_metrics_D6.txt > MDUP_Metrics.txt

awk ' FNR== 7 ' output.Alignment_Summary_A1.txt > ALIGSUM_Header.txt

awk ' FNR== 8 ' output.Alignment_Summary_A1.txt output.Alignment_Summary_A2.txt \
output.Alignment_Summary_A3.txt output.Alignment_Summary_A4.txt \
output.Alignment_Summary_A5.txt output.Alignment_Summary_A6.txt \
output.Alignment_Summary_B1.txt output.Alignment_Summary_B2.txt \
output.Alignment_Summary_B3.txt output.Alignment_Summary_B4.txt \
output.Alignment_Summary_B5.txt output.Alignment_Summary_B6.txt \
output.Alignment_Summary_C1.txt output.Alignment_Summary_C2.txt \
output.Alignment_Summary_C3.txt output.Alignment_Summary_C4.txt \
output.Alignment_Summary_C5.txt output.Alignment_Summary_C6.txt \
output.Alignment_Summary_D1.txt output.Alignment_Summary_D2.txt \
output.Alignment_Summary_D3.txt output.Alignment_Summary_D4.txt \
output.Alignment_Summary_D5.txt output.Alignment_Summary_D6.txt > AlIGSUM_Metrics.txt

paste GC_Header.txt RNA_Header.txt MDUP_Header.txt ALIGSUM_Header.txt > Header.txt

paste GC_Metrics.txt RNA_Metrics.txt MDUP_Metrics.txt AlIGSUM_Metrics.txt > Metrics.txt

cat Header.txt Metrics.txt > Full_Metrics.txt

awk -F, '{$1=++i FS $1;}1' OFS=, Full_Metics.txt
</code></pre>
