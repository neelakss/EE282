## Homework 4  
### Summarize partitions of a genome assembly
We will be revisiting the Drosophila melanogaster genome. As with Homework 3, start at flybase.org. Go to the most current download genomes section and download the gzipped fasta file for all chromosomes.
> To download the data we can use wget or rsync:
<pre><code>mkdir ./HW4
cd ./HW4 
wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-chromosome-r6.24.fasta.gz 
gunzip *.gz 
</code></pre> 
Hint: The partitioning can be accomplished in many different ways. In my opinion, the easiest way is by using bioawk and faSize. The bioawk tool can be found in the module jje/jjeutils and the fa* utilities can be found in the module jje/kent.
#### Calculate the following for all sequences ≤ 100kb and all sequences > 100kb: 
Because the calculations will be for two genome partitions, there will be 6 total responses.
> To filter out the data awk was used and based on the size of the sequences two different .fa files were created:
<pre><code>zcat *.fasta* | awk 'BEGIN{RS=">";ORS=""}length($0)<=100000{print ">"$0}' > fa_lte100.fa
zcat *.fasta* | awk 'BEGIN{RS=">";ORS=""}length($0)>100000{print ">"$0}' > fa_mte100.fa
</code></pre>  
>1. Total number of nucleotides:
<pre><code>grep -v "^>" fa_mt100.fa | tr -d -C 'A\T\G\C\N' | wc -m
</code></pre>  
<pre><code>137547960
</code></pre> 
<pre><code>grep -v "^>" fa_lte100.fa | tr -d -C 'A\T\G\C\N' | wc -m
</code></pre>
<pre><code>6178042
</code></pre> 

>2. Total number of Ns:
<pre><code>grep -v "^>" fa_lte100.fa | tr -d -C 'N' | wc -m
</code></pre> 
<pre><code>662593
</code></pre> 
<pre><code>grep -v "^>" fa_mt100.fa | tr -d -C 'N' | wc -m
</code></pre> 
<pre><code>490385
</code></pre> 

>3. Total number of sequences:
<pre><code>grep -c "^>" fa_lte100.fa
</code></pre>
<pre><code>1863
</code></pre>
<pre><code>grep -c "^>" fa_mt100.fa
</code></pre>
<pre><code>7
</code></pre>

#### Plots of the following for the whole genome, for all sequences ≤ 100kb, and all sequences > 100kb:
Hint: bioawk has a function called gc(). Don't forget about the CDF plotting utility we used in class.Because the calculations will be for the whole genome and two genome partitions, there will be 9 total plots.
> In order to carry out the mentioned commands we download the modules loaded:
<pre><code>module load jje/jjeutils
module load emboss
module load rstudio/0.99.9.9
</code></pre>
>1. Sequence length distribution
<pre><code>infoseq -auto -nocolumns -delimiter ',' -only -noheading -name -length fa_lte100.fa > hist_lte.txt
infoseq -auto -nocolumns -delimiter ',' -only -noheading -name -length fa_mt100.fa > hist_mt.txt
infoseq -auto -nocolumns -delimiter ',' -only -noheading -name -length dmel-all-chromosome-r6.24.fasta > hist_wg.txt
</code></pre>
![Histogram of SLD for the sequence sizes less than and equal to 100kb](hist_lte.jpg)
![Histogram of SLD for the sequence sizes more than 100kb](hist_mt.jpg)
![Histogram of SLD for the whole genome](hist_wg.jpg)
>2. Sequence GC% distribution
<pre><code>bioawk -c fastx '{ print $name, gc($seq) }' fa_lte100.fa > gc_lte100.txt
bioawk -c fastx '{ print $name, gc($seq) }' fa_mt100.fa > gc_mt100.txt
bioawk -c fastx '{ print $name, gc($seq) }' dmel-all-chromosome-r6.24.fasta > gc_wg.txt
</code></pre>
![Histogram of GC% for the sequence sizes less than and equal to 100kb](gc_lte.jpg)
![Histogram of GC% for the sequence sizes more than 100kb](gc_mt.jpg)
![Histogram of GC% for the whole genome](gc_wg.jpg)
>3. Cumulative genome size sorted from largest to smallest sequences
<pre><code>bioawk -c fastx ' { print length($seq) } ' fa_lte100.fa | sort -rn | awk ' BEGIN { print "Assembly\tLength\nseq_length\t0" } { print "seq_length\t" $1 } ' > len_lte100.length
plotCDF2  len_lte100.length len_lte100.png
bioawk -c fastx ' { print length($seq) } ' fa_mt100.fa | sort -rn | awk ' BEGIN { print "Assembly\tLength\nseq_length\t0" } { print "seq_length\t" $1 } ' > len_mt100.length
plotCDF2  len_mt100.length len_mt100.png
bioawk -c fastx ' { print length($seq) } ' dmel-all-chromosome-r6.24.fasta | sort -rn | awk ' BEGIN { print "Assembly\tLength\nseq_length\t0" } { print "seq_length\t" $1 } ' > len_wg.length
plotCDF2  len_wg.length len_wg.png
</code></pre>
![Cummulative genome sizes for the sequence sizes less than and equal to 100kb](len_lte100.png)
![Cummulative genome sizes for the sequence sizes more than 100kb](len_mt100.png)
![Cummulative genome sizes for the whole genome](len_wg.png)
The code for the Rscript used is given below:
<pre><code>module load rstudio
rstudio
reads<- read.csv("hist_lte.txt",sep=",",header=F)
reads <- matrix(as.numeric(unlist(reads)))
jpeg('hist_lte.jpg')
hist(reads, breaks = 12, col = "cyan", border = "pink",main = "Histogram for Sequence Length Distribution (Size =< 100kb)")
reads<- read.csv("hist_mt.txt",sep=",",header=F)
reads <- matrix(as.numeric(unlist(reads)))
jpeg('hist_mt.jpg')
hist(reads, breaks = 12, col = "cyan", border = "pink",main = "Histogram for Sequence Length Distribution (Size > 100kb)")
reads<- read.csv("hist_wg.txt",sep=",",header=F)
reads <- matrix(as.numeric(unlist(reads)))
jpeg('hist_wg.jpg')
hist(reads, breaks = 12, col = "cyan", border = "pink",main = "Histogram for Sequence Length Distribution (Whole Genome)")
reads <- read.csv("gc_lte100.txt",sep=",",header=F)
req <- seq(0,nrow(reads),2)
reads <- reads[req,]
reads <- as.numeric(levels(reads)[reads])
jpeg('gc_lte.jpg')
hist(reads, breaks = 12, col = "cyan", border = "pink",main = "Histogram for GC% (Size =< 100kb)")
reads <- read.csv("gc_mt100.txt",sep=",",header=F)
req <- seq(0,nrow(reads),2)
reads <- reads[req,]
reads <- as.numeric(levels(reads)[reads])
jpeg('gc_mt.jpg')
hist(reads, breaks = 12, col = "cyan", border = "pink",main = "Histogram for GC% (Size > 100kb)")
reads <- read.csv("gc_wg.txt",sep=",",header=F)
req <- seq(0,nrow(reads),2)
reads <- reads[req,]
reads <- as.numeric(levels(reads)[reads])
jpeg('gc_wg.jpg')
hist(reads, breaks = 12, col = "cyan", border = "pink",main = "Histogram for GC% (Whole Genome)")
</code></pre>
### Genome assembly
#### Assemble a genome from MinION reads.  
Hint: Read up on miniasm here. We're using one of the simplest assembly approaches possible. This assembly can literally be accomplished with three lines of code. This will literally take only 3 command lines.

>1. Download the reads
<pre><code>wget https://hpc.oit.uci.edu/~solarese/ee282/iso1_onp_a2_1kb.fastq.gz
gunzip *.gz
ln -sf iso1_onp_a2_1kb.fastq reads.fq
</code></pre>
>2. Use minimap to overlap reads
<pre><code>qrsh -q epyc,abio128,free88i,free72i -pe openmp 32      
minimap -t 32 -Sw5 -L100 -m0 reads.fq{,} | gzip -1 > onp.paf.gz
</code></pre>
>3. Use miniasm to construct an assembly
<pre><code>miniasm -f reads.fq onp.paf.gz > reads.gfa
</code></pre>

#### Assembly assessment
Hint: For MUMmer, you should run nucmer, delta-filter, and mummerplot.
>1. Calculate the N50 of your assembly (this can be done with only faSize+awk+sort or with bioawk+awk+sort) and compare it to the Drosophila community reference's contig N50
<pre><code>awk ' $0 ~/^S/ { print ">" $2" \n" $3 } ' reads.gfa \
| tee >(n50 /dev/stdin > n50.txt) \
| fold -w 60 \
> unitigs.fa
n50 () {
  bioawk -c fastx ' { print length($seq); n=n+length($seq); } END { print n; } ' unitigs.fa \
  | sort -rn \
  | gawk ' NR == 1 { n = $1 }; NR > 1 { ni = $1 + ni; } ni/n > 0.5 { print $1; exit; } '
}
</code></pre>
<pre><code>4494246
</code></pre>
>2. Compare your assembly to the contig assembly (not the scaffold assembly!) from Drosophila melanogaster on FlyBase using a dotplot constructed with MUMmer (Hint: use faSplitByN as demonstrated in class)
<pre><code>faSplitByN dmel-all-chromosome-r6.24.fasta wg_contigassembly.fasta 10
source /pub/jje/ee282/bin/.qmbashrc
module load gnuplot/4.6.0
REF="wg_contigassembly.fasta"
PREFIX="flybase"
SGE_TASK_ID=1
QRY=$(ls u*.fa | head -n $SGE_TASK_ID | tail -n 1)
PREFIX=${PREFIX}_$(basename ${QRY} .fa)
nucmer -l 100 -c 150 -d 10 -banded -D 5 -prefix ${PREFIX} ${REF} ${QRY}
mummerplot --fat --layout --filter -p ${PREFIX} ${PREFIX}.delta -R ${REF} -Q ${QRY} --png
</code></pre>
![MUMmer plot for comparison](flybase_unitigs.png)
>3. Compare your assembly to both the contig assembly and the scaffold assembly from the Drosophila melanogaster on FlyBase using a contiguity plot (Hint: use plotCDF2 as demonstrated in class)
<pre><code>bioawk -c fastx '{ print length($seq) }' unitigs.fa | sort -rn | awk ' BEGIN { print "Assembly\tLength\nMini_Unitigs\t0" } { print "Mini_Unitigs\t" $1 } ' > mini_unitigs
bioawk -c fastx '{ print length($seq) }' wg_contigassembly.fasta | sort -rn | awk ' BEGIN { print "Assembly\tLength\nWG_Contigass\t0" } { print "WG_Contigass\t" $1 } ' > wg_contigass
bioawk -c fastx '{ print length($seq) }' dmel-all-chromosome-r6.24.fasta | sort -rn | awk ' BEGIN { print "Assembly\tLength\nWG_Scaffold\t0" } {print "WG_Scaffold\t" $1 } ' > wg_scaffold
plotCDF2 mini_unitigs wg_contigass wg_scaffold CDF_compare.png
</code></pre>
![CDF plot for comparison](CDF_compare.png)
>4. Calculate BUSCO scores of both assemblies and compare them
<pre><code>module load augustus/3.2.1
module load blast/2.2.31 
module load hmmer/3.1b2 
module load boost/1.54.0
source /pub/jje/ee282/bin/.buscorc
 BUSCO.py -c 32 -i unitigs.fa -m geno -o Busco_unitigs -l /pub/jje/ee282/bin/busco/lineages/diptera_odb9
</code></pre>
<pre><code>Results:
INFO    C:0.5%[S:0.5%,D:0.0%],F:1.1%,M:98.4%,n:2799
INFO    13      Complete BUSCOs (C)
INFO    13      Complete and single-copy BUSCOs (S)
INFO    0       Complete and duplicated BUSCOs (D)
INFO    32      Fragmented BUSCOs (F)
INFO    2754    Missing BUSCOs (M)
INFO    2799    Total BUSCO groups searched
</code></pre>
<pre><code>module load augustus/3.2.1
module load blast/2.2.31 
module load hmmer/3.1b2 
module load boost/1.54.0
source /pub/jje/ee282/bin/.buscorc
BUSCO.py -c 32 -i wg_contigassembly.fasta -m geno -o Busco_contig_assembly -l /pub/jje/ee282/bin/busco/lineages/diptera_odb9
</code></pre>
<pre><code>Results:
INFO    C:98.7%[S:98.2%,D:0.5%],F:0.8%,M:0.5%,n:2799
INFO    2763    Complete BUSCOs (C)
INFO    2749    Complete and single-copy BUSCOs (S)
INFO    14      Complete and duplicated BUSCOs (D)
INFO    21      Fragmented BUSCOs (F)
INFO    15      Missing BUSCOs (M)
INFO    2799    Total BUSCO groups searched
</code></pre>
