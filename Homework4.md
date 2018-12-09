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

>2. Sequence GC% distribution
<pre><code>bioawk -c fastx '{ print $name, gc($seq) }' fa_lte100.fa > gc_lte100.txt
bioawk -c fastx '{ print $name, gc($seq) }' fa_mt100.fa > gc_mt100.txt
bioawk -c fastx '{ print $name, gc($seq) }' dmel-all-chromosome-r6.24.fasta > gc_wg.txt
</code></pre>

>3. Cumulative genome size sorted from largest to smallest sequences
<pre><code>bioawk -c fastx ' { print length($seq) } ' fa_lte100.fa | sort -rn | awk ' BEGIN { print "Assembly\tLength\nseq_length\t0" } { print "seq_length\t" $1 } ' > len_lte100.length
bioawk -c fastx ' { print length($seq) } ' fa_mt100.fa | sort -rn | awk ' BEGIN { print "Assembly\tLength\nseq_length\t0" } { print "seq_length\t" $1 } ' > len_mt100.length
bioawk -c fastx ' { print length($seq) } ' dmel-all-chromosome-r6.24.fasta | sort -rn | awk ' BEGIN { print "Assembly\tLength\nseq_length\t0" } { print "seq_length\t" $1 } ' > len_wg.length
</code></pre>
