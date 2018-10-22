## Homework 3  
### Summarize a genome assembly
#### Getting the data
We will be working with the Drosophila melanogaster genome. You can start at flybase.org. Go to the most current download genomes section and download the gzipped fasta file for all chromosomes.
>In order to download the data from the given website, we can use wget or curl. 
<pre><code>
$ wget -r -np -N ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.24_FB2018_05/fasta/md5sum.txt & wget -r -np -N ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.24_FB2018_05/fasta/dmel-all-chromosome-r6.24.fasta.gz
</code></pre> 
<pre><code>
--2018-10-22 13:06:53--  ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.24_FB2018_05/fasta/dmel-all-chromosome-r6.24.fasta.gz
           => "ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.24_FB2018_05/fasta/dmel-all-chromosome-r6.24.fasta.gz"
==> CWD not required.
==> PASV ... unlink: No such file or directory
--2018-10-22 13:06:53--  ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.24_FB2018_05/fasta/md5sum.txt
           => "ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.24_FB2018_05/fasta/md5sum.txt"
==> CWD not required.
==> PASV ... done.    ==> RETR dmel-all-chromosome-r6.24.fasta.gz ... done.    ==> RETR md5sum.txt ... done.
Length: 42465980 (40M)

 0% [                                   ] 0           --.-K/s              done.
Length: 1489 (1.5K)

100%[==================================>] 1,489       --.-K/s   in 0s      

2018-10-22 13:06:53 (241 MB/s) - "ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.24_FB2018_05/fasta/md5sum.txt" saved [1489]
</code></pre>
#### File integrity
> Verify the file integrity of the gzipped fasta file using a checksum
<pre><code>
$ md5sum -c md5sum.txt
</code></pre>
<pre><code>
md5sum: dmel-all-aligned-r6.24.fasta.gz: No such file or directory
dmel-all-aligned-r6.24.fasta.gz: FAILED open or read
md5sum: dmel-all-CDS-r6.24.fasta.gz: No such file or directory
dmel-all-CDS-r6.24.fasta.gz: FAILED open or read
dmel-all-chromosome-r6.24.fasta.gz: OK
md5sum: dmel-all-clones-r6.24.fasta.gz: No such file or directory
dmel-all-clones-r6.24.fasta.gz: FAILED open or read
md5sum: dmel-all-exon-r6.24.fasta.gz: No such file or directory
dmel-all-exon-r6.24.fasta.gz: FAILED open or read
md5sum: dmel-all-five_prime_UTR-r6.24.fasta.gz: No such file or directory
dmel-all-five_prime_UTR-r6.24.fasta.gz: FAILED open or read
md5sum: dmel-all-gene_extended2000-r6.24.fasta.gz: No such file or directory
dmel-all-gene_extended2000-r6.24.fasta.gz: FAILED open or read
md5sum: dmel-all-gene-r6.24.fasta.gz: No such file or directory
dmel-all-gene-r6.24.fasta.gz: FAILED open or read
md5sum: dmel-all-intergenic-r6.24.fasta.gz: No such file or directory
dmel-all-intergenic-r6.24.fasta.gz: FAILED open or read
md5sum: dmel-all-intron-r6.24.fasta.gz: No such file or directory
dmel-all-intron-r6.24.fasta.gz: FAILED open or read
md5sum: dmel-all-miRNA-r6.24.fasta.gz: No such file or directory
dmel-all-miRNA-r6.24.fasta.gz: FAILED open or read
md5sum: dmel-all-miscRNA-r6.24.fasta.gz: No such file or directory
dmel-all-miscRNA-r6.24.fasta.gz: FAILED open or read
md5sum: dmel-all-ncRNA-r6.24.fasta.gz: No such file or directory
dmel-all-ncRNA-r6.24.fasta.gz: FAILED open or read
md5sum: dmel-all-predicted-r6.24.fasta.gz: No such file or directory
dmel-all-predicted-r6.24.fasta.gz: FAILED open or read
md5sum: dmel-all-pseudogene-r6.24.fasta.gz: No such file or directory
dmel-all-pseudogene-r6.24.fasta.gz: FAILED open or read
md5sum: dmel-all-sequence_features-r6.24.fasta.gz: No such file or directory
dmel-all-sequence_features-r6.24.fasta.gz: FAILED open or read
md5sum: dmel-all-synteny-r6.24.fasta.gz: No such file or directory
dmel-all-synteny-r6.24.fasta.gz: FAILED open or read
md5sum: dmel-all-three_prime_UTR-r6.24.fasta.gz: No such file or directory
dmel-all-three_prime_UTR-r6.24.fasta.gz: FAILED open or read
md5sum: dmel-all-transcript-r6.24.fasta.gz: No such file or directory
dmel-all-transcript-r6.24.fasta.gz: FAILED open or read
md5sum: dmel-all-translation-r6.24.fasta.gz: No such file or directory
dmel-all-translation-r6.24.fasta.gz: FAILED open or read
md5sum: dmel-all-transposon-r6.24.fasta.gz: No such file or directory
dmel-all-transposon-r6.24.fasta.gz: FAILED open or read
md5sum: dmel-all-tRNA-r6.24.fasta.gz: No such file or directory
dmel-all-tRNA-r6.24.fasta.gz: FAILED open or read
md5sum: WARNING: 21 of 22 listed files could not be read
</code></pre>

#### Calculate the following for the whole genome
> 1. Total number of nucleotides
<pre><code>
$ grep -c 'A\|T\|G\|C\|N' *.gz
</code></pre> 
<pre><code>
130659
</code></pre> 
> 2. Total number of Ns
<pre><code>
$ grep -c 'N' *.gz
</code></pre> 
<pre><code>
81114
</code></pre> 
