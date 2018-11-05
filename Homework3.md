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
$ zgrep -v  "^>" *.fasta*| tr -d -C 'A\T\G\C\N' | wc -m
</code></pre> 
<pre><code>
143726002
</code></pre> 
> 2. Total number of Ns
<pre><code>
$ zgrep -v  "^>" *.fasta* | tr -d -C 'N' | wc -m
</code></pre> 
<pre><code>
1152978
</code></pre> 
> 3. Total number of Sequences
<pre><code>
$zgrep -c "^>" *.fasta*
</code></pre> 
<pre><code>
1870
</code></pre>

### Summarize an annotation file
#### Getting the data
Go to the most current download genomes section at flybase.org and download the gzipped gtf annotation file for D. melanogaster.
>In order to download the data from the given website, we can use wget or curl. 
<pre><code>
$ wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.24_FB2018_05/gtf/*
</code></pre> 
<pre><code>
--2018-11-02 13:09:55--  ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.24_FB2018_05/gtf/*
           => ".listing"
Resolving ftp.flybase.net... 52.23.126.124
Connecting to ftp.flybase.net|52.23.126.124|:21... connected.
Logging in as anonymous ... Logged in!
==> SYST ... done.    ==> PWD ... done.
==> TYPE I ... done.  ==> CWD (1) /genomes/Drosophila_melanogaster/dmel_r6.24_FB2018_05/gtf ... done.
==> PASV ... done.    ==> LIST ... done.

    [ <=>                                   ] 266         --.-K/s   in 0s      

2018-11-02 13:09:56 (45.4 MB/s) - ".listing" saved [266]

Removed ".listing".
--2018-11-02 13:09:56--  ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.24_FB2018_05/gtf/dmel-all-r6.24.gtf.gz
           => "dmel-all-r6.24.gtf.gz"
==> CWD not required.
==> PASV ... done.    ==> RETR dmel-all-r6.24.gtf.gz ... done.
Length: 3897419 (3.7M)

100%[======================================>] 3,897,419   2.66M/s   in 1.4s    

2018-11-02 13:09:58 (2.66 MB/s) - "dmel-all-r6.24.gtf.gz" saved [3897419]

--2018-11-02 13:09:58--  ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.24_FB2018_05/gtf/md5sum.txt
           => "md5sum.txt.1"
==> CWD not required.
==> PASV ... done.    ==> RETR md5sum.txt ... done.
Length: 56

100%[======================================>] 56          --.-K/s   in 0s      

2018-11-02 13:09:58 (8.34 MB/s) - "md5sum.txt.1" saved [56]
</code></pre>
#### File integrity
> Verify the file integrity of the gzipped gtf annotation using a checksum
<pre><code>
$ md5sum -c md5sum.txt.1
</code></pre>
<pre><code>
dmel-all-r6.24.gtf.gz: OK
</code></pre>

#### Print a summary report with the following information:
> 1. Total number of features of each type, sorted from the most common to the least common
<pre><code>
$ zcat *.gtf* | cut -f3 | sort | uniq -c | sort -nr | nl
</code></pre> 
<pre><code>
     1	 187315 exon
     2	 161014 CDS
     3	  46339 5UTR
     4	  33358 3UTR
     5	  30591 start_codon
     6	  30533 stop_codon
     7	  30507 mRNA
     8	  17772 gene
     9	   2961 ncRNA
    10	    485 miRNA
    11	    334 pseudogene
    12	    312 tRNA
    13	    299 snoRNA
    14	    262 pre_miRNA
    15	    115 rRNA
    16	     32 snRNA
</code></pre> 
> 2. Total number of genes per chromosome arm (X, Y, 2L, 2R, 3L, 3R, 4)
<pre><code>
$ zcat *.gtf* | grep -P '^\S+\s\S+\s+gene\b'| cut -f1 | grep -Eo 'X|Y|2L|2R|3L|3R|4' | sort |uniq -c | sort -nr | nl
</code></pre> 
<pre><code>
     1	   4202 3R
     2	   3628 2R
     3	   3501 2L
     4	   3464 3L
     5	   2676 X
     6	    122 4
     7	    113 Y
</code></pre> 
