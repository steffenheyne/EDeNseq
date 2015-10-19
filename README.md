# EDeNseq

EDeNseq is a tool for fast clustering and classification of (metagenomics) reads.

# 1. Installation

EDeNseq requires C++11 support! GCC 4.8 should be sufficient. 

1. clone or download EDeNseq

	`git clone https://github.com/steffenheyne/EDeNseq.git`
	
2. Then change into source folder: `cd EDeNseq/src` 

3. Build it with: `make`  - that's it!

# 2. Example: classify metagenomic reads against bacterial genomes

The folder [`test_data/`](https://github.com/steffenheyne/EDeNseq/blob/master/test_data/) 
contains the test files for the following test run (please adapt path names):

```EDeNseq -a CLASSIFY -i test_data/test.reads.fna.gz 
--index_seqs test_data/test.genomes.fa.gz  --index_bed test_data/test.small.bed  
--numThreads 4  -b 30 -F 5 --num_hash_shingles 3 --num_repeat_hash_functions 2
-r 4 -d 7 --min_radius 4 --min_distance 7 --seq_window 70 --index_seq_shift 0.15 
--seq_shift 0.13 --pure_approximate_sim 0```

- [**`test_data/test.reads.fna.gz`**](https://github.com/steffenheyne/EDeNseq/blob/master/test_data/test.reads.fna.gz) contains 231.578 sequences (100nt) sampled from 29 bacterial chromosomes
- [**`test_data/test.genomes.fa.gz`**](https://github.com/steffenheyne/EDeNseq/blob/master/test_data/test.genomes.fa.gz) contains 62 bacterial chromosomes
- [**`test_data/test.small.bed`**](https://github.com/steffenheyne/EDeNseq/blob/master/test_data/test.small.bed) 
indicates regions from the genomes that will be indexed and against that reads are classified. I
n the given example BED file, column 4 contains the taxonomic id (NCBI) at the organism level. 
Although there are 62 different entries/sequences given, they map to only 29 organisms. 
Hence, reads are classified between these 29 organisms. The output file [**`test.reads.fna.gz.classified.tab.gz`**](https://github.com/steffenheyne/EDeNseq/blob/master/test_data/test.reads.fna.gz.classified.tab.txt) 
gives a mapping table in the header. 

- The command above produces the index file **`test.small.bed.bhi`** and a results file **`test.reads.fna.gz.classified.tab.gz`**. 
An example output is shown in 
[**`test_data/test.reads.fna.gz.classified.tab.txt`**](https://github.com/steffenheyne/EDeNseq/blob/master/test_data/test.reads.fna.gz.classified.tab.txt)
When re-running the same command the existing index file (*.bhi) is used and not created again. 
All lines in the result file starting with "#" are header lines.
The tags "#HIST_IDX" give the mapping between the internal ID and the 
feature ID provided in the BED file. 

- [**`test_data/test.log`**](https://github.com/steffenheyne/EDeNseq/blob/master/test_data/test.log) shows the terminal output of the above example

# 3. Sequence Classification

EDeNseq classifies sequences (e.g. reads) against an index of genomes or 
any desired set of target sequences.

Sequence classifciation mode is started with 

```EdENseq -a CLASSIFY```

Please see all available paramters with   

```EdENseq -a CLASSIFY -h``` or ```EdENseq -a CLASSIFY -H```


## 3.1 Build the Index

### 3.1.1 Index BED File

In order to build the index, EDeNseq needs a BED file and fasta file. 
The BED file is used to define regions on sequences that are
indexed under a certain label. 

The labeling makes the whole classification procedure very powerful and 
general, as this can be the sequence name (BED column 1) or a certain feature 
identifier (BED file column 4). This allows to index different 
sequences under the same label. For example, genomes of different 
strains or related organisms, homologue transcripts or sequences of a certain 
RNA family (Rfam).

The BED file has the following format:

`<SEQNAME>	<START>	<END>	<LABEL>	.	0	<SEQ_DESC>	<LABEL_DESC>`

The first 4 cols are required. Columns 7 and 8 are optional, but then also 
columns 5 and 6 have to be provided.

`<LABEL>` can be any string or number. Internally each unique label is mapped to
a number. 

`<LABEL_DESC>` is used to annotate the results. 

During indexing, EDeNseq scans through the fasta file and for each found seq, 
it selects all provided regions in BED file for that sequence and index them. 
This means also that nothing is indexed for a seq without BED entry as well as 
if there is a BED entry without sequence in the FASTA file. 

### 3.1.2 Create BED file for large fasta file

One can use samtools for this purpose: `samtools faidx my_genomes.fa.gz` 

Then best use `awk` or a simlar tool to get the columns in the right order! 


# 4. Sequence Clustering

EDeNseq can be used to cluster large-scale sequence dataset.

In clustering mode, no BED file is necessary and all 
provided sequences are used!

**NOTE:** Clustering is under development!

# Notes

## Download bacterial reference genomes from NCBI

Here is an example about creating the index file for bacterial reference genomes. 

1. Get table with reference genomes from NCBI ftp server (prok_reference_genomes.txt)
2. Download latest version for each refseq entry via eutils in fasta format, zip it and save it in $outdir 

`wget ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prok_reference_genomes.txt`

`outdir="refseq_bact_genomes"; mkdir -p $outdir;` 

`for i in $(cat prok_reference_genomes.txt | awk '{FS="\t";OFS="\t";if (NR==1) next; split($4,CHR,","); for (i in CHR){gsub(" ","_",$3);print CHR[i]}}'); do
	get "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="$i"&rettype=fasta&retmode=text" -q -O - | gzip > $outdir/$i.fa.gz; 
done`


