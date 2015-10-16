# EDeNseq

EDeNseq is a tool for fast clustering and classification of (metagenomics) reads.

# Read Classification

EDeNseq classifies sequences (e.g. reads) against an index of e.g. genomes or 
any desired set of target sequences.

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

Note: Currently only in clustering mode, no BED file is necessary and all 
provided sequences are used!

## Example 
`EDeNseq -a CLASSIFY -f FASTA -i test_data/test.reads.fna.gz --index_seqs test_data/test.genomes.fa.gz  --num_repeat_hash_functions 2 --num_hash_shingles 3 --numThreads 4 --index_bed test_data/test.small.bed --seq_window 70 --index_seq_shift 0.15 --seq_shift 0.13 -b 30 -F 5 -r 4 -d 7 --min_radius 4 --min_distance 7 --pure_approximate_sim 0`


## Download bacterial reference genomes from NCBI

Here is an example about creating the index file for bacterial reference genomes. 

1. Get table with reference genomes from NCBI ftp server (prok_reference_genomes.txt)
2. Download latest version for each refseq entry via eutils in fasta format, zip it and save it in $outdir 

`wget ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prok_reference_genomes.txt`

`outdir="refseq_bact_genomes"; mkdir -p $outdir;` 

`for i in $(cat prok_reference_genomes.txt | awk '{FS="\t";OFS="\t";if (NR==1) next; split($4,CHR,","); for (i in CHR){gsub(" ","_",$3);print CHR[i]}}'); do
	get "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="$i"&rettype=fasta&retmode=text" -q -O - | gzip > $outdir/$i.fa.gz; 
done`


