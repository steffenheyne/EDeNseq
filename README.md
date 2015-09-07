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

## Download bacterial reference genomes from NCBI

Here is an example about creating the index file for bacterial reference genomes. 

1. Get table with reference genomes from NCBI ftp server (prok_reference_genomes.txt)
2. Download latest version for each refseq entry via eutils in fasta format, zip it and save it in $outdir 
3. create EDeNseq index table (refseq_bact_genomes.index_list)

`wget ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prok_reference_genomes.txt`

`outdir="refseq_bact_genomes"; mkdir -p $outdir;` 

`for i in $(cat prok_reference_genomes.txt | awk '{FS="\t";OFS="\t";if (NR==1) next; split($4,CHR,","); for (i in CHR){gsub(" ","_",$3);print CHR[i]}}'); do
	get "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="$i"&rettype=fasta&retmode=text" -q -O - | gzip > $outdir/$i.fa.gz; 
done`

`cat prok_reference_genomes.txt | awk '{FS="\t";OFS="\t";if (NR==1) next;idx++;split($4,CHR,","); for (i in CHR){gsub(" ","_",$3);print CHR[i],idx,$3}}' |  awk -v PATH=$(pwd)/$outdir '{OFS="\t";fasta=PATH"/"$1".fa.gz"; if ((getline _ < fasta)>=0){ print $2,fasta,$1"|"$3 }}' > refseq_bact_genomes.index_list`


## Example 
`EDeNseq -a CLASSIFY -f FASTA -i test_10k_redsea.eden -b 30 --num_hash_shingles 2 --numThreads 4 --index_data_file refseq_bact_genomes.bed --seq_window 80 --seq_shift 0.3 -r 6 -d 14 --min_radius 6 --min_distance 14  --num_repeat_hash_functions 0  -F 10`