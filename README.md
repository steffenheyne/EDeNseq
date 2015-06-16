# EDeNseq

EDeNseq is a tool for fast clustering and classification of (metagenomics) reads.

# Read Classification

EDeNseq classifies reads against an index of e.g. genomes. In order to build the index,
EDeNseq needs a tab-delimited file with three columns:

`<idx>	<fa-file-path>	<description>`

The first column, `<idx>`, is a running number (>=1) that denotes the bin/slot 
in the index. 
Second column, `<fa-file-path>`, is the path to a fasta file used for that idx. 
The last column is some description of the file/slot. 
This is used only for output purposes. 

Note: there can multiple lines with the same `<idx>`! This allows to put several 
fasta files (genomes) in the same index slot! This is useful e.g. for bacteria 
with serveral chromosomes or in order to put very similar genomes together. 

## Download bacterial reference genomes from NCBI

Here is an example about creating the index file for bacterial reference genomes. 

1. Get table with reference genomes from NCBI ftp server (prok_reference_genomes.txt)
2. Download latest version for each refseq entry via eutils in fasta format, zip it and save it in $outdir 
3. create EDeNseq index table (refseq_bact_genomes.index_list)

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prok_reference_genomes.txt

outdir="refseq_bact_genomes"; mkdir -p $outdir; 

for i in $(cat prok_reference_genomes.txt | awk '{FS="\t";OFS="\t";if (NR==1) next; split($4,CHR,","); for (i in CHR){gsub(" ","_",$3);print CHR[i]}}'); do
	get "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id="$i"&rettype=fasta&retmode=text" -q -O - | gzip > $outdir/$i.fa.gz; 
done

cat prok_reference_genomes.txt | awk '{FS="\t";OFS="\t";if (NR==1) next;idx++;split($4,CHR,","); for (i in CHR){gsub(" ","_",$3);print CHR[i],idx,$3}}' |  awk -v PATH=$(pwd)/$outdir '{OFS="\t";fasta=PATH"/"$1".fa.gz"; if ((getline _ < fasta)>=0){ print $2,fasta,$1"|"$3 }}' > refseq_bact_genomes.index_list

```