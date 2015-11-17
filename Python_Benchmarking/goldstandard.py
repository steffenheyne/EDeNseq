""" This module is designed for generating goldstandards for (meta)genomic data analysis tool comparison.

@note: The functionality of the module strongly depends on the file structure of the NCBI ftp server. Changes within the
file structure can lead to errors while iterating through the file system of the ftp server or failure of fetching data.
Therefore, correctness highly depends on the version of the module.

@author: Katharina Moos
@version: 0.1.1 19.10.2015
"""


# IMPORTS ##############################################################################################################


import sys
import argparse
import subprocess as sub

import gzip
import os
import shutil
from ftplib import FTP

import random
import time

import numpy


# METHODS ##############################################################################################################


# Getting Reference Sequences ------------------------------------------------------------------------------------------


def ncbi_ftp_tree(ofile=""):
    """ This function accesses the NCBI ftp server for genomes and either build a hash or writes a file containing the
    server's directory tree.

    The hash assigns an organism group a nested list with the subdirectories of all organisms belonging to that group.
    Similarly, the file is segmented into listings of directories belonging to one organism group.
    Organism groups are: eucaryotes (eu) (not including fungi), fungi (fu), bacteria (ba) and viruses (vi).

    hash structure: {'eu': [[ directory of organism 1, [directory of chromosome 1, files of chromosome 1],
                                                       [directory of chromosome 2, ...], ...],
                            [directory of organism 2, ...],
                            [directory of organism 3, ...],
                            ...],
                    fu, bac or vi: [[subdirectory of organism 1, files of organism 1],
                                    [subdirectory of organism 2, files of organism 2],
                                    ...]
                    }

    file structure: #eu\n
                    org:dir_org1\n
                    chr:dir_chr1\ttab-delimited list of files
                    #fu, bac or vi\n
                    org:dir_org1\ttab-delimited list of files

    @note: Only nucleotide data (.fna files) are fetched.

    @param ofile:   in case the hash should be saved as file, this is the output file name; else, give an empty string
                    DEFAULT: empty string

    @return: a hash representing the NCBI ftp server directory tree structure, or a corresponding file
    """

    # initializing hash or outputfile

    if ofile:
        ftp_tree = None
        ofile = open(ofile, 'w')
        ofile.write('#eu\n')
    else:
        ftp_tree = {'eu': [], 'fu': [], 'ba': [], 'vi': []}

    # sets storing ftp directories and files containing unrelevant data

    eu_unallowed = frozenset(['ASSEMBLY_BACTERIA', 'ASSEMBLY_REPORTS', 'Bacteria', 'Bacteria_DRAFT', 'CLUSTERS',
                              'Chloroplasts', 'Fungi', 'GENOME_REPORTS', 'HUMAN_MICROBIOM', 'IDS', 'INFLUENZA',
                              'MITOCHONDRIA', 'MapView', 'PLANTS', 'Plasmids', 'Protozoa', 'README-09.14.2011.old',
                              'README.txt', 'README_AGP_FORMAT_CHANGE', 'TARGET', 'TOOLS', 'Viruses', 'all', 'genbank',
                              'old_genomeID2nucGI', 'refseq'])

    fung_unallowed = frozenset(['README.txt'])

    bac_unallowed = frozenset(
        ['AFLN00000000.1.11117.log', 'AMHJ00000000.1.11222.log', 'CLUSTERS', 'ERR', 'NZ_AAEK01000176', 'ReadMe.txt',
         'SameSpecies.gi', 'all.GeneMark.tar.gz', 'all.Glimmer3.tar.gz', 'all.Prodigal.tar.gz', 'all.asn.tar.gz',
         'all.faa.tar.gz', 'all.ffn.tar.gz', 'all.fna.tar.gz', 'all.frn.tar.gz', 'all.gbk.tar.gz', 'all.gff.tar.gz',
         'all.ptt.tar.gz', 'all.rnt.tar.gz', 'all.rpt.tar.gz', 'all.val.tar.gz', 'summary.txt'])

    vir_unallowed = frozenset(
        ['Viroids_RefSeq_and_neighbors_genome_data.tab', 'Viruses_RefSeq_and_neighbors_genome_data.tab',
         'all.asn.tar.gz', 'all.faa.tar.gz', 'all.ffn.tar.gz', 'all.fna.tar.gz', 'all.frn.tar.gz', 'all.gbk.tar.gz',
         'all.gff.tar.gz', 'all.ptt.tar.gz', 'all.rnt.tar.gz', 'all.rpt.tar.gz', 'all.val.tar.gz', 'SameSpecies.gi'])

    # initializing ftp

    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    ftp.cwd('genomes')

    # browsing through directories (need to go one directory deeper for eucaryotic data)

    # EUCARYOTIC DATA

    print('storing data tree: Eucaryotes ... ')

    for element in ftp.nlst():

        if element not in eu_unallowed:  # iterating through valid subdirectories (containing animal genomic data)

            org_entry = []  # entry storing an organism's subdirectory tree

            print(element + '... ', end='')
            ftp.cwd(element)            # CWD (down): organism's directory
            org_entry.append(element)   # storing organism's directory

            for subdir in ftp.nlst():   # iterating through organism's subdirectories (containing chromosome subdirs)

                if subdir.startswith('CHR'):
                    chrom_entry = []    # entry storing a chromosome's files

                    ftp.cwd(subdir)                         # CWD (down): chromosome's subdirectory
                    chrom_entry.append(subdir)              # storing chromosome's directory

                    for file in ftp.nlst():                 # iterating through chromosome's files

                        if file.endswith('.fna') or file.endswith('fa.gz'):  # checking for nucleotide sequence data
                            chrom_entry.append(file)                         # appending to chromosome's files

                    if len(chrom_entry) > 1:                # checking if chromosomal nucleotide sequence data exists,
                        org_entry.append(chrom_entry)       # appending to organism's list

                    ftp.cwd('..')       # CWD (up): organism's directory

            if len(org_entry) > 1:      # checking if chromosomal nucleotide data exists,
                                        # appending organism tree to hash
                if ofile:
                    ofile.write('org:' + org_entry.pop(0))
                    for chrom in org_entry:
                        ofile.write('\nchr:' + chrom.pop(0))
                        for chrom_file in chrom:
                            ofile.write('\t' + chrom_file)
                else:
                    ftp_tree['eu'].append(org_entry)
                print('stored')
            else:
                print('not stored')

            ftp.cwd('..')       # CWD (up): superdirectory (containing organisms)

    # FUNGI, BACTERIAL AND VIRAL DATA

    for organism in [['Fungi', 'fu', fung_unallowed], ['Bacteria', 'ba', bac_unallowed],
                     ['Viruses', 'vi', vir_unallowed]]:

        print('storing data tree: ' + organism[0] + ' ... ')

        key = organism[1]
        if ofile:
            ofile.write('\n#' + key)

        ftp.cwd(organism[0])                # CWD (down): organism group directory
        for element in ftp.nlst():          # iterating through valid subdirectories (containing genomic data)
            if element not in organism[2]:

                print(element)
                entry = [element]           # entry storing an organism's files
                ftp.cwd(element)                    # CWD (down): organism's directory

                for file in ftp.nlst():             # iterating through organism's files

                    if file.endswith('.fna')\
                            or file.endswith('fa.gz'):  # checking for nucleotide sequence data,
                        entry.append(file)              # appending to chromosome's files

                if len(entry) > 1:                  # checking if chromosomal nucleotide data exists,
                                                    # appending organism tree to hash
                    if ofile:
                        ofile.write('\norg:' + entry.pop(0))
                        for chrom_file in entry:
                            ofile.write('\t' + chrom_file)
                    else:
                        ftp_tree[key].append(entry)

                ftp.cwd('..')               # CWD (up): superdirectory (containing organisms)
        ftp.cwd('..')                   # CWD (up): suerdirectory (containing all organism groups)

    # finally, closing ftp and returning

    ftp.quit()

    if ofile:
        ofile.close()
        return
    else:
        return ftp_tree


def wget_fetch_file(ftp_file, outfile):
    """ This function fetches a file from a given ftp server.

    @param ftp_file:    ftp file path of file to be fetched
    @param outfile:     output file path to store ftp file to

    @return: tuple containing wget standard output and standard error
    """

    # printing call to stdout
    print("wget --output-document " + outfile + " " + ftp_file)

    # calling subprocess
    wget_process = sub.Popen("wget --output-document " + outfile + " " + ftp_file, shell=True, stdout=sub.PIPE,
                             stderr=sub.PIPE)

    # waiting for process to terminate and catching stdout and stderr streams
    output, errors = (s.decode() for s in wget_process.communicate())

    return output, errors


def rand_sequence(outfile, out_format, description='goldstandard random sequence', length=1000, rand_seed=0):
    """ This function generates a random nucleotide sequence.

    Output file format: fa or fa.gz

    An even distribution of nucleotides is assumed.

    @param outfile:     output file path to store sequence (format suffix .fa or .fa.gz is appended automatically)
    @param out_format:  output file format (valid formats: fa, fa.gz)
    @param description: description line to use for the random sequence within the fasta file ('>' is appended
                        automatically)
                        DEFAULT: goldstandard random sequence
    @param length:      length of the sequence
                        DEFAULT: 1000
    @param rand_seed:   random seed to be used
                        DEFAULT: current time

    @return: void
    """

    # initializing variables

    alphabet = (b'A', b'T', b'G', b'C')                 # nucleotides to use for the random sequence
    if rand_seed:                                       # initializing random number generator
        random.seed(rand_seed)
    else:
        random.seed(time.clock())

    # generating sequence and writing to file

    if out_format == 'fa':                   # opening uncompressed output file
        o_file = open(outfile + '.fa', 'w')
    elif out_format == 'fa.gz':              # opening compressed output file
        o_file = gzip.open(outfile + '.fa.gz', 'wb')
    else:                                   # raising error if unknown file format is used
        raise ValueError("Output file format '" + out_format + "' as specified by 'outformat' is invalid; valid formats"
                         " are: 'fa', 'fa.gz'.\n")

    o_file.write(b'>' + bytes(description, 'UTF-8'))    # writing sequence description line
    for nt_pos in range(length):
        if nt_pos % 100 == 0:                           # every 100 (99: 0-based) characters, a line-break is inserted
            o_file.write(b'\n')
        o_file.write(alphabet[random.randint(0, 3)])    # randomly choosing a nucleotide, writing to file directly)

    # closing output file and returning

    o_file.close()
    return


# Getting Reads --------------------------------------------------------------------------------------------------------


def get_reads(infile, in_format, outfile, out_format='fa.gz', read_length=100, step_min=20, step_max=120, append=False,
              rand_seed=0):
    """ This function generates error-free reads from a reference sequence.

    Input file format: fa or fa.gz
    Output file format: fa.gz or fq.gz

    Reads are generated sliding over the reference sequence, using a randomly drawn step size from a given step size
    interval for each new read. Step size is defined as the distance between the starting positions of adjacent reads
    and must be larger than. Coverage can be influenced by the interval's mean step size but cannot be given
    explicitly.

    @note: As error-free reads are generated, the fastq quality string is set to a sequence of '~'.

    @param infile:      reference sequence file to generate reads from (fa or fa.gz)
    @param in_format:   reference sequence input file format (valid formats: fa or fa.gz)
    @param outfile:     output file path to store reads generated (format suffix .fa.gz or .fq.gz is appended
                        automatically)
    @param out_format:  output file format (valid formats: fa.gz or fq.gz)
                        DEFAULT: fa.gz
    @param read_length: length of the reads
                        DEFAULT: 100
    @param step_min:    minimum step size to be taken for the next read
                        DEFAULT: 10
    @param step_max:    maximum step size to be taken for the next read; this parameter is set to step_min if its value
                        is less than the value specified by step_min
                        DEFAULT: 120
    @param append:      if output file 'outfile' already exists, append reads to this file
                        DEFAULT: False
    @param rand_seed:   random seed to be taken
                        DEFAULT: current time

    @return: void
    """

    # checking step_min and step_max values validity, initializing random number generator, opening input file,
    # opening output file, initializing additional read information for fastq format

    if step_min < 1:            # checking step_min
        raise ValueError("Minimum step size '" + str(step_min) + "' as defined by 'step_min' is invalid; minimum step"
                         " size must be larger than 0.\n")
    if step_max < step_min:     # checking step_max
        step_max = step_min

    if rand_seed:               # initializing random seed
        random.seed(rand_seed)
    else:
        random.seed(time.clock())

    if in_format == 'fa':               # opening uncompressed input file
        ref_seq = open(infile, 'rb')
    elif in_format == 'fa.gz':          # opening compressed input file
        ref_seq = gzip.open(infile, 'rb')
    else:                               # raising error message if file format is invalid
        raise ValueError("Input file format '" + in_format + "' as specified by 'informat' is invalid; valid formats "
                         "are: 'fa', 'fa.gz'.\n")

    if os.path.exists(outfile + '.' + out_format):
        if append:
            o_file = gzip.open(outfile + '.' + out_format, 'ab')    # opening output file (for appending)
        else:
            raise IOError("Output file '" + outfile + "' as given by 'outfile' already exists.")
    else:
        o_file = gzip.open(outfile + '.' + out_format, 'wb')    # opening output file

    if out_format == 'fq.gz':           # initializing read information for fastq format
        descriptor = b'@'                   # leading character of fastq description line
        fastq_appending = b'\n+\n'
        for i in range(read_length):
            fastq_appending += b'~'
    elif out_format == 'fa.gz':         # initializing read information for fasta format
        descriptor = b'>'                   # leading character of fasta description line
        fastq_appending = b''               # leaving additional read information empty in case of fasta format
    else:                               # raising error message if file format is invalid
        raise ValueError("Output file format '" + out_format + "' as specified by 'informat' is invalid; valid formats "
                         "are: 'fa.gz', 'fq.gz'.\n")

    # iterating through input file, writing to output file

    ref_seq_description = next(ref_seq).strip()[1:]     # storing and skipping description line of reference sequence
    # (needed for read description)
    read_counter = 0                                    # counter for number of reads (needed for read description)
    current_step = random.randint(step_min, step_max)   # initializing current step size

    current_seq = next(ref_seq).strip()[current_step:]  # initializing current subsequence of reference sequence to
    # get reads from, shifting its beginning to read's starting position as defined by the current step size

    for line in ref_seq:        # iterating through input file
        current_seq += line.strip()     # appending next input file line to current reference sequence subsequence

        # as soon as current subsequence is long enough, a read is generated
        if len(current_seq) >= read_length:
            o_file.write(descriptor + ref_seq_description + b'_' + b'read' + bytes(str(read_counter), 'UTF-8') + b'\n' +
                         current_seq[0:read_length] + fastq_appending + b'\n')    # writing read to outfile
            read_counter += 1       # incrementing read counter

        # updating current step size and current subsequence
        current_step = random.randint(step_min, step_max)
        current_seq = current_seq[current_step:]

    # closing file streams and returning

    ref_seq.close()
    o_file.close()
    return


def art_illumina(infile, outfile, mode, read_length, coverage, frag_mean, frag_std,
                 ins1=0, ins2=0, del1=0, del2=0, id_tag=None, rand_seed=0):
    """ This function executes the ART simulation tool art_illumina to generate ngs read data.

    Output file format: fq.gz

    @param infile:      reference sequence input file to be used (fa)
    @param outfile:     output file directory
    @param mode:        type of read simulation; allowed values are 'single' for single-end, 'paired' for paired-end,
                        and 'mate' for mate-pair reads
    @param read_length: length of reads
    @param coverage:    coverage of reads to be simulated
    @param frag_mean:   mean fragment size for paired-end or mate-pair simulation (in case of mode 'single' this
                        parameter is ignored, in case of mode 'paired' this parameter cannot exceed a value of 2000 and
                        will automatically be scaled down if being above 2000)
    @param frag_std:    standard deviation of fragment size (in case of mode 'single' this parameter is ignored)

    @param ins1:        insertion rate for first read to be used
                        DEFAULT: 0
    @param ins2:        insertion rate for first read to be used (will be ignored in case of mode 'single')
                        DEFAULT: 0
    @param del1:        deletion rate for first read to be used
                        DEFAULT: 0
    @param del2:        deletion rate for first read to be used (will be ignored in case of mode 'single')
                        DEFAULT: 0
    @param id_tag:      a prefix tag for the read ID
                        DEFAULT: None
    @param rand_seed:   seed to be used for random number generators; using the same seed will yield the same
                        simulations
                        DEFAULT: 0 (random seed is chosen by art_illumina)

    @return: error messages raised during illumina call
    """

    # GENERATING ART COMMAND

    # command to use for subprocess calling
    call = "art_illumina -i " + infile + " -l " + str(read_length) + " -f " + str(coverage)

    # checking mode and mode-dependent parameters
    if mode == 'paired':
        if frag_mean > 2000:
            frag_mean = 2000
        call += " -p -m " + str(frag_mean) + "-s " + str(frag_std) + " -ir2 " + str(ins2) + " -dr2 " + str(del2)
    elif mode == 'mate':
        call += " -mp -m " + str(frag_mean) + "-s " + str(frag_std) + " -ir2 " + str(ins2) + " -dr2 " + str(del2)

    # else (mode == single) frag_mean and frag_std are left out

    call += " -ir " + str(ins1) + " -dr " + str(del1)
    if id_tag:
        call += " -d " + str(id_tag)
    if rand_seed:
        call += " -rs " + str(rand_seed)

    # finally, adding outfile directory
    call += " -o " + outfile

    # CALLING SUBPROCESS

    sim_process = sub.Popen(call, shell=True, stdout=sub.PIPE, stderr=sub.PIPE)  # calling subprocess
    output, errors = (s.decode() for s in sim_process.communicate())        # waiting for process to terminate and
    # catching stdout and stderr streams

    # RETURNING ERRORS

    return errors


# Manipulating Reads ---------------------------------------------------------------------------------------------------


def mutate_reads(infile, in_format, outfile, out_format, snp_min=0.1, snp_max=0, snp_const=True,
                 in_min=0.1, in_max=0, in_const=True, del_min=0.1, del_max=0, del_const=True,
                 append=False, rand_seed=0):
    """ This function mutates nucleotide reads from a given input file and stores them in a new file.

    Input file format: fa, fa.gz
    Output file format: fa.gz, fq.gz

    It is possible to define to use either mutation rates (per read) or mutation probabilities (per nucleotide).
    Furthermore, mutation rate/probability can be constant for all reads or randomly drawn from an probability interval
    for each single read. Both using a rate or probability and the corresponding values can be specified distinctly
    for SNPs (short nucleotide polymorphisms), insertions and deletions.

    @note: Insertions, SNPs and deletions are introduced in this order.
    @note: Below a maximum SNP rate of 0.5, single nucleotide positions are not mutated multiple times. If the maximum
           SNP rate defined exceeds 0.5, multiple mutations per nucleotide position are enabled.
    @note: Silent SNP mutations are not introduced. Though, it is possible that silent mutations occur due to multiple
           mutations or insertions and deletions canceling out, and that SNPs are inserted or removed by insertions and
           deletions.

    @note: In case of fasta input and fastq output, the fastq quality string is set to a sequence of '~'.
    @note: This method can also be used to convert fasta to fastq and vice verca, adjusting all relevant parameters.

    @param infile:      input file containing reads (fa, fa.gz, fq, fq.gz)
    @param in_format:   input file format (valid formats: fa, fa.gz, fq, fq.gz)
    @param outfile:     output file path to store mutated reads (format suffix fa.gz or fq.gz is appended automatically)
    @param out_format:  output file format (valid formats are: fa.gz, fq.gz)

    @param snp_min:     minimum snp rate/probability (per nucleotide) to be used
                        DEFAULT: 0.1
    @param snp_max:     maximum snp probability (per nucleotide) to be used; set to 0 to use a constant probability as
                        specified by snp_min; if snp_const is set to True, this parameter will be ignored
                        DEFAULT: 0
    @param snp_const:   set to True to use a mutation rate, set to False to use a mutation probability
                        DEFAULT: False
    @param in_min:      minimum insertion rate/probability (per nucleotide) to be used
                        DEFAULT: 0
    @param in_max:      maximum insertion probability (per nucleotide) to be used; set to 0 to use a constant
                        probability as specified by in_min
                        DEFAULT: 0
    @param in_const:    set to True to use a insertion rate, set to False to use a insertion probability
                        DEFAULT: True
    @param del_min:     minimum deletion rate/probability (per nucleotide) to be used
                        DEFAULT: 0
    @param del_max:     maximum deletion probability (per nucleotide) to be used; set to 0 to use a constant probability
                        as specified by del_min
                        DEFAULT: 0
    @param del_const:   set to True to use a deletion rate, set to False to use a deletion probability
                        DEFAULT: True

    @param append:      If output file 'outfile' already exists, append read to this file
                        DEFAULT: False
    @param rand_seed:   random seed to be taken
                        DEFAULT: current time

    @return: void
    """

    # checking parameter values ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if snp_min < 0:
        raise ValueError("Minimum SNP rate '" + str(snp_min) + "' as defined by 'snp_min' is invalid; minimum snp rate"
                         " must at least be 0.\n")
    elif snp_min > 1:
        raise ValueError("Minimum SNP rate '" + str(snp_min) + "' as defined by 'snp_min' is invalid; snp rate cannot"
                         " exceed a value of 1.\n")
    if snp_max < snp_min:
        snp_max = snp_min
    if snp_max > 0.5:
        multiple = True
        print("NOTE: As maximum SNP rate '" + str(snp_max) + "' as defined by 'snp_max' exceeds 50%, multiple mutations"
              " per nucleotide position are enabled.\n")
    else:
        multiple = False

    if in_min < 0:
        raise ValueError("Minimum insertion rate '" + str(in_min) + "' as defined by 'in_min' is invalid; minimum "
                         " insertion rate must at least be 0.\n")
    if in_max < in_min:
        in_max = in_min

    if del_min < 0:
        raise ValueError("Minimum deletion rate '" + str(del_min) + "' as defined by 'del_min' is invalid; minimum "
                         " deletion rate must at least be 0.\n")
    if del_max < del_min:
        del_max = del_min

    if out_format not in ('fa.gz', 'fq.gz'):
        raise ValueError("Output file format '" + out_format + "' as specified by 'outformat' is invalid; valid formats"
                         " are: 'fa.gz', 'fq.gz'.\n")

    # initializations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # initializing random number generator

    if rand_seed:
        random.seed(rand_seed)
    else:
        random.seed(time.clock())

    # opening input and output file

    if in_format in ('fa', 'fq'):            # opening uncompressed input file
        i_file = open(infile, 'rb')
    elif in_format in ('fa.gz', 'fq.gz'):    # opening compressed input file
        i_file = gzip.open(infile, 'rb')
    else:                                   # raising error message in case of invalid input file format
        raise ValueError("Input file format '" + in_format + "' as specified by 'informat' is invalid; valid formats "
                         "are: 'fa', 'fa.gz', 'fq', 'fq.gz'.\n")

    if os.path.exists(outfile):
        if append:
            o_file = gzip.open(outfile + '.' + out_format, 'ab')    # opening output file (for appending)
        else:
            raise IOError("Output file '" + outfile + "' as given by 'outfile' already exists.")
    else:
        o_file = gzip.open(outfile + '.' + out_format, 'wb')    # opening output file

    # implementing mutation approaches ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # implementing snps: method 'snp'

    if snp_min == 0 and snp_max == 0:   # no SNPs are to be introduced, return original read
        def snp(read, snp_rate, multiple, alphabet=b'ATGC'):
            return read

    elif snp_const:                     # SNP rate

        def snp(read, snp_rate, multiple, alphabet=b'ATGC'):

            read_len = len(read)                        # read length
            read_nts = list(read)                       # converting read string to list of single nucleotides
            snp_num = int(round(read_len * snp_rate))   # total number of SNPs to be introduced
            mutated = {}                                # hash storing positions already mutated

            for i in range(snp_num):                    # introducing as many mutations as defined by snp_num

                pos = random.randint(0, read_len-1)     # drafting position to be mutated
                if pos in mutated:                          # checking if position is already mutated
                    if not multiple:                        # if multiple mutations are disabled ...
                        while pos in mutated:                   # ... drafting until a not mutated position is selected
                            pos = random.randint(0, read_len-1)
                        mutated[pos] = b''                      # ... and storing this position within hash 'mutated'

                old_nt = read_nts[pos]                      # storing read nucleotide to be mutated

                new_nt = random.choice(alphabet)            # randomly choosing a new nucleotide ...
                while new_nt == old_nt:                     # ... until this new nucleotide differs from the old one
                    new_nt = random.choice(alphabet)

                read_nts[pos] = new_nt                      # introducing mutation

            return bytes(read_nts)   # returning mutated read

    else:                               # SNP probability

        def snp(read, snp_rate, multiple, alphabet=b'ATGC'):
            read_nts = list(read)   # converting read string to list of single nucleotides
            for idx, nts in enumerate(read_nts):        # iterating through read nucleotides
                if random.random() <= snp_rate:             # with probability 'snp_rate', current nucleotide is mutated
                    new_nts = random.choice(alphabet)       # randomly drafting a new nucleotide ...
                    while new_nts == nts:                   # ... until it differs from the old one
                        new_nts = random.choice(alphabet)
                    read_nts[idx] = new_nts                 # introducing mutation
            return bytes(read_nts)  # returning mutated read

    # implementing insertions: method 'insertion'

    if in_min == 0 and in_max == 0:     # no insertions are to be introduced, return original read
        def insertion(read, in_rate, alphabet=b'ATGC'):
            return read

    elif in_const:                      # insertion rate

        def insertion(read, in_rate, alphabet=b'ATGC'):
            read_len = len(read)    # storing read length
            read_nts = {}           # converting read string to hash of single nucleotides
            for i in range(read_len):
                read_nts[i] = list(read[i:i+1])
            ins_num = int(round(in_rate * read_len, 0))     # total number of insertions to be introduced

            for i in range(ins_num):            # introducing as many insertions as defined by ins_num
                read_nts[random.randint(0, read_len-1)].append(random.choice(alphabet))     # inserting random
                #  nucleotide at random position

            mutated_read = b''                      # assembling mutated read from hash elements
            for i in range(read_len):
                mutated_read += bytes(read_nts[i])

            return mutated_read     # returning mutated read

    else:                               # insertion probability

        def insertion(read, in_rate, alphabet=b'ATGC'):
            read_len = len(read)    # storing read length
            read_nts = {}           # converting read string to hash of single nucleotides
            for i in range(read_len):
                read_nts[i] = list(read[i:i+1])

            for pos in range(read_len):             # iterating through nucleotide positions to introduce insertions
                if random.random() <= in_rate:      # introduce insertions with probability 'in_rate'
                    read_nts[random.randint(0, read_len-1)].append(random.choice(alphabet))     # inserting random
                    # nucleotide at random position

            mutated_read = b''                      # assembling mutated read from hash elements
            for i in range(read_len):
                mutated_read += bytes(read_nts[i])

            return mutated_read     # returning mutated read

    # implementing deletions: method 'deletion'

    if del_min == 0 and del_max == 0:   # no deletions are to be introduced, return original read
        def deletion(read, del_rate):
            return read

    elif del_const:                     # deletion rate

        def deletion(read, del_rate):

            read_nts = list(read)   # converting read string to list of single nucleotides
            del_num = int(round(del_rate * len(read), 0))   # total number of deletions to be introduced

            for i in range(del_num):    # introducing as many deletions as defined by ins_num
                read_nts.pop(random.randint(0, len(read_nts)-1))    # removing nucleotide at random position

            return bytes(read_nts)  # returning mutated read

    else:                               # deletion probability
        def deletion(read, del_rate):
            read_nts = list(read)   # converting read string to list of single nucleotides
            del_counter = 0         # counting number of deletions for correct list index via indices
            for idx in range(len(read_nts)):    # iterating through nucleotide positions to introduce deletions
                if random.random() <= del_rate:     # introduce deletions with probability 'del_rate'
                    read_nts.pop(idx - del_counter)     # deleting nucleotide at current position
                    del_counter += 1                    # updating deletion counter
            return bytes(read_nts)  # returning mutated read

    # implementing wrapper

    def snp_indels(read, snp_rate, in_rate, del_rate, multiple, alphabet=b'ATGC'):
        read = insertion(read, in_rate, alphabet)       # introducing insertions
        read = snp(read, snp_rate, multiple, alphabet)  # introducing SNPs
        read = deletion(read, del_rate)                 # introducing deletions
        return read

    # implementing file format conversion, including mutation of read sequence,
    # defining a function 'process_read' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if out_format == 'fa.gz':

        # both file formats are fasta / fasta.gz

        def process_read(filestream):
            # first line (read description)
            first_line = next(filestream)                   # getting read description (no need to modify anything)
            # second line (read sequence)
            read = next(filestream).strip()                 # getting read sequence removing line break
            snp_rate = random.uniform(snp_min, snp_max)     # randomly choosing snp_rate from given interval
            in_rate = random.uniform(in_min, in_max)        # randomly choosing in_rate from given interval
            del_rate = random.uniform(del_min, del_max)     # randomly choosing del_rate from given interval
            read = snp_indels(read, snp_rate, in_rate, del_rate, multiple)  # mutating read
            # returning
            return first_line + read + b'\n'

    else:

        # input is fasta / fasta.gz and output is fastq.gz

        def process_read(filestream):
            # first line (read description)
            first_line = next(filestream)                   # getting read description
            first_line = b'@' + first_line.split(b'>')[1]   # adapting read description to fastq format
            # second line (read sequence)
            read = next(filestream).strip()                 # getting read sequence removing line break
            snp_rate = random.uniform(snp_min, snp_max)     # randomly choosing snp_rate from given interval
            in_rate = random.uniform(in_min, in_max)        # randomly choosing in_rate from given interval
            del_rate = random.uniform(del_min, del_max)     # randomly choosing del_rate from given interval
            read = snp_indels(read, snp_rate, in_rate, del_rate, multiple)  # mutating read
            # third line ('+')
            third_fourth_line = b'\n+\n'                    # fastq third read entry line: '+'
            # fourth line (quality string)
            for i in range(len(read)):
                third_fourth_line += b'~'                   # fastq fourth read entry line: quality score
            # returning
            return first_line + read + b'\n' + third_fourth_line + b'\n'

    # iterating through input file,
    # writing to output file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    try:                        # generating read entries from input file read and writing to output file,
    # until whole input file is processed (exception StopIteration is raised)
        read_entry = process_read(i_file)
        while read_entry:
            o_file.write(read_entry)
            read_entry = process_read(i_file)
    except StopIteration:       # as soon as entire input file is processed
        pass
    finally:                    # in case any exception occurred, file streams are closed
        i_file.close()
        o_file.close()

    return


def shuffle_reads(infile, in_format, outfile, append=False, rand_seed=0):
    """ This function shuffles reads from a given file and stores them in a new file.

    Input file format: fa, fa.gz, fq, fq.gz
    Output file format: gzipped fasta format (as quality strings don't make any sense for shuffled reads)

    @param infile:      input file with reads to be shuffled (fa, fa.gz, fq, fq.gz)
    @param in_format:   input file format (valid formats:fa, fa.gz, fq, fq.gz)
    @param outfile:     output file path to store shuffled reads (format suffix .fa.gz or .fq.gz is appended
                        automatically)
    @param append:      If output file 'outfile' already exists, append read to this file
                        DEFAULT: False
    @param rand_seed:   random seed to be taken
                        DEFAULT: current time

    @return: void
    """

    # initializing random seed

    if rand_seed:
        random.seed(rand_seed)
    else:
        random.seed(time.clock())

    # opening input and output files

    if in_format in ('fa.gz', 'fq.gz'):
        i_file = gzip.open(infile, 'rb')
    elif in_format in ('fa', 'fq'):
        i_file = open(infile, 'rb')
    else:
        raise ValueError("Input file format '" + in_format + "' as specified by 'in_format' is invalid; valid formats"
                         " are: 'fa', 'fa.gz', 'fq', 'fq.gz'")

    if os.path.exists(outfile):
        if append:
            o_file = gzip.open(outfile + '.fa.gz', 'ab')    # opening output file (for appending)
        else:
            raise IOError("Output file '" + outfile + "' as given by 'outfile' already exists.")
    else:
        o_file = gzip.open(outfile + '.fa.gz', 'wb')    # opening output file

    # iterating through input file reads and writing to output file

    for line in i_file:
        o_file.write(line)                      # writing read description to output file
        read = (list(next(i_file).strip()))     # getting read sequence as list
        random.shuffle(read)                    # shuffling read
        o_file.write(bytes(read) + b'\n')       # writing shuffled read to output file
        if in_format in ('fq', 'fq.gz'):
            next(i_file)
            next(i_file)

    # closing file streams and returning

    i_file.close()
    o_file.close()
    return


def read_filter(infile, in_format, prop, rand_seed=0):
    """ This function randomly filters a given input read data file.

    Input file format: fq or fq.gz
    Output file format: fq.gz

    Filtered read file is stored in the same directory as the original read file with the appending '_filtered.fq.gz'.

    @param infile:      input file to be filtered (fq or fq.gz)
    @param in_format:   input file format (valid formats: fq, fq.gz)
    @param prop:        proportion of reads TO KEEP (number between 0 and 1, excluding 0 and 1)
    @param rand_seed:   seed to use for initialization of random number generator
                        DEFAULT: current time

    @return: the number of reads written
    """

    if rand_seed:   # initializing random number generator
        random.seed(rand_seed)
    else: random.seed(time.clock())
    o_file = gzip.open(os.path.join(infile + '_filtered.fq.gz'), 'wb')  # opening output file
    reads_written = 0   # initializing number of reads written

    # iterating through lines of input fastq file

    if in_format == 'fq':            # HANDLING OF UNCOMPRESSED DATA
        i_file = open(infile, 'rb')
    elif in_format == 'fq.gz':       # HANDLING OF COMPRESSED DATA
        i_file = gzip.open(infile, 'rb')
    else:
        raise ValueError("Input file format '" + in_format + "' as specified by 'informat' is invalid; valid formats"
                         " are: 'fq', 'fq.gz'.\n")

    for line in i_file:
        if random.random() < prop:  # if random number is below proportion value, read is taken
            o_file.write(line)
            reads_written += 1
            for i in range(3):
                o_file.write(next(i_file))
        else:  # else, read is skipped
            for i in range(3):
                next(i_file)

    # finally, closing files and returning

    i_file.close()
    o_file.close()
    return reads_written


# Wrappers to Get Goldstandards ----------------------------------------------------------------------------------------


def gold_random(out_dir, refseq_num=10, refseq_len=1000, read_num=1, read_len=100, read_merge=True, step_min=[20],
                step_max=[120], snp_min=[0.1], snp_max=[0], snp_const=[1], in_min=[0], in_max=[0], in_const=[1],
                del_min=[0], del_max=[0], del_const=[1], shuffle=[0], rand_seed=0):
    """ This function generates a set of random reference sequences and one or multiple sets of corresponding reads.

    Reference sequences output file format: fa.gz
    Reads output file format: fa.gz

    Reference sequences to be generated are either stored in separate or joined into one fasta or fasta.gz file.
    Afterwards, reads are derived. One read set is defined as one read derivation of all reference sequences and stored
    within either one file ore separate files (one read file per reference sequence). As many such sets are derived as
    defined by 'read_num', differing in read composition, coverage and mutation of reads by giving the parameters
    multiple values (one parameter value for each set of reads). Read files are named including the set number and, in
    case of multiple read files per set, the underlying reference sequence.
    For more information on how mutated reads are generated, please see the function 'goldstandard.mutate_reads'.

    @note:  Parameters for read generation: Give one value to use for all read set generations, or as many values as
            read sets to be created ('read_num').
    @note:  For each set to be generated, an additional set of corresponding non-mutated raw reads is produced. If both
            parameters 'step_min' and 'step_max' only have one parameter value assigned, only one set of non-mutated
            raw reads is generated and serves as template for all read sets specified. If all parameter values defining
            read mutations are 0, only the non-mutated raw read set will be created.
    @note:  Read lengths may differ from the read length as defined by 'read_len', if insertions and / or deletions are
            to be introduced (internally, first non-mutated reads of length 'read_len' are generated and then further
            processed to insert mutations).
    @note:  As the random seed is given to subprocesses producing the reads, same reads will be generated if a
            subprocess is called with the same parameter values and same random seed.

    @param out_dir:         directory to store gold_random results in (must be existent)
    @param refseq_num:      number of reference sequences to be generated
                            DEFAULT: 10
    @param refseq_len:      length of reference sequences to be generated
                            DEFAULT: 1000
    @param read_num:        number of read sets to be generated
                            DEFAULT: 1
    @param read_len:        length of reads to be generated; if insertions and / or deletions are to be inserted into
                            reads, the final read length may differ from the read length defined here
                            DEFAULT: 100
    @param read_merge:      set to 'True' to create one read file per read set, set to 'False' to create one read file
                            per reference sequence per read set
                            DEFAULT: True
    @param rand_seed:       random seed to be used
                            DEFAULT: current time

    @param step_min:    parameter for read generation: minimum step size to be taken for the next read
                        DEFAULT: 10
    @param step_max:    parameter for read generation: maximum step size to be taken for the next read; this parameter
                        is set to step_min if its value is less than the value specified by step_min
                        DEFAULT: 120
    @param snp_min:     minimum snp rate/probability (per nucleotide) to be used
                        DEFAULT: 0.1
    @param snp_max:     maximum snp probability (per nucleotide) to be used; set to 0 to use a constant probability as
                        specified by snp_min; if snp_const is set to True, this parameter will be ignored
                        DEFAULT: 0
    @param snp_const:   set to any non-zero integer to use a mutation rate, set to 0 to use a mutation probability
                        DEFAULT: 1
    @param in_min:      minimum insertion rate/probability (per nucleotide) to be used
                        DEFAULT: 0
    @param in_max:      maximum insertion probability (per nucleotide) to be used; set to 0 to use a constant
                        probability as specified by in_min
                        DEFAULT: 0
    @param in_const:    set to any non-zero integer to use a insertion rate, set to 0 to use a insertion probability
                        DEFAULT: 1
    @param del_min:     minimum deletion rate/probability (per nucleotide) to be used
                        DEFAULT: 0
    @param del_max:     maximum deletion probability (per nucleotide) to be used; set to 0 to use a constant probability
                        as specified by del_min
                        DEFAULT: 0
    @param del_const:   set to any non-zero integer to use a deletion rate, set to 0 to use a deletion probability
                        DEFAULT: 1
    @param shuffle:     set to any non-zero integer to generate a read file with all read sequences shuffled (shuffling
                        is only done for 'mutated' reads, not for non-mutated template reads).
                        DEFAULT: 0

    @return: void
    """

    # checking and adapting parameter values ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    rgp =   [step_min, step_max, snp_min, snp_max, list(bool(value) for value in snp_const), in_min, in_max,
             list(bool(value) for value in in_const), del_min, del_max, list(bool(value) for value in del_const),
             list(bool(value) for value in shuffle)]  # 'read_gen_params',
             # storing all parameters needed for read generation within one list, converting integer values for
             # parameters which are to be interpreted as booleans to such
    raw_reads_unique = False    # storing whether or not only one raw reads file is to be generated
    raw_reads_only = False      # storing whether or not only raw read (no additional read sets) are to be generated

    # initializing random seed

    if rand_seed:
        random.seed(rand_seed)
    else:
        random.seed(time.clock())

    # checking output directories

    if not os.path.exists(out_dir):
        raise ValueError("Directory to store results as defined by 'out_dir' must be existent (directory given: '"
                         + out_dir + "').")     # checking for existent directory

    refseq_dir = os.path.join(out_dir, 'gold_random_refseqs')   # defining subdirectory to store reference sequences
    read_dir = os.path.join(out_dir, 'gold_random_reads')       # defining subdirectory to store read sets
    os.mkdir(refseq_dir)        # create refseq directory
    os.mkdir(read_dir)          # create read sets directory

    # checking and adapting number of mutation relevant parameter values

    for idx, param in enumerate(rgp):   # iterating through parameters and checking for their number of values
        value_num = len(param)
        if value_num == 1:                  # if one value is given, this value is multiplied to the number of read sets
            # to be generated ('read_num')
            rgp[idx] = numpy.full(read_num, param[0])
        elif value_num != read_num:         # if multiple values are given, but not as many as read sets to be
            # generated, an exception is raised
            raise ValueError("The number of parameter values given for read generation parameters is not consistently"
                             "1 or 'read_num' = " + str(read_num) + " (stopped after first error occurred for parameter"
                             " no.: " + str(idx+1) + ").\n")

    # checking: one or multiple raw read files (raw_reads_unique)

    if len(rgp[0]) == 1 and len(rgp[1]) == 1:   # if only one parameter value is specified for both 'step_min' and
        # 'step_max', it is stored that only one raw read set is to be produced.
        raw_reads_unique = True
    raw_read_unique_exists = False  # in case a unique raw reads set is to be produced, this flag is set as soon as the
    # raw reads file was generated

    # checking mutation parameters, checking: read sets (raw_reads_only)

    mutation_values = 0                         # summing up all values given for parameters defining read mutations,
    # if they sum up to 0 this means that no mutations are to be inserted and only the raw reads file is produced

    for mutation_parameter in rgp[2:]:          # iterating to all relevant parameters
        for value in mutation_parameter:            # getting all values of a parameter
            if not type(value) == bool:                 # if it is not a boolean parameter, ...
                if value >= 0:                              # ... check if value is not below zero and add it
                    mutation_values += value
                else:                                       # ... if value is below zero, raise error
                    raise ValueError("At least one invalid parameter value was found: mutation rates / probabilities "
                                     "must at least be 0.\n")
            else:                                       # if parameter is boolean, skip and continue with next parameter
                break

    if mutation_values == 0:                    # if all mutation parameter values are 0, it is stored that only the
        # raw read set is to be produced
        raw_reads_only = True

    # generating reference sequences ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    for seq in range(refseq_num):   # generating as many sequences as defined by refseq_num

        rand_sequence(os.path.join(refseq_dir, 'refseq' + str(seq)), 'fa.gz',
                      description='reference sequence ' + str(seq), length=refseq_len, rand_seed=rand_seed)

    # generating reads ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    reads_raw = ''          # raw reads file name
    reads_set = ''          # read set file name
    reads_shuffle = ''      # shuffled reads file name

    for set_number in range(read_num):  # generating as many read sets as defined

        for seq in range(refseq_num):       # for each reference sequence, generating the corresponding read sets

            # declaring file names (distinguishing unique/non-unique raw read files and merged/not-merged read files)

            if read_merge:
                if raw_reads_unique:
                    reads_raw = os.path.join(read_dir, 'RAW')
                else:
                    reads_raw = os.path.join(read_dir, 'RAW_' + str(set_number))
                reads_set = os.path.join(read_dir, 'SET_' + str(set_number))
                reads_shuffle = os.path.join(read_dir, 'SHUFFLE_' + str(set_number))
            else:
                if raw_reads_unique:
                    reads_raw = os.path.join(read_dir, 'refseq' + str(seq) + '_RAW')
                else:
                    reads_raw = os.path.join(read_dir, 'refseq' + str(seq) + '_RAW_' + str(set_number))

                reads_set = os.path.join(read_dir, 'refseq' + str(seq) + '_SET_' + str(set_number))
                reads_shuffle = os.path.join(read_dir, 'refseq' + str(seq) + '_SHUFFLE_' + str(set_number))

            # generating raw reads (only if there is not a unique raw reads file already existing)

            if not raw_read_unique_exists:
                get_reads(os.path.join(refseq_dir , 'refseq' + str(seq)) + '.fa.gz', 'fa.gz', reads_raw, 'fa.gz',
                          read_length=read_len, step_min=rgp[0][set_number], step_max=rgp[1][set_number],
                          append=read_merge, rand_seed=rand_seed)  # generating raw reads

            # generating mutated and shuffled reads (only if not only raw reads are to be produced) ...

                # ... in case 'not merged' (within reference sequence iteration, because one raw read file is produced
                # for each reference sequence)
            if (not read_merge) and (not raw_reads_only):

                mutate_reads(reads_raw + '.fa.gz', 'fa.gz', reads_set, 'fa.gz', snp_min=rgp[2][set_number],
                             snp_max=rgp[3][set_number], snp_const=rgp[4][set_number], in_min=rgp[5][set_number],
                             in_max=rgp[6][set_number], in_const=rgp[7][set_number], del_min=rgp[8][set_number],
                             del_max=rgp[9][set_number], del_const=rgp[10][set_number], rand_seed=False)

                if rgp[11][set_number]:
                    shuffle_reads(reads_set + '.fa.gz', 'fa.gz', reads_shuffle)

                # ... in case 'merged' (outside of reference sequence iteration, because one raw read file is
                # produced for each set)

        if (read_merge) and (not raw_reads_only):
            mutate_reads(reads_raw + '.fa.gz', 'fa.gz', reads_set, 'fa.gz', snp_min=rgp[2][set_number],
                         snp_max=rgp[3][set_number], snp_const=rgp[4][set_number], in_min=rgp[5][set_number],
                         in_max=rgp[6][set_number], in_const=rgp[7][set_number], del_min=rgp[8][set_number],
                         del_max=rgp[9][set_number], del_const=rgp[10][set_number], rand_seed=False)

            if rgp[11][set_number]:
                shuffle_reads(reads_set + '.fa.gz', 'fa.gz', reads_shuffle)

        # finally, as soon as one reads set was generated, and only one unique raw reads file is to be generated,
        # 'raw_read_unique_exists' is set to true

        if raw_reads_unique:
            raw_read_unique_exists = True

    # returning

    return


def gold_mgen(outfile, tmpdir, num, props,
              ftp_tree_file='', mode='single', read_length=100, coverage=30, frag_mean=None, frag_std=None,
              ins1=0, ins2=0, del1=0, del2=0, id_tag=None, rand_seed=0):
    """ This function generates a metagenomic data goldstandard.

    Output file format: fq.gz

    It randomly fetches genomic data (nucleotide sequences only) from the NCBI ftp server and uses the ART simulation
    tool to generate reads from these reference sequences.
    First, genomes to build the data on are fetched from the NCBI ftp server (randomly, but w.r.t. the organism
    proportions specified).
    Second, for each genome ngs read data is simulated using the ART simulation tools (choosing the coverage randomly).
    Third, those reads are filtered randomly to create a non-uniform read distribution as usually seen for sequencing
    data.
    Finally, ngs read data of single genomes is written into a common file (compressed fastq).

    @note: It is also possible to generate single-organism genomic data by giving an artificial ftp_tree_file containing
    the directory of only one organism and adjusting all relevant parameter values.

    @param outfile:     path to store goldstandard output file (format suffix fq.gz is appended automatically)
    @param tmpdir:      directory to temporarily store fetched genome data
    @param num:         total number of genomes to be fetched; as long as this number is not too high, multiple fetching
                        of a genome is disabled
    @param props:       list containing the proportions of eucaryotes (not including fungi), fungi, bacteria and viruses
                        (in this order) to be contained within the data as decimal numbers; according to fetched
                        genomes, not to number of reads within metagenome data; due to rounding, the total number of
                        genomes fetched may differ up to a value of 2 from the value of parameter num specified

    @param ftp_tree_file:   in case the ftp tree was already saved as file, this is the file name; else, set to False
                            (don't set this flag); file format must be equivalent to format specified within method
                            'ncbi_ftp_tree'
                            DEFAULT: False (flag not set)
    @param mode:            ART simulation parameter; type of read simulation; allowed values are 'single' for
                            single-end, 'paired' for paired-end, and 'mate' for mate-pair reads
                            DEFAULT: single
    @param read_length:     ART simulation parameter; length of reads
                            DEFAULT: 100
    @param coverage:        ART simulation parameter; coverage of reads to be simulated
                            DEFAULT: 30
    @param frag_mean:       ART simulation parameter; mean fragment size for paired-end or mate-pair simulation; in case
                            of mode 'single' this parameter is ignored, in case of mode 'paired' this parameter cannot
                            exceed a value of 2000 and will automatically be scaled down if above 2000
                            DEFAULT: None
    @param frag_std:        ART simulation parameter; standard deviation of fragment size; in case of mode 'single' this
                            parameter is ignored
                            DEFAULT: None
    @param ins1:            ART simulation parameter; insertion rate for first read to be used
                            DEFAULT: 0
    @param ins2:            ART simulation parameter; insertion rate for first read to be used (will be ignored in case
                            of mode 'single')
                            DEFAULT: 0
    @param del1:            ART simulation parameter; deletion rate for first read to be used
                            DEFAULT: 0
    @param del2:            ART simulation parameter; deletion rate for first read to be used (will be ignored in case
                            of mode 'single')
                            DEFAULT: 0
    @param id_tag:          ART simulation parameter; a prefix tag for the read ID
                            DEFAULT: None
    @param rand_seed:       random number generator seed to be used for executions using random numbers; default:
                            current time; using the same seed will yield the same goldstandard (this enables
                            reproducible goldstandards)
                            DEFAULT: 0 (random seed is chosen by each process itself)

    @return: a list containing the number of genomes fetched for each organism group
    """

    # generating temporary directory strings

    gold_tmpdir = os.path.join(tmpdir, 'gold_mgen')
    gold_tmpdir_refseq = os.path.join(gold_tmpdir, 'gold_mgen_refseq')
    gold_tmpdir_art = os.path.join(gold_tmpdir, 'gold_mgen_art')

    try:    # try ... finally statement guarantees to delete all temporary files in case of unexpected execution abort

        # INITIALIZING TEMPORARY FILE DIRECTORIES

        os.mkdir(gold_tmpdir)
        os.mkdir(gold_tmpdir_refseq)
        os.mkdir(gold_tmpdir_art)

        # INITIALIZING RANDOM NUMBER GENERATOR -------------------------------------------------------------------------

        if rand_seed:
            random.seed(rand_seed)
        else:
            random.seed(time.clock())

        # FETCHING GENOMIC DATA ----------------------------------------------------------------------------------------

        # preparing fetching

        # computing number or genomes to be fetched for each organism group

        eu_num = round(num * props[0])
        fung_num = round(num * props[1])
        bac_num = round(num * props[2])
        vir_num = round(num * props[3])

        # setting list of organism numbers
        org_nums = [eu_num, fung_num, bac_num, vir_num]

        # getting ncbi ftp tree (as hash)

        # hash structure: {'eu': [[ directory of organism 1, [directory of chromosome 1, files of chromosome 1],
        #                                                    [directory of chromosome 2, ...], ...],
        #                         [directory of organism 2, ...],
        #                         [directory of organism 3, ...],
        #                        ...],
        #                 fu, bac or vi: [[subdirectory of organism 1, files of organism 1],
        #                                 [subdirectory of organism 2, files of organism 2],
        #                                 ...]
        #                 }

        # file structure: #eu\n
        #                 org:dir_org1\n
        #                 chr:dir_chr1\ttab-delimited list of files
        #                 #fu, bac or vi\n
        #                 org:dir_org1\ttab-delimited list of files

        print('\n--- STORING NCBI FTP TREE ---\n')

        if ftp_tree_file:  # CASE: building hash from file

            ftp_tree = {'eu': [], 'fu': [], 'ba': [], 'vi': []}     # initializing hash
            hash_file = open(ftp_tree_file, 'r')                    # opening ftp tree file

            # iteration for eucaryotes

            if not next(hash_file).startswith('#eu'):  # skipping first line, checking for correct file format
                raise ValueError(
                    "File does not fulfill the format specifications. For more information, see the function "
                    "'goldstandard.ncbi_ftp_tree'.\n")

            line = next(hash_file)  # initializing first line

            while not line.startswith('#'):  # iterating through tree for eucaryotes
                org_dir = [line.strip().split(':')[1]]  # getting organism directory, saving it within a list
                line = next(hash_file)

                while line.startswith('chr'):  # iterating through lines holding information on current organism
                    print(line)
                    line_parts = line.strip().split('\t')  # splitting tab-delimited line into single parts
                    line_parts[0] = line_parts[0].split(':')[1]  # getting chromosome directory
                    org_dir.append(line_parts)
                    line = next(hash_file)

                ftp_tree['eu'].append(org_dir)

            # iteration for non-eucaryotes

            org = line.strip().split('#')[1]  # initializing first key

            for line in hash_file:  # iterating through remaining file
                if not line.startswith('#'):  # building hash for current organism
                    line_parts = line.strip().split('\t')
                    line_parts[0] = line_parts[0].split(':')[1]
                    ftp_tree[org].append(line_parts)
                else:  # initializing key for next organism
                    org = line.strip().split('#')[1]

            hash_file.close()  # finally closing ftp tree file

        else:
            ftp_tree = ncbi_ftp_tree()  # CASE: buidling hash from ftp server using method 'ncbi_ftp_tree'

            # fetching genomic data for single organism groups
            #   choosing organisms randomly
            #   fetching a genomic data file via wget

        print('\n--- FETCHING FILES ---\n\n')
        for org in ['eu', 'fung', 'bac', 'vir']:  # ITERATING THROUGH ORGANISM GROUPS

            counter = 0  # counter for progess report to stdout

            # depending on organism group,
            # getting number of genomes to be fetched, getting organism list from hash, checking for multiple fetching
            # checking if number of genomes to be fetched doesn't exceed size of subdirectory list (genomes available),
            # initializing ftp path (for fetching with wget), setting filename prefix

            print('fetching ' + org + ' data')

            if org == 'eu':
                org_num = org_nums[0]
                org_list = ftp_tree['eu']
                if len(org_list) < 0.95 * org_num:
                    orgs_fetched = False
                else:
                    orgs_fetched = set('filling_value')
                ftp_path = 'ftp://ftp.ncbi.nlm.nih.gov/genomes'
                filename_pref = 'eu'

            elif org == 'fung':
                org_num = org_nums[1]
                org_list = ftp_tree['fu']
                if len(org_list) < 0.95 * org_num:
                    orgs_fetched = False
                else:
                    orgs_fetched = set('filling_value')
                ftp_path = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/Fungi'
                filename_pref = 'fung'

            elif org == 'bac':
                org_num = org_nums[2]
                org_list = ftp_tree['ba']
                if len(org_list) < 0.95 * org_num:
                    orgs_fetched = False
                else:
                    orgs_fetched = set('filling_value')
                ftp_path = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria'
                filename_pref = 'bac'

            else:
                org_num = org_nums[3]
                org_list = ftp_tree['vi']
                if len(org_list) < 0.95 * org_num:
                    orgs_fetched = False
                else:
                    orgs_fetched = set('filling_value')
                ftp_path = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses'
                filename_pref = 'vir'

            for i in range(org_num):  # ITERATION PER ORGANISM GROUP, selecting as many organisms as defined by org_num

                # CHOOSING ORGANISM

                last_idx = len(org_list) - 1  # getting last (highest) index of organism list
                rand_idx = random.randint(0,
                                          last_idx)  # choosing random index from possible indices of subdirectory list

                if orgs_fetched:  # checking if multiple fetching is disabled,
                    while rand_idx in orgs_fetched:  # choose new index as long as index chosen has already been used
                        rand_idx = random.randint(0, last_idx)

                organism = org_list[rand_idx]

                ftp_path = os.path.join(ftp_path, organism.pop(0))  # ftp_path (down): organism directory

                # FETCHING ORGANISM'S FILES

                if org == 'eu':     # EUCARYOTIC DATA (need to go one more subdirectory deeper to get files)

                    for chromosome in organism:  # iterating through chromosomes

                        ftp_path = os.path.join(ftp_path, chromosome.pop(0))  # ftp_path (down): chromosome directory
                        for file in chromosome:  # iterating through files within chromosome subdirectory

                            stdout, stderr = wget_fetch_file(os.path.join(ftp_path, file),
                                                             os.path.join(gold_tmpdir_refseq,
                                                                          filename_pref + str(i) + '_' + file))
                            print(stdout)
                            print(stderr)
                            # add chosen index to set orgs_fetched
                            orgs_fetched.add(rand_idx)

                        ftp_path = os.path.split(ftp_path)[0]  # ftp_path (up): organism directory

                else:               # FUNGI, BACTERIAL and VIRAL DATA (already at file level)

                    for file in organism:  # iterating through files within chromosome subdirectory

                        stdout, stderr = wget_fetch_file(os.path.join(ftp_path, file),
                                                         os.path.join(gold_tmpdir_refseq,
                                                                      filename_pref + str(i) + '_' + file))
                        print(stdout)
                        print(stderr)
                        # add chosen index to set orgs_fetched
                        orgs_fetched.add(rand_idx)

                counter += 1
                if counter % 10 == 0:
                    print(str(counter) + ' ... ')

                ftp_path = os.path.split(ftp_path)[
                    0]  # ftp_path (up): superdirectory (containing all organisms)

    # EXECUTING ART (READ SIMULATION) ----------------------------------------------------------------------------------

        print('\n--- EXECUTING ART SIMULATION ---\n')
        print(os.listdir(gold_tmpdir_refseq))       # printing reference sequence files to stdout

        # iterating through fetched files, executing art (illumina, single-end)

        for idx, file in enumerate(os.listdir(gold_tmpdir_refseq)):
            art_illumina(os.path.join(gold_tmpdir_refseq, file), os.path.join(gold_tmpdir_art, str(idx)),
                         mode, read_length, coverage, frag_mean, frag_std, ins1, ins2, del1, del2, id_tag, rand_seed)

    # FILTERING READS AND MERGING TO OUTPUT FILE -----------------------------------------------------------------------

        print('\n--- FILTERING READS ---\n')

        # iterating through genomic data files, filtering reads

        print(os.listdir(gold_tmpdir_art))          # printing files generated by art to stdout
        for file in os.listdir(gold_tmpdir_art):
            if file.endswith('fq'):

                # doing filtering
                prop = numpy.random.gumbel(0.02, 0.02)  # randomly draw proportion of reads to be fetched (according to
                # a Gumbel distribution)
                while prop > 1 or prop < 0:
                    prop = numpy.random.gumbel(0.02, 0.02)

                read_filter(os.path.join(gold_tmpdir_art, file), 'fq', prop, rand_seed)  # filtering

        # finally, merging single filtered files

        o_file = gzip.open(outfile + '.fq.gz', 'wb')  # opening output file

        for file in os.listdir(gold_tmpdir_art):    # iterating through read files,
            if file.endswith('filtered'):           # only considering files with filtered reads
                reads_filtered = gzip.open(os.path.join(gold_tmpdir_art, file), 'rb')
                for line in reads_filtered:
                    o_file.write(line)
                reads_filtered.close()

        o_file.close()  # closing output file again

    # DELETING TEMPORARY FILES -----------------------------------------------------------------------------------------

    # finally, temporary files are deleted:

    finally:

        print('\n--- REMOVING TEMPORARY FILES; RETURNING ---')
        if os.path.exists(os.path.join(gold_tmpdir)):
            shutil.rmtree(gold_tmpdir)

    return org_nums


# Main Routine #########################################################################################################


if __name__ == '__main__':

    # Main routine handling goldstandard calls

    # defining function to generate parameter help messages from docstrings

    def help_messages(docstring):

        help_strings = docstring.split('@param')[1:]  # getting descriptions as defined by '@param' tags
        help_strings[-1] = help_strings[-1].split('@return')[0]     # splitting the end of the docstring away, too

        for idx, message in enumerate(help_strings):    # iterating through '@param' descriptions
            help_mss = ''
            message = message.split(':', maxsplit=1)[1]         # getting rid of parameter name
            for line in message.split('\n'):                    # iterating through single lines
                help_mss += line.strip() + '\n'                     # removing leading/trailing whitespaces
            help_mss.rstrip()                                   # finally removing trailing newline
            help_strings[idx] = help_mss                        # storing processed help message
        return help_strings

    # initializing argument parser

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('gold_method', type=str,
                        help="goldstandard method to be executed: 'gold_mgen' or 'gold_random'")

# gold_mgen ------------------------------------------------------------------------------------------------------------

    if sys.argv[1] == 'gold_mgen':

        # printing user input

        print("\n------ GOLDSTANDARD: GOLD_MGEN ------\n")

        # getting help messages from gold_mgen docstring

        help_strings = help_messages(gold_mgen.__doc__)

        # adding arguments

        parser.add_argument('outfile', type=str,
                            help=help_strings[0])
        parser.add_argument('tmpdir', type=str, help=help_strings[1])
        parser.add_argument('genome_number', type=int, help=help_strings[2])
        parser.add_argument('org_ratios', type=float, nargs=4,
                            help=help_strings[3])

        parser.add_argument('--ftp_tree', type=str, default='', help=help_strings[4])
        parser.add_argument('--mode', type=str, default='single', help=help_strings[5])
        parser.add_argument('--read_length', type=int, default=100, help=help_strings[6])
        parser.add_argument('--coverage', type=int, default=30, help=help_strings[7])
        parser.add_argument('--frag_mean', type=float, default=None, help=help_strings[8])
        parser.add_argument('--frag_std', type=float, default=None, help=help_strings[9])
        parser.add_argument('--ins1', type=float, default=0, help=help_strings[10])
        parser.add_argument('--ins2', type=float, default=0, help=help_strings[11])
        parser.add_argument('--del1', type=float, default=0, help=help_strings[12])
        parser.add_argument('--del2', type=float, default=0, help=help_strings[13])
        parser.add_argument('--id_tag', type=str, default=None, help=help_strings[14])
        parser.add_argument('--rand_seed', type=float, default=0, help=help_strings[15])

        # parsing arguments and calling gold_mgen

        goldstandard_call = sys.stdin.realine().strip().split(' ', maxsplit=1)[1]   # from the goldstandard call,
        # removing trailing linebreak and leading method specification (not needed for ArgumentParser)
        args = parser.parse_args(goldstandard_call.split())
        gold_mgen(args.outfile, args.tmpdir, args.rand_seed, args.genome_number, args.org_ratios,
                  ftp_tree_file=args.ftp_tree, mode=args.mode, read_length=args.read_length, coverage=args.coverage,
                  frag_mean=args.frag_mean, frag_std=args.frag_std, ins1=args.ins1, ins2=args.ins2, del1=args.del1,
                  del2=args.del2, id_tag=args.id_tag, rand_seed=args.rand_seed)

# gold_random ----------------------------------------------------------------------------------------------------------

    elif sys.argv[1] == 'gold_random':

        # printing user input

        print("\n------ GOLDSTANDARD: GOLD_RANDOM ------\n")

        # getting help messages from gold_random docstring

        help_strings = help_messages(gold_random.__doc__)

        # adding arguments

        parser.add_argument('out_dir', type=str, help=help_strings[0])
        parser.add_argument('--refseq_num', type=int, default=10, help=help_strings[1])
        parser.add_argument('--refseq_len', type=int, default=1000, help=help_strings[2])
        parser.add_argument('--read_num', type=int, default=3, help=help_strings[3])
        parser.add_argument('--read_len', type=int, default=100, help=help_strings[4])
        parser.add_argument('--read_merge', action='store_true', help=help_strings[5])

        parser.add_argument('--step_min', type=int, nargs='+', default=[20], help=help_strings[6])
        parser.add_argument('--step_max', type=int, nargs='+', default=[120], help=help_strings[7])
        parser.add_argument('--snp_min', type=float, nargs='+', default=[0.1], help=help_strings[8])
        parser.add_argument('--snp_max', type=float, nargs='+', default=[0], help=help_strings[9])
        parser.add_argument('--snp_const', type=int, nargs='+', default=[1], help=help_strings[10])
        parser.add_argument('--in_min', type=float, nargs='+', default=[0], help=help_strings[11])
        parser.add_argument('--in_max', type=float, nargs='+', default=[0], help=help_strings[12])
        parser.add_argument('--in_const', type=int, nargs='+', default=[1], help=help_strings[13])
        parser.add_argument('--del_min', type=float, nargs='+', default=[0], help=help_strings[14])
        parser.add_argument('--del_max', type=float, nargs='+', default=[0], help=help_strings[15])
        parser.add_argument('--del_const', type=int, nargs='+', default=[1], help=help_strings[16])
        parser.add_argument('--shuffle', type=int, nargs='+', default=[0], help=help_strings[17])
        parser.add_argument('--rand_seed', type=int, default=0, help=help_strings[18])

        # parsing arguments, writing a log file into 'read_dir' to report parameter values used, and calling gold_random

        args = parser.parse_args()  # parsing arguments
        logfile = open(os.path.join(args.out_dir, "logfile.txt"), 'w')  # opening log file
        logfile.write(str(args))    # writing to log file
        logfile.close()     # closing log file

        gold_random(args.out_dir, refseq_num=args.refseq_num, refseq_len=args.refseq_len,
                    read_num=args.read_num, read_len=args.read_len, read_merge=args.read_merge,
                    step_min=args.step_min, step_max=args.step_max,
                    snp_min=args.snp_min, snp_max=args.snp_max, snp_const=args.snp_const,
                    in_min=args.in_min, in_max=args.in_min, in_const=args.in_const,
                    del_min=args.del_min, del_max=args.del_max, del_const=args.del_const,
                    shuffle=args.shuffle, rand_seed=args.rand_seed)

# invalid user input ---------------------------------------------------------------------------------------------------

    else:
        sys.stdout.write("\nInvalid Goldstandard method parsed. Returning.\n")
