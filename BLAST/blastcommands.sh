#!/bin/bash

# how to blast multiple immune genes, identifying the immune genes across multiple species

# load blast
module load blast

# setwd for output
cd /project/berglandlab/madi/blast_outputs

# blast list of proteins against our genome, query is genes, subject is the genome you are blasting to:
blastn |
-query /project/berglandlab/madi/blast_files/genes/gnbps.fa \
-subject /project/berglandlab/madi/blast_files/genomes/nadpulex.fa \
-out nadpulex.txt

blastn |
-query /project/berglandlab/madi/blast_files/genes/gnbps.fa \
-subject /project/berglandlab/madi/blast_files/genomes/dma.fa \
-out gnbpgff.txt 

# protein - protein blast
blastp |
-query /project/berglandlab/madi/blast_files/proteins/gnbp.fa \
-subject /project/berglandlab/madi/blast_files/genomes/protein_seq/magna.fa \
-outfmt 5 \
-out magnaprotein.xml

# add before "-out" if wanting xml
-outfmt 5
