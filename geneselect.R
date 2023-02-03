# Get candidate gene list and convert to 
# 12.5.2022
# ijob -c 1 --mem=50G -p largemem -A berglandlab
# module load gcc/9.2.0 openmpi/3.1.6 R/4.2.1; R

# log in to terminal
# ssh med7xdv@rivanna.hpc.virginia.edu

# load R
# module load gcc/7.1.0
# module load openmpi/3.1.4
# module load intel/18.0 
# module load intelmpi/18.0
# module load R/4.0.3

# Libraries
library(data.table)
library(tidyverse)
library(readxl)

# Working directory
setwd("/project/berglandlab/madi/")

# Read in gene annotations
panth <- data.table(read_excel("/project/berglandlab/daphnia_ref/Daphnia_annotation_PANTHER.xls"))
colnames(panth)[36:38] <- c("bio_func", "mol_func", "cell_func")

# GTF file
gene.gtf <- data.table(fread("/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gtf"))
colnames(gene.gtf)[1:5] <- c("chrom", "file", "sec", "start", "stop")

# Get gene info
extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- strsplit(gtf_attributes, "; ")
  att <- gsub("\"","",unlist(att))
  if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
    return( unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
  }else{
    return(NA)}
}
gene.gtf$gene <- unlist(str_remove_all(lapply(gene.gtf$V9, extract_attributes, "transcript_id"), pattern = ";"))

# Adding gene splice
gene.gtf <- data.table(gene.gtf %>% 
                         mutate(gene.id=tstrsplit(gene, "-")[[1]],
                                splice=tstrsplit(gene, "-")[[2]]))

# Gene list
out.ld <- unique(panth[bio_func %like% "mmune" |
                         bio_func %like% "athogen" |
                         bio_func %like% "efens"]$qseqid)

# Output list
write.table(out.ld, 
            file = "/project/berglandlab/madi/immune_genes_list_names",
            quote = F, 
            row.names = FALSE, 
            col.names = FALSE, 
            sep = "\t")


# Output list of genes - bed file
write.table(gene.gtf[gene %in% out.ld][sec=="transcript"] %>% 
              select(-c(file, sec, V6:V9)), 
            file = "/project/berglandlab/madi/immune_genes.bed",
            quote = F, 
            row.names = FALSE, 
            col.names = FALSE, 
            sep = "\t")