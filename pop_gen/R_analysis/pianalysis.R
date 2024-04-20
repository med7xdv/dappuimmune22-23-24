# Read in exon output and visualize
# Connor Murray 10.13.2023 and madi 
# ijob -c 10 --mem=10G -p standard -A berglandlab
# module load goolf/7.1.0_3.1.4 R/4.0.3; module load gdal geos proj; R

.libPaths(c("/scratch/med7xdv/R-packages/", .libPaths()))
# Packages
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(tibble)
library(readr)
library(stringr)
library(foreach)
library(readxl)
library(cowplot)
library(patchwork)
library(RColorBrewer)
library(ggsignif)
require(scales)

# Working directory
setwd('/scratch/med7xdv/')

# Sample metadata
meta <- data.table(fread("samples.9.8.22.csv"))

# Read in gtf file
gene.gtf <- data.table(fread("/project/berglandlab/daphnia_ref/Daphnia.aed.0.6.gtf"))
colnames(gene.gtf)[1:5] <- c("chrom", "file", "sec", "start", "stop")

extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- strsplit(gtf_attributes, "; ")
  att <- gsub("\"","",unlist(att))
  if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
    return( unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
  }else{
    return(NA)}
}
gene.gtf$gene <- unlist(str_remove_all(lapply(gene.gtf$V9, extract_attributes, "transcript_id"), pattern = ";"))

# Read in gene annotations
panth <- data.table(read_excel("/project/berglandlab/daphnia_ref/Daphnia_annotation_PANTHER.xls"))
colnames(panth)[36:38] <- c("bio_func", "mol_func", "cell_func")

panthfunc <- panth %>%
  select(qseqid, bio_func, mol_func, cell_func) %>%
  rename(gene = qseqid)

# PI
files <- system("ls -f -R /scratch/med7xdv/popgen/allpiout/*", intern = TRUE)
listi <- lapply(files, fread)
setattr(listi, 'names', list.files(path = "/scratch/med7xdv/popgen/allpiout/", pattern = "pi_nomiss"))

# Bind list 
pi.dt <- data.table(rbindlist(listi, use.names = T, idcol = T) %>% 
                      select(-c(".id")))
colnames(pi.dt)[1:8] <- c("chrom", "start", "stop", "num.sites", "pi", 
                          "exon.start", "exon.stop", "gene")


###################################################### my code ######################################################

# subset columns to genes and other helpful columns
merged <- merge(pi.dt, panthfunc, by = "gene", all.x = TRUE)

### grab only RA genes and classify genes as immune or non immune; THIS IS THE REAL SET OF GENES
filter1 <- merged %>%
  filter(grepl("-RA$", gene)) %>%
  mutate(exon.size = exon.stop - exon.start) %>%
  mutate(classification = ifelse(gene == "Daphnia11727-RA" | 
                                   gene == "Daphnia00294-RA" | 
                                   gene == "Daphnia00295-RA" | 
                                   gene == "Daphnia00298-RA" | 
                                   gene == "Daphnia06247-RA" | 
                                   gene == "Daphnia01717-RA" | 
                                   gene == "Daphnia11353-RA" | 
                                   gene == "Daphnia11354-RA" |
                                   gene == "Daphnia11812-RA" |
                                   gene == "Daphnia05615-RA" |
                                   gene == "Daphnia13369-RA" |
                                   gene == "Daphnia13364-RA" |
                                   gene == "Daphnia12722-RA" | 
                                   gene == "Daphnia04857-RA" |
                                   gene == "Daphnia06045-RA" |
                                   gene == "Daphnia08032-RA" |
                                   gene == "Daphnia10030-RA" |
                                   gene == "Daphnia02997-RA" |
                                   gene == "Daphnia07409-RA" |
                                   gene == "Daphnia04063-RA" |
                                   gene == "Daphnia04695-RA" |
                                   gene == "Daphnia10460-RA" |
                                   gene == "Daphnia03330-RA" |
                                   gene == "Daphnia12703-RA" |
                                   gene == "Daphnia06369-RA" |
                                   gene == "Daphnia12616-RA" |
                                   gene == "Daphnia10736-RA" |
                                   gene == "Daphnia03756-RA" |
                                   gene == "Daphnia08434-RA" |
                                   gene == "Daphnia08343-RA" |
                                   gene == "Daphnia09231-RA" |
                                   gene == "Daphnia08340-RA" |
                                   gene == "Daphnia12428-RA" |
                                   gene == "Daphnia05542-RA", "immune", "non-immune"))


# filter to genes that are present in each pond, number of exons, number of sites, and find avg pi per pond 
avgpiperpond <- filter1 %>%
  group_by(gene) %>%
  summarize(avgpi = mean(pi, na.rm = TRUE))

filter2 <- merge(avgpiperpond, filter1, by = "gene", all.x = TRUE)

# name exons and filter to distinct exons
pgdt <- as.data.table(filter2)
pgdt <- pgdt[,exon_id:=paste(gene, start, stop, sep = "_")]
pigenesunique <- distinct(pgdt) 


### test ###

onlyimmune <- pigenesunique %>%
  filter(classification == "immune")

# get random non-immune exons
sizefilter <- pigenesunique %>%
  filter(classification == "non-immune", exon.size > 47, exon.size < 4272)

#### immune vs genome avg ####
boot <- foreach(i = 1:1000, .combine = "rbind" )%dopar% {
   one_exon_per_gene <- onlyimmune %>%
    group_by(gene) %>%
    sample_n(size = 1, replace = FALSE) %>%
    ungroup()
  
  iexons <- one_exon_per_gene %>%
    sample_n(size = 10, replace = FALSE)
  
  one_niexon_per_gene <- sizefilter %>%
    group_by(gene) %>%
    sample_n(size = 1, replace = FALSE) %>%
    ungroup()
  
  niexons <- one_niexon_per_gene %>%
    sample_n(size = 10, replace = FALSE)
  
  ianiexons <- rbind(niexons, iexons)
  
  t_test_result <- kruskal.test(pi ~ classification, data = ianiexons)
  
  #pvals <-
  boot <- data.table(pval = t_test_result$p.value, t = t_test_result$statistic, iteration = i, NonImmune_Average_Pi = mean(niexons$pi), Immune_Average_Pi = mean(iexons$pi))
  
}

boot <- boot %>%
  pivot_longer(cols = c(NonImmune_Average_Pi, Immune_Average_Pi))

p <- ggplot(boot, aes(x = t, fill = "plum1")) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = 0) +
  labs(x = "t Statistic", y = "Count") +
  ggtitle("t Statistics for 1000 Bootstraps") +
  guides(fill = "none")


p1 <- ggplot(boot, aes(x = pval, fill = "peachpuff")) +
  geom_histogram(bins = 30) +
  labs(x = "P-Value", y = "Count") +
  ggtitle("P-Value for 1000 Bootstraps") +
  guides(fill = "none")

p2 <- ggplot(boot, aes(x = name, y = value, color = name)) +
  geom_boxplot() +
  scale_color_manual(values = c("Immune_Average_Pi"="royalblue1", "NonImmune_Average_Pi"="lightpink1")) +
  labs(x = "Classification", y = "Average Pi") +
  ggtitle("Average Immune versus Non-Immune Diversity") +
  guides(fill = "none") +
  theme_light()


layout <- "
AACC
AACC
BBCC
BBCC"

p3 <- p +
  p1 +
  p2 +
  plot_layout(design = layout, guides = "collect") +
  plot_annotation(tag_levels = "A", title = "Bootstrap Test for Immune versus Non-Immune Genes")
