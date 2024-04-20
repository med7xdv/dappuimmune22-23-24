library(ggtree)
library(ggplot2)
library(tidyverse)
library(treeio)
library(RColorBrewer)

setwd("/Users/madidoceti/Desktop/lab/CurrentResearch/Trees+Alignments")

myd88tree <- read.tree("myd88.treefile")

patterns <- c("XP_0466", "XP_0464", "Daphnia", "KAI95", "XP_0573", "XP_032", "CAH0", "XP_045", "YP_0091")
species <- c("pulicaria", "NAMpulex", "Epulex", "sinesis", "carinata", "magna", "galeata", "magna", "carinata")
palette <- brewer.pal(n = 7, name = "YlOrRd")

tip_labels <- myd88tree$tip.label

speciesdf <- data.frame(labels = tip_labels, species = NA)


for (i in seq_along(patterns)) {
  matching_indices <- grepl(patterns[i], tip_labels)
  speciesdf$species[matching_indices] <- species[i]
}

psamp <- ggtree(myd88tree)

# Plot the tree with circular layout
p <- psamp %<+% speciesdf + 
  geom_tiplab(aes(fill = species), size = 3,
              offset = 0.5, 
              hjust = 1.2, 
              align = TRUE,
              color = "black", # color for label font
              geom = "label",  # labels not text
              label.padding = unit(0.15, "lines"), # amount of padding around the labels
              label.size = 0) + # size of label border
  scale_fill_manual(values = palette) + # Use the RColorBrewer palette
  theme(legend.title = element_text(size = 12), # Set legend title font size
        legend.text = element_text(size = 10), # Set legend text font size
        legend.key = element_blank())

msaplot(p, "OG0006048.trans.fa", offset = .5, color = c("khaki1", "lightgreen", "azure", "azure", "steelblue1", "azure", "azure", "azure", "azure", "azure", "lightpink1")) +
  ggtitle("MyD88") +
  guides(fill = FALSE)


