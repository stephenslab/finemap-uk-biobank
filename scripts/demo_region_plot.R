# TO DO: Briefly explain here what this script is for, and how to use
# it.
library(readr)
library(ggplot2)
library(cowplot)

# Load the gene data.
genes        <- read_delim("../data/seq_gene.md.gz",delim = "\t",quote = "")
class(genes) <- "data.frame"
genes        <- subset(genes,
                       group_label == "GRCh37.p5-Primary Assembly" &
                       feature_type == "GENE")

# Load the SNP data and GeneATLAS summary statistics.
geneatlas        <- read_csv("../data/sample_locuszoom_data.csv.gz")
class(geneatlas) <- "data.frame"

# Get the genes overlapping with the region.
start.pos <- min(geneatlas$pos)
stop.pos  <- max(geneatlas$pos)
plot.genes <- subset(genes,
                     chromosome == 3 &
                     ((chr_start > start.pos & chr_start < stop.pos) |
                      (chr_stop > start.pos & chr_start < stop.pos)))
n               <- nrow(plot.genes)
plot.genes$vpos <- rep(1:12,length.out = n)

# Plot the association p-values for all SNPs, and compare against the
# gene coding regions.
geneatlas <- transform(geneatlas,
                       pos    = pos/1e6,
                       pvalue = -log10(pvalue + 1e-300))
plot.genes <- transform(plot.genes,
                        chr_start = chr_start/1e6,
                        chr_stop  = chr_stop/1e6)
p <- ggplot(geneatlas,aes(x = pos,y = pvalue)) +
  geom_point(color = "darkblue",size = 1) +
  scale_x_continuous(breaks = plot.genes$chr_start,
                     labels = plot.genes$feature_name) +
  labs(x = "base-pair position on chromosome 3 (Mb)",
       y = "-log10 p-value") + 
  theme_cowplot(font_size = 10) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        axis.line   = element_blank()) 
