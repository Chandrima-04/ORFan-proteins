rm(list=ls())
library(dplyr)
library(stringr)
library(ggplot2)
library(Unicode)

#1 Load cluster stats
cluster_stats <- read.csv("1_out/cluster_stats.csv",stringsAsFactors = FALSE)
colnames(cluster_stats)[2:3] <- c("n_seqs","n_raw_seqs")
colnames(cluster_stats)[6:7] <- c("RNAcode_p_value","ORF_length")
colnames(cluster_stats)

#Replace non-ascii hyphen character with a normal one for parsing to double
strange_hyphen <- cluster_stats$RNAcode_p_value[[1]] %>% substring(7,7)

cluster_stats$RNAcode_p_value <- str_replace_all(cluster_stats$RNAcode_p_value,strange_hyphen,"-") %>% as.numeric()

#dev.off()
plot_color <- "#76bf72"
dir.create("2_plots")

plt <- ggplot(cluster_stats,aes(x="Clusters",y=n_seqs)) + geom_boxplot(fill=plot_color) +
   xlab("") + ylab("Sequences in cluster") + theme_bw()
ggsave("2_plots/orfan_seqs_x_cluster_dist.pdf",width=3.5,height=7,plt,dpi=300)

plt <- ggplot(cluster_stats,aes(x="Clusters",y=RNAcode_p_value)) + geom_boxplot(fill=plot_color) +
  xlab("") + ylab("RNAcode p-values") + theme_bw()
ggsave("2_plots/orfan_rnacode_pval_dist.pdf",width=3.5,height=7,plt,dpi=300)

plt <- ggplot(cluster_stats,aes(x="Clusters",y=ORF_length)) + geom_boxplot(fill=plot_color) +
  xlab("") + ylab("ORF length(aa)") + theme_bw()
ggsave("2_plots/orfan_orf_len_dist.pdf",width=3.5,height=7,plt,dpi=300)
