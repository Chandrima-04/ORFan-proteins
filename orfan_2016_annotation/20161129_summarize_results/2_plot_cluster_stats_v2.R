rm(list=ls())
library(tidyverse)
library(stringr)
library(Unicode)
library(RColorBrewer)
library(ggrepel)

# setwd("~/Box Sync/ORFan_sonnhammer/2016.03_protein_family_db_search/20161129_summarize_results")


theme_set(theme_bw())

#****************************************
# 1 Load cluster stats
#****************************************

cluster_stats <- read.csv("1_out/cluster_stats.csv",stringsAsFactors = FALSE)
colnames(cluster_stats)[2:3] <- c("n_seqs","n_raw_seqs")
colnames(cluster_stats)[6:7] <- c("RNAcode_p_value","ORF_length")
colnames(cluster_stats)

#Replace non-ascii hyphen character with a normal one for parsing to double
strange_hyphen <- cluster_stats$RNAcode_p_value[[1]] %>% substring(7,7)

cluster_stats$RNAcode_p_value <- str_replace_all(cluster_stats$RNAcode_p_value,strange_hyphen,"-") %>% as.numeric()

#****************************************
# Load cluster hits
#****************************************
nr_hmp_hits <- read_delim("2_out/nr_hmp_hits.csv",delim=",")

all_hit_summary <- read_csv("2_out/dbhits_x_cluster.csv")

nr_hmp_hits.filt <- nr_hmp_hits %>% filter(tool=="hmmsearch") %>%
  mutate(species=( str_extract(s_description,'\\[.+?\\]') %>% str_replace_all('\\[|\\]',''))) %>%
  mutate(species=str_replace(species,"Chlamydia trachomatis","TTV-like virus")) %>% select(cluster,species)

metagenomic_hits <- all_hit_summary %>%
  select(Cluster,dbs) %>%
  mutate(metahit=grepl("metahit_2014",dbs),ncbi_env=grepl("env_",dbs),no_annot=is.na(dbs)) %>%
  mutate(in_metagenome=metahit | ncbi_env) %>%
  mutate(annotation=ifelse(no_annot,"Unannotated",ifelse(in_metagenome,"Public metagenomes","Other dbs(nr,Pfam)"))) %>%
  select(-dbs)

#****************************************
# Plot
#****************************************
brewer.pal(3,"Set1")

cluster.info <- cluster_stats %>%
  inner_join(metagenomic_hits,by="Cluster") %>%
  left_join(nr_hmp_hits.filt,by=c("Cluster"="cluster"))

ggplot(cluster.info,aes(x=ORF_length,y=RNAcode_p_value)) +
  geom_point(aes(size=n_seqs,color=annotation),alpha=0.7) +
  geom_text_repel(aes(label=species), point.padding = unit(1, 'lines'), Segment.color = '#cccccc', segment.size = 0.5) +
  scale_color_manual(values=brewer.pal(3,"Set1")) +
  scale_size_area(limits=c(0,12)) +
  xlab("ORF length") +
  ylab("RNAcode p-value")
ggsave("orfan_proteins_summary.pdf",width=7,height=5,units="in")
