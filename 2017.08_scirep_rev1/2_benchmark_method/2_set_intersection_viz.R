rm(list=ls())

library(tidyverse)
library(stringr)
library(UpSetR)
library(RColorBrewer)

setwd("~/Box Sync/ORFan_sonnhammer/2017.07_scireports_rev1/2_benchmark_method")

tools.df <- read_delim("./1_merged_viruses.tsv",delim="\t")

tools.upset <- tools.df %>% mutate(present = 1L) %>%
  select(-abundance) %>%
  #mutate(clade=factor(clade)) %>%
  spread(tool,present,fill=0L) %>% data.frame

#text.scale:
#c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars).
tools <- c("metaphlan2","kraken","kaiju","ORFan")
#tools <- rev(tools)

my.colors <- brewer.pal(5,"Set1")


q1 <- list(query=intersects,
           color=my.colors[2],
           params=list("ORFan"),
           active=TRUE)

q2 <- list(query=intersects,
           color=my.colors[2],
           params=list("ORFan","kaiju"),
           active=TRUE)

q3 <- list(query=intersects,
           color=my.colors[2],
           params=list("ORFan","kraken","kaiju"),
           active=TRUE)

queries <- list(q1,q2,q3)

#NEEDS illustrator tweaking to put bars in proper order
# Ideally, should be, ORFan first, non-orfan second
# and within each group, ordered by degree (most shared to less shared)
pdf(file="method_validation_tools.pdf",width=7,height=5,useDingbats=FALSE)
upset(tools.upset,sets=tools,
      keep.order = T,
      order.by = "degree",
      point.size = 4, line.size = 1.5,
      mainbar.y.label = "Viruses in common", sets.x.label = "Detected viruses",
      text.scale =1.5,
      queries=queries)
dev.off()

