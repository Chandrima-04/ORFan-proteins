rm(list=ls())

library(tidyverse)
library(stringr)
library(UpSetR)
library(RColorBrewer)

setwd("~/Box Sync/ORFan_sonnhammer/2017.07_scireports_rev1/2_benchmark_method")

tools.df <- read_delim("./1_merged_viruses.tsv",delim="\t")
tools.df.withzeros <- tools.df %>% spread(tool,abundance,fill=0) %>% gather(key="tool",value="abundance",-clade)

theme_set(theme_bw())

#Calculate relative abundance
tools.rel_ab <- tools.df.withzeros %>%
  group_by(tool) %>% summarise(total.viral.reads=sum(abundance)) %>%
  inner_join(tools.df.withzeros,"tool") %>% mutate(rel.ab=abundance/total.viral.reads)

#Filter out viruses that are only called by kaiju and nothing else
tools.rel_ab.filt <- tools.rel_ab %>%
  group_by(clade) %>% filter(rel.ab > 0) %>% count() %>%
  inner_join(tools.rel_ab,"clade") %>%
  filter(tool != "kaiju" || n > 1)

#Be able to split
tools.rel_ab.filt2 <- tools.rel_ab.filt %>% filter(tool=="ORFan") %>% filter(abundance > 0 ) %>% select(clade) %>%
  mutate(in_ORFan=TRUE) %>% right_join(tools.rel_ab.filt,"clade") %>% mutate(in_ORFan=!is.na(in_ORFan)) %>% filter(tool != "ORFan")

# FAcets and common-coord plot
# tools.rel_ab.filt2 %>% ggplot(aes(rel.ab,tool)) + geom_point(aes(size=rel.ab)) + facet_wrap(~clade)

#Use facets & barplots
tools.rel_ab.filt2 %>% filter(!grepl("POOphage",clade)) %>% filter(in_ORFan) %>%
  ggplot(aes(tool,rel.ab)) +
  geom_bar(stat="identity",aes(fill=tool)) +
  geom_text(aes(label=paste("(",abundance,")",sep="")),hjust=-0.5,size=2) +
  facet_wrap(~clade) + coord_flip(ylim=c(0,1)) +
  scale_fill_manual(drop=FALSE,values=brewer.pal(3,"Set2")) +
  scale_y_continuous(labels=scales::percent) +
  xlab("") + ylab("") + guides(fill=FALSE) +
  theme(axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=6),
        strip.text.x = element_text(size = 7))
ggsave("virus_in_orfan.pdf",width=6,height=1.5,units="in")

tools.rel_ab.filt2 %>% filter(!grepl("POOphage",clade)) %>% filter(!in_ORFan) %>%
  ggplot(aes(tool,rel.ab)) +
  geom_bar(stat="identity",aes(fill=tool)) +
  geom_text(aes(label=paste("(",abundance,")",sep="")),hjust=-0.5,size=2) +
  facet_wrap(~clade) + coord_flip(ylim=c(0,1)) +
  scale_fill_manual(drop=FALSE,values=brewer.pal(3,"Set2")) +
  scale_y_continuous(labels=scales::percent) +
  xlab("") + ylab("") + guides(fill=FALSE) +
  theme(axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.title.y=element_text(size=6),
        strip.text.x = element_text(size = 7))
ggsave("virus_not_in_orfan.pdf",width=6,height=2.4,units="in")


# # all viruses in same plot
# tools.rel_ab.filt2 %>% filter(!grepl("POOphage",clade)) %>% filter(in_ORFan) %>%
#   ggplot(aes(clade,rel.ab)) + geom_bar(stat="identity",position="dodge",aes(fill=tool))
#
# tools.rel_ab.filt2 %>% filter(!grepl("POOphage",clade)) %>% filter(!in_ORFan) %>%
#   ggplot(aes(clade,rel.ab)) + geom_bar(stat="identity",position="dodge",aes(fill=tool))
