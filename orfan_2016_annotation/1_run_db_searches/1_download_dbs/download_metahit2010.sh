#!/bin/bash

ORIGINAL_PWD=$PWD

mkdir -p fasta/metahit

cd fasta/metahit
wget -N ftp://public.genomics.org.cn/BGI/gutmeta/UniSet/UniGene.cds.gz
wget -N ftp://public.genomics.org.cn/BGI/gutmeta/UniSet/UniGene.pep.gz
wget -N ftp://public.genomics.org.cn/BGI/gutmeta/UniSet/Annotation/*.gz
wget -N ftp://public.genomics.org.cn/BGI/gutmeta/Unique_gene_set/Gene2Sample.list.gz

wget -N ftp://public.genomics.org.cn/BGI/gutmeta/UniGene_Taxonomy/UniqGene_NR.tax.catalog
wget -N ftp://public.genomics.org.cn/BGI/gutmeta/NOG/UniGene.nog.catelog
wget -N ftp://public.genomics.org.cn/BGI/gutmeta/KEGG/UniGene.kegg.map.catalog
wget -N ftp://public.genomics.org.cn/BGI/gutmeta/KEGG/UniGene.kegg.ko.catalog
