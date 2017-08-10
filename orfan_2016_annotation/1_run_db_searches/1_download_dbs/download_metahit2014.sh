#!/bin/bash

ORIGINAL_PWD=$PWD

mkdir -p fasta/metahit_2014

cd fasta/metahit_2014
wget -N ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/1.GeneCatalogs/IGC.fa.gz
wget -N ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/1.GeneCatalogs/IGC.pep.gz

wget -N ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/3.IGC.AnnotationInfo/IGC.annotation_OF.summary.gz
wget -N ftp://climb.genomics.cn/pub/10.5524/100001_101000/100064/6.PreviousGeneCatalog/PGC.gene.KO.list.gz
