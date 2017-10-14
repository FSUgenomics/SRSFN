################################################################################
library(data.table) #Loading library for fast and memory efficien file loading:
library(fdasrvf)
library(parallel) 
rm(list=ls(all=TRUE));

TSSregions.bed.name="data/chr1.TSSgencode.bed";
filenames=c("data/chr1.GH1_cap_R1_clip.fastq_q20.bed", "data/chr1.GH2_cap_R1_clip.fastq_q20.bed",
            "data/chr1.IH1_cap_R1_clip.fastq_q20.bed", "data/chr1.IH2_cap_R1_clip.fastq_q20.bed")

results<-SRSFN_wrapper(filenames,TSSregions)
