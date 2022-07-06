#!/usr/bin/env Rscript

#clear the object from memory
rm(list=ls())

# install required libraries before running this script

# load libraries
library(dada2)
library(Biostrings)
library(ShortRead)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
library(phyloseq)
library(xlsx)
library(vegan)
library(tibble)
library(dplyr)
library(DECIPHER)
library(seqRFLP)

#import paths
paths <- read.delim("paths.txt", header=F)
source <- paths[1,1]
taxa_database <- paths[3,1] #fungi db

#set working directory to source file location
setwd(source)

#generate output files
output <- T
#assign taxonomy
taxonomy <- T

#list of fungi runs
runs <- list.files("../cmd/primer_cutted/fungi")

#folder to store interim results
temp <- "fungi_interim/"
dir.create(temp)


# ------------------------------------------------------------------------
# Loop over all fungi runs
# ------------------------------------------------------------------------
for (run in runs) {
  print(paste(run, "started"))
  
  # ------------------------------------------------------------------------
  # get the samples
  # ------------------------------------------------------------------------
  
  #Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
  path <- paste("../cmd/primer_cutted/fungi/", run, sep="")
  fnFs <- sort(list.files(path, pattern="r1_cutted.fastq", full.names = TRUE))
  fnRs <- sort(list.files(path, pattern="r2_cutted.fastq", full.names = TRUE))
  
  #Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
  sample.names <- sapply(strsplit(basename(fnFs), "-r._", fixed=F), `[`, 1)
  
  #Place filtered files in filtered/ subdirectory
  filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
  filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names
  
  #give them a unique name
  assign(paste(run, "_fnFs", sep=""), fnFs)
  assign(paste(run, "_fnRs", sep=""), fnRs)
  assign(paste(run, "_filtFs", sep=""), filtFs)
  assign(paste(run, "_filtRs", sep=""), filtRs)
  
  # ------------------------------------------------------------------------
  #Quality filtering and trimming
  # ------------------------------------------------------------------------
  
  #filter
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,170),
                       maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                       compress=TRUE, multithread=TRUE)
  
  #unique name and save
  assign(paste(run, "_out", sep=""), out)
  saveRDS(out, paste(temp, run, "_out.RDS", sep=""))
  
  #Avoids to raise an error message if no reads from a sample pass the filter
  missing_fltFs <- c()
  missing_fltRs <- c()
  for (file in filtFs) {if(!file.exists(file)) {missing_fltFs <- c(missing_fltFs, file)}} #collect names of the missing files
  for (file in filtRs) {if(!file.exists(file)) {missing_fltRs <- c(missing_fltRs, file)}}
  filtFs <- setdiff(filtFs, missing_fltFs)
  filtRs <- setdiff(filtRs, missing_fltRs)

  # ------------------------------------------------------------------------
  # Learn the error rate
  # ------------------------------------------------------------------------
  
  #error rate
  errF <- learnErrors(filtFs, multithread=TRUE)
  errR <- learnErrors(filtRs, multithread=TRUE)
  
  #unique name and save
  assign(paste(run, "_errF", sep=""), errF)
  assign(paste(run, "_errR", sep=""), errR)
  saveRDS(errF, paste(temp, run, "_errF.RDS", sep=""))
  saveRDS(errR, paste(temp, run, "_errR.RDS", sep=""))
  
  # ------------------------------------------------------------------------
  # Dereplicate the files
  # ------------------------------------------------------------------------
  
  #dereplicate
  derepFs <- derepFastq(filtFs, verbose=TRUE)
  derepRs <- derepFastq(filtRs, verbose=TRUE)
  
  #unique name and save
  assign(paste(run, "_derepFs", sep=""), derepFs)
  assign(paste(run, "_derepRs", sep=""), derepRs)
  saveRDS(derepFs, paste(temp, run, "_derepFs.RDS", sep=""))
  saveRDS(derepRs, paste(temp, run, "_derepRs.RDS", sep=""))
  
  # ------------------------------------------------------------------------
  #Sample inference algorithm
  # ------------------------------------------------------------------------
  
  # inference algortihm
  dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
  dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
  
  #unique name and save
  assign(paste(run, "_dadaFs", sep=""), dadaFs)
  assign(paste(run, "_dadaRs", sep=""), dadaRs)
  saveRDS(dadaFs, paste(temp, run, "_dadaFs.RDS", sep=""))
  saveRDS(dadaRs, paste(temp, run, "_dadaRs.RDS", sep=""))
  
}

# ------------------------------------------------------------------------
#combine all fungi runs
# ------------------------------------------------------------------------

#combine and save
fnFs <- vector()
for(i in runs){fnFs <- c(fnFs, get(paste(i, "_fnFs", sep="")))}
fnRs <- vector()
for(i in runs){fnRs <- c(fnRs, get(paste(i, "_fnRs", sep="")))}
saveRDS(fnFs, paste0(temp, "fnFs.RDS"))
saveRDS(fnRs, paste0(temp, "fnRs.RDS"))

filtFs <- vector()
for(i in runs){filtFs <- c(filtFs, get(paste(i, "_filtFs", sep="")))}
filtRs <- vector()
for(i in runs){filtRs <- c(filtRs, get(paste(i, "_filtRs", sep="")))}
saveRDS(filtFs, paste0(temp, "filtFs.RDS"))
saveRDS(filtRs, paste0(temp, "filtRs.RDS"))

out <- get(paste(runs[1], "_out", sep=""))
if(length(runs)>1){for(i in runs[2:length(runs)]){out <- rbind(out, get(paste(i, "_out", sep="")))}}
saveRDS(out, paste0(temp, "out.RDS"))

derepFs <- vector()
for(i in runs){derepFs <- c(derepFs, get(paste(i, "_derepFs", sep="")))}
derepRs <- vector()
for(i in runs){derepRs <- c(derepRs, get(paste(i, "_derepRs", sep="")))}
saveRDS(derepFs, paste0(temp, "derepFs.RDS"))
saveRDS(derepRs, paste0(temp, "derepRs.RDS"))

dadaFs <- vector()
for(i in runs){dadaFs <- c(dadaFs, get(paste(i, "_dadaFs", sep="")))}
dadaRs <- vector()
for(i in runs){dadaRs <- c(dadaRs, get(paste(i, "_dadaRs", sep="")))}
saveRDS(dadaFs, paste0(temp, "dadaFs.RDS"))
saveRDS(dadaRs, paste0(temp, "dadaRs.RDS"))

#names from all fungi runs
sample.names <- sapply(strsplit(basename(fnFs), ".", fixed=T), `[`, 1)

# ------------------------------------------------------------------------
# Merge the forward and reverse reads
# ------------------------------------------------------------------------

#merge
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
saveRDS(mergers, paste0(temp, "mergers.RDS"))

# ------------------------------------------------------------------------
# ASV table
# ------------------------------------------------------------------------

#create table and remove bimera
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab, paste0(temp, "seqtab.RDS"))
saveRDS(seqtab.nochim, paste0(temp, "seqtab_nochim.RDS"))

#rename rownames
rownames(seqtab.nochim) <- gsub("-ITS2_(F|R)_filt.fastq","",rownames(seqtab.nochim))
rownames(seqtab.nochim) <- gsub("-ITS1F-runBE[0-9]{2}","",rownames(seqtab.nochim))

# ------------------------------------------------------------------------
# track table
# ------------------------------------------------------------------------

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
saveRDS(track, paste0(temp, "track.RDS"))

# ------------------------------------------------------------------------
# Assign taxanomy
# ------------------------------------------------------------------------
if (taxonomy) {
  taxa <- assignTaxonomy(seqtab.nochim, taxa_database, multithread=TRUE)
  saveRDS(taxa, paste(temp, "taxa.RDS", sep=""))
}

# ------------------------------------------------------------------------
# Header names
# ------------------------------------------------------------------------

#change headers from the long sequence to simple ASV numbers
seq <- colnames(seqtab.nochim)
ASV <- c()
for (i in 1:length(seq)) {ASV[i] <- paste("ASV",i , sep="") }
ASV_seq <- data.frame(ASV=ASV, seq=seq)
#ASV_seq[,1] <- paste(">", ASV_seq[,1], sep="")
saveRDS(ASV_seq, paste0(temp, "sequences.RDS"))

#replace sequence headers by ASV numbers in seqtab.nochim
colnames(seqtab.nochim) <- ASV

if (taxonomy) {
  #replace the sequence headers by the correct ASV number in taxa
  for (i in 1:nrow(taxa)){
    for (j in 1:nrow(ASV_seq)){
      if(rownames(taxa)[i] == ASV_seq[j,2]){
        rownames(taxa)[i] <- ASV_seq[j,1]
      }
    }
  }
}

# ------------------------------------------------------------------------
# Clustering
# ------------------------------------------------------------------------

nproc <- 16

#change ASV names to corresponding sequences
seqtab.nochim_ASV <- seqtab.nochim
colnames(seqtab.nochim_ASV) <- ASV_seq$seq

#align sequences and compute distance matrix
ASV_seqs <- Biostrings::DNAStringSet(ASV_seq$seq) #convert to DNA string set
ASV_seqs_aln <- DECIPHER::AlignSeqs(ASV_seqs, processors = nproc)
ASV_dist <- DECIPHER::DistanceMatrix(ASV_seqs_aln, processors = nproc)

#OTUs
ASV97_cluster <- DECIPHER::IdClusters(ASV_dist, method = "complete", processors = nproc, cutoff = 0.03) #cutoff = 0.03 -> 97% OTU 
ASV98_cluster <- DECIPHER::IdClusters(ASV_dist, method = "complete", processors = nproc, cutoff = 0.02) #cutoff = 0.02 -> 98% OTU 
ASV99_cluster <- DECIPHER::IdClusters(ASV_dist, method = "complete", processors = nproc, cutoff = 0.01) #cutoff = 0.01 -> 99% OTU

#loop over the 3 thresholds
for (n in c("97", "98", "99")) {
  
  #get the data
  ASV_cluster <- get(paste0("ASV", n, "_cluster"))
  ASV_cluster$ASV <- gsub(">","",ASV_seq$ASV) #add ASV name
  ASV_cluster$ASV_seq <- ASV_seq$seq #add sequences
  
  # Number of sequences per ASV
  ASV_cluster$ASV_abu <- NA
  for (seq in ASV_cluster$ASV_seq) {
    abu <- sum(seqtab.nochim_ASV[,colnames(seqtab.nochim_ASV) == seq])
    ASV_cluster$ASV_abu[ASV_cluster$ASV_seq == seq] <- abu
  }
  
  #most abundant sequences for each OTU
  ASV_cluster$OTU_seq <- NA
  ASV_cluster$OTU_abu <- NA
  
  for (OTU in ASV_cluster$cluster) {
    ASV_cluster_subset <- ASV_cluster[ASV_cluster$cluster == OTU,] #all entries for current OTU
    n_seq <- sum(ASV_cluster_subset$ASV_abu)
    most_abu_seq <- ASV_cluster_subset$ASV_seq[ASV_cluster_subset$ASV_abu == max(ASV_cluster_subset$ASV_abu)]# OTU(s) with most seqs
    most_abu_seq <- most_abu_seq[1] #take the first entry in case multiple OTUs have max number of seqs
    ASV_cluster$OTU_seq[ASV_cluster$cluster == OTU] <- most_abu_seq #add seq to table
    ASV_cluster$OTU_abu[ASV_cluster$cluster == OTU] <- n_seq #add num to table
  }
  
  #sort by OTU abundance
  ASV_cluster_sort <- ASV_cluster[order(ASV_cluster$OTU_abu, decreasing = T),] #sort by OTU name
  
  #rename OTU by abundance
  clusters_sort <- unique(ASV_cluster_sort$cluster)
  ASV_cluster$OTU <- NA
  for(i in seq(length(clusters_sort))){
    ASV_cluster$OTU[ASV_cluster$cluster == clusters_sort[i]] <- paste0("OTU", i) #rename
  }
  
  ### FINAL OTU FILES
  
  #OTU ASV_seq
  OTU_ASV_seq_all <- ASV_cluster[order(ASV_cluster$OTU_abu, decreasing = T),]
  OTU_ASV_seq_all <- data.frame(OTU = OTU_ASV_seq_all$OTU, seq = OTU_ASV_seq_all$OTU_seq)
  OTU_ASV_seq <- data.frame(OTU = unique(OTU_ASV_seq_all$OTU), seq = NA)
  for (OTU in OTU_ASV_seq$OTU) {
    seq <- OTU_ASV_seq_all$seq[OTU_ASV_seq_all$OTU == OTU][1] #get the sequence
    OTU_ASV_seq$seq[OTU_ASV_seq$OTU == OTU] <- seq #add to seqtab.nochima frame
  }
  
  
  #OTU seqtab.nochim
  OTU_seqtab.nochim <- seqtab.nochim_ASV %>% t %>% rowsum(ASV_cluster$cluster) %>% t #cluster
  for(i in seq(length(clusters_sort))){
    colnames(OTU_seqtab.nochim)[colnames(OTU_seqtab.nochim) == clusters_sort[i]] <- paste0("OTU", i) #rename
  }
  OTU_seqtab.nochim <- OTU_seqtab.nochim[,OTU_ASV_seq$OTU] #sort
  
  
  #OTU taxa
  OTU_taxa <- data.frame(matrix(nrow = nrow(OTU_ASV_seq), ncol = ncol(taxa)))
  rownames(OTU_taxa) <- OTU_ASV_seq$OTU
  colnames(OTU_taxa) <- colnames(taxa)
  
  for(OTU in OTU_ASV_seq$OTU){
    seq <- OTU_ASV_seq$seq[OTU_ASV_seq$OTU == OTU] #OTU seq
    ASV <- ASV_cluster$ASV[ASV_cluster$ASV_seq == seq] #ASV belong to OTU seq
    
    taxa_ASV <- taxa[rownames(taxa) == ASV[1], ] #taxa belong to the first ASV
    OTU_taxa[rownames(OTU_taxa) == OTU, ] <- taxa_ASV #transfer to OTU taxa
  }
  
  #OTU_ASV
  OTU_ASV <- data.frame(OTU = ASV_cluster$OTU, ASV = ASV_cluster$ASV, 
                        OTU_seq = ASV_cluster$OTU_seq, ASV_seq = ASV_cluster$ASV_seq,
                        OTU_abu = ASV_cluster$OTU_abu, ASV_abu = ASV_cluster$ASV_abu)
  
  OTU_ASV <- OTU_ASV[order(ASV_cluster$OTU_abu, decreasing = T),] #sort
  
  #store clustered data
  assign(paste0("OTU_ASV_seq", n), OTU_ASV_seq)
  assign(paste0("OTU_seqtab.nochim", n), OTU_seqtab.nochim)
  assign(paste0("OTU_taxa", n), OTU_taxa)
  assign(paste0("OTU_ASV", n), OTU_ASV)
  
}

# ------------------------------------------------------------------------
# Output
# ------------------------------------------------------------------------

if(output){
  
  #QUALITY CHECK
  
  setwd(source)
  dir.create("../../4_output/fungi_ASV")
  dir.create("../../4_output/fungi_ASV/quality_check")
  
  path_out <- c("../../4_output/fungi_ASV/")
  path_out_qc <- c("../../4_output/fungi_ASV/quality_check/")
  
  #quality profiles unfiltered and untrimmed
  pdf(paste(path_out_qc,"fungi_reads_quality_unfilt_untrim.pdf",sep=""))
  for (i in 1:length(sample.names)) {
    try(figure <- plotQualityProfile(c(fnFs[i],fnRs[i])))
    try(print(figure))
    try(rm(figure))
  }
  dev.off()
  
  #quality profiles filtered and trimmed
  pdf(paste(path_out_qc,"fungi_reads_quality_filt_trim.pdf",sep=""))
  for (i in 1:length(sample.names)) {
    try(figure <- plotQualityProfile(c(filtFs[i],filtRs[i])))
    try(print(figure))
    try(rm(figure))
  }
  dev.off()
  
  #error rate learning
  runs_errF <- c (paste(runs, "_errF", sep=""))
  runs_errR <- c (paste(runs, "_errR", sep=""))
  
  pdf(paste(path_out_qc,"fungi_error_rate.pdf",sep=""))
  for (i in 1:length(runs)) {
    plot(plotErrors(get(runs_errF[i]), nominalQ=TRUE))
    grid.text(paste(runs[i], "forward"),hjust=-1.2, vjust = -27.5, rot = 90)
    plot(plotErrors(get(runs_errR[i]), nominalQ=TRUE))
    grid.text(paste(runs[i], "reverse"),hjust=-1.2, vjust = -27.5, rot = 90)
  }
  dev.off()
  
  #rarefaction plots
  pdf(paste(path_out_qc,"fungi_rarefaction.pdf",sep=""))
  for (run in runs) {
    DAT <- seqtab.nochim[grep(run, rownames(seqtab.nochim)),]
    rarecurve(DAT, step = 20, label=F, main = paste("Rarefaction", run), ylab="ASVs", xlab="Sequencing depth")
  }
  dev.off()
  
  #number of reads
  write.xlsx(track, paste(path_out_qc, "fungi_track.xlsx", sep=""))
  
  # ASV TABLES
  
  #97
  dir.create(paste0(path_out,"ASV97"))
  write.table(OTU_seqtab.nochim97, paste0(path_out, "ASV97/fungi_DAT97.tab"), sep="\t")
  write.table(OTU_taxa97, paste0(path_out, "ASV97/fungi_TAXA97.tab"), sep="\t")
  write.table(OTU_ASV_seq97, paste0(path_out, "ASV97/fungi_SEQ97.tab"), sep="\t", row.names = F)
  dataframe2fas(OTU_ASV_seq97, paste0(path_out, "ASV97/fungi_SEQ97.fasta")) #library seqRFLP
  write.table(OTU_ASV97, paste0(path_out, "ASV97/fungi_ABU97.tab"), sep="\t", row.names = F)
  
  #98
  dir.create(paste0(path_out,"ASV98"))
  write.table(OTU_seqtab.nochim98, paste0(path_out, "ASV98/fungi_DAT98.tab"), sep="\t")
  write.table(OTU_taxa98, paste0(path_out, "ASV98/fungi_TAXA98.tab"), sep="\t")
  write.table(OTU_ASV_seq98, paste0(path_out, "ASV98/fungi_SEQ98.tab"), sep="\t", row.names = F)
  dataframe2fas(OTU_ASV_seq98, paste0(path_out, "ASV98/fungi_SEQ98.fasta")) #library seqRFLP
  write.table(OTU_ASV98, paste0(path_out, "ASV98/fungi_ABU98.tab"), sep="\t", row.names = F)
  
  #99
  dir.create(paste0(path_out,"ASV99"))
  write.table(OTU_seqtab.nochim99, paste0(path_out, "ASV99/fungi_DAT99.tab"), sep="\t")
  write.table(OTU_taxa99, paste0(path_out, "ASV99/fungi_TAXA99.tab"), sep="\t")
  write.table(OTU_ASV_seq99, paste0(path_out, "ASV99/fungi_SEQ99.tab"), sep="\t", row.names = F)
  dataframe2fas(OTU_ASV_seq99, paste0(path_out, "ASV99/fungi_SEQ99.fasta")) #library seqRFLP
  write.table(OTU_ASV99, paste0(path_out, "ASV99/fungi_ABU99.tab"), sep="\t", row.names = F)
  
  #100
  dir.create(paste0(path_out,"ASV100"))
  write.table(seqtab.nochim, paste0(path_out, "ASV100/fungi_DAT100.tab"), sep="\t")
  write.table(taxa, paste0(path_out, "ASV100/fungi_TAXA100.tab"), sep="\t")
  #write.table(taxa_conf, paste0(path_out, "ASV100/fungi_TAXA100_conf.tab"), sep="\t")
  write.table(ASV_seq, paste0(path_out, "ASV100/fungi_SEQ100.tab"), sep="\t")
  dataframe2fas(ASV_seq, paste0(path_out, "ASV100/fungi_SEQ100.fasta")) #library seqRFLP
  
}
