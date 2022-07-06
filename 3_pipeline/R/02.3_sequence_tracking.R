#!/usr/bin/env Rscript

print("sequences tracking started")

#clear the object from memory
rm(list=ls())

# install required libraries before running this script

#libaries
library(tidyverse)
library(xlsx)

#paths
paths <- read.delim("paths.txt", header=F)
source <- paths[1,1]

#set working directory to source file location
setwd(source)


# ------------------------------------------------------------------------
# Import
# ------------------------------------------------------------------------

seqs <- read.table("../../2_data/seqs.txt", header = F) #import from dem.
taxa <- list.dirs("../cmd/primer_cutted",full.names = F)[-1] #get taxa

if ("bacteria" %in% taxa) {
  runs_bac <- list.files("../cmd/primer_cutted/bacteria")
  track_bac <- readRDS("bacteria_interim/track.RDS")
} 

if ("fungi" %in% taxa) {
  runs_fun <- list.files("../cmd/primer_cutted/fungi")
  track_fun <- readRDS("fungi_interim/track.RDS")
}

#run data.frame
runs_all <- unique(c(runs_bac, runs_fun))
runs <- data.frame(runs=runs_all, bac=F, fun=F)
try(runs$bac <- runs_bac %in% runs$runs, silent = T)
try(runs$fun <- runs_fun %in% runs$runs, silent = T)

##############
# track_bac <- readRDS("track_bac.RDS")
# track_fun <- track_bac / 20
# seqs <- read.table("seqs.txt", header = F)
# run="run04"
##############

# ------------------------------------------------------------------------
# Loop over all bacteria runs
# ------------------------------------------------------------------------

#track summary table
for (run in runs$runs) {
  
  # ------------------------------------------------------------------------
  # sequences tracking summary table
  # ------------------------------------------------------------------------
  
  #summary of read tracking
  run_seqs <- seqs[,3][seqs$V1 == run]
  run_dem_seqs_bac <- seqs[,5][seqs$V1 == run]
  run_dem_seqs_fun <- seqs[,7][seqs$V1 == run]
  ifelse(runs$bac[runs$runs==run], 
         run_track_bac <- colSums(track_bac[grepl(run,rownames(track_bac)),]), 
         run_track_bac <- rep(NA, 6))
  ifelse(runs$fun[runs$runs==run], 
         run_track_fun <- colSums(track_fun[grepl(run,rownames(track_fun)),]), 
         run_track_fun <- rep(NA, 6))

  run_track_sum <- data.frame(run=run,
                              step=c("raw", "demultiplexed", "input", "filtered", "denoised", "merged", "nonchim"),
                              bacteria=c(NA, run_dem_seqs_bac, run_track_bac[1:2], mean(run_track_bac[3:4]), run_track_bac[5:6]),
                              fungi=c(NA, run_dem_seqs_fun, run_track_fun[1:2], mean(run_track_fun[3:4]), run_track_fun[5:6])) #mean of denoising

  run_track_sum$step <- factor(run_track_sum$step, levels = run_track_sum$step)
  run_track_sum$total <- rowSums(run_track_sum[,3:4], na.rm = T)
  run_track_sum$total[1] <- run_seqs #add raw reads
  run_track_sum$percent <- paste(round(run_track_sum$total / run_seqs * 100, 1), "%")
  
  # ------------------------------------------------------------------------
  # Sequences tracking plot
  # ------------------------------------------------------------------------
  
  #long format for plot
  run_track_sum_long <- pivot_longer(run_track_sum, c(total, bacteria, fungi), names_to = "taxa", values_to = "seqs")
  run_track_sum_long$percent[run_track_sum_long$taxa != "total"] <- NA
  run_track_sum_long$taxa <- factor(run_track_sum_long$taxa, levels = c("total", "bacteria", "fungi"))

  #plot
  run_track_plot <- ggplot(run_track_sum_long,aes(step, seqs, group=taxa))+
                      geom_point(position=position_dodge(0.2))+
                      geom_line(aes(col=taxa),position=position_dodge(0.2))+
                      #geom_line(aes(step, bacteria), col="red")+
                      geom_text(aes(label=percent), hjust=0.2, vjust=-.5)+
                      scale_y_continuous(trans="log2",limits = c(min(run_track_sum$fungi)*0.5,max(run_track_sum$total)*1.05))+
                      theme_bw()+
                      ggtitle(paste("Sequences tracking:", run))
  
  
  #store for each run
  assign(paste0(run, "_track_sum"), run_track_sum)
  assign(paste0(run, "_track_plot"), run_track_plot)
  
}

# ------------------------------------------------------------------------
# Output
# ------------------------------------------------------------------------

setwd(source)
dir.create("../../4_output/sequences_tracking")
path_out <- "../../4_output/sequences_tracking/"


#tables
for (run in runs$runs) {
  if(run == runs[1]){track_sum <- get(paste0(run,"_track_sum"))}
  else{track_sum <- rbind(track_sum, get(paste0(run,"_track_sum")))}
}
write.xlsx(track_sum, paste0(path_out, "sequences_tracking.xlsx"))

#plots
pdf(paste0(path_out,"sequences_tracking.pdf"))
for (run in runs$runs) {
  plot(get(paste0(run, "_track_plot")))
}
dev.off()

print("sequences tracking completed")

