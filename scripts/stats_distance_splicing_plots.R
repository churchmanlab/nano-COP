install.packages("rrcov")
library(rrcov)
library(dplyr)
library(tidyr)


# This script takes as input datasets produced in the scripts from Figures 2, 3, S5 and S7

Dir="/path/to/datasets/from_figures_2-3-S5-S7"


################

### Human vs. Drosophila (Fig. 2B)
human_splice<-read.table(paste(Dir,"human_distance_transcribed_med_plot_df.txt",sep=""),stringsAsFactors=F,header=T)
droso_splice<-read.table(paste(Dir,"drosophila_distance_transcribed_med_plot_df.txt",sep=""),stringsAsFactors=F,header=T)

human_splice$Group <- "human"
droso_splice$Group <- "drosophila"

fields <- c("range", "name", "percentSpliced", "Group")

test_df <- rbind(human_splice,droso_splice) %>% select(fields)

res.aov2 <- aov(percentSpliced ~ Group + range, data = test_df)
summary(res.aov2)
#<2e-16


#####################


### Drosophila short vs. long introns (Fig. 2D)
droso_short_splice<-read.table(paste(Dir,"drosophila_distance_transcribed_short_introns_med_plot_df.txt",sep=""),stringsAsFactors=F,header=T)
droso_long_splice<-read.table(paste(Dir,"drosophila_distance_transcribed_long_introns_med_plot_df.txt",sep=""),stringsAsFactors=F,header=T)

droso_short_splice$Group <- "short"
droso_long_splice$Group <- "long"

test_df <- rbind(droso_short_splice,droso_long_splice) %>% select(fields)

res.aov2 <- aov(percentSpliced ~ Group + range, data = test_df)
summary(res.aov2)
#3.05e-14

### Human short vs. long introns (Fig. 2D)
human_short_splice<-read.table(paste(Dir,"human_distance_transcribed_short_introns_med_plot_df.txt",sep=""),stringsAsFactors=F,header=T)
human_long_splice<-read.table(paste(Dir,"human_distance_transcribed_long_introns_med_plot_df.txt",sep=""),stringsAsFactors=F,header=T)

human_short_splice$Group <- "short"
human_long_splice$Group <- "long"

test_df <- rbind(human_short_splice,human_long_splice) %>% select(fields)

res.aov2 <- aov(percentSpliced ~ Group + range, data = test_df)
summary(res.aov2)
#0.20024

#####################


### Drosophila const vs. alt introns (Fig. 2E)
droso_const_splice<-read.table(paste(Dir,"drosophila_distance_transcribed_const_introns_med_plot_df.txt",sep=""),stringsAsFactors=F,header=T)
droso_alt_splice<-read.table(paste(Dir,"drosophila_distance_transcribed_alt_introns_med_plot_df.txt",sep=""),stringsAsFactors=F,header=T)

droso_const_splice$Group <- "const"
droso_alt_splice$Group <- "alt"

test_df <- rbind(droso_const_splice,droso_alt_splice) %>% select(fields)

res.aov2 <- aov(percentSpliced ~ Group + range, data = test_df)
summary(res.aov2)
#1.15e-06


### Human const vs. alt introns (Fig. 2E)
human_const_splice<-read.table(paste(Dir,"human_distance_transcribed_const_introns_med_plot_df.txt",sep=""),stringsAsFactors=F,header=T)
human_alt_splice<-read.table(paste(Dir,"human_distance_transcribed_alt_introns_med_plot_df.txt",sep=""),stringsAsFactors=F,header=T)

human_const_splice$Group <- "const"
human_alt_splice$Group <- "alt"

test_df <- rbind(human_const_splice,human_alt_splice) %>% select(fields)

res.aov2 <- aov(percentSpliced ~ Group + range, data = test_df)
summary(res.aov2)
#0.0297

#####################

#### Human minimap2 vs. GMAP (Fig. S5C)
human_minimap_splice<-read.table(paste(Dir,"human_distance_transcribed_minimap2_plot_df",sep=""),stringsAsFactors=F,header=T)
human_gmap_splice<-read.table(paste(Dir,"human_distance_transcribed_gmap_plot_df",sep=""),stringsAsFactors=F,header=T)

human_minimap_splice$Group <- "minimap"
human_gmap_splice$Group <- "gmap"

test_df <- rbind(human_minimap_splice,human_gmap_splice) %>% select(fields)

res.aov2 <- aov(percentSpliced ~ Group + range, data = test_df)
summary(res.aov2)
#0.572

#### Drosophila minimap2 vs. GMAP (Fig. S5C)
drosophila_minimap_splice<-read.table(paste(Dir,"drosophila_distance_transcribed_minimap2_plot_df",sep=""),stringsAsFactors=F,header=T)
drosophila_gmap_splice<-read.table(paste(Dir,"drosophila_distance_transcribed_gmap_plot_df",sep=""),stringsAsFactors=F,header=T)

drosophila_minimap_splice$Group <- "minimap"
drosophila_gmap_splice$Group <- "gmap"

test_df <- rbind(drosophila_minimap_splice,drosophila_gmap_splice) %>% select(fields)

res.aov2 <- aov(percentSpliced ~ Group + range, data = test_df)
summary(res.aov2)
#0.715


##################

### Human stringency (Fig. S5D)
human_no_splice<-read.table(paste(Dir,"human_distance_transcribed_noStringency_plot_df.txt",sep=""),stringsAsFactors=F,header=T)
human_low_splice<-read.table(paste(Dir,"human_distance_transcribed_lowStringency_plot_df.txt",sep=""),stringsAsFactors=F,header=T)
human_med_splice<-read.table(paste(Dir,"human_distance_transcribed_mediumStringency_plot_df.txt",sep=""),stringsAsFactors=F,header=T)
human_high_splice<-read.table(paste(Dir,"human_distance_transcribed_highStringency_plot_df.txt",sep=""),stringsAsFactors=F,header=T)

human_no_splice$Group <- "no"
human_low_splice$Group <- "low"
human_med_splice$Group <- "med"
human_high_splice$Group <- "high"

test_df <- rbind(human_no_splice,human_low_splice,human_med_splice,human_high_splice) %>% select(fields)

res.aov2 <- aov(percentSpliced ~ Group + range, data = test_df)
summary(res.aov2)
#0.00522

### Drosophila stringency (Fig. S5D)
drosophila_no_splice<-read.table(paste(Dir,"drosophila_distance_transcribed_noStringency_plot_df.txt",sep=""),stringsAsFactors=F,header=T)
drosophila_low_splice<-read.table(paste(Dir,"drosophila_distance_transcribed_lowStringency_plot_df.txt",sep=""),stringsAsFactors=F,header=T)
drosophila_med_splice<-read.table(paste(Dir,"drosophila_distance_transcribed_mediumStringency_plot_df.txt",sep=""),stringsAsFactors=F,header=T)
drosophila_high_splice<-read.table(paste(Dir,"drosophila_distance_transcribed_highStringency_plot_df.txt",sep=""),stringsAsFactors=F,header=T)

drosophila_no_splice$Group <- "no"
drosophila_low_splice$Group <- "low"
drosophila_med_splice$Group <- "med"
drosophila_high_splice$Group <- "high"

test_df <- rbind(drosophila_no_splice,drosophila_low_splice,drosophila_med_splice,drosophila_high_splice) %>% select(fields)

res.aov2 <- aov(percentSpliced ~ Group + range, data = test_df)
summary(res.aov2)
#<2e-16


#########################

### Human high vs. medium RPKM (Fig. S5H)
human_low_splice<-read.table(paste(Dir,"human_distance_transcribed_lowRPKM_plot_df.txt",sep=""),stringsAsFactors=F,header=T)
human_medium_splice<-read.table(paste(Dir,"human_distance_transcribed_medRPKM_plot_df.txt",sep=""),stringsAsFactors=F,header=T)
human_high_splice<-read.table(paste(Dir,"human_distance_transcribed_highRPKM_plot_df.txt",sep=""),stringsAsFactors=F,header=T)

human_low_splice$Group <- "low"
human_med_splice$Group <- "med"
human_high_splice$Group <- "high"

test_df <- rbind(human_low_splice,human_med_splice,human_high_splice) %>% select(fields)

res.aov2 <- aov(percentSpliced ~ Group + range, data = test_df)
summary(res.aov2)
#4.38e-05

### Drosophila high vs. medium RPKM (Fig. S5H)
drosophila_low_splice<-read.table(paste(Dir,"drosophila_distance_transcribed_lowRPKM_plot_df.txt",sep=""),stringsAsFactors=F,header=T)
drosophila_medium_splice<-read.table(paste(Dir,"drosophila_distance_transcribed_medRPKM_plot_df.txt",sep=""),stringsAsFactors=F,header=T)
drosophila_high_splice<-read.table(paste(Dir,"drosophila_distance_transcribed_highRPKM_plot_df.txt",sep=""),stringsAsFactors=F,header=T)

drosophila_low_splice$Group <- "low"
drosophila_med_splice$Group <- "med"
drosophila_high_splice$Group <- "high"

test_df <- rbind(drosophila_low_splice,drosophila_med_splice,drosophila_high_splice) %>% select(fields)

res.aov2 <- aov(percentSpliced ~ Group + range, data = test_df)
summary(res.aov2)
#5.56e-15


##############################

### Drosophila vs. Human without first and last introns (Fig. S5I)
human_splice<-read.table(paste(Dir,"human_distance_transcribed_med_noFirstLast_plot_df.txt",sep=""),stringsAsFactors=F,header=T)
droso_splice<-read.table(paste(Dir,"drosophila_distance_transcribed_med_noFirstLast_plot_df.txt",sep=""),stringsAsFactors=F,header=T)

human_splice$Group <- "human"
droso_splice$Group <- "drosophila"

test_df <- rbind(human_splice,droso_splice) %>% select(fields)

res.aov2 <- aov(percentSpliced ~ Group + range, data = test_df)
summary(res.aov2)
#<2e-16

################################

### Human high vs. low splicing yield (Fig. S7G)
human_low_splice<-read.table(paste(Dir,"human_distance_transcribed_lowYield_plot_df.txt",sep=""),stringsAsFactors=F,header=T)
human_med_splice<-read.table(paste(Dir,"human_distance_transcribed_medYield_plot_df.txt",sep=""),stringsAsFactors=F,header=T)
human_high_splice<-read.table(paste(Dir,"human_distance_transcribed_highYield_plot_df.txt",sep=""),stringsAsFactors=F,header=T)

human_low_splice$Group <- "low"
human_med_splice$Group <- "med"
human_high_splice$Group <- "high"

test_df <- rbind(human_low_splice,human_med_splice,human_high_splice) %>% select(fields)

res.aov2 <- aov(percentSpliced ~ Group + range, data = test_df)
summary(res.aov2)
#6.35e-09

#################################

# Human DMSO vs. PlaB (Fig. 3)
human_dmso_splice<-read.table(paste(Dir,"Human_DMSO_distance_transcribed_plot_df.txt",sep=""),stringsAsFactors=F,header=T)
human_plaB_splice<-read.table(paste(Dir,"Human_PlaB_distance_transcribed_plot_df.txt",sep=""),stringsAsFactors=F,header=T)

human_dmso_splice$Group <- "dmso"
human_plaB_splice$Group <- "plaB"

test_df <- rbind(human_dmso_splice,human_plaB_splice) %>% select(fields)

res.aov2 <- aov(percentSpliced ~ Group + range, data = test_df)
summary(res.aov2)
#5.77e-06

# Drosophila DMSO vs. PlaB (Fig. 3)
droso_dmso_splice<-read.table(paste(Dir,"Drosophila_DMSO_distance_transcribed_plot_df.txt",sep=""),stringsAsFactors=F,header=T)
droso_plaB_splice<-read.table(paste(Dir,"Drosophila_PlaB_distance_transcribed_plot_df.txt",sep=""),stringsAsFactors=F,header=T)

droso_dmso_splice$Group <- "dmso"
droso_plaB_splice$Group <- "plaB"

test_df <- rbind(droso_dmso_splice,droso_plaB_splice) %>% select(fields)

res.aov2 <- aov(percentSpliced ~ Group + range, data = test_df)
summary(res.aov2)
#5.21e-06



