#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

### SUB-SAMPLING ANALYSIS
# ---------------------
# this script will produce "guide" files, each of which will specify a different subsampling regime. There are 3 such regimes:
# 1) Random downsampling of all populations
# 2) Uneven downsampling of *some* but not all populations
# 3) Even downsampling of *all* populations
# Run this script as part of WGS_subsampling_wrapper.sh, which will also include code for running SMC++ and stairway plot on all subsets.

# importing population data
setwd(args[1])
inbu_pops <- read.delim("INBU_pop_map.txt", header = F, col.names = c("sample", "population")) ; inbu_pops$sample <- as.character(inbu_pops$sample) ; inbu_pops$population <- as.character(inbu_pops$population)
lazb_pops <- read.delim("LAZB_pop_map.txt", header = F, col.names = c("sample", "population")) ; lazb_pops$sample <- as.character(lazb_pops$sample) ; lazb_pops$population <- as.character(lazb_pops$population)
setwd(args[2])

pops <- rbind(inbu_pops, lazb_pops)
pops$sample <- as.character(pops$sample) ; pops$population <- as.character(pops$population)

# 1) Random subsampling
set.seed(14)
rows_to_keep <- sample(row.names(pops), size = round(nrow(pops) * 0.75))
sample1 <- pops[rows_to_keep,]
colnames(sample1) <- NULL
write.table(sample1, "./random_sample.txt", row.names = F, sep = "\t", quote = F)

# 2) Uneven downsampling
set.seed(14)
uneven_down_lazb <- sample(unique(lazb_pops$population), size = 1)
lazb_unchanged <- lazb_pops[lazb_pops$population!=uneven_down_lazb[1],]
sub <- lazb_pops[lazb_pops$population==uneven_down_lazb[1],]
rows_to_keep <- sample(row.names(sub), size = round(nrow(sub) * 0.75))
lazb_down <- sub[rows_to_keep,]

set.seed(14)
uneven_down_inbu <- sample(unique(inbu_pops$population), size = 2) # use ls(pat = "inbu_down_") to check what the states are this run, and then change line 46 accordingly
inbu_sub1 <- inbu_pops[inbu_pops$population!=uneven_down_inbu[1],]
inbu_unchanged <- inbu_sub1[inbu_sub1$population!=uneven_down_inbu[2],] ; remove(inbu_sub1)
sub1 <- inbu_pops[inbu_pops$population==uneven_down_inbu[1],]
sub2 <- inbu_pops[inbu_pops$population==uneven_down_inbu[2],]
inbu_down <- rbind(sub1, sub2) ; remove(sub2) ; remove(sub1)
for (i in uneven_down_inbu){
	sub <- inbu_down[inbu_down$population==i,]
	rows_to_keep <- sample(row.names(sub), size = round(nrow(sub) * 0.75))
	samp <- sub[rows_to_keep,]
	nam <- paste("inbu_down", i, sep = "_")
	assign(nam, samp)
}
if ("kansas" %in% pops$population){
  inbu_down <- rbind(inbu_down_kansas, inbu_down_new_york)
} else {
  inbu_down <- rbind(inbu_down_louisiana, inbu_down_new_york)
}
sample2 <- rbind(lazb_unchanged, lazb_down, inbu_unchanged, inbu_down)
colnames(sample2) <- NULL
write.table(sample2, "./uneven_sample.txt", row.names = F, sep = "\t", quote = F)

# 3) Even downsampling
for (i in unique(pops$population)){
	sub <- pops[pops$population==i,]
	rows_to_keep <- sample(row.names(sub), size = round(nrow(sub) * 0.75))
	samp <- sub[rows_to_keep,]
	nam <- paste("sample", i, sep = "_")
	assign(nam, samp)
}

if ("washington" %in% pops$population){
  sample3 <- rbind(sample_new_york, sample_louisiana, sample_kansas, sample_washington, sample_wyoming)
} else {
  sample3 <- rbind(sample_new_york, sample_louisiana, sample_nevada)
}

colnames(sample3) <- NULL
write.table(sample3, "./even_sample.txt", row.names = F, sep = "\t", quote = F)

# Separating into the two species
sample1_lazb <- sample1[!grepl("INBU", sample1[,1]),]
sample1_inbu <- sample1[!grepl("LAZB", sample1[,1]),]
write.table(sample1_lazb, "./random_sample_lazb.txt", row.names = F, sep = "\t", quote = F)
write.table(sample1_inbu, "./random_sample_inbu.txt", row.names = F, sep = "\t", quote = F)

sample2_lazb <- sample2[!grepl("INBU", sample2[,1]),]
sample2_inbu <- sample2[!grepl("LAZB", sample2[,1]),]
write.table(sample2_lazb, "./uneven_sample_lazb.txt", row.names = F, sep = "\t", quote = F)
write.table(sample2_inbu, "./uneven_sample_inbu.txt", row.names = F, sep = "\t", quote = F)

sample3_lazb <- sample3[!grepl("INBU", sample3[,1]),]
sample3_inbu <- sample3[!grepl("LAZB", sample3[,1]),]
write.table(sample3_lazb, "./even_sample_lazb.txt", row.names = F, sep = "\t", quote = F)
write.table(sample3_inbu, "./even_sample_inbu.txt", row.names = F, sep = "\t", quote = F)
