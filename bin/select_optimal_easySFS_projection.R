#!/usr/bin/env Rscript

# import projection options saved from easySFS
options <- read.delim("projection_options.txt", header = F, sep = ";")
skip <- which(grepl("best and then rerun easySFS with the", options$V1))+4
options <- read.delim("projection_options.txt", header = F, , sep = ";",skip = skip)

# Reformat those options into a data frame
options_df <- data.frame("population" = NA,
						 "sample_size" = as.integer(NA),
						 "SNP_count" = as.integer(NA))
for (i in which(!grepl(",", options$V1))){
  
  option_vector <- options$V1[i+1]
  option_vector <- unlist(strsplit(option_vector, " "))
  option_vector <- unlist(strsplit(option_vector, "(", fixed = T))
  option_vector <- unlist(strsplit(option_vector, ")", fixed = T))
  option_vector <- unlist(strsplit(option_vector, ","))
  option_vector <- as.integer(option_vector[option_vector!=""])
  
  for (j in seq(1,length(option_vector), 2)){
	
	new_row <- c(options$V1[i], option_vector[j], option_vector[j+1])
	options_df <- rbind(options_df, new_row)
	
  }
  
}
options_df <- options_df[2:nrow(options_df),] ; rownames(options_df) <- NULL
options_df$sample_size <- as.integer(options_df$sample_size)
options_df$SNP_count <- as.integer(options_df$SNP_count)

# Identify the projections for each population that provide the most SNPs
projections <- ""
for (i in unique(options_df$population)){
  
  sub <- options_df[options_df$population==i,]
  max <- max(sub$SNP_count)
  proj <- sub[sub$SNP_count==max, "sample_size"]
  
  if (length(proj) > 1){
	
	proj <- proj[length(proj)]
  }
  
  if (projections == ""){
	projections <- proj
  } else {
	projections <- paste(projections, proj, sep = ",")
  }
  
}

write.table(projections, "projection.txt",
			quote = F, sep = "", col.names = F, row.names = F)  
