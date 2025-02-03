options(warn = -1)
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(ape))
suppressMessages(library(reshape2))
suppressMessages(library(tidyverse))
suppressMessages(library(parallel))
suppressMessages(library(pheatmap))
suppressMessages(library(dplyr))
args <- commandArgs(trailingOnly=TRUE)

wd <- args[1]
matrix_file <- args[2]
treelistfile <- args[3]
decimals <- as.numeric(args[4])
selection <- args[5]
verbose <- FALSE


setwd(wd)
cat(paste("Working directory: ", getwd(), "\n"))

cat("Will plot quartet similarity heatmap based on the provided simmilarity matrix file\n")
data <- read.csv(matrix_file)
data$X <- NULL
#data <- data[1:10,1:10]
rownames(data) <- colnames(data)
data <- as.matrix(t(data))
mode(data) <- "numeric"

# for better formatting fo decimal places
cat(paste0("Will truncate values to ", decimals, " decimal places.\n"))
formatstring <- paste0("%.",decimals, "f") 
div <- 10**decimals
data <- trunc(data*div)/div


#data <- melt(data)
#colnames(data) <- c("First", "Second", "Similarity")
#data$First <- factor(data$First, levels=levels(data$Second))


roundb <- 2 # rounding bounds

#reduce_data <- function(data) {
#	selection <- strsplit(selection, ",")
#	data2 <- subset(data, grepl(selection))
#	return(data2)
#}

if (treelistfile == "none") {
  cat("Plotting without detailed tree name information.\n")
  if (selection != "none") {
	  allnames <- colnames(data)
	  sel <- as.vector(strsplit(selection, ",")[[1]])
          cat( "Will select only trees: ", paste(sel, " "), "\n")
	  data <- data[,sel]
	  data <- data[sel,]
  }
  if (nrow(data) == 0 || ncol(data) == 0) {
	  cat("Nothing left to plot. Check your filters / treelist.\n")
	  quit()
  }
  if (nrow(data) <= 10) { # more sensitive plot size, this is still a bit hit or miss...
  	pwidth <- 1 * nrow(data)
  	pheight <- 0.75 * nrow(data)
  } else {
  	pwidth <- 0.38 * nrow(data)
  	pheight <- 0.25 * nrow(data)
  }

  pdf(file=paste0("quartet-similarity-heatmap-", nrow(data), "-trees.pdf"), width=pwidth, height=pheight)
  #colnames(data) <- rownames(data)
  p <- pheatmap(data, display_numbers=TRUE, number_format = formatstring, treeheight_col=0, treeheight_row=0, main=paste0("Similarity (percentage of identical quartets of tips) of pairs of ", nrow(data), " trees"), angle_col=0)
  #p <- ggplot(data, aes(First, Second, fill=Similarity)) + geom_tile() + scale_fill_gradient(low = "red", high = "white") + ggtitle("% similarity of pairs of trees") + geom_text(aes(label = format(round(Similarity, roundb), nsmall=2)))
  garbage <- dev.off()
} else {
  cat("Using detailed tree names and reduce to tree list while plotting...\n")
  treelist <- read.csv(treelistfile, sep="\t", header=F)
  print(treelist$V2)
  #remove according to selection
  if (selection != "none") {
	sel <- as.vector(strsplit(selection, ",")[[1]])
	cat( "Will select only trees containing: ", paste(sel, " "), "\n")
	treelist2 <- filter(treelist, grepl(paste(sel, collapse="|"), V2))
  }
  else { treelist2 <- treelist }
  data <- data[treelist2$V1,]
  data <- data[,treelist2$V1]
  # get new tree names
  newtreenames <- c() 
  for (names in strsplit(treelist2$V2,"/")){
    newtreenames <-c(newtreenames, paste0(names[3], "-", names[5], "-", strsplit(names[4], "-")[[1]][3]))
  }
  # rename trees in dataframe
  treelist2$V2 <- newtreenames
  colnames(data) <- treelist2$V2[match(colnames(data), treelist2$V1)]
  rownames(data) <- treelist2$V2[match(rownames(data), treelist2$V1)]

  if (nrow(data) <= 10) { # more sensitive plot size, recalculate
  	pwidth <- 1 * nrow(data)
  	pheight <- 0.75 * nrow(data)
  } else {
  	pwidth <- 0.38 * nrow(data)
  	pheight <- 0.25 * nrow(data)
  }
  if (nrow(data) == 0 || ncol(data) == 0) {
	  cat("Nothing left to plot. Check your filters / treelist.\n")
	  quit()
  }
  #p <- ggplot(data, aes(First, Second, fill=Similarity)) + geom_tile() + scale_fill_gradient(low = "red", high = "white") + ggtitle("% similarity of pairs of trees") + geom_text(aes(label = format(round(Similarity, roundb), nsmall=2)))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  p <- pheatmap(data, display_numbers=TRUE, number_format = formatstring, treeheight_col=0, treeheight_row=0, main=paste0("Similarity (percentage of identical quartets of tips) of pairs of ", nrow(data), " trees"))
  pdf(file=paste0("quartet-similarity-heatmap-",nrow(data),"-trees.pdf"), width=pwidth, height=pheight)
  print(p)
  garbage <- dev.off()
  
}

cat("Similarity plotting is done.\n")
cat(paste0("Output: quartet-similarity-heatmap-", nrow(data), "-trees.pdf - Heatmap of % similarity based on quartets of tips.\n"))



