options(warn = -1)
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(ape))
suppressMessages(library(reshape2))
suppressMessages(library(tidyverse))
suppressMessages(library(parallel))
suppressMessages(library(pheatmap))
args <- commandArgs(trailingOnly=TRUE)

wd <- args[1]
matrix_file <- args[2]
treelistfile <- args[3]
decimals <- as.numeric(args[4])

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


if (nrow(data) <= 10) { # more sensitive plot size, this is still a bit hit or miss...
	pwidth <- 1 * nrow(data)
	pheight <- 0.75 * nrow(data)
} else {
	pwidth <- 0.38 * nrow(data)
	pheight <- 0.25 * nrow(data)
}
roundb <- 2 # rounding bounds
if (treelistfile == "none") {
  cat("Plotting without detailed tree name information.\n")
  pdf(file=paste0("quartet-similarity-heatmap-", nrow(data), "-trees.pdf"), width=pwidth, height=pheight)
  #colnames(data) <- rownames(data)
  p <- pheatmap(data, display_numbers=TRUE, number_format = formatstring, treeheight_col=0, treeheight_row=0, main=paste0("Similarity (percentage of identical quartets of tips) of pairs of ", nrow(data), " trees"), angle_col=0)
  #p <- ggplot(data, aes(First, Second, fill=Similarity)) + geom_tile() + scale_fill_gradient(low = "red", high = "white") + ggtitle("% similarity of pairs of trees") + geom_text(aes(label = format(round(Similarity, roundb), nsmall=2)))
  print(p)
  garbage <- dev.off()
} else {
  cat("Using detailed tree names and reduce to tree list while plotting...\n")
  treelist <- read.csv(treelistfile, sep="\t", header=F)
  newtreenames <- c() 
  for (names in strsplit(treelist$V2,"/")){
    newtreenames <-c(newtreenames, paste0(names[3], "-", names[5], "-", strsplit(names[4], "-")[[1]][3]))
  }
  #remove rows and columns of tree which are not in the treelist file
  data <- data[treelist$V1,]
  data <- data[,treelist$V1]
  treelist$V2 <- newtreenames
  colnames(data) <- treelist$V2[match(colnames(data), treelist$V1)]
  rownames(data) <- treelist$V2[match(rownames(data), treelist$V1)]

  #p <- ggplot(data, aes(First, Second, fill=Similarity)) + geom_tile() + scale_fill_gradient(low = "red", high = "white") + ggtitle("% similarity of pairs of trees") + geom_text(aes(label = format(round(Similarity, roundb), nsmall=2)))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  p <- pheatmap(data, display_numbers=TRUE, number_format = formatstring, treeheight_col=0, treeheight_row=0, main=paste0("Similarity (percentage of identical quartets of tips) of pairs of ", nrow(data), " trees"))
  pdf(file=paste0("quartet-similarity-heatmap-",nrow(data),"-trees.pdf"), width=pwidth, height=pheight)
  print(p)
  garbage <- dev.off()
  
}

cat("Similarity plotting is done.\n")
cat(paste0("Output: quartet-similarity-heatmap-", nrow(data), "-trees.pdf - Heatmap of % similarity based on quartets of tips.\n"))



