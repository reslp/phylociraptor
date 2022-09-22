options(warn = -1)
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(ape))
suppressMessages(library(reshape2))
suppressMessages(library(tidyverse))
suppressMessages(library(parallel))

args <- commandArgs(trailingOnly=TRUE)

wd <- args[1]
matrix_file <- args[2]
treelistfile <- args[3]

verbose <- FALSE

setwd(wd)
cat(paste("Working directory: ", getwd(), "\n"))

cat("Will plot quartet similarity heatmap based on the provided simmilarity matrix file\n")
data <- read.csv(matrix_file)
data <- melt(data)
colnames(data) <- c("First", "Second", "Similarity")
data$First <- factor(data$First, levels=levels(data$Second))
pwidth <- 30
pheight <- 14
if (treelistfile == "none") {
  cat("Plotting without detailed tree name information.\n")
  pdf(file=paste0("quartet-similarity-heatmap-", length(unique(data$First)), "-trees.pdf"), width=pwidth, height=pheight)
  p <- ggplot(data, aes(First, Second, fill=Similarity)) + geom_tile() + scale_fill_gradient(low = "red", high = "white") + ggtitle("% similarity of pairs of trees") + geom_text(aes(label = format(round(Similarity, 4), nsmall=2)))
  print(p)
  dev.off()
} else {
  cat("Using detailed tree name information for plotting.\n")
  treelist <- read.csv(treelistfile, sep="\t", header=F)
  newtreenames <- c() 
  for (names in strsplit(treelist$V2,"/")){
    newtreenames <-c(newtreenames, paste0(names[3], "-", names[4], "-", strsplit(names[2], "-")[[1]][2]))
  }
  treelist$V2 <- newtreenames
  data$First <- treelist$V2[match(data$First, treelist$V1)]
  data$Second <- treelist$V2[match(data$Second, treelist$V1)]

  p <- ggplot(data, aes(First, Second, fill=Similarity)) + geom_tile() + scale_fill_gradient(low = "red", high = "white") + ggtitle("% similarity of pairs of trees") + geom_text(aes(label = format(round(Similarity, 4), nsmall=2)))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  pdf(file=paste0("quartet-similarity-heatmap-",length(unique(data$First)),"-trees.pdf"), width=pwidth, height=pheight)
  print(p)
  dev.off()
  
}

cat("Similarity plotting is done.\n")
cat('Output: quartet-similarity-heatmap-*trees.pdf - Heatmap of % similarity based on quartets of tips. * is the number of trees in comparison.\n')



