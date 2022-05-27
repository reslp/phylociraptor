library(ggplot2)
library(reshape2)
library(ape)
library(reshape2)
library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)

wd <- args[1]
matrix_file <- args[2]
treelistfile <- args[3]
verbose <- TRUE

setwd(wd)
cat(paste("Working directory: ", getwd(), "\n"))
data <- read.csv(matrix_file)

data <- melt(data)
colnames(data) <- c("First", "Second", "Similarity")
data$First <- factor(data$First, levels=levels(data$Second))

pdf(file="similarity-matrix.pdf", width=20, height=10)
ggplot(data, aes(First, Second, fill=Similarity)) + geom_tile() + ggtitle("% similarity of pairs of trees") + geom_text(aes(label = format(round(Similarity, 4), nsmall=2)))
dev.off()
cat("Plot file: similarity-matrix.pdf created successfully.\n")


if (treelistfile != "none") {
  ### taken from: https://github.com/andersgs/harrietr/blob/master/R/melt_distance.R
  melt_dist <- function(dist, order = NULL, dist_name = 'dist') {
    if(!is.null(order)){
      dist <- dist[order, order]
    } else {
      order <- row.names(dist)
    }
    diag(dist) <- NA
    dist[upper.tri(dist)] <- NA
    dist_df <- as.data.frame(dist)
    dist_df$iso1 <- row.names(dist)
    dist_df <- dist_df %>%
      tidyr::gather_(key = "iso2", value = lazyeval::interp("dist_name", dist_name = as.name(dist_name)), order, na.rm = T)
    return(dist_df)
  }
  ###
  
  
  
  treelist <- read.csv(treelistfile, header=F, check.names=FALSE, sep="\t")
  colnames(treelist) <- c("tree", "path")
  
  cat("Reading trees from treelist file:\n")
  trees <- list()
  i <- 1
  for (treepath in treelist$path) {
    
    trees[[i]]<-read.tree(treepath)
    i <- i + 1
  }
  
  all_tips <- c()
  
  cat("Extracting all tips from tree...")
  for (i in 1:length(trees)) {
    all_tips <- c(all_tips, trees[[i]]$tip.label)
  }
  all_tips <- unique(all_tips)
  cat("done\n")
  
  dists <- list()
  cat("Now calculating distances, please be patient...\n")
  for (i in 1:length(trees)) {
    cat(paste0("Calculating pairwise distances for tree ", i, "...\n"))
    dist <- cophenetic.phylo(trees[[i]])
    melted_dist <- melt_dist(dist)
    visited <- c()
    remaining_tips <- all_tips
    names <- c()
    distances <- c()
    for (sp1 in all_tips) {
      visited <- c(visited, sp1)
      remaining_tips <- remaining_tips[!remaining_tips %in% visited]
      for (sp2 in remaining_tips) {
        if (verbose == TRUE) {
          cat(paste0("Tree ",i, ": ", sp1, " - ", sp2, "\n"))
        }
        if (sp1 == sp2) {next}
        distance <- subset(melted_dist, (iso1==sp1 & iso2==sp2))$dist
        if (length(distance) == 0) {
          distance <- subset(melted_dist, (iso1==sp2 & iso2==sp1))$dist
        } 
        if (length(distance) == 0) {distance <- NA}
        distances <- c(distances, distance)
        names <- c(names, paste0(sp1, "-",sp2))
      }
    }
    df <- data.frame(distances)
    rownames(df) <- names
    colnames(df) <- paste0("tree", i)
    dists[[i]] <- df
  }
  complete_df <- do.call("cbind", dists)
  complete_df <- t(complete_df) 
  write.csv(file="pairwise-distance-matrix.csv", complete_df, sep=",")
}

