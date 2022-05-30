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
if (treelistfile == "none") {
  cat("Plotting without detailed tree name information.\n")
  pdf(file="similarity-matrix.pdf", width=20, height=10)
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
  pdf(file="similarity-matrix.pdf", width=20, height=10)
  p <- ggplot(data, aes(First, Second, fill=Similarity)) + geom_tile() + scale_fill_gradient(low = "red", high = "white") + ggtitle("% similarity of pairs of trees") + geom_text(aes(label = format(round(Similarity, 4), nsmall=2)))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  print(p)
  dev.off()
  
}
cat("Plot file: similarity-matrix.pdf created successfully.\n")
quit()




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
system.time(
  for (i in 1:1) {
    cat(paste0("Calculating pairwise distances for tree ", i, "...\n"))
    dist <- cophenetic.phylo(trees[[i]])
    melted_dist <- melt_dist(dist)
    print(nrow(melted_dist))
    visited <- c()
    remaining_tips <- all_tips
    names <- c()
    distances <- c()
    all_combinations <- t(combn(all_tips, 2))
    #nrow(all_combinations)
    for (j in 1:100) {
      sp1 <- all_combinations[j, 1]
      sp2 <- all_combinations[j, 2]
      if (verbose == TRUE) {
        cat(paste0("Tree ",i, ": ", sp1, " - ", sp2, "\n"))
      }
      if (sp1 %in% melted_dist$iso1 && sp2 %in% melted_dist$iso2) {
        distance <- subset(melted_dist, (iso1==sp1 & iso2==sp2))$dist
      } else if (sp2 %in% melted_dist$iso1 && sp1 %in% melted_dist$iso2) {
        distance <- subset(melted_dist, (iso1==sp2 & iso2==sp1))$dist
      } else {
        distance <- NA
      }
      distances <- c(distances, distance)
      names <- c(names, paste0(sp1, "-",sp2))
    }
    df <- data.frame(distances)
    rownames(df) <- names
    colnames(df) <- paste0("tree", i)
    dists[[i]] <- df
  }
)


sp2 %in% melted_dist$iso1
sp1 %in% melted_dist$iso1
df
    # for (sp1 in all_tips[1:10]) {
    #   visited <- c(visited, sp1)
    #   remaining_tips <- remaining_tips[!remaining_tips %in% visited]
    #   for (sp2 in remaining_tips[1:10]) {
    #     if (verbose == TRUE) {
    #       cat(paste0("Tree ",i, ": ", sp1, " - ", sp2, "\n"))
    #     }
    #     if (sp1 == sp2) {next}
    #     distance <- subset(melted_dist, (iso1==sp1 & iso2==sp2))$dist
    #     if (length(distance) == 0) {
    #       distance <- subset(melted_dist, (iso1==sp2 & iso2==sp1))$dist
    #     } 
    #     if (length(distance) == 0) {distance <- NA}
    #     distances <- c(distances, distance)
    #     names <- c(names, paste0(sp1, "-",sp2))
    #   }
    # }
    # df <- data.frame(distances)
    # rownames(df) <- names
    # colnames(df) <- paste0("tree", i)
    # dists[[i]] <- df
  }
  dists
  complete_df <- do.call("cbind", dists)
  complete_df <- t(complete_df) 
  complete_df
  write.csv(file="pairwise-distance-matrix.csv", complete_df, sep=",")
}

