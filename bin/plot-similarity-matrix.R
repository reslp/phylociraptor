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
threads <- strtoi(args[4])
ndistances <- args[5]
randomseed <- args[6]

verbose <- FALSE
#if (quiet == "true") {
#  verbose <- FALSE
#}
if (randomseed != "random"){
  cat(paste0("Random seed: ", randomseed, "\n"))
  set.seed(strtoi(randomseed))
}

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

  p <- ggplot(data, aes(First, Second, fill=Similarity)) + geom_tile() + scale_fill_gradient(low = "red", high = "white") + ggtitle("% similarity of pairs of trees") + geom_text(aes(label = format(round(Similarity, 4), nsmall=2)))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  pdf(file=paste0("quartet-similarity-heatmap-",length(unique(data$First)),"-trees.pdf"), width=20, height=10)
  print(p)
  dev.off()
  
}

cat("\nA treelist file was specified. Will therefore also calculate tip 2 tip distances as a measure of tree similarity.\n")
if (treelistfile != "none") {
  ### the melt_dist function is from: https://github.com/andersgs/harrietr/blob/master/R/melt_distance.R
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
  
  cat("Reading trees from treelist file.\n")
  trees <- list()
  i <- 1
  for (treepath in treelist$path) {
    
    trees[[i]]<-read.tree(treepath)
    i <- i + 1
  }
  cat(paste0("Found ", length(trees), " trees.\n"))
  
  all_tips <- c()
  
  cat("Extracting all tips from all trees...\n")
  for (i in 1:length(trees)) {
    all_tips <- c(all_tips, trees[[i]]$tip.label)
  }
  all_tips <- unique(all_tips)
  cat("Creating pairwise combinations of tips...\n")
  all_combinations <- t(combn(all_tips, 2))
  
  for (k in 1:nrow(all_combinations)){ #sort species pairs alphabetically
    sp1 <- all_combinations[k, 1]
    sp2 <- all_combinations[k, 2]
    all_combinations[k, 1] <- paste0(sort(c(sp1, sp2))[1], "-", sort(c(sp1, sp2))[2])
    #all_combinations[k, 2] <- sort(c(sp1, sp2))[2]
  }
  if (ndistances != "all") {
    cat("Will choose a random sample of tip-combinations.\n")
    all_combinations <- all_combinations[sample(nrow(all_combinations), strtoi(ndistances)),]
  } 
  cat(paste0("Number of tip-pairs used to calculate distances: ", length(all_combinations[,1]), "\n"))

  

  cat("Now calculating tip2tip distances, please be patient this could take a few minutes per tree...\n")
  cat(paste0("Will be using ", threads, " threads for this step.\n"))

  dist_single_tree <- function(i) {
    cat(paste0(Sys.time(), " Calculating pairwise distances for tree ", i, "...\n"))
    
    trees[[i]]$edge.length[trees[[i]]$edge.length == "NaN"] <- 0
    dist <- cophenetic.phylo(trees[[i]])
    
    melted_dist <- melt(dist)
    melted_dist <- as.matrix(melted_dist) #convert to matrix; this is faster than df
    combined <- c()
    combined[1:nrow(melted_dist)] <- 1:nrow(melted_dist)
    for (k in 1:nrow(melted_dist)){ #sort species pairs alphabetically
      sp1 <- melted_dist[k, 1] # matrix is faster than df
      sp2 <- melted_dist[k, 2]
      combined[k] <- paste(sort(c(sp1, sp2)), collapse="-")
    }
    melted_dist <- as.data.frame(melted_dist)
    melted_dist$combined <- combined
    #remove duplicated rows for faster processing
    melted_dist$Var1 <- NULL
    melted_dist$Var2 <- NULL
    melted_dist <- unique(melted_dist)
    
    cat(paste0(Sys.time()," Checking tip2tip distances from tree ", i, " against all combinations from all trees...\n"))
    names <- c()
    distances <- c()
    names[1:nrow(all_combinations)] <- 1:nrow(all_combinations)
    distances[1:nrow(all_combinations)] <- 1:nrow(all_combinations)
    
    for (j in 1:nrow(all_combinations)) {
      comb <- all_combinations[j, 1]
      if (verbose == TRUE) {
        cat(paste0("Tree ", i, " Combination: ", j, "/", nrow(all_combinations)))
      }
      distance <- melted_dist$value[melted_dist$combined == comb]
      if (length(distance) == 0) {
        distance <- NA
      }
      distances[j] <- distance
      names[j] <- comb
    }
    cat(paste0(Sys.time(), " Distances extracted for tree ", i, ": ", length(distances), "\n"))
    df <- data.frame(as.numeric(distances))
    rownames(df) <- names
    colnames(df) <- paste0("tree", i)
    return(df)
  }
  Sys.time()
  dists <- list()
  dists <- mclapply(1:length(trees), dist_single_tree, mc.cores=threads)
  Sys.time()

  complete_df <- do.call("cbind", dists)
  complete_df <- t(complete_df)
  
  #complete_df <- t(complete_df) 
#  for (i in 1:12) {
#    complete_df[,i] <- as.numeric(complete_df[,i])
#  }
  
  cat("Calculating PCA of tip2tip distances now.\n")
  pca <-prcomp(complete_df,na.action = na.omit)
  sumpca <- summary(pca)
  PC1Variance <- sumpca$importance[2,1]*100
  PC2Variance <- sumpca$importance[2,2]*100
  
  plotdf <- as.data.frame(pca$x)
  tnames <- rownames(plotdf)
  plotdf["treenames"] <- tnames
  
  treebuilders <- c()
  for (x in strsplit(newtreenames, "-")){
    treebuilders <- c(treebuilders, x[1])
    
  }
  alitrim <- c()
  for (x in strsplit(newtreenames, "-")){
    alitrim <- c(alitrim, paste0(x[2],"-",x[3]))
    
  }
  bootstrapcutoffs <- c()
  for (x in strsplit(newtreenames, "-")){
    bootstrapcutoffs <- c(bootstrapcutoffs, x[4])
    
  }
  plotdf$treebuilder <- treebuilders
  plotdf$bscutoff <- bootstrapcutoffs
  plotdf$alitrim <- alitrim
  
  p <- ggplot(data = plotdf, aes(x=PC1, y=PC2, label=treenames, color=alitrim, shape=treebuilder, size=bscutoff)) + geom_point() +xlab(paste0("PC1 (", PC1Variance, "%)")) +ylab(paste0("PC2 (", PC2Variance,"%)")) +ggtitle(paste0("Similarity of trees based on tip to tip distance matrices.\nBased on ", ndistances, " distances between tips."))+ scale_shape_manual(values = c(15, 16, 17, 18))+ scale_color_brewer(palette="Set1")
  
  pdf(file=paste0("tip2tip-PCA-",ndistances,".pdf"))
  print(p)
  dev.off()
  write.csv(file=paste0("pairwise-distance-matrix-",length(trees),"-trees.csv"), complete_df, sep=",", quote=FALSE)
 
} # if statement end

cat("Similarity plotting is done.\n")
cat('Output: quartet-similarity-heatmap-*trees.pdf - Heatmap of % similarity based on quartets of tips. * is the number of trees in comparison.\n')
if (treelistfile != "none") {
  cat('Output: tip2tip-PCA-*.pdf - PCA of tip to tip distances in trees. * will be the number of distances specified with --ndistances\n')
  cat('Output: pairwise-distance-matrix-*-trees.csv - tip2tip distance matrix used to create the PCA. * will be the number of trees in comparison.\n')
}



