options(warn = -1)
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(ape))
suppressMessages(library(reshape2))
suppressMessages(library(tidyverse))
suppressMessages(library(parallel))
suppressMessages(library(RColorBrewer))
suppressMessages(library(dplyr))
args <- commandArgs(trailingOnly=TRUE)

wd <- args[1]
treelistfile <- args[2]
threads <- strtoi(args[3])
ndistances <- args[4]
randomseed <- args[5]
select <- args[6]
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
cat("\nWill calculate tip 2 tip distances as a measure of tree similarity.\n")
treelist <- read.csv(treelistfile, sep="\t", header=F)
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
if (select != "none") {
	cat("Will only select trees based on specifications in --select:", select, "\n")
	sel <- str_split(select, ",")[[1]]
	treelist <- filter(treelist, grepl(paste(sel, collapse="|"), path))
}

trees <- list()
i <- 1

for (treepath in treelist$path) {  
  trees[[i]]<-read.tree(treepath)
  i <- i + 1
}

cat(paste0("Found ", length(trees), " trees.\n"))
newtreenames <- c() 
for (names in strsplit(treelist$path,"/")){
  if (grepl("unpartitioned", names[3], fixed=TRUE) == TRUE) { # check if some trees come from unpartitioned analyses
	  treebuilder <- gsub("-", "_", names[3])
  } else { treebuilder <- names[3] } 
  newtreenames <-c(newtreenames, paste0(treebuilder, "-", strsplit(names[5], ".", fixed=T)[[1]][1], "-", strsplit(names[4], "-")[[1]][3]))
}
treelist$path <- newtreenames
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
  #cat(paste0(Sys.time(), " Calculating pairwise distances for tree ", i, "...\n"))
  
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
  #cat(paste0(Sys.time(), " Distances extracted for tree ", i, ": ", length(distances), "\n"))
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

# make sure everyting is saved as numeric values
for (i in 1:ncol(complete_df)) {
   complete_df[,i] <- as.numeric(complete_df[,i])
 }

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
plotdf$bootstrap_cutoff <- bootstrapcutoffs
plotdf$aligner_trimmer <- alitrim

cols <- colorRampPalette(brewer.pal(length(alitrim), "Set1"))(length(alitrim))
names(cols) <- sort(alitrim)

p <- ggplot(data = plotdf, aes(x=PC1, y=PC2, label=treenames, color=aligner_trimmer, shape=treebuilder, size=bootstrap_cutoff)) + geom_point() +xlab(paste0("PC1 (", PC1Variance, "%)")) +ylab(paste0("PC2 (", PC2Variance,"%)")) +ggtitle(paste0("Similarity of trees based on tip to tip distance matrices.\nBased on ", ndistances, " distances between tips."))+ scale_shape_manual(values = c(15, 16, 17, 18))+ scale_color_manual(values=cols)

if (select != "none") {
	selection <- paste0("-select-", select)
} else { selection <- "" }

pdf(file=paste0("tip2tip-PCA-",ndistances,selection,".pdf"))
print(p)
garbage <- dev.off()
write.csv(file=paste0("pairwise-tip2tip-distance-matrix-",length(trees),"-trees.csv"), complete_df, sep=",", quote=FALSE)
cat("PCA plotting is done.\n")
cat('Output: tip2tip-PCA-*.pdf - PCA of tip to tip distances in trees. * will be the number of distances specified with --ndistances and info from --select (if specified)\n')
cat('Output: pairwise-tip2tip-distance-matrix-*-trees.csv - tip2tip distance matrix used to create the PCA. * will be the number of trees in comparison.\n')
