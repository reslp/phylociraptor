options(warn=-1)
suppressMessages(library(phytools))
suppressMessages(library(ape))
suppressMessages(library(ggtree))
suppressMessages(library(magrittr))
suppressMessages(library(tidytree))
suppressMessages(library(patchwork))
suppressMessages(library(cowplot))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(colorspace))

args <- commandArgs(trailingOnly=TRUE)

# variables which need to be passed to the script from the command line
wd <- args[1]
setwd(wd)

cat("This script plots conflicts between two trees based on quartets of tips ...\n")
cat("Plotting can take some time. Please be patient ...\n")
seed <- args[2]
treenames <- args[3]
outgroup <- args[4]
lineage_file <- args[5]
level <- args[6]
conflictfile <- args[7]
treelistfile <- args[8]

# reformat commandline argument:
outgroup <- strsplit(outgroup,",")[[1]]
treenames <- strsplit(treenames,",")[[1]]

#set seed if specified
if (seed != "random") {
  set.seed(seed)
}

#load and check quartelistt and treelist
if (conflictfile != "none") {
  cat(paste0("Loading quartet file:", conflictfile, "..."))
  conflicts <- read.csv(conflictfile, header=T, check.names=FALSE)
  rownames(conflicts) <- conflicts$quartet
  conflicts$quartet <- NULL
  conflicts <- t(conflicts)
  cat("DONE\n")
  cat(paste0("Loading treelist file:", treelistfile, "..."))
  treelist <- read.csv(treelistfile, header=F, check.names=FALSE, sep="\t")
  colnames(treelist) <- c("tree", "path")
  cat("DONE\n")
  #print(conflicts)
} else {
  cat("Conflicts file not found. Will stop.\n")
  quit(1)
}
  
reroot_my_tree <- function (tree, outgroup){
	tryCatch(
		{	
			rootnode <- getMRCA(tree, outgroup)
			tree <- root(tree, node=rootnode, resolve.root = TRUE)
			return(tree)
 		},
		error = function(e) {
			cat("There was an error while rerooting the tree. Will plot the tree as is. Please check manually\n")
			return(tree)
		})
}

# this function adds the missing tips (the ones for which no lineage info was downloaded) to the plotting dataframe so the tips get colors
add_missing_tips <- function(simpdf, tree) {
      for (tip in tree$tip.label) {
	if ( tip %in% simpdf$name) {
		next
	} else {
		df <- data.frame("name"=tip, "lineage"="missing")
		simpdf <- rbind(simpdf, df)
		}
      }
      return(simpdf)
}

# load lineage information file and fill missing values
if (lineage_file == "none") {
  cat("No lineage file specified. Will not plot lineage information.\n")
  level <- "none"
} else {
  lineages <- read.csv(lineage_file, header=T, sep="\t", na=c(""))
  lineages[is.na(lineages)] <- "missing"
}

get_conflicts_and_support <- function(tree, conflict_quartets) {
  if (length(conflict_quartets) == 0) {cat("WARNING: NO CONFLICTS IN SELECTION!\n")}
  edge_thickness <- rep(1, length(tree$edge.length)+1)
  i <- 1
  for (quat in names(conflict_quartets)) {
    if (quat == "") {next} # not sure why this happens, but apparently it does sometimes
    if (length(strsplit(quat, ",")[[1]]) < 3) {
      cat(paste0("Not a proper Quartet. Will skip: ", quat, "\n"))
      next
      }
    cat(paste0("\rQuartet ", i, "/", length(names(conflict_quartets)), ": ", quat))
    left <- strsplit(quat, "-")[[1]][1]
    right <- strsplit(quat, "-")[[1]][2]
    left1 <- strsplit(left, ",")[[1]][1]
    left2 <- strsplit(left, ",")[[1]][2]
    right1 <- strsplit(right, ",")[[1]][1]
    right2 <- strsplit(right, ",")[[1]][2]
    #identify nodes between tips
    mcanodel <- getMRCA(tree, c(which(tree$tip.label == left1), which(tree$tip.label == left2)))
    mcanoder <- getMRCA(tree, c(which(tree$tip.label == right1), which(tree$tip.label == right2)))
    pathl1 <- nodepath(tree, from=mcanodel, to=which(tree$tip.label == left1))
    pathl2 <- nodepath(tree, from=mcanodel, to=which(tree$tip.label == left2))
    pathr1 <- nodepath(tree, from=mcanoder, to=which(tree$tip.label == right1))
    pathr2 <- nodepath(tree, from=mcanoder, to=which(tree$tip.label == right2))
    npancestor <- nodepath(tree, from=mcanodel, to=mcanoder)
    all_nodes_to_highlight <- c(pathl1, pathl2, pathr1, pathr2)
    
    # now select which edges in the tree
    p <- ggtree(tr = tree, ladderize = FALSE) 
    edgeorder <- data.frame(parent=p$data$parent, node=p$data$node)
    df <- data.frame(one=edgeorder$parent %in% all_nodes_to_highlight, two=edgeorder$node %in% all_nodes_to_highlight)
    
    #df <- data.frame(one=tree$edge[,2] %in% all_nodes_to_highlight, two=tree$edge[,1] %in% all_nodes_to_highlight)
    boolCols <- sapply(df, is.logical)
    selected_rows <- rowSums(df[,boolCols]) == sum(boolCols)
    #selected_rows <- rev(selected_rows)
    edges_with_conflict <- which(selected_rows == TRUE)
    if (conflict_quartets[quat] == 0 ) {
      for (edge in edges_with_conflict) {
        edge_thickness[edge] = edge_thickness[edge] + 1
      }
    } 
    
    i <- i + 1 
  }
  thickness <- data.frame(edge=1:length(edge_thickness), conflict=edge_thickness)
  thickness$logthick <- log(thickness$conflict+1)
  return(thickness)
}


generate_colors <- function(ncols) {
  color_palette <- qualitative_hcl(ncols, alpha=c(0.6,1), l=c(30,100), c=200) # the parameters used here have been working well in a number of scenarios
  return(color_palette)
}
 
get_pdf_height <- function(tree) {
	size_per_tip = 0.2 # in inches
	return(round(length(tree$tip.label) * size_per_tip, digits=0))
}

cat("\nPlotting conflicts...")
treepath1 <- treelist$path[treelist["tree"] == treenames[1]]
bs_cutoff <- strsplit(strsplit(treepath1,"/")[[1]][2],"-")[[1]][2]
algorithm <- strsplit(treepath1,"/")[[1]][3]
alitrim <- strsplit(treepath1,"/")[[1]][4]
prefix1 <- paste( algorithm, alitrim, bs_cutoff, sep="-")

treepath2 <- treelist$path[treelist["tree"] == treenames[2]]
bs_cutoff <- strsplit(strsplit(treepath2,"/")[[1]][2],"-")[[1]][2]
algorithm <- strsplit(treepath2,"/")[[1]][3]
alitrim <- strsplit(treepath2,"/")[[1]][4]
prefix2 <- paste( algorithm, alitrim, bs_cutoff, sep="-")

tree1 <- read.tree(treepath1)
tree2 <- read.tree(treepath2)

if (outgroup != "none") { #reroot tree in case an outgroup was specified
  tree1 <- reroot_my_tree(tree1,outgroup)
  tree2 <- reroot_my_tree(tree2,outgroup)  
} else {cat("No outgroup was set. The plotted tree comparison may look weird.\n")}

conflicts_t <- conflicts[paste0(treenames[1], "-", treenames[2]),]
if (length(names(conflicts_t[conflicts_t == 0])) == 0) {
  print("No conflicts found. Nothing to plot.")
  quit()
}

cat(paste0("Plot will be based on ", as.character(length(names(conflicts_t[conflicts_t == 0]))), " conflicting quartets.\n"))

cat("Extracting conflicts for first tree:\n")
conflicts_info1 <- get_conflicts_and_support(tree1, conflicts_t[conflicts_t == 0])
cat("\nExtracting conflicts for second tree:\n")
conflicts_info2 <- get_conflicts_and_support(tree2, conflicts_t[conflicts_t == 0])

cat("\nPlot conflitcs between trees...\n")
if (lineage_file != "none") {
  cat("Will add lineage information...\n")
  simpdf <- lineages[c("name",level)]
  colnames(simpdf) <- c("name", "lineage")
  simpdf <- add_missing_tips(simpdf, tree2)  
      
  # generate some colors for lineage info:
  cols <- generate_colors(length(na.omit(unique(lineages[,level]))))
  names(cols) <- na.omit(unique(lineages[,level]))
  cols["missing"] <- "black"

  t1 <- ggtree(tree1, branch.length='none', aes(size=conflicts_info1$conflict)) %<+% simpdf + geom_tiplab(aes(color = factor(lineage)),size=4, hjust=0, geom="text") + scale_color_manual(values=cols) + theme(legend.position = c("none")) + scale_size_continuous(range = c(0.2, 5))
  minx <- ggplot_build(t1)$layout$panel_params[[1]]$x.range[1]
  maxx <- ggplot_build(t1)$layout$panel_params[[1]]$x.range[2]
  t1 <- t1+xlim(minx, maxx+40) 
 
  t2 <- ggtree(tree2, branch.length='none', aes(size=conflicts_info2$conflict)) %<+% simpdf + geom_tiplab(aes(color = factor(lineage)),size=4, offset=-40, geom="text") + scale_color_manual(values=cols) + theme(legend.position = c("none")) + scale_size_continuous(range = c(0.2, 5))
  #reverse coordinates and create space for labels
  minx <- ggplot_build(t2)$layout$panel_params[[1]]$x.range[1]
  maxx <- ggplot_build(t2)$layout$panel_params[[1]]$x.range[2]
  t2 <- t2+xlim(maxx+40, minx) 
} else {
  cat("Plotting without lineage information...\n")
  t1 <- ggtree(tree1, branch.length='none', aes(size=conflicts_info1$conflict)) +theme(legend.position = c("none")) + scale_size_continuous(range = c(0.2, 5))
  minx <- ggplot_build(t1)$layout$panel_params[[1]]$x.range[1]
  maxx <- ggplot_build(t1)$layout$panel_params[[1]]$x.range[2]
  t1 <- t1+xlim(minx, maxx+40) +geom_tiplab(size=4, hjust=0)
  
  t2 <- ggtree(tree2, branch.length='none', aes(size=conflicts_info2$conflict)) +theme(legend.position = c("none")) + scale_size_continuous(range = c(0.2, 5))
  #reverse coordinates and create space for labels
  minx <- ggplot_build(t2)$layout$panel_params[[1]]$x.range[1]
  maxx <- ggplot_build(t2)$layout$panel_params[[1]]$x.range[2]
  t2 <- t2+xlim(maxx+40, minx) +geom_tiplab(size=4, offset=-40)
  
}
layout <- "
AABB
"
cat(paste0("Output PDF: conflicts-", treenames[1], "-", treenames[2], "-", as.character(length(names(conflicts_t[conflicts_t == 0]))), "-quartets.pdf\n")) 
pdf(file=paste0("conflicts-", treenames[1],"-",treenames[2],"-", as.character(length(names(conflicts_t[conflicts_t == 0]))), "-quartets.pdf"), width=10, height=get_pdf_height(tree1))
  print(t1 + t2 + theme(legend.position="none")  + plot_layout(design = layout))#+ plot_layout(guides = 'none')# & theme(legend.position='bottom')
garbage <- dev.off()
  

