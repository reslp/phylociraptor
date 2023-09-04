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
suppressMessages(library(ggnewscale))
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
outgroups <- outgroup #for information on PDF
outgroup <- strsplit(outgroup,",")[[1]]
treenames <- strsplit(treenames,",")[[1]]

# some treename checks to avoid problem when first tree has larger number than second tree:
if (strtoi(strsplit(treenames[1], "T")[[1]][2]) > strtoi(strsplit(treenames[2], "T")[[1]][2])){
   tn2 <- treenames[2]
   treenames[2] <- treenames[1]
   treenames[1] <- tn2
}
print(treenames[1])
print(treenames[2])

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
  cat(paste0("Loading treelist file: ", treelistfile, "..."))
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
    if (quat == "") {next} # not sure why this happens, but apparently it sometimes does
    if (length(strsplit(quat, ",")[[1]]) < 3) {
      cat(paste0("Not a proper Quartet. Will skip: ", quat, "\n"))
      next
      }
    cat(paste0("\rConflict quartet ", i, "/", length(names(conflict_quartets)), ": ", quat))
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
  thickness$scaledconflict <- scale_numbers(thickness$conflict)
  return(thickness)
}

scale_numbers <- function(range) {
 minr <- min(range)
 maxr <- max(range)
 normalized_range <- c()
 for (value in range) {
   normalized_value <- (10-1)*((value-minr)/(maxr-minr))+1
   normalized_range <- c(normalized_range, normalized_value)
 }
 return(normalized_range)
}

get_node_names_for_bars <- function(tree) {
    node_names <- c()
    node_names_support <- c()
    nodes_to_collapse <- c()
    node_supports <- c()
    nodes_singletons <- c()
    node_names_singletons <- c()
    ntips <- length(tree$tip.label)
    for (name in unique(lineages[,level])) {
      which_tips <- lineages$name[lineages[level] == name][lineages$name[lineages[level] == name] %in% tree$tip.label]
      node <- getMRCA(tree, c(which_tips))
      if (length(node) != 0)  {
        #check if this taxon level is monophyletic
        descendants <- tree$tip.label[getDescendants(tree=tree, node=node)]
        descendants <- descendants[!is.na(descendants)]
        tiptax <- lineages[level][,1][lineages$name %in% descendants]
        
        # the next check is true if all tips have the same taxonomic level
        if (length(unique(tiptax)) == 1){
          #cat(paste0("    OK ", name, "\n"))
          #implement check if all descendents have the same label
          nodes_to_collapse <- c(nodes_to_collapse, node) 
          node_supports <- c(node_supports, is_node_supported(as.double(tree$node.label[node-ntips])))
          node_names <- c(node_names, name)
          node_names_support <- c(node_names_support, name)
        } else {
      	if (name != "missing") { # do not print this if taxononmy level is "missing"
          	cat(paste0("    ", name, " appear to be PARAPHYLETIC. Tips in the tree will be colored but group will not be collapsed or otherwise highlighted.\n"))
          }
          node_supports <-c(node_supports, "notmono")
          node_names_support <- c(node_names_support, name)
        }
      } else if (length(node) == 0) {
      	if (length(which_tips) == 1) {
        		nodes_singletons <- c(nodes_singletons, match(which_tips, tree$tip.label))
        		node_names_singletons <- c(node_names_singletons, name)
        		#cat(paste0("    SINGLETON ", name, "\n"))
      	} else {
      		#debug code:
      		#cat(paste0("  ", name, " is present in lineage file but not in the tree.\n"))
      		next
      		}
      } 
    }

    names(nodes_to_collapse) <- node_names
    names(nodes_singletons) <- node_names_singletons
    clade_label_df <- as.data.frame(c(nodes_to_collapse, nodes_singletons))
    clade_label_df$name <- rownames(clade_label_df)
    colnames(clade_label_df) <- c("node", "name")
    return(clade_label_df)
}

is_node_supported <- function(support) {
  bs_support <- 90
  pb_support <- 0.95
  if (length(support) == 0) {return("no")}
  if (support <= 1) { #we are dealing with posterior probabilities
    if (support <= pb_support){
      return("no") # no support
    } else {return("yes")}
  } else { # we are dealing with bootstrap values
    if (support <= bs_support){
      return("no") #no support
    } else {return("yes")}
    
  }
}

generate_colors <- function(ncols) {
  color_palette <- qualitative_hcl(ncols, alpha=c(0.6,1), l=c(30,100), c=200) # the parameters used here have been working well in a number of scenarios
  return(color_palette)
}
 
get_pdf_height <- function(tree) {
	size_per_tip = 0.2 # in inches
	return(round(length(tree$tip.label) * size_per_tip, digits=0))
}

cat("Plotting conflicts...\n")
treepath1 <- treelist$path[treelist["tree"] == treenames[1]]
bs_cutoff1 <- strsplit(strsplit(treepath1,"/")[[1]][2],"-")[[1]][2]
algorithm1 <- strsplit(treepath1,"/")[[1]][3]
alitrim1 <- strsplit(treepath1,"/")[[1]][4]
prefix1 <- paste( algorithm1, alitrim1, bs_cutoff1, sep="-")

treepath2 <- treelist$path[treelist["tree"] == treenames[2]]
bs_cutoff2 <- strsplit(strsplit(treepath2,"/")[[1]][2],"-")[[1]][2]
algorithm2 <- strsplit(treepath2,"/")[[1]][3]
alitrim2 <- strsplit(treepath2,"/")[[1]][4]
prefix2 <- paste( algorithm2, alitrim2, bs_cutoff2, sep="-")

tree1 <- read.tree(treepath1)
tree2 <- read.tree(treepath2)

if (outgroup != "none") { #reroot tree in case an outgroup was specified
  tree1 <- reroot_my_tree(tree1,outgroup)
  tree2 <- reroot_my_tree(tree2,outgroup)  
} else {cat("No outgroup was set. The plotted tree comparison may look weird.\n")}

conflicts_t <- conflicts[paste0(treenames[1], "-", treenames[2]),]
if (length(names(conflicts_t[conflicts_t == 0])) == 0) {
  cat("No conflicts found. Nothing to plot.\n")
  quit()
}

cat(paste0("Plot will be based on ", as.character(length(names(conflicts_t[conflicts_t == 0]))), " conflicting quartets.\n"))

cat("\nExtracting conflicts for first tree:")
conflicts_info1 <- get_conflicts_and_support(tree1, conflicts_t[conflicts_t == 0])
cat("\n\nExtracting conflicts for second tree:")
conflicts_info2 <- get_conflicts_and_support(tree2, conflicts_t[conflicts_t == 0])

cat("\n\nPreparing tree plots...\n")
if (lineage_file != "none") {
  cat("Will add lineage information...\n")
  simpdf <- lineages[c("name",level)]
  colnames(simpdf) <- c("name", "lineage")
  simpdf <- add_missing_tips(simpdf, tree2)  
      
  # generate some colors for lineage info:
  cols <- generate_colors(length(na.omit(unique(lineages[,level]))))
  names(cols) <- na.omit(unique(lineages[,level]))
  cols["missing"] <- "black"
  #cols <- cols[names(cols) != "missing"]
  t1_clade_df <- get_node_names_for_bars(tree1)
  # sort colors so they match in the plot:
  cols2 <- c()
  for (n in t1_clade_df$name) {
     cols2 <- c(cols2, cols[n])
  }
  t1_clade_df$cols <- cols2
  # first tree:
  cat(paste0("Tree 1: ", treenames[1], "-", prefix1, "\n"))
  t1 <- ggtree(tree1, branch.length='none', aes(color=conflicts_info1$scaledconflict), size=1) + scale_color_continuous(low="black", high="red")
  t1 <- t1 + geom_cladelab(data=t1_clade_df, node=t1_clade_df$node, label=t1_clade_df$name, textcolor=t1_clade_df$cols, barcolor=t1_clade_df$cols, fontsize = 2, offset=70, offset.text=0.3)
  t1 <- t1 + new_scale("color")
  t1 <- t1 %<+% simpdf + geom_tiplab(aes(color=factor(lineage)), size=4, hjust=0, geom="text") +scale_color_manual(values=cols) 
  minx <- ggplot_build(t1)$layout$panel_params[[1]]$x.range[1]
  maxx <- ggplot_build(t1)$layout$panel_params[[1]]$x.range[2]
  t1 <- t1+xlim(minx, maxx+40) 
  t1 <- t1 + theme(legend.position="none") + ggtitle(prefix1) + theme(plot.title = element_text(hjust=0.5))
  
  # second tree:
  cat(paste0("Tree 2: ", treenames[2], "-", prefix2, "\n"))
  t2 <- ggtree(tree2, branch.length='none', aes(color=conflicts_info2$scaledconflict), size=1) +scale_color_continuous(low="black", high="red") 
  t2 <- t2 + new_scale("color") 
  t2 <- t2 %<+% simpdf + geom_tiplab(aes(color = factor(lineage)),size=4, offset=-40, geom="text") + scale_color_manual(values=cols) + theme(legend.position = c("none")) + ggtitle(prefix2) + theme(plot.title = element_text(hjust=0.5))
  #reverse coordinates and create space for labels
  minx <- ggplot_build(t2)$layout$panel_params[[1]]$x.range[1]
  maxx <- ggplot_build(t2)$layout$panel_params[[1]]$x.range[2]
  t2 <- t2+xlim(maxx+40, minx+20) 
} else {
  cat("Plotting without lineage information...\n")
  cat(paste0("Tree 1: ", treenames[1], "-", prefix1, "\n"))
  t1 <- ggtree(tree1, branch.length='none', aes(color=conflicts_info1$conflict), size=1) +theme(legend.position = c("none")) + scale_color_continuous(low="black", high="red") + ggtitle(prefix1) + theme(plot.title = element_text(hjust=0.5))
  minx <- ggplot_build(t1)$layout$panel_params[[1]]$x.range[1]
  maxx <- ggplot_build(t1)$layout$panel_params[[1]]$x.range[2]
  t1 <- t1 + new_scale("color")
  t1 <- t1+xlim(minx, maxx+40) +geom_tiplab(size=4, hjust=0)
  
  cat(paste0("Tree 2: ", treenames[2], "-", prefix2, "\n"))
  t2 <- ggtree(tree2, branch.length='none', aes(color=conflicts_info2$conflict), size=1) +theme(legend.position = c("none")) + scale_color_continuous(low="black", high="red") + ggtitle(prefix2) + theme(plot.title = element_text(hjust=0.5))
  #reverse coordinates and create space for labels
  minx <- ggplot_build(t2)$layout$panel_params[[1]]$x.range[1]
  maxx <- ggplot_build(t2)$layout$panel_params[[1]]$x.range[2]
  t2 <- t2 + new_scale("color")
  t2 <- t2+xlim(maxx+40, minx) +geom_tiplab(size=4, offset=-40)
  
}
layout <- "
AABB
"
# get some numbers for plot title:
nconflicts <- as.character(length(names(conflicts_t[conflicts_t == 0])))
nquartets <- as.character(length(names(conflicts_t)))
cat(paste0("Output PDF: conflicts-", treenames[1], "-", treenames[2], "-", nconflicts, "-quartets.pdf\n")) 
pdf(file=paste0("conflicts-", treenames[1],"-",treenames[2],"-", nconflicts, "-quartets.pdf"), width=10, height=get_pdf_height(tree1))
  print(t1 + t2 + theme(legend.position="none")  + plot_layout(design = layout) + plot_annotation(title = paste0("An estimation of topological conflict between two trees based on\n", nconflicts, " conficts found in ", nquartets, " analyzed quartets of tips"), caption=paste0("Taxonomic level: ", level,". Random seed: ", seed,". Trees: ",  treenames[1], "(",prefix1,"), ", treenames[2], "(", prefix2, ")\nOutgroup: ", outgroups),  theme = theme(plot.title=element_text(hjust=0.5, size=16))))#+ plot_layout(guides = 'none')# & theme(legend.position='bottom')
garbage <- dev.off()
  

