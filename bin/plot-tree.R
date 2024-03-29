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
cat("This script will plot trees...\n")
seed <- args[2]
treenames <- args[3]
outgroup <- args[4]
lineage_file <- args[5]
level <- args[6]
outprefix <- args[7]
single <- args[8]
  
# reformat commandline argument:
cat(paste0("Used outgroups: ", outgroup, "\n"))
outgroups <- outgroup # for plot annotation
outgroup <- strsplit(outgroup,",")[[1]]
cat(paste0("Trees to plot: ", treenames, "\n")) 
treenames <- strsplit(treenames,",")[[1]]

#set seed if specified
if (seed != "random") {
  set.seed(seed)
  cat(paste0("Random seed: ", seed,"\n"))
} else {
  seed <- sample(0:100000, 1)
  cat(paste0("Random seed: ", seed, "\n"))
}

if (single == "yes") {
  cat("Will only plot one tree instead of two.\nWARNING: With this option plots are hard to standardize, so the quality could be suboptimal.\n")
}

# load lineage information file and fill missing values
if (lineage_file == "none") {
  cat("No lineage file specified. Will not plot lineage information.\n")
  level <- "none"
} else {
  lineages <- read.csv(lineage_file, header=T, sep="\t", na=c(""))
  lineages[is.na(lineages)] <- "missing"
}

bs_support <- 90
pb_support <- 0.95

# the code to draw triangles comes from here:
# https://gist.github.com/jeanmanguy/fe6b1ee46476f29149455124e2a14529

get_offsprings <- function(node_to_collapse, phylo) {
  # todo: assert that node is scalar
  phylo %>%
    tidytree::as_tibble() %>%
    tidytree::offspring(.node = node_to_collapse) %>%
    dplyr::pull(node)
}

get_offspring_tips <- function(phylo, node_to_collapse) {
  # todo: assert that node is scalar
  phylo %>%
    ggtree::fortify() %>%
    tidytree::offspring(.node = node_to_collapse) %>%
    dplyr::filter(isTip) %>%
    dplyr::pull(node)
}

remove_collapsed_nodes <- function(phylo, nodes_to_collapse) {
  nodes <- purrr::map(nodes_to_collapse, get_offsprings, phylo = phylo) %>% unlist()
  phylo %>%
    ggtree::fortify() %>%
    tibble::as_tibble() %>%
    dplyr::filter(!node %in% nodes)
}

get_collapsed_offspring_nodes_coordinates <- function(phylo, nodes) {
  phylo %>%
    ggtree::fortify() %>%
    tibble::as_tibble() %>%
    dplyr::filter(node %in% nodes) %>%
    dplyr::summarise(xmax = max(x), xmin = min(x), ymax = max(y), ymin = min(y))
}

get_collapsed_node_coordinates <- function(phylo, node_to_collapse) {
  # todo: assert that node is scalar
  phylo %>%
    ggtree::fortify() %>%
    tibble::as_tibble() %>%
    dplyr::filter(node == node_to_collapse) %>%
    dplyr::select(x, y)
}

get_triangle_coordinates_ <- function(node, phylo, mode = c("max", "min", "mixed")) {
  mode <- match.arg(mode)
  # for one
  tips_to_collapse <- get_offspring_tips(phylo, node)
  node_to_collapse_xy <- get_collapsed_node_coordinates(phylo, node)
  tips_to_collapse_xy <- get_collapsed_offspring_nodes_coordinates(phylo, tips_to_collapse)
  
  triange_df <- mode %>% switch(
    max = dplyr::data_frame(
      x = c(node_to_collapse_xy$x, tips_to_collapse_xy$xmax, tips_to_collapse_xy$xmax),
      y = c(node_to_collapse_xy$y, tips_to_collapse_xy$ymin, tips_to_collapse_xy$ymax)
    ),
    min = data_frame(
      x = c(node_to_collapse_xy$x, tips_to_collapse_xy$xmin, tips_to_collapse_xy$xmin),
      y = c(node_to_collapse_xy$y, tips_to_collapse_xy$ymin, tips_to_collapse_xy$ymax)
    ),
    mixed = data_frame(
      x = c(node_to_collapse_xy$x, tips_to_collapse_xy$xmin, tips_to_collapse_xy$xmax),
      y = c(node_to_collapse_xy$y, tips_to_collapse_xy$ymin, tips_to_collapse_xy$ymax)
    )
  )
  return(triange_df)
}

get_triangle_coordinates <- function(phylo, nodes, mode = c("max", "min", "mixed")) {
  mode <- match.arg(mode)
  # todo: make sure there is no conflict between nodes (nesting...)
  purrr::map(nodes, get_triangle_coordinates_, phylo = phylo, mode = mode) %>%
    dplyr::bind_rows(.id = "node_collapsed")
}

# till here

is_node_supported <- function(support) {
  if (is.na(support) == TRUE) {return("no")} #treat nodes without support values as unsupported
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
	space_for_legend <- 1.2
	size_per_tip = 0.2 # in inches 
	return(round(space_for_legend + (length(tree$tip.label) * size_per_tip), digits=0))
}
 
get_conflicts <- function(tree, conflict_quartets) {
  if (length(conflict_quartets) == 0) {cat("WARNING: NO CONFLICTS IN SELECTION!\n")}
  edge_thickness <- rep(1, length(tree$edge.length)+1)
  for (quat in conflict_quartets) {
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
    all_nodes_to_highlight <- c(pathl1, pathl2, pathr1, pathr2, npancestor)
    #all_nodes_to_highlight <- c(npancestor)
    
    # now select which edges in the tree
    p <- ggtree(tr = tree, ladderize = FALSE) 
    edgeorder <- data.frame(parent=p$data$parent, node=p$data$node)
    df <- data.frame(one=edgeorder$parent %in% all_nodes_to_highlight, two=edgeorder$node %in% all_nodes_to_highlight)
    
    #df <- data.frame(one=tree$edge[,2] %in% all_nodes_to_highlight, two=tree$edge[,1] %in% all_nodes_to_highlight)
    boolCols <- sapply(df, is.logical)
    selected_rows <- rowSums(df[,boolCols]) == sum(boolCols)
    #selected_rows <- rev(selected_rows)
    edges_with_conflict <- which(selected_rows == TRUE)
    for (edge in edges_with_conflict) {
      edge_thickness[edge] = edge_thickness[edge] + 1
    }
    
  }
  thickness <- data.frame(edge=1:length(edge_thickness), conflict=edge_thickness)
  thickness$logthick <- log(thickness$conflict+1)
  return(thickness)
  
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

reroot_my_tree <- function (tree, outgroup){
	tryCatch(
		{	
			rootnode <- getMRCA(tree, outgroup)
			tree <- root(tree, node=rootnode, resolve.root = TRUE)
			return(c(tree, 1))
 		},
		error = function(e) {
			cat("There was an error while rerooting the tree. Will plot the tree as is. Please check manually\n")
			return(c(tree, 0))
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

if (level == "none" || lineage_file == "none") {
  cat("Plotting tree(s) without lineage information.\n")
  all_supports_list <- list()
  all_names <- c()
  singlet <- FALSE
  for (ntree in 1:length(treenames)) {
    #extract filename information:
    treename <- treenames[ntree]

    bs_cutoff <- strsplit(strsplit(treename,"cutoff-")[[1]][2],"/")[[1]][1]
    algorithm <- strsplit(treename,"/")[[1]][3]
    alitrim <- strsplit(strsplit(treename,".",fixed=T)[[1]][1], "/")[[1]][5]
    hash <- strsplit(strsplit(treename,".",fixed=T)[[1]][2], "/")[[1]][1]
    prefix <- paste( algorithm, alitrim, bs_cutoff, sep="-")

    #reroot tree
    tree <- read.tree(treename)
    if (outgroup != "none") { #reroot tree in case an outgroup was specified
      ret <- reroot_my_tree(tree,outgroup)
      if (ret[[2]] == 0) { outgroups <- "none"}
      tree <- ret[[1]]
    }
    #plot(tree)
    ntips <- length(tree$tip.label)
    
    # this is where we decide how to plot (conflicts or not). Maybe this will be refactored later...
    cat(paste0("Plot tree: ", prefix, "\n"))
    t2 <- ggtree(tree) + theme(legend.position = c("none")) +geom_tiplab()
    t2 <- t2 + coord_cartesian(clip = 'off')
    minx <- ggplot_build(t2)$layout$panel_params[[1]]$x.range[1]
    maxx <- ggplot_build(t2)$layout$panel_params[[1]]$x.range[2]
    t2 <- t2+xlim(minx, maxx+2) # to create space for the labels

    # extract legend then remove it
    cat(paste0("    Write PDF: ",prefix,"-",level,"-tree.pdf\n"))
    pdf(file = paste0(prefix,"-",level,"-tree.pdf"), width = 10, height=get_pdf_height(tree))
    print(t2 + theme(legend.position="none") +  plot_annotation(title = prefix, caption=paste0("Taxonomy: ", level,". Random seed: ", seed,".\nOutgroup: ", outgroups, "\nHash: ", hash), theme=theme(plot.title=element_text(hjust=0.5, size=16))))#+ plot_layout(guides = 'none')# & theme(legend.position='bottom')
    garbage <- dev.off()
    if (singlet == TRUE) {break}
  }
} else { # plot tree when lineage information is available.
  cols <- generate_colors(length(na.omit(unique(lineages[,level]))))
  # keep the older color code below for reference:
  #if (length(na.omit(lineages[,level])) <= 11) {
    # prettier colors when there are not too many different models
    #cols <- brewer.pal(length(na.omit(lineages[,level])), "BrBG")	
  #} else {
    # we need more colors

    #color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    #if (length(na.omit(unique(lineages[,level]))) > length(color)){
    #  # if there are more taxa than available colors sample colors with replacement, this could result in some groups having the same color
    #  cols <- sample(color, length(na.omit(unique(lineages[,level]))), replace=T)
    #} else {cols <- sample(color, length(na.omit(unique(lineages[,level]))))}
  #}
  names(cols) <- na.omit(unique(lineages[,level]))
  cols["missing"] <- "black"

  cat("The awesome code for plotting subtrees as triangles comes from here: https://jean.manguy.eu/post/subtrees-as-triangles-with-ggtree/\n")
  all_supports_list <- list()
  all_names <- c()
  singlet <- FALSE
  for (ntree in 1:length(treenames)) {
    #extract filename information:
    treename <- treenames[ntree]
    bs_cutoff <- strsplit(strsplit(treename,"cutoff-")[[1]][2],"/")[[1]][1]
    algorithm <- strsplit(treename,"/")[[1]][3]
    alitrim <- strsplit(strsplit(treename,".",fixed=T)[[1]][1], "/")[[1]][5]
    hash <- strsplit(strsplit(treename,".",fixed=T)[[1]][2], "/")[[1]][1]
    prefix <- paste( algorithm, alitrim, bs_cutoff, sep="-")
    cat(paste0("Plot tree: ", prefix, "\n"))
    #reroot tree
    tree <- read.tree(treename)
    if (outgroup != "none") { #reroot tree in case an outgroup was specified
      ret <- reroot_my_tree(tree,outgroup)
      if (ret[[2]] == 0) { outgroups <- "none"}
      tree <- ret[[1]]
    }
    #plot(tree)
    ntips <- length(tree$tip.label)
    
    node_names <- c()
    node_names_support <- c()
    nodes_to_collapse <- c()
    node_supports <- c()
    nodes_singletons <- c()
    node_names_singletons <- c()
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
    
    #gather support values
    cat("    Gather support values...\n")
    names(node_supports) <- node_names_support
    node_supports <- node_supports[!is.na(names(node_supports))] # get rid of missing values from NAs when taxon level is missing
    node_supports <- node_supports[names(node_supports) != "missing"]
    all_supports_list[[ntree]] <- node_supports
    all_names <- c(all_names, prefix)
 
    cat("    Collapse tree...\n")
    names(nodes_to_collapse) <- node_names

    if ("missing" %in% names(nodes_to_collapse)) { # nodes whith missing taxonomy info do not need to be collapsed, even if the form a monophyletic group
      nodes_to_collapse <- nodes_to_collapse[names(nodes_to_collapse) != "missing"]
    }

    collapsed_tree_df <- tree %>%
      remove_collapsed_nodes(nodes = nodes_to_collapse)
    
    triangles_df <- tree %>%
      get_triangle_coordinates(nodes_to_collapse)
    
    cat("    Plot right (collapsed) tree...\n")
    t1 <- ggtree::ggtree(collapsed_tree_df) +
      geom_polygon(
        data = triangles_df,
        mapping = aes(group = node_collapsed, fill = node_collapsed),
        color = "#333333"
      ) +
      #scale_fill_brewer(palette = "Set1") +
      scale_fill_manual(values = cols) +
      theme(
        strip.background = element_blank()
      )
    t1 <- t1 + scale_x_reverse()
    
    #need to create a dataframe with just the right taxonomic level to be used for coloring the tree...
    simpdf <- lineages[c("name",level)]
    colnames(simpdf) <- c("name", "lineage")
    simpdf <- add_missing_tips(simpdf, tree)
    cat("    Plot left (uncollapsed) tree...\n")
    if (single == "yes") {
        t2 <- ggtree(tree) %<+% simpdf + geom_tiplab(aes(color = factor(lineage)), size=2, geom="text") +scale_color_manual(values = cols) +theme(legend.position = c("none"))
    } else {
        t2 <- ggtree(tree, branch.length='none') %<+% simpdf + geom_tiplab(aes(color = factor(lineage)), size=2, align=TRUE, geom="text") +scale_color_manual(values = cols) +theme(legend.position = c("none"))
    } 
    #t2 <- ggtree(tree, branch.length='none') %<+% simpdf + geom_tiplab(aes(color = factor(lineage)), size=2, align=TRUE, geom="text") + theme(legend.position = c("none"))
    t2 <- t2 + coord_cartesian(clip = 'off')
    minx <- ggplot_build(t2)$layout$panel_params[[1]]$x.range[1]
    maxx <- ggplot_build(t2)$layout$panel_params[[1]]$x.range[2]
    if (single == "yes") {
       	t2 <- t2+xlim(minx, maxx+2) # to create space for the labels
    } else {
       	t2 <- t2+xlim(minx, maxx+40) # to create space for the labels
    }
    # now we create clade labels for the tree; this now also includes singleton nodes
    names(nodes_singletons) <- node_names_singletons
    clade_label_df <- as.data.frame(c(nodes_to_collapse, nodes_singletons))
    clade_label_df$name <- rownames(clade_label_df)
    colnames(clade_label_df) <- c("node", "name")
    if (single == "yes") {
       t2 <- t2 + geom_cladelab(data = clade_label_df, mapping = aes(node = node, label = name, color = name), fontsize = 2, offset=1, offset.text=0.1)
    } else {
       t2 <- t2 + geom_cladelab(data = clade_label_df, mapping = aes(node = node, label = name, color = name), fontsize = 2, offset=27, offset.text=0.3)
    }
    # extract legend then remove it
    legend <- get_legend(t1)
    t1 <- t1 + theme(legend.position="none")
  
    layout <- "
    AABB
    "

    cat(paste0("    Write PDF: ",prefix,"-",level,"-tree.pdf\n"))
    pdf(file = paste0(prefix,"-",level,"-tree.pdf"), width = 10, height=get_pdf_height(tree))
    if (single == "yes") {
        print(t2 + theme(legend.position="none") + plot_annotation(title = prefix, caption=paste0("Taxonomy: ", level,". Random seed:", seed,".\nOutgroup: ", outgroups), theme=theme(plot.title=element_text(hjust=0.5, size=16))))#+ plot_layout(guides = 'none')# & theme(legend.position='bottom')
    } else {
        print(t2 + t1 + theme(legend.position="none") +  plot_annotation(tag_levels = 'A', title = prefix, caption=paste0("A) Topology without branch lengths. B) Topology with branch lengths and collapsed taxonomic groups. Taxonomy: ", level,". Random seed: ", seed,".\nOutgroup: ", outgroups, "\nHash: ", hash), theme=theme(plot.title=element_text(hjust=0.5, size=16))) + plot_layout(design = layout))#+ plot_layout(guides = 'none')# & theme(legend.position='bottom')
    } 
    garbage <- dev.off()
    if (single == TRUE) {break}
  }
# this last parts creates an overview plot indicating (un)supported groups based on the specified lineage level.
  if (length(all_names) > 1) {  # only necessary if there is more than 1 tree
    cat("Generating overview support plot now...\n")
    support_df <- as.data.frame(do.call(cbind, all_supports_list))
    colnames(support_df) <- all_names
    support_df[all_names] <- lapply(support_df[all_names] , factor)
  
    support_cols <- c("#d53e4f", "#ffffbf","#3288bd")
    names(support_cols) <- c("no", "notmono", "yes")
    sdf <- melt(t(support_df))
    colnames(sdf) <- c("tree", "group", "supported")
    mywidth <- round(((0.1*length(sdf$tree)) + 5)/2, digits=0)
    pdf(file=paste0("compare-", level, ".pdf"), width=mywidth, height=10)
      print(ggplot(sdf, aes(x = tree, y = group)) + geom_tile(aes(fill=supported),colour = "white") + scale_fill_manual(values=support_cols) +theme_minimal()+theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(angle = 45, vjust = 0, hjust=0))+scale_x_discrete(position = "top"))
    garbage <- dev.off()
    cat(paste0("Overview support plot saved to: compare-", level, ".pdf\n")) 
  }
}




