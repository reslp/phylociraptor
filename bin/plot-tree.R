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

args <- commandArgs(trailingOnly=TRUE)

# variables which need to be passed to the script from the command line
runmode <- args[1] # different run modes should be: plot, conflicts and
wd <- args[2]
setwd(wd)
if (runmode == "plot") {
  cat("Will plot trees...\n")
  seed <- args[3]
  treenames <- args[4]
  outgroup <- args[5]
  lineage_file <- args[6]
  level <- args[7]
  outprefix <- args[8]
  
} else if (runmode == "conflicts") {
  cat("Will plot conflicts...\n")
  seed <- args[3]
  treenames <- args[4]
  outgroup <- args[5]
  lineage_file <- args[6]
  level <- args[7]
  conflictfile <- args[8]
  treelistfile <- args[9]
} else {
  cat("Runmode not recognized...\n")
  quit(1)
}


# reformat commandline argument:
outgroup <- strsplit(outgroup,",")[[1]]
treenames <- strsplit(treenames,",")[[1]]

#print(outgroup)
#print(treenames)
#print(lineage_file)

#set seed if specified
if (seed != "random") {
  set.seed(seed)
}

# load lineage information file and fill missing values
if (lineage_file == "none") {
  cat("No lineage file specified. Will not plot lineage information.\n")
} else {
  lineages <- read.csv(lineage_file, header=T, sep="\t", na=c(""))
  lineages[is.na(lineages)] <- "missing"
}


#load and format quartet conflicts file:
if (runmode == "conflicts") {
  if (conflictfile != "none") {
    conflicts <- read.csv(conflictfile, header=T, check.names=FALSE)
    rownames(conflicts) <- conflicts$quartet
    conflicts$quartet <- NULL
    conflicts <- t(conflicts)
    treelist <- read.csv(treelistfile, header=F, check.names=FALSE, sep="\t")
    colnames(treelist) <- c("tree", "path")
    #print(conflicts)
  } else {
    cat("Conflicts file not found. Will stop.\n")
    quit(1)
  }
  
}




# select comparison for visualization: this is debug code!
#conflicts_t <- t(conflicts[13,])
#table(conflicts_t)
#rownames(conflicts_t)
#conflicts_t["Drosophila_azteca,Drosophila_rufa-Musca_domestica,Phortica_variegata",]
#conflict_quartets <- conflicts_t
#conflict_quartets <- names(conflicts_t[conflicts_t == 0,])


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

generate_color <- function(label) {
  colors <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  
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

# generate some colors
if (runmode == "plot") {
  if (level == "none") {
    cat("Plotting tree(s) without lineage information.\n")
    all_supports_list <- list()
    all_names <- c()
    single <- FALSE
    for (ntree in 1:length(treenames)) {
      #extract filename information:
      treename <- treenames[ntree]
      bs_cutoff <- strsplit(strsplit(treename,"/")[[1]][2],"-")[[1]][2]
      algorithm <- strsplit(treename,"/")[[1]][3]
      alitrim <- strsplit(treename,"/")[[1]][4]
      prefix <- paste( algorithm, alitrim, bs_cutoff, sep="-")
      
      #reroot tree
      tree <- read.tree(treename)
      if (outgroup != "none") { #reroot tree in case an outgroup was specified
        rootnode <- getMRCA(tree, outgroup)
        tree <- root(tree, node=rootnode, resolve.root = TRUE)
      }
      #plot(tree)
      ntips <- length(tree$tip.label)
      
      
      # this is where we decide how to plot (conflicts or not). Maybe this will be refactored later...
      cat(paste0("Plot tree: ", prefix, "\n"))
      if (runmode == "conflicts") {
        cat("Extracting Conflicts:\n")
        conflicts_info <- get_conflicts_and_support(tree, conflict_quartets)
        
        t2 <- ggtree(tree, branch.length='none', aes(size=conflicts_info$conflict)) +theme(legend.position = c("none")) +geom_tiplab()
        
      } else {
        cat("    plotting without conflicts...\n")
        t2 <- ggtree(tree, branch.length='none') + theme(legend.position = c("none")) +geom_tiplab()
        
      }
      t2 <- t2 + coord_cartesian(clip = 'off')
      minx <- ggplot_build(t2)$layout$panel_params[[1]]$x.range[1]
      maxx <- ggplot_build(t2)$layout$panel_params[[1]]$x.range[2]
      t2 <- t2+xlim(minx, maxx+40) # to create space for the labels
      # now we create clade labels for the tree
      #clade_label_df <- as.data.frame(nodes_to_collapse)
      #clade_label_df$name <- rownames(clade_label_df)
      #colnames(clade_label_df) <- c("node", "name")
      
      #t2 <- t2 + geom_cladelab(data = clade_label_df, mapping = aes(node = node, label = name, color = name), fontsize = 2, offset=27, offset.text=0.3)
      #ggplot_build(t2)
      #factor(lineages["order"][,1])
      
      # extract legend then remove it
      
      
      
      cat(paste0("    write PDF: ",prefix,"-",level,"-tree.pdf\n"))
      pdf(file = paste0(prefix,"-",level,"-tree.pdf"), width = 10, height=70)
      print(t2 + theme(legend.position="none"))#+ plot_layout(guides = 'none')# & theme(legend.position='bottom')
      #plot.new()
      #print(legend)
      dev.off()
      if (single == TRUE) {break}
    }
  } else {
    if (length(na.omit(lineages[,level])) <= 11) {
      # prettier colors when there are not too many different models
      cols <- brewer.pal(length(na.omit(lineages[,level])), "BrBG")	
    } else {
      # we need more colors
      color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
      if (length(na.omit(unique(lineages[,level]))) > length(color)){
        # if there are more taxa than available colors sample colors with replacement, this could result in some groups having the same color
        cols <- sample(color, length(na.omit(unique(lineages[,level]))), replace=T)
      } else {cols <- sample(color, length(na.omit(unique(lineages[,level]))))}
    }
    
    names(cols) <- na.omit(unique(lineages[,level]))
    
    cat("The awesome code for plotting subtrees as triangles comes from here: https://jean.manguy.eu/post/subtrees-as-triangles-with-ggtree/\n")
    all_supports_list <- list()
    all_names <- c()
    single <- FALSE
    for (ntree in 1:length(treenames)) {
      #extract filename information:
      treename <- treenames[ntree]
      bs_cutoff <- strsplit(strsplit(treename,"/")[[1]][2],"-")[[1]][2]
      algorithm <- strsplit(treename,"/")[[1]][3]
      alitrim <- strsplit(treename,"/")[[1]][4]
      prefix <- paste( algorithm, alitrim, bs_cutoff, sep="-")
      cat(paste0("Plot tree: ", prefix, "\n"))
      #reroot tree
      tree <- read.tree(treename)
      if (outgroup != "none") { #reroot tree in case an outgroup was specified
        rootnode <- getMRCA(tree, outgroup)
        tree <- root(tree, node=rootnode, resolve.root = TRUE)
      }
      #plot(tree)
      ntips <- length(tree$tip.label)
      
      node_names <- c()
      node_names_support <- c()
      nodes_to_collapse <- c()
      node_supports <- c()
      for (name in unique(lineages[,level])) {
        which_tips <- lineages$name[lineages[level] == name][lineages$name[lineages[level] == name] %in% tree$tip.label]
        #print(which_tips)
        node <- getMRCA(tree, c(which_tips))
        #print(node)
        if (length(node) != 0)  {
          #check if this taxon level is monophyletic
          descendants <- tree$tip.label[getDescendants(tree=tree, node=node)]
          descendants <- descendants[!is.na(descendants)]
          tiptax <- lineages[level][,1][lineages$name %in% descendants]
          
          # the next check is true if all tips have the same taxonomic level
          if (length(unique(tiptax)) == 1){
            cat(paste0("    OK ", name, "\n"))
            #implement check if all descendents have the same label
            nodes_to_collapse <- c(nodes_to_collapse, node) 
            node_supports <- c(node_supports, is_node_supported(as.double(tree$node.label[node-ntips])))
            node_names <- c(node_names, name)
            node_names_support <- c(node_names_support, name)
          } else {cat(paste0("    PARAPHYLETIC ", name, "\n"))
            node_supports <-c(node_supports, "notmono")
            node_names_support <- c(node_names_support, name)
          }
        } else {
          cat(paste0("    SINGLETON ", name, "\n"))
        }
      }
      
      #gather support values
      cat("    Gather support values...\n")
      names(node_supports) <- node_names
      node_supports <- node_supports[!is.na(names(node_supports))] # get rid of missing values from NAs when taxon level is missing
      all_supports_list[[ntree]] <- node_supports
      all_names <- c(all_names, prefix)
      
      cat("    Collapse tree...\n")
      names(nodes_to_collapse) <- node_names
      collapsed_tree_df <- tree %>%
        remove_collapsed_nodes(nodes = nodes_to_collapse)
      
      triangles_df <- tree %>%
        get_triangle_coordinates(nodes_to_collapse)
      
      cat("    Plot collapsed tree...\n")
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
      

      cat("    Plot second tree without conflicts...\n")
      t2 <- ggtree(tree, branch.length='none') %<+% simpdf + geom_tiplab(aes(color = factor(lineage)),size=2, align=TRUE, geom="text") +scale_color_manual(values=cols) +theme(legend.position = c("none"))

      t2 <- t2 + coord_cartesian(clip = 'off')
      minx <- ggplot_build(t2)$layout$panel_params[[1]]$x.range[1]
      maxx <- ggplot_build(t2)$layout$panel_params[[1]]$x.range[2]
      t2 <- t2+xlim(minx, maxx+40) # to create space for the labels
      # now we create clade labels for the tree
      clade_label_df <- as.data.frame(nodes_to_collapse)
      clade_label_df$name <- rownames(clade_label_df)
      colnames(clade_label_df) <- c("node", "name")
      
      t2 <- t2 + geom_cladelab(data = clade_label_df, mapping = aes(node = node, label = name, color = name), fontsize = 2, offset=27, offset.text=0.3)
      
      # extract legend then remove it
      legend <- get_legend(t1)
      t1 <- t1 + theme(legend.position="none")
    
      layout <- "
      AABB
      "

      cat(paste0("    Write PDF...",prefix,"-",level,"-tree.pdf\n"))
      pdf(file = paste0(prefix,"-",level,"-tree.pdf"), width = 10, height=70)
      print(t2 + t1 + theme(legend.position="none")  + plot_layout(design = layout))#+ plot_layout(guides = 'none')# & theme(legend.position='bottom')
      dev.off()
      if (single == TRUE) {break}
    }
    
    cat("Generating overview support plot now...\n")
    support_df <- as.data.frame(do.call(cbind, all_supports_list))
    support_df
    colnames(support_df) <- all_names
    support_df
    support_df[all_names] <- lapply(support_df[all_names] , factor)
    
    support_cols <- c("#d53e4f", "#ffffbf","#3288bd")
    names(support_cols) <- c("no", "notmono", "yes")
    sdf <- melt(t(support_df))
    colnames(sdf) <- c("tree", "group", "supported")
    pdf(file=paste0("compare-", level, ".pdf"), width=10, height=10)
    print(ggplot(sdf, aes(x = tree, y = group)) + geom_tile(aes(fill=supported),colour = "white") + scale_fill_manual(values=support_cols) +theme_minimal()+theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x = element_text(angle = 45, vjust = 0, hjust=0))+scale_x_discrete(position = "top"))
    dev.off()
  }
}else if (runmode == "conflicts") {
  cat("Plotting conflicts...")
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
    rootnode <- getMRCA(tree1, outgroup)
    tree1 <- root(tree1, node=rootnode, resolve.root = TRUE)
    rootnode <- getMRCA(tree2, outgroup)
    tree2 <- root(tree2, node=rootnode, resolve.root = TRUE)
  } else {cat("No outgroup was set. The plotted tree comparison may look weird.\n")}
  
  conflicts_t <- conflicts[paste0(treenames[1], "-", treenames[2]),]
  if (length(names(conflicts_t[conflicts_t == 0])) == 0) {
    print("No conflicts found. Nothing to plot.")
    quit()
  }

  cat(paste0("Plot will be based on ", as.character(length(names(conflicts_t[conflicts_t == 0]))), " conflicting quartets.\n"))
 
  cat("Extracting Conflicts for first tree:\n")
  conflicts_info1 <- get_conflicts_and_support(tree1, conflicts_t[conflicts_t == 0])
  cat("Extracting Conflicts for second tree:\n")
  conflicts_info2 <- get_conflicts_and_support(tree2, conflicts_t[conflicts_t == 0])

  cat("Plot conflitcs between trees...\n")
  if (lineage_file != "none") {
    cat("Will add lineage information...\n")
    simpdf <- lineages[c("name",level)]
    colnames(simpdf) <- c("name", "lineage")
    t1 <- ggtree(tree1, branch.length='none', aes(size=conflicts_info1$conflict)) %<+% simpdf + geom_tiplab(aes(color = factor(lineage)),size=4, hjust=0, geom="text")  +theme(legend.position = c("none")) + scale_size_continuous(range = c(0.2, 5))
    minx <- ggplot_build(t1)$layout$panel_params[[1]]$x.range[1]
    maxx <- ggplot_build(t1)$layout$panel_params[[1]]$x.range[2]
    t1 <- t1+xlim(minx, maxx+40) 
    
    t2 <- ggtree(tree2, branch.length='none', aes(size=conflicts_info2$conflict)) %<+% simpdf + geom_tiplab(aes(color = factor(lineage)),size=4, offset=-40, geom="text")+theme(legend.position = c("none")) + scale_size_continuous(range = c(0.2, 5))
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
  pdf(file=paste0("conflicts-", treenames[1],"-",treenames[2], ".pdf"), width=10, height=70)
    print(t1 + t2 + theme(legend.position="none")  + plot_layout(design = layout))#+ plot_layout(guides = 'none')# & theme(legend.position='bottom')
  dev.off()
  

}
  





