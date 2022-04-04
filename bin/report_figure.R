args <- commandArgs(trailingOnly=TRUE)

wd <- getwd()
setwd(paste0(wd,"/bin"))

library(yaml)
config_data <- read_yaml(paste0("../",args[1]))

library(patchwork)
library(RColorBrewer)
library(ggpubr)
library(stringr)
library(reshape2)

downloaded_genomes_statistics_file <- "../results/statistics/downloaded_genomes_statistics.txt"
failed_genome_downloads_file <- "../results/downloaded_genomes/not_downloaded.txt"
successfull_genome_downloads_file <- "../results/downloaded_genomes/successfully_downloaded.txt"
local_species_file <- "../results/statistics/local_species.txt"

pars_sites <-config_data$filtering$min_parsimony_sites
pars_sites <- strtoi(pars_sites)

# check if seed was specified or not:
if (is.null(config_data$seed) || config_data$seed == "") {
	config_data$seed <- "random"
} else {
	set.seed(strtoi(config_data$seed))
	} 

## functions:

read_my_overview <-function(dat) {
  return(read.csv(dat,header=F,sep="\t"))
}

read_my_csv <-function(dat) {
  return(read.csv(dat,header=T,sep="\t"))
}

multi_merge <- function(d) {
  df <- d[[1]]
  if (length(d) >= 2) {
    for (i in 2:length(dat)) {
      df <- merge(df, d[[i]], by="alignment", all=T, suffixes=c("",""))
      #df <- full_join(df, d[[i]], by="alignment", suffix = c("", ""))
    }
  }
  return(df)
}

multi_merge2 <- function(d) {
  df <- d[[1]]
  if (length(d) >= 2) {
    for (i in 2:length(dat)) {
      df <- merge(df, d[[i]], by="gene", all=T, suffixes=c("",""))
      #df <- full_join(df, d[[i]], by="gene", suffix = c("", ""))
    }
  }
  return(df)
}


## functions end
if (file.exists(failed_genome_downloads_file))
{
  if (file.size(failed_genome_downloads_file) > 0) {
    failed_sp <- read.table(failed_genome_downloads_file, header=F)
    failed <- paste(failed_sp$V1, sep="\n")
  } else { failed <- numeric()}
} else {
  failed <- numeric()
}
if (file.exists(successfull_genome_downloads_file)){
  #print(file.info(successfull_genome_downloads_file)$size)
  if (file.size(successfull_genome_downloads_file) > 0) {
    success_sp <- read.table(successfull_genome_downloads_file, header=F)
    success <- paste(success_sp$V1, sep="\n")
  } else {success <- numeric()}
} else {
  success <- numeric()
}
if (file.exists(local_species_file)){
  if (file.size(local_species_file) > 0) {
    local_sp <- read.table(local_species_file, header=F)
    local <- paste(local_sp$V1, sep="\n")
  } else {local <- numeric()}
} else {
  local <- numeric()
}

### definitions:

colors <- c(total="grey69", pass="#2166ac", success="#2166ac", fail="#b2182b", failed="#b2182b", local="#7fbc41")
total <- length(failed) + length(success) + length(local)

genome_information <- data.frame(total=c(total), success=c(length(success)), failed=c(length(failed)), local=c(length(local)))
genome_information <- melt(genome_information)
colnames(genome_information) <- c("category", "no. of genomes")
setup_plot <- ggplot(data=genome_information, aes(x=category, y=`no. of genomes`, fill=category)) +geom_bar(stat="identity") + scale_fill_manual(values=colors) 
setup_plot <- setup_plot + ggtitle("Download overview") + theme_minimal() + theme(legend.position="bottom", plot.title = element_text(size = 9))
setup_plot <- setup_plot + theme(axis.title.x = element_blank()) + theme(legend.key.size = unit(8, "pt")) + theme(legend.margin=margin(0,0,0,0, "cm"), legend.title = element_blank())
setup_plot <- setup_plot + theme(legend.position = "none") # remove legend altogether, it is actually not needed
layout_setup <- 
"AA
BB
BB
BB
"

#plotsetup <- setup_plot_annotations + setup_plot + plot_layout(design=layout_setup) 
#plotsetup <- plotsetup + theme(text = element_text(size = 8))
plotsetup <- setup_plot +plot_annotation(title="setup", theme = theme(plot.title = element_text(hjust = 0.5)))
#plotsetup

#### ORTHOLOGY

busco_set_dir <- list.dirs("../results/orthology/busco/busco_set/",recursive=F)
busco_overview_file <- paste0(busco_set_dir,"/dataset.cfg")
print(busco_overview_file)
busco_summary_file <- "../results/statistics/busco_summary.txt"
orthology_filtering_genomes_file <- "../results/statistics/orthology_filtering_genomes_statistics.txt"
orthology_filtering_genes_file <- "../results/statistics/orthology_filtering_gene_statistics.txt"


if (file.exists(busco_overview_file)){
  busco_overview <- read.table(busco_overview_file, sep="=", header=F)
  busco_set <- busco_overview$V2[busco_overview$V1 == "name"]
  nbuscos <- busco_overview$V2[busco_overview$V1 == "number_of_BUSCOs"]
} else {
  busco_set <- "unknown (maybe BUSCO set was not downloaded yet)"
  nbuscos <- "unknown (maybe BUSCO set was not downloaded yet)"
}

if (file.exists(orthology_filtering_genomes_file)) {
  ortho_filter_genomes <- read.table(orthology_filtering_genomes_file, sep="\t", header=F)
  total_genomes <- nrow(ortho_filter_genomes)
  colnames(ortho_filter_genomes) <- c("sample", "status", "completeness", "cutoff")
    if (length(ortho_filter_genomes$sample[ortho_filter_genomes$status == "FAILED"]) > 0 ){
      failed_genomes <-length(ortho_filter_genomes$sample[ortho_filter_genomes$status == "FAILED"])
    } else {failed_genomes <- 0}
  successfull_genomes <- total_genomes - failed_genomes
  cutoff <- ortho_filter_genomes$cutoff[1]
  orthology_genome_data <- data.frame(total=c(total), pass=c(successfull_genomes), fail=c(failed_genomes))
  orthology_genome_data <- melt(orthology_genome_data)
  orthology_genome_data
  orthology_genome_data$variable <- gsub("_", "\n", orthology_genome_data$variable)
  colnames(orthology_genome_data) <- c("BUSCO analysis status", "no. of samples")
  orthology_genome_data <- within(orthology_genome_data , `BUSCO analysis status` <- factor(`BUSCO analysis status`, levels=c("total", "pass", "fail")))
  
  orthology_plot <- ggplot(data=orthology_genome_data, aes(x=`BUSCO analysis status`, y=`no. of samples`, fill=`BUSCO analysis status`)) +geom_bar(stat="identity") + scale_fill_manual(values=colors)
  orthology_plot <- orthology_plot + ggtitle("samples") + theme_minimal() + theme(legend.position="bottom",legend.title = element_blank())
  orthology_plot <- orthology_plot + theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))+ theme(text = element_text(size = 8)) +theme(legend.key.size = unit(0.2, "cm"))
  orthology_plot <- orthology_plot + theme(legend.position = "none")
  } else {
  orthology_plot <- ggplot()
  # do nothing
}#


if (file.exists(orthology_filtering_genes_file)) {
  ortho_filter_genes <- read.table(orthology_filtering_genes_file, sep="\t", header=F)
  total_genes <- nrow(ortho_filter_genes)
  total_genes
  colnames(ortho_filter_genes) <- c("gene", "status", "nseqs", "cutoff")
  if (length(ortho_filter_genes$gene[ortho_filter_genes$status == "FAILED"]) > 0 ){
    failed_genes <-length(ortho_filter_genes$gene[ortho_filter_genes$status == "FAILED"])
    failed_genes
  } else {failed_genes <- 0}
  failed_genes
  successful_genes <- total_genes - failed_genes
  successful_genes
  cutoff_gene <- ortho_filter_genes$cutoff[1]
  orthology_gene_data <- data.frame(total=c(total_genes), pass=c(successful_genes), fail=c(failed_genes))
  orthology_gene_data <- melt(orthology_gene_data)
  orthology_gene_data$variable <- gsub("_", "\n", orthology_gene_data$variable)
  orthology_gene_data
  
  colnames(orthology_gene_data) <- c("BUSCO analysis status", "no. of genes")
  orthology_gene_data <- within(orthology_gene_data , `BUSCO analysis status` <- factor(`BUSCO analysis status`, levels=c("total", "pass", "fail")))
  
  orthology_geneplot <- ggplot(data=orthology_gene_data, aes(x=`BUSCO analysis status`, y=`no. of genes`, fill=`BUSCO analysis status`)) +geom_bar(stat="identity") + scale_fill_manual(values=colors)
  orthology_geneplot <- orthology_geneplot + ggtitle("genes") + theme_minimal() + theme(legend.position="bottom",legend.title = element_blank())
  
  orthology_geneplot <- orthology_geneplot + theme(axis.title.x = element_blank(),plot.title = element_text(hjust = 0.5))+ theme(text = element_text(size = 8)) +theme(legend.key.size = unit(0.2, "cm"))
  orthology_geneplot <- orthology_geneplot + theme(legend.position = "none")
} else {
  # do nothing
  orthology_geneplot <- ggplot()
}

plotortho <- orthology_plot + orthology_geneplot + plot_layout(ncol=2) #& theme(legend.position="bottom",legend.title = element_blank())
plotortho <- plotortho &  theme(legend.margin = margin(0,0,0,0, "cm"))
plotortho <- plotortho + plot_annotation(title="orthology", theme = theme(plot.title = element_text(hjust = 0.5)))
plotortho
########## alignments


dirs <- list.dirs("../results/statistics/", recursive=F)
dirs <- grep("align-*", dirs, value = TRUE)

if (length(dirs) != 0) {
  dat <- list()
  i <- 1
  for (dir in dirs) { # this code needs to be expanded and the dataframe integrated into the one created in the next chunk
    files <- list.files(path=dir, pattern="*overview*", full.names=T)
    data <- do.call(rbind,lapply(files,read_my_overview))	
    dat[[i]] <- data
    i <- i + 1
  }
  
  alignment_data_combined_overview <- do.call(cbind, dat)
}


if (length(dirs)==0) {
  cat("<br><b> Alignment statistic files not found. Did you run phylociraptor align?</b>\n")
} else {
  #cat("<br><b>Used aligner(s) and settings:</b>")
  #pars_sites <- tail(strsplit(alignment_data_combined_overview[1,1], " ")[[1]], 1)
  #pars_sites <- strtoi(pars_sites)
  dat <- list()
  i <- 1
  alignment_statistics <- data.frame(aligner=c(), total=numeric(), pass=numeric(),fail=numeric())
  aligner_names <- c()
  for (dir in dirs) {
    print(dir)
    files <- list.files(path=dir, pattern="*statistics*", full.names=T)
    data <- do.call(rbind,lapply(files,read_my_csv))	
    aligner_name <- strsplit(dir,split="align-")[[1]][2]
    aligner_names <- c(aligner_names, aligner_name)
    #cat(paste0("<br><b>Number of alignments for ", strsplit(dir,split="align-")[[1]][2],": </b>", nrow(data), "\n"))
    #data[,6] <- cell_spec(data[,6], background = ifelse(data[,6] < pars_sites, "red", ""))
    if (i ==1) {
      data <- data[,c(1,4:8)]
    } else{
      data <- data[,c(4:8)]
    }
    below_rows <- length(which(data$nparsimony < pars_sites))
    above_rows <- nrow(data) - below_rows
    alignment_statistics <- rbind(alignment_statistics, c(aligner_name, nrow(data), above_rows, below_rows))
  }
}

colnames(alignment_statistics) <- c("aligner", "total", "pass", "fail")
maxtotal <- strtoi(max(alignment_statistics$total))
#rownames(alignment_statistics) <- alignment_statistics$aligner

#alignment_statistics$total <- NULL
alignment_statistics_melted <- melt(alignment_statistics, id="aligner")
alignment_statistics_melted$value <- as.numeric(alignment_statistics_melted$value)
colnames(alignment_statistics_melted) <- c("aligner", "status", "no. of alignments")
alignment_statistics_melted <- within(alignment_statistics_melted, `status` <- factor(`status`, levels=c("total", "pass", "fail")))

full_alignments_plot <- ggplot(alignment_statistics_melted, aes(fill=status, y=`no. of alignments`, x=aligner)) + 
  geom_bar(position="dodge", stat="identity")+ ggtitle("full") + theme_minimal() + theme(legend.position="bottom",legend.title = element_blank(),plot.title = element_text(size = 8))  + scale_fill_manual(values=colors)
full_alignments_plot <- full_alignments_plot + theme(axis.title.x = element_blank(), legend.margin=margin(0,0,0,0, "cm"), legend.key.size = unit(9, "pt"), text=element_text(size=7), axis.text.x = element_text(angle=60, vjust=1, hjust=1), plot.title = element_text(hjust = 0.5, size=8))
full_alignments_plot <- full_alignments_plot + scale_y_continuous(limits = c(0, maxtotal))
alignment_legend <- get_legend(full_alignments_plot) #extract legend
full_alignments_plot <- full_alignments_plot+ theme(legend.position = "none") #now remove
full_alignments_plot
#trimmed
dirs <- list.dirs("../results/statistics/", recursive=F)
dirs <- grep("/trim-*", dirs, value = TRUE)
dirs

#if (length(dirs) != 0) {
#  dat <- list()
#  i <- 1
#  for (dir in dirs) { # this code needs to be expanded and the dataframe integrated into the one created in the next chunk
#    files <- list.files(path=dir, pattern="*overview*", full.names=T)
#    data <- do.call(rbind,lapply(files,read_my_overview))	
#    dat[[i]] <- data
#    i <- i + 1
#  }
#  
#  alignment_data_combined_overview <- do.call(cbind, dat)
#}


if (length(dirs)==0) {
  cat("<br><b> Alignment statistic files not found. Did you run phylociraptor align?</b>\n")
} else {
  read_my_csv <-function(dat) {
    return(read.csv(dat,header=T,sep="\t"))
  }
  #cat("<br><b>Used aligner(s) and settings:</b>")
  #pars_sites <- tail(strsplit(alignment_data_combined_overview[1,1], " ")[[1]], 1)
  #pars_sites <- strtoi(pars_sites)
  #first_aligner <- paste(head(strsplit(alignment_data_combined_overview[1,1], " ")[[1]], length(strsplit(alignment_data_combined_overview[1,1], " ")[[1]])-1), collapse= " ")
  #first_aligner <- gsub("None", "", first_aligner)
  #cat(paste0("<br><b>1. </b>", first_aligner,"<br>"))	
  #if (length(alignment_data_combined_overview[1,]) == 2) {
  #  second_aligner <- paste(head(strsplit(alignment_data_combined_overview[1,2], " ")[[1]], length(strsplit(alignment_data_combined_overview[1,2], " ")[[1]])-1), collapse= " ")
  #  second_aligner <- gsub("None", "", second_aligner)
  #  cat(paste0("<b>2. </b>", second_aligner,"<br>"))	
  #}
  #cat(paste0("<b>Parsimony informative sites cutoff: </b>", pars_sites, "<br>"))
  dat <- list()
  i <- 1
  trimmed_alignment_statistics <- data.frame(aligner=c(),trimmer=c(), total=numeric(), pass=numeric(), fail=numeric())
  trimmer_names <- c()
  for (dir in dirs) {
    print(dir)
    files <- list.files(path=dir, pattern="*statistics*", full.names=T)
    data <- do.call(rbind,lapply(files,read_my_csv))	
    combination <- strsplit(dir,split="trim-")[[1]][2]
    #aligner_name <- strsplit(combination,split="-")[[1]][1]
    trimmer_names <- c(trimmer_names, strsplit(combination,split="-")[[1]][2])
    #cat(paste0("<br><b>Number of alignments for ", strsplit(dir,split="align-")[[1]][2],": </b>", nrow(data), "\n"))
    #data[,6] <- cell_spec(data[,6], background = ifelse(data[,6] < pars_sites, "red", ""))
    if (i ==1) {
      data <- data[,c(1,4:8)]
    } else{
      data <- data[,c(4:8)]
    }
    below_rows <- length(which(data$nparsimony < pars_sites))
    above_rows <- nrow(data) - below_rows
    trimmed_alignment_statistics <- rbind(trimmed_alignment_statistics, c(combination, nrow(data), above_rows, below_rows))
  }
}

colnames(trimmed_alignment_statistics) <- c("combination", "total", "pass", "fail")
trimmed_alignment_statistics$total <- NULL
trimmed_alignment_statistics_melted <- melt(trimmed_alignment_statistics, id="combination")
trimmed_alignment_statistics_melted$combination <- gsub("-", "\n",trimmed_alignment_statistics_melted$combination)

trimmed_alignment_statistics_melted$value <- as.numeric(trimmed_alignment_statistics_melted$value)


colnames(trimmed_alignment_statistics_melted) <- c("combination", "status", "no. of alignments")
trimmed_alignment_statistics_melted <- within(trimmed_alignment_statistics_melted, `status` <- factor(`status`, levels=c("total", "pass", "fail")))

trimmed_alignments_plot <- 
  ggplot(trimmed_alignment_statistics_melted, aes(fill=status, y=`no. of alignments`, x=combination)) + 
  geom_bar(position="dodge", stat="identity")+ ggtitle("trimmed") + theme_minimal() + theme(legend.position="bottom",legend.title = element_blank())  + scale_fill_manual(values=colors)
trimmed_alignments_plot <- trimmed_alignments_plot + theme(axis.title.x = element_blank(), legend.margin=margin(0,0,0,0, "cm"), legend.key.size = unit(4, "pt"), text=element_text(size=7), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.text.x = element_text(angle=60, vjust=1, hjust=1),plot.title = element_text(hjust = 0.5, size=8))
trimmed_alignments_plot <- trimmed_alignments_plot + scale_y_continuous(limits = c(0, maxtotal)) + theme(legend.position="none")
trimmed_alignments_plot <- trimmed_alignments_plot + theme(legend.position = "none")

#filtered
dirs <- list.dirs("../results/statistics/", recursive=F)
dirs <- grep("/filter-*", dirs, value = TRUE)
dirs

if (length(dirs)==0) {
  cat("<br><b> Alignment statistic files not found. Did you run phylociraptor align?</b>\n")
} else {
  read_my_csv <-function(dat) {
    return(read.csv(dat,header=T,sep="\t"))
  }
  #cat("<br><b>Used aligner(s) and settings:</b>")
  #pars_sites <- tail(strsplit(alignment_data_combined_overview[1,1], " ")[[1]], 1)
  #pars_sites <- strtoi(pars_sites)
  #first_aligner <- paste(head(strsplit(alignment_data_combined_overview[1,1], " ")[[1]], length(strsplit(alignment_data_combined_overview[1,1], " ")[[1]])-1), collapse= " ")
  #first_aligner <- gsub("None", "", first_aligner)
  #cat(paste0("<br><b>1. </b>", first_aligner,"<br>"))	
  #if (length(alignment_data_combined_overview[1,]) == 2) {
  #  second_aligner <- paste(head(strsplit(alignment_data_combined_overview[1,2], " ")[[1]], length(strsplit(alignment_data_combined_overview[1,2], " ")[[1]])-1), collapse= " ")
  #  second_aligner <- gsub("None", "", second_aligner)
  #  cat(paste0("<b>2. </b>", second_aligner,"<br>"))	
  #}
  #cat(paste0("<b>Parsimony informative sites cutoff: </b>", pars_sites, "<br>"))
  dat <- list()
  i <- 1
  filtered_alignment_statistics <- data.frame(aligner=c(),trimmer=c(), total=numeric(), pass=numeric(), fail=numeric())
  for (dir in dirs) {
    print(dir)
    files <- list.files(path=dir, pattern="*statistics*", full.names=T)
    data <- do.call(rbind,lapply(files,read_my_csv))	
    combination <- strsplit(dir,split="filter-")[[1]][2]
    #aligner_name <- strsplit(combination,split="-")[[1]][1]
    #trimmer_name <- strsplit(combination,split="-")[[1]][2]
    #cat(paste0("<br><b>Number of alignments for ", strsplit(dir,split="align-")[[1]][2],": </b>", nrow(data), "\n"))
    #data[,6] <- cell_spec(data[,6], background = ifelse(data[,6] < pars_sites, "red", ""))
    if (i ==1) {
      data <- data[,c(1,4:8)]
    } else{
      data <- data[,c(4:8)]
    }
    below_rows <- length(which(data$nparsimony < pars_sites))
    above_rows <- nrow(data) - below_rows
    filtered_alignment_statistics <- rbind(filtered_alignment_statistics, c(combination, nrow(data), above_rows, below_rows))
  }
}

colnames(filtered_alignment_statistics) <- c("combination", "total", "pass", "fail")
#filtered_alignment_statistics
filtered_alignment_statistics$total <- NULL
filtered_alignment_statistics_melted <- melt(filtered_alignment_statistics, id="combination")
filtered_alignment_statistics_melted$combination <- gsub("-", "\n",filtered_alignment_statistics_melted$combination)

filtered_alignment_statistics_melted$value <- as.numeric(filtered_alignment_statistics_melted$value)
colnames(filtered_alignment_statistics_melted) <- c("combination", "status", "no. of alignments")
filtered_alignment_statistics_melted <- within(filtered_alignment_statistics_melted, `status` <- factor(`status`, levels=c("total", "pass", "fail")))

filtered_alignments_plot <- 
  ggplot(filtered_alignment_statistics_melted, aes(fill=status, y=`no. of alignments`, x=combination)) + 
  geom_bar(position="dodge", stat="identity")+ ggtitle("filtered") + theme_minimal() + theme(legend.position="bottom",legend.title = element_blank())  + scale_fill_manual(values=colors)
filtered_alignments_plot <- filtered_alignments_plot + theme(axis.title.x = element_blank(),legend.margin=margin(0,0,0,0, "cm"), legend.key.size = unit(4, "pt"), text=element_text(size=7), axis.title.y=element_blank(), axis.text.y=element_blank(),axis.text.x = element_text(angle=60, vjust=1, hjust=1),plot.title = element_text(hjust = 0.5))
filtered_alignments_plot <- filtered_alignments_plot + scale_y_continuous(limits = c(0, maxtotal)) + theme(legend.position="none")
filtered_alignments_plot<- filtered_alignments_plot + theme(legend.position = "none")

layout <- "
ABBBCCC
ABBBCCC
ABBBCCC
ABBBCCC
ABBBCCC
ABBBCCC
ABBBCCC
ABBBCCC
ABBBCCC
ABBBCCC
ABBBCCC
ABBBCCC
ABBBCCC
ABBBCCC
ABBBCCC
ABBBCCC
ABBBCCC
ABBBCCC
ABBBCCC
ABBBCCC
ABBBCCC
ABBBCCC
ABBBCCC
DDDDDDD
"
combined_alignment_plot <- full_alignments_plot + trimmed_alignments_plot + filtered_alignments_plot + alignment_legend + plot_layout(design=layout)# + theme(legend.position="bottom",legend.title = element_blank())
combined_alignment_plot 
#combined_alignment_plot <- combined_alignment_plot & theme(legend.position="bottom",legend.title = element_blank())
combined_alignment_plot <- combined_alignment_plot + plot_annotation(title="alignments", theme = theme(plot.title = element_text(hjust = 0.5)))
#combined_alignment_plot


############### Modeltest
files <- list.files(path="../results/modeltest/", pattern ="best_models_*", recursive=F)
if (length(files) > 0) {
  i <- 1
  dat <- list()
  all_genes <- c()
  names <- c("alignment")
  titles <- c()
  plots <- list()
  models <- c()
  for (file in files) {
    # get reusable and meaningful name
    title <- str_split(file, "_")[[1]][3:4]
    # use other stringsplit function which supports fixed=T
    title[2] <- strsplit(title[2], ".", fixed=T)[[1]][1]
    title <- paste(title, collapse="-")
    names <- c(names, title)
    titles <- c(titles, title)
    dat[[i]] <- read.table(paste0("../results/modeltest/", file), sep="\t", header=F)
    colnames(dat[[i]]) <- c("alignment", "model")
    all_genes <- c(all_genes,dat[[i]]$alignment)
    models <- c(models, unique(dat[[i]]$model))
    i <- i + 1
  }
  # get unqiue number of models in dataset and create colors for them:
  models <- unique(models)
  if (length(models) <= 11) {
    # prettier colors when there are not too many different models
    cols <- brewer.pal(length(models), "BrBG")	
  } else {
    # we need more colors
    color <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    cols <- sample(color, length(models))
  }
  names(cols) <- models
  for (i in 1:length(dat)) {
    dd <- as.data.frame(table(dat[[i]]$model))
    colnames(dd) <- c("model", "count")
    p <- ggplot(dd, aes(x="", y=count, fill=model)) + geom_bar(stat="identity", width=1) + coord_polar("y", start=0) + theme_void() + ggtitle(titles[i]) + scale_fill_manual(limits=models, values=cols) + theme(text = element_text(size=7), legend.margin=margin(5,5,5,5, "pt"), legend.key.size = unit(4, "pt")) + 
      theme(plot.title = element_text(hjust = 0.5))
    plots[[i]] <- p
  }
  if (length(plots) == 6) { # the numbers below are hardcoded atm which is not very elegant
    modeltest_plots <- ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]], common.legend = TRUE, legend="bottom")
  }
  if (length(plots) == 5) {
    modeltest_plots <- ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], common.legend = TRUE, legend="bottom")
  }
  if (length(plots) == 4) {
    modeltest_plots <- ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], common.legend = TRUE, legend="bottom")
  }
  if (length(plots) == 3) {	
    modeltest_plots <- ggarrange(plots[[1]], plots[[2]], plots[[3]], common.legend = TRUE, legend="bottom")
  }
  if (length(plots) == 2) {	
    modeltest_plots <- ggarrange(plots[[1]], plots[[2]], common.legend = TRUE, legend="bottom")
  }
  if (length(plots) == 1) {	
    modeltest_plots <- ggarrange(plots[[1]], common.legend = TRUE, legend="right")
  }
  #data_combined <- do.call(multi_merge, list(dat))
  #colnames(data_combined) <- names
  #data_combined %>% kbl(escape=F) %>% kable_paper("hover", full_width = F)%>% scroll_box(width = "100%", height = "600px")	
}
modeltest_plots <- annotate_figure(modeltest_plots, top="modeltest")
modeltest_plots <- modeltest_plots + theme(legend.title = element_blank())
modeltest_plots
##### Gene Trees
files <- list.files(path="../results/modeltest/", pattern ="genetree_filter_*", recursive=F)
files
if (length(files) > 0) {
  i <- 1
  dat <- list()
  all_genes <- c()
  names <- c("gene", "bootstrap_mean", "check")
  all_colnames <- c("gene", "aligner", "trimmer", "bootstrap_sum", "bootstrap_n", "bootstrap_mean")
  titles <- c()
  plots <- list()
  models <- c()
  thresh <- 0
  genes_above_thresh_outstring <- c()
  for (file in files) {
    # get reusable and meaningful name
    title <- str_split(file, "_")[[1]][3:4]
    # use other stringsplit function which supports fixed=T
    title[2] <- strsplit(title[2], ".", fixed=T)[[1]][1]
    title <- paste(title, collapse="-")
    names <- c(names,  c("gene", "bootstrap_mean"))
    titles <- c(titles, title)
    dat[[i]] <- read.table(paste0("../results/modeltest/", file), sep="\t", header=F)
    colnames(dat[[i]]) <- all_colnames 
    #thresh <- dat[[i]]$bootstrap_cutoff[1]
    dat[[i]] <- dat[[i]][,c("gene", "bootstrap_mean")]
    
    #genes_above_thresh_outstring <- c(genes_above_thresh_outstring, paste0("<b>Number of ", title, " trees above bootstrap threshold (",thresh,"): ", nrow(dat[[i]][dat[[i]]$check == "OK", ]), "</b><br>\n"))
    #all_genes <- c(all_genes,dat[[i]]$gene)
    i <- i + 1
    
  }
  genetree_data_combined <- do.call(multi_merge2, list(dat))
  #plot_data <- genetree_data_combined[, grepl("bootstrap_mean", names(data_combined))]
  #plot_data <- as.data.frame(plot_data) #make sure data structure is dataframe, this is necessary because this will be a vector for a single parameter combination
  #colnames(plot_data) <- titles
  #p <- p %>% layout(hovermode = "y unified") #unfortunately this does not work together with custom text.
  #print(htmltools::tagList(list(ggplotly(p)))) # this works only when the plot is being shown outside the if statement below! weird!!
  colnames(genetree_data_combined) <- c("gene", titles)		
  #options(knitr.kable.NA = '') # do not display NA values
  #data_combined %>% kbl(escape=F) %>% kable_paper("hover", full_width = F)%>% scroll_box(width = "100%", height = "600px")	
}

genetree_data_combined_melted <- melt(genetree_data_combined)
colnames(genetree_data_combined_melted) <- c("gene", "combination", "mean bootstrap support per tree") 
genetree_data_combined_melted$combination <- gsub("-", "\n", genetree_data_combined_melted$combination)
genetree_plot <- ggplot(data=genetree_data_combined_melted, aes(x=combination, y=`mean bootstrap support per tree`, fill=combination)) + geom_violin() + theme_minimal() + theme(legend.position="none")  + scale_fill_brewer(palette="Spectral")
genetree_plot <- genetree_plot + theme(axis.title.x = element_blank(), text=element_text(size=7))
genetree_plot <- annotate_figure(genetree_plot, top="genetrees")
#genetree_plot
#if (p != "") {
# p <- p  %>% layout(autosize = F, width = 10, height = 10, showlegend=FALSE) %>% config(displayModeBar=FALSE)# this works only in combination with the print statement above
#} else { p <- "gene tree result statistics not found" }
#p


#additional information plot

trimmers_outstring <- (if (length(trimmer_names) > 1) {paste(unique(trimmer_names), collapse=", ")} else {unique(trimmer_names)})

get_aligner_settings_string <- function(){
  first_aligner <- paste(head(strsplit(alignment_data_combined_overview[1,1], " ")[[1]], length(strsplit(alignment_data_combined_overview[1,1], " ")[[1]])-1), collapse= " ")
  first_aligner <- gsub("None", "", first_aligner)
  outstring <- c(paste0("   ", first_aligner, "\n"))
  if (length(alignment_data_combined_overview[1,]) == 2) {
    second_aligner <- paste(head(strsplit(alignment_data_combined_overview[1,2], " ")[[1]], length(strsplit(alignment_data_combined_overview[1,2], " ")[[1]])-1), collapse= " ")
    second_aligner <- gsub("None", "", second_aligner)
    #cat(paste0("<b>2. </b>", second_aligner,"<br>"))
    outstring <- c(outstring, paste0("   ", second_aligner, "\n"))
  }
  return(paste(outstring, collapse=""))
}


setup_text <- paste0(
"Total samples: ", toString(total),
"\nSuccessfully downloaded: ", toString(length(success)),
"\nFailed to download: ", toString(length(failed)),
"\nLocally provided: ", toString(length(local))
)
orthology_text <- paste0(
  "\n\nMethod: BUSCO ", toString(config_data$busco$version),
  "\nBUSCO set: ", busco_set,
  "\nNo. of BUSCO genes: ", toString(nbuscos),
  "\nMinimum BUSCO completeness: ", config_data$filtering$cutoff
)
config_data
align_text <- paste0(
  "\n\nAligners and settings:\n", get_aligner_settings_string(),
  "Sequence type: ", config_data$filtering$seq_type,
  "\nParsimony informative sites cut-off: ", toString(pars_sites),
  "\nAlignment trimmers:\n   ", trimmers_outstring,
  "\nFiltering duplicated sequences: ", config_data$filtering$dupseq,
  "\nMinimum number of sequences per alignment: ", config_data$filtering$minsp  
)

reprod_text <- paste0(
  "\n\nAnalysis seed: ", config_data$seed
)

additional_info <- ggplot() +theme_void() +
  annotate("text", x = 0, y=10, label="setup:", hjust=0, size=4, fontface="bold.italic")+
  annotate("text", x = 0.2, y = 8.7, label = setup_text, hjust = 0, size=3)+
  annotate("text", x = 0, y = 7.5, label="orthology:", hjust=0, size=4, fontface="bold.italic")+
  annotate("text", x = 0.2, y = 6.7, label = orthology_text, hjust = 0, size=3)+
  annotate("text", x = 0, y = 5, label="alignments:", hjust=0, size=4, fontface="bold.italic")+
  annotate("text", x = 0.2, y = 3.1, label = align_text, hjust = 0, size=3)+
  annotate("text", x = 0.2, y = 0.8, label = reprod_text, hjust = 0, size=3, fontface="bold")+
  coord_cartesian(ylim = c(0, 10), xlim=c(0,10), clip = "off") #+ theme(text = element_text(size = 5)) 
additional_info <- annotate_figure(additional_info, top="overview")
additional_info
#additional_info <- ggparagraph(size=7, text=paste0("Total samples: ",toString(total),"\nSuccessfully downloaded: ", toString(length(success)),"\nFailed download: ", toString(length(failed)),"\nLocally provided: ", toString(length(local))))



pdf(file="report-figure.pdf", width=11.3, height=8.7)
p <- ggarrange(plotsetup, plotortho,combined_alignment_plot,modeltest_plots,genetree_plot,additional_info, ncol=3, nrow=2,labels="AUTO")
p
dev.off()
