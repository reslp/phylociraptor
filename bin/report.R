library(ggplot2)
library(reshape2)
library(viridis)
args <- commandArgs(trailingOnly=TRUE)

busco_set <- args[1]
mywidth <- as.numeric(args[2])
table_file <- args[3]
summary_file <- args[4]
outfile1 <- args[5]
outfile2 <- args[6]


data2 <- read.csv(summary_file, sep="\t", header=T, stringsAsFactors = F, row.names = 1)
total <- data2$total[1]

data <- read.csv(table_file, sep="\t", header=T, stringsAsFactors = F)
rownames(data) <- data$species
data$species <- NULL
data$percent_complete <- NULL

mdata <- melt(as.matrix(data))
base_size <- 7
p <- ggplot(mdata, aes(Var2, Var1, fill=value)) +geom_tile() + scale_fill_gradient(low="white", high="orange") +theme(axis.ticks = element_blank(),legend.position="None",axis.text.y = element_text(size=base_size*0.7), axis.text.x = element_text(size = base_size *0.8,angle=-90, hjust = 0,vjust=0.5, colour = "grey50"))+  xlab("Species") + ylab("BUSCO gene") + ggtitle(paste("Overview of complete single-copy BUSCOs per species (total genes=", as.character(total),") BUSCO set: ",busco_set, sep=""))


pdf(file=outfile1, width=mywidth)
print(p)
dev.off()

as.character(total)
data2$total <- NULL
data2$complete <- NULL
mdata2 <- melt(as.matrix(data2))
base_size <- 12
colnames(mdata2) <- c("species", "BUSCO category", "no. of genes")
p2 <- ggplot(mdata2, aes_string(fill="`BUSCO category`", y="`no. of genes`", x="`species`")) +geom_bar(position="stack", stat="identity")+theme(axis.ticks.x = element_blank(), panel.grid=element_blank(), panel.background = element_blank(), axis.text.x = element_text(size = base_size *0.8,angle=-45, hjust = 0,vjust=1, colour = "grey50"))+scale_fill_manual(values=c("#6198B8", "#E37332", "#C6E88B", "#BCB2CF"))+ggtitle(paste("Overview of BUSCO results (total genes=", as.character(total),") BUSCO set: ",busco_set, sep=""))

pdf(file=outfile2, width=mywidth/5)
print(p2)
dev.off()


