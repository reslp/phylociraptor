library(biomartr)
args = commandArgs(trailingOnly=TRUE)
species <- args[1]
wd <- args[2]
species1 <- gsub("_", " ", species)
print(file.path(wd, "results", "downloaded_genomes"))

all_species <- unlist(strsplit(species, ","))

for (sp in all_species) {
  if (file.exists(file.path(wd,"results","downloaded_genomes",paste(sp,"_genomic_genbank.fna",sep="")))) {
	print(sp)
	print("File exists, will skip")
	next
	}
  sp <- gsub("_", " ", sp)
  print(sp)
  getGenome(db="genbank", organism=sp, path=file.path(wd, "results","downloaded_genomes"), gunzip=T)
  #getGFF(db="genbank", organism=sp, path=file.path(wd, "results","downloaded_genomes"))
}
#getGenome(db="genbank", organism=species1, path=file.path(wd, "results","downloaded_genomes"))
#is.genome.available(db="genbank", organism="Ascobolus immersus")
#?getGenome()
