library(biomartr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
species <- args[1]
wd <- args[2]
species1 <- gsub("_", " ", species)
print(file.path(wd, "results", "downloaded_genomes"))

all_species <- unlist(strsplit(species, ","))

messagereader <- function (fun, ...) { # this function will capture the output from the genome download function
  a <- textConnection("me", "w", local=TRUE) # create a text connection, capturing all text. the first variable (me) will be containing the messages captured
  sink(a, type="message")
  res <- tryCatch (
    expr = {
      fun(...) # evaluate the genome download call
      message("Species download finished")
    },
    error = function(cond){
      message("An error occurred during download!") # message when an error occurs
      return(cond)
    },
    finally = {
      message("trycatch done")
    }
  )
  sink(type="message")
  close(a)
  list(res, messages = me)
}
successful <- c()
failed <- c()
warnings <- c()
species_data <- data.frame(species=character(0), assembly= character(0), status= character(0), stringsAsFactors = F)

for (sp in all_species) {
  if (file.exists(file.path(wd,"results","downloaded_genomes",paste(sp,"_genomic_genbank.fna",sep="")))) {
	print(sp)
	print("File exists, will skip")
	successful <- c(sp, successful)
	next
	}
  sp <- gsub("_", " ", sp)
  print(sp)
  
  mess <- messagereader(getGenome ,db="genbank", organism=sp, path=file.path(wd, "results","downloaded_genomes"), gunzip=T)
  
  my_message <- paste(mess$messages, collapse="\n")
  print(my_message)
  sp <- gsub(" ", "_", sp)
  if (grepl("Unzipping downloaded file", my_message, fixed=T)) {
      path <- file.path(wd, "results","downloaded_genomes",paste(sp,"_genomic_genbank.fna",sep=""))
      species_data[nrow(species_data)+1, ] <- c(sp, path, "successful")
      successful <- c(sp, successful)
      next
    } else if (grepl("no genome file could be found", my_message, fixed=T) == TRUE) {
    species_data[nrow(species_data)+1, ] <- c(sp,"not_downloaded", "failed")
    failed <- c(sp, failed)
    next
    } else {
        print("Some other error occurred!")
	species_data[nrow(species_data)+1, ] <- c(sp,"not_downloaded", "failed")
	failed <- c(sp, failed)
    }

}
#species_data$assembly[species_data$status=="bli"]
#print(species_data)
#message("Completed species: ")
#message(successful)
#message("Failed species:")
#message(failed)
#message("There were warnings for:")
#message(warnings)
if (length(successful)>0) {
fwrite(list(species_data$species[species_data$status=="successful"]), file = file.path(wd, "results","downloaded_genomes", "successfully_downloaded.txt"))
}
if (length(failed)>0) {
  fwrite(list(species_data$species[species_data$status=="failed"]), file = file.path(wd, "results","downloaded_genomes", "not_downloaded.txt"))
}
write.csv(species_data, file.path(wd, "results","downloaded_genomes", "download_overview.txt"), row.names = F, quote=F)

