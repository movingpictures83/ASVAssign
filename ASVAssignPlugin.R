#First must have dada2 and Bioconductor installed: https://benjjneb.github.io/dada2/dada-installation.html 
#Must also have phyloseq installed: http://joey711.github.io/phyloseq/install.html
dyn.load(paste("RPluMA", .Platform$dynlib.ext, sep=""))
source("RPluMA.R")

library(dada2); 
packageVersion("dada2")

set.seed(1234)

input <- function(inputfile) {
  pfix = prefix()
  if (length(pfix) != 0) {
     pfix <- paste(pfix, "/", sep="")
  }

  parameters <- read.table(inputfile, as.is=T);
  rownames(parameters) <- parameters[,1]
  nochim <<- paste(pfix, toString(parameters["nochim", 2]), sep="")
  DB <<- paste(pfix, toString(parameters["DB", 2]), sep="")
  specDB <<- paste(pfix, toString(parameters["specDB", 2]), sep="")
  seqtab.nochim <<- readRDS(nochim)

}

run <- function() {
dim(seqtab.nochim) #how many ASVs you have

# Inspect distribution of sequence lengths

table(nchar(getSequences(seqtab.nochim)))

#sum(seqtab.nochim)/sum(seqtab) #percentage of merged sequences that are NOT chimeras

#Track reads through pipeline---------------------------

#getN <- function(x) sum(getUniques(x))
#track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(merger1.nochim, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
#colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
#rownames(track) <- sample.names
#head(track)
#mean(track[,5]/track[,4]) #percentage of denoised sequences that have merged

#Assign taxonomy----------------------------------------
#print("ASSIGNING...")
#print(DB)
#print(specDB)
taxa <<- assignTaxonomy(seqtab.nochim, DB, multithread=TRUE)
#print("SPECIES...")
taxa <<- addSpecies(taxa, specDB)
#print("DONE")
taxa.print <<- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#Extract info--------------------------------------------

asv_seqs <<- colnames(seqtab.nochim)
asv_headers <<- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <<- paste(">ASV", i, sep="_")
}
#print(asv_seqs)
#print(asv_headers)
#print("DONE RUN")
}

output <- function(outputfile) {
# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, paste(outputfile, ".fa", sep=""))
asv_fasta

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
#print(colnames(asv_tab))
colnames(asv_tab) <- gsub('_F_filt.fastq.gz','',colnames(asv_tab))
#print(colnames(asv_tab))
write.table(asv_tab, paste(outputfile, ".counts.tsv", sep=""), sep="\t", quote=F, col.names=NA)
asv_tab

asv_tab1 = as.data.frame(asv_tab)
#print(colnames(asv_tab1))
colnames(asv_tab1) <- gsub('_F_filt.fastq.gz','',colnames(asv_tab1))
write.csv(asv_tab1, paste(outputfile, ".tab.csv", sep=""))


# tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, paste(outputfile, ".taxonomy.tsv", sep=""), sep="\t", quote=F, col.names=NA)
asv_tax

#for export of the tax table as a csv

#asv_tax1 = as.data.frame(asv_tax) 
#write.csv(asv_tax1, file = '~/Desktop/asv_tax.csv')
}


#Bonus: Handoff to phyloseq------------------------------

#library(phyloseq); packageVersion("phyloseq")
#library(Biostrings); packageVersion("Biostrings")
#library(ggplot2); packageVersion("ggplot2")

#Creating abundance bar plots------------------------------------

#TAB = otu_table(asv_tab, taxa_are_rows=TRUE)
#TAX = tax_table(asv_tax)

#physeq = phyloseq(TAB, TAX)

#plot_bar(physeq, fill = "Phylum")
#plot_bar(physeq, fill = "Family")
#plot_bar(physeq, fill = "Genus")

#nsamples(physeq)
#sample_names(physeq)
#sample_variables(physeq)

