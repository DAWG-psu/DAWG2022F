#Install Packages####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("dada2", "Biostrings", "ShortRead"))
library(dada2)
library(Biostrings)
library(ShortRead)
  
## Below command should be run under "Terminmal"

## Install cutadapt - primer trimming tool

# pip3 install --user --upgrade cutadapt
# cutadapt path - "/storage/home/t/tuc289/.local/bin"

## Download today's files

# cd ~/scratch    # move to your own scratch directory
# wget --no-check-certificate "https://pennstateoffice365-my.sharepoint.com/:u:/g/personal/tuc289_psu_edu/ESbdgUo93X1HgdMFsDNpo-sB0YiZ-UApG4RXrAtJOSYSDw?e=zRrTF5&download=1" -O DAWG.zip    
# unzip DAWG.zip    # decompress .zip file
# cd DAWG     # move to the directory with files
# ls     # check if you can see the list of files
# pwd   # copy the output of this command - your file directory

## Processing/filtering 16S rRNA Illumina amplicon dataset

#1-3. Set Working directory and Sort Files ####
setwd("~/where your downloaded fastq files are")
path <- setwd("~/where your downloaded fastq files are")

fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz"))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz"))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`,1)
sample.names

#4. Check for primers ####
FWD <- "GTGYCAGCMGCCGCGGTAA"
REV <- "GGACTACNVGGGTWTCTAAT"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) # Filter all "N"s 

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))

# 5. Remove primers using cutadapt
cutadapt <- "/change to wherever cutadapt is located/cutadapt" 
system2(cutadapt, args = "--version")

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
        "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
        fnFs.filtN[i], fnRs.filtN[i])) # input files
}
# sanity check!
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# 6. Check Complexity and Quality ####
plot(seqComplexity(getSequences(fnFs.cut[1])))
plot(seqComplexity(getSequences(fnRs.cut[1])))

plotQualityProfile(fnFs.cut[1:4])
plotQualityProfile(fnRs.cut[1:4])

#7. Filter/Trim reads based on the quality metrics
filt_path <- file.path(path, "filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Filter
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs.cut[i], fnRs.cut[i]), c(filtFs[i], filtRs[i]),
                    truncLen=c(240,160), 
                    maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, verbose=TRUE)
}

#8. Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#9. Learn the error rates
dadaFs.lrn <- dada(derepFs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errF <- dadaFs.lrn[[1]]$err_out
dadaRs.lrn <- dada(derepRs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errR <- dadaRs.lrn[[1]]$err_out

dadaFs.lrn[1]
dadaRs.lrn[1]

plotErrors(dadaFs.lrn[[1]], nominalQ=TRUE)

#10. Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]

#11. Merge Forward/Reverse Reads ####
merged <- mergePairs(dadaFs.lrn, derepFs, dadaRs.lrn, derepRs, verbose = TRUE)
# Inspect the merger data.frame from the first sample
head(merged[[1]])

#12. Constructing the sequence table - sequence counting table
seqtab <- makeSequenceTable(merged)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#Sequences that are much longer or shorter than expected may be the result of non-specific priming, 
#and may be worth removing
seqtab1 <- seqtab[,nchar(colnames(seqtab)) %in% seq(250, 256)]

#13. Remove chimeric sequences ####
seqtab.nochim <- removeBimeraDenovo(seqtab1, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab1)

#13-2. How many Seqs were lost at each step? ####
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(dadaFs.lrn, getN), sapply(merged, getN), rowSums(seqtab1), rowSums(seqtab.nochim))
colnames(track) <- c("denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
track
write.table(track, "Processed-16S-Sequencing-Counts.txt")

##14. Assign taxonomy ####
#Below lines should be run uner terminal
#wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz?download=1 -O silva_nr99_v138.1_train_set.fa.gz
#wget https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz?download=1 -O silva_species_assignment_v138.1.fa.gz

set.seed(128)
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", multithread = TRUE)
taxa.species <- addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz", verbose=TRUE)

##15.Combine everything as one phyloseq object and save it in your computer

# Install phyloseq package
BiocManager::install("phyloseq")
library(phyloseq)

# Making phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), tax_table(taxa.species))
ps

# Saving files

# Save an object to a file
saveRDS(ps, file = "phyloseq.rds")
# Restore the object
readRDS(file = "phyloseq.rds")

# Save entire environment to a file
save.image("dada2.Rdata")
load(file = "dada2.Rdata")
