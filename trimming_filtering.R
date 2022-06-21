# code to trim, merge, filter and make an OTU table.

# DADA2 for trimming off primers and merging the reads

library(tidyr)
# library(microseq)
# library(R.utils)
library(ShortRead)
library(dada2)
library(Biostrings)
library(stats)
library(devtools)
# devtools::install_github("shaunpwilkinson/insect")
library(insect)
#FW: 23
#REV: 33
#trimLeft = c(FWD_PRIMER_LEN, REV_PRIMER_LEN)

path <- file.path("miseq_data/hotel_pilot_1")
# countFastq(path)

# list of files
list.files(path)
fnFs <- sort(list.files(path, pattern = "R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2_001.fastq", full.names = TRUE))
fnFs
fnRs

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
# filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)
filterAndTrim(fwd = fnFs, filt = fnFs.filtN, 
              rev = fnRs, filt.rev = fnRs.filtN, 
              trimLeft = c(23, 33))
# but they are not readable/legible/viewabale on command line with less, mojibake ??
path_cut <- file.path("miseq_data/hotel_pilot_1/filtN")
list.files(path_cut)

fnFs <- sort(list.files(path_cut, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path_cut, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names
# plotQualityProfile(fnFs[1:2]) # no specific quality drop, trim 10 from end
# plotQualityProfile(fnRs[1:2]) # no specific quality drop, trim 10 from end
# trimRight = 10

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
filtFs
filtRs

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(23, 33), trimRight = 10,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, # remove PhiX reads too
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

errF <- learnErrors(filtFs)
# 109054285 total bases in 1156718 reads from 10 samples will be used for learning the error rates.
plotErrors(errF, nominalQ=TRUE)

errR <- learnErrors(filtRs)
# 101349979 total bases in 1354312 reads from 12 samples will be used for learning the error rates.
plotErrors(errR, nominalQ=TRUE)

# sample inference
dadaFs <- dada(filtFs, err=errF)
# Sample 1 - 187475 reads in 8511 unique sequences.
# Sample 2 - 202196 reads in 7767 unique sequences.
# Sample 3 - 140741 reads in 3084 unique sequences.
# Sample 4 - 1775 reads in 464 unique sequences.
# Sample 5 - 131490 reads in 3220 unique sequences.
# Sample 6 - 86314 reads in 2647 unique sequences.
# Sample 7 - 93308 reads in 1427 unique sequences.
# Sample 8 - 130820 reads in 3648 unique sequences.
# Sample 9 - 1142 reads in 372 unique sequences.
# Sample 10 - 181457 reads in 5534 unique sequences.
# Sample 11 - 87291 reads in 2560 unique sequences.
# Sample 12 - 110303 reads in 3456 unique sequences.
# Sample 13 - 2080 reads in 350 unique sequences.
# Sample 14 - 76370 reads in 1835 unique sequences.
# Sample 15 - 55486 reads in 11130 unique sequences.
dadaRs <- dada(filtRs, err=errR)
# Sample 1 - 187475 reads in 28534 unique sequences.
# Sample 2 - 202196 reads in 14157 unique sequences.
# Sample 3 - 140741 reads in 6042 unique sequences.
# Sample 4 - 1775 reads in 900 unique sequences.
# Sample 5 - 131490 reads in 5472 unique sequences.
# Sample 6 - 86314 reads in 5067 unique sequences.
# Sample 7 - 93308 reads in 5380 unique sequences.
# Sample 8 - 130820 reads in 7761 unique sequences.
# Sample 9 - 1142 reads in 888 unique sequences.
# Sample 10 - 181457 reads in 11022 unique sequences.
# Sample 11 - 87291 reads in 5826 unique sequences.
# Sample 12 - 110303 reads in 13019 unique sequences.
# Sample 13 - 2080 reads in 767 unique sequences.
# Sample 14 - 76370 reads in 10884 unique sequences.
# Sample 15 - 55486 reads in 22215 unique sequences.

# pool = TRUE?

dadaFs[[3]]

# Merge the paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

# Construct the ASV table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
write.csv(seqtab, "seqtab.csv")
# head(seqtab)
# [1]  15 250

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# 98 100 106 114 115 116 117 118 119 120 121 122 126 
# 1   1   1   1  21  37  53  94  22  10   7   1   1 

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# Identified 125 bimeras out of 250 input sequences.
dim(seqtab.nochim)
# [1]  15 125
sum(seqtab.nochim)/sum(seqtab)
# [1] 0.9885216 # less than 2% of the merged reads are chimeras
head(seqtab.nochim)
write.csv(seqtab.nochim, "seqtab_nochim.csv")

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, "filtering_results.csv")

# Assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "~/tax/silva_nr_v132_train_set.fa.gz", multithread=TRUE)

'''
assignTaxonomy(...) expects a training fasta file (or compressed fasta file) in which the taxonomy corresponding to each sequence is encoded in the id line in the following fashion (the second sequence is not assigned down to level 6):
  
  >Level1;Level2;Level3;Level4;Level5;Level6;
ACCTAGAAAGTCGTAGATCGAAGTTGAAGCATCGCCCGATGATCGTCTGAAGCTGTAGCATGAGTCGATTTTCACATTCAGGGATACCATAGGATAC
>Level1;Level2;Level3;Level4;Level5;
CGCTAGAAAGTCGTAGAAGGCTCGGAGGTTTGAAGCATCGCCCGATGGGATCTCGTTGCTGTAGCATGAGTACGGACATTCAGGGATCATAGGATAC
assignSpecies(...) expects the training data to be provided in the form of a fasta file (or compressed fasta file), with the id line formatted as follows:
  
  >ID Genus species
ACCTAGAAAGTCGTAGATCGAAGTTGAAGCATCGCCCGATGATCGTCTGAAGCTGTAGCATGAGTCGATTTTCACATTCAGGGATACCATAGGATAC
>ID Genus species
CGCTAGAAAGTCGTAGAAGGCTCGGAGGTTTGAAGCATCGCCCGATGGGATCTCGTTGCTGTAGCATGAGTACGGACATTCAGGGATCATAGGATAC
'''
# fish 12S database compatible with dada2? 







