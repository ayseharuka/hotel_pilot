
# Training the classifier and assigning taxonomy

library(dplyr)
library(tidyr)
# library(microseq)
# library(R.utils)
library(ShortRead)
library(dada2)
library(Biostrings)
library(stats)
library(insect)

# Load the dataset/ASV table (chimeras removed)
# test <- data(samoa)
seqtab.nochim <- read.csv("seqtab_nochim.csv", header = TRUE, row.names = 1)
seqtab.nochim
## read in the example seqtab.nochim ASV table
# data(seqtab.nochim)
## get sequences from table column names
x <- char2dna(colnames(seqtab.nochim))
x
## name the sequences sequentially
names(x) <- paste0("ASV", seq_along(x))
## optionally remove column names that can flood the console when printed
colnames(seqtab.nochim) <- NULL 
head(seqtab.nochim)

# Load the classifier
classifier <- readRDS("mifish_classifier.rds")
classifier
names(attributes(classifier))

longDF <- classify(x, classifier, threshold = 0.7)
longDF <- cbind(longDF, t(seqtab.nochim))
View(longDF)

