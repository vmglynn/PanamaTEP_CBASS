
library("knitr")
library(ggplot2)
library(gridExtra)
library(phyloseq)
library(dada2)
library(DECIPHER)
library(phangorn)
library(tidyverse)
library(Biostrings)
library(ShortRead)


knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
set.seed(100)




## -----------------------------------------------------------------------------
#Path will change for each run
path <- "/home/vmglynn/projects/def-barrett/vmglynn/CBASS/16S/cutadapt"
list.files(path)






## -----------------------------------------------------------------------------

# Sort ensures forward/reverse reads are in same order; these are cutadapt files already
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)





plotQualityProfile(fnFs[4:5])

#Trim at 200, as ITS-2 not shorter than ~200 bp or longer than ~300 bp length (Hume et al. 2018)  and reads are exceptionally good quality. DADA2 is quite stringent, so better keep more than less reads at this step.

#Yet, as the ITS2 primer sets are not the same length, will need to trim F and R reads differently to a fixed length.  




plotQualityProfile(fnRs[4:5])

#Trim at 160, as ITS-2 not shorter than ~200 bp or longer than ~300 bp length (Hume et al. 2018). DADA2 is quite stringent, so better keep more than less reads at this step.

## The reverse reads are of significantly worse quality, especially at the end, which is common in Illumina sequencing. This isn't too worrisome, as DADA2 incorporates quality information into its error model which makes the algorithm robust to lower quality sequence, but trimming as the average qualities crash will improve the algorithm's sensitivity to rare sequence variants. Based on these profiles, we will truncate the reverse reads at position 160 where the quality distribution crashes.



## ----filter and trim----------------------------------------------------------

# Place filtered files in filtered/ subdirectory

filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))

any(duplicated(c(fnFs, fnRs)))
any(duplicated(c(filtFs, filtRs)))

length(fnFs)
length(fnRs)




out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)







## ----estimate error-----------------------------------------------------------

#Estimate error from first 40 samples
## The learnErrors method learns from a parametric error model from the data, by alternating estimation of the error rates and inference of sample composition until they converge on a jointly consistent solution. 
## As in many machine-learning problems, the algorithm must begin with an initial guess, for which the maximum possible error rates in this data are used (the error rates if only the most abundant sequence is correct and all the rest are errors).

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)





## ----plot errors--------------------------------------------------------------

# Plot estimated error rates, sanity check

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

## nominalQ, if TRUE, plot the expected error rates (red line) if quality scores exactly matched their nominal definition: Q = -10 log10(p_err).

## The black line shows the estimated error rates after convergence of the machine-learning algorithm. The red line shows the error rates expected under the nominal definition of the Q-score. 




## ----dereplication step-------------------------------------------------------

## Dereplication combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance” equal to the number of reads with that unique sequence. Dereplication substantially reduces computation time by eliminating redundant comparisons.

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

## verbose=TRUE, results in the service generating more output (will show you both WARNING and INFO log levels), normally you will only see WARNING or higher (ERROR for example).

# Name the derep-class objects by the sample names

names(derepFs) <- sample.names
names(derepRs) <- sample.names







# Sample inference with pseudo-pooling

dadaFs <- dada(derepFs, errF, pool="pseudo", multithread = TRUE)
dadaRs <- dada(derepRs, errR, pool="pseudo", multithread = TRUE)

dadaFs[[1]]
dadaRs[[1]]

## Pseudo-pooling is where samples are processed independently after sharing information between samples, approximating pooled sample inference in linear time.

## To pool or not to pool? Rarity biological-relevance vs PCR-artifact.

## OMEGA_A:  The key sensitivity parameters, controls the p-value threshold at which to call new ASVs.

## OMEGA_C: The error-correction threshold. One alternative is to turn off error-correction.







## ----merge paired reads-------------------------------------------------------
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# Inspect merger data.frame. from first sample

head(mergers[[1]])
mergers[[1]]

## Merging is performed by aligning the denoised forward reads with the reverse-complement of the corresponding denoised reverse reads, and then constructing the merged “contig” sequences. By default, merged sequences are only output if the forward and reverse reads overlap by at least 12 bases, and are identical to each other in the overlap region.

## Most of your reads should successfully merge. If that is not the case upstream parameters may need to be revisited: Did you trim away the overlap between your reads?





## ----sequence table-----------------------------------------------------------

#Construct Sequence Table 
seqtab <- makeSequenceTable(mergers)

#Explore size of seqs
dim(seqtab)

## Inspect distribution of sequence lengths

table(nchar(getSequences(seqtab)))

hist(nchar(getSequences(seqtab)))

hist(nchar(getSequences(seqtab)), 
     main="Length distribution of seqs", 
     xlab="ITS-2 Size", 
     xlim=c(200,450),
     las=1, 
     breaks=50)




## ----remove chimeras----------------------------------------------------------

## Chimeric sequences are identified if they can be exactly reconstructed by combining a left-segment and a right-segment from two more abundant “parent” sequences.

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)


## ----seq length filter--------------------------------------------------------

# Sequence table length filter

MINLEN <- 150
MAXLEN <- 350

# MIN/MAXLEN can change based on histogram plot 

seqlens <- nchar(getSequences(seqtab.nochim)) 

st_len_filt <- seqtab.nochim[,seqlens >= MINLEN & seqlens <= MAXLEN]

## st_len_filt <- seqtab.nochim(seqlens = 150) 

## st_len_filt <- seqtab.nochim(seqlens = 350) 

hist(nchar(getSequences(st_len_filt)), 
     main="Length Distribution of seqs", 
     xlab="ITS-2 Size", 
     xlim=c(150,350),
     las=1, 
     breaks=50)




## ----seq abundance filter-----------------------------------------------------

# Sequence table abundance filter

MINABUND <- 1 # set to 1 perhaps, can change
abundances <- colSums(st_len_filt)
st_len_abu_filt <- st_len_filt[,abundances >= MINABUND]

hist(nchar(getSequences(st_len_abu_filt)), 
     main="Abundance distribution of seqs", 
     xlab="Abundance", 
     xlim=c(150, 400),
     las=1, 
     breaks=50)




## ----further exploration------------------------------------------------------

seqtab.nochim2 <- as_tibble(seqtab.nochim)
glimpse(seqtab.nochim2)
str(seqtab.nochim2)

#Return sequences with col 1 with no heading but number and column 2 heading = x

seqs.nochim <- getSequences(seqtab.nochim)
head(seqs.nochim)
str(seqs.nochim)

## Compactly display the internal structure of an R object, a diagnostic function and an alternative to summary (and to some extent, dput). Ideally, only one line for each ‘basic’ structure is displayed. It is especially well suited to compactly display the (abbreviated) contents of (possibly nested) lists; str(object, …) = any R object about which you want to have some information.

## head() = return first part of vector, matrix, table, df, or function. 

## remove: colnames(seqs.nochim) <- c("sequence")
### Received "Error in `colnames<-`(`*tmp*`, value = "sequence") : attempt to set 'colnames' on an object with less than two dimensions"



# Write seq table

write.csv(seqtab.nochim, file = "seqtab.nochim.csv")
write.csv(st_len_filt, file = "st_len_filt.csv")
write.csv(st_len_abu_filt, file = "st_len_abu_filt.csv")


## ----compare filter steps 1---------------------------------------------------
hist(rowSums(seqtab.nochim))
hist(rowSums(st_len_filt))
hist(rowSums(st_len_abu_filt))

row_seqtab = tibble(rowSums(seqtab.nochim))
row_slaf = tibble(rowSums(st_len_abu_filt))

qplot(row_seqtab$`rowSums(seqtab.nochim)`, geom="histogram", binwidth=5000)
qplot(row_slaf$`rowSums(st_len_abu_filt)`, geom="histogram", binwidth=5000)

hist(colSums(st_len_filt))
hist(colSums(st_len_abu_filt))
## REMINDER: abundances <- colSums(st_len_filt) | st_len_abu_filt <- st_len_filt[,abundances >= MINABUND]



## ----compare filter steps 2 -- not working------------------------------------
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


## ----tibbles tracking reads through pipeline----------------------------------

track2 <- as_tibble(track)
glimpse(rownames(track2))

track2 %>% 
  summarize_all(median)

# Write track
write.csv(track, file = "track.csv")





## ----seqs to FASTA------------------------------------------------------------
getSequences(seqtab.nochim) 

#export seqs as fasta
uniquesToFasta(getUniques(seqtab.nochim), fout="./uniqueSeqs_ITS_14June2022.fasta", ids=paste0("Seq", seq(length(getUniques(seqtab.nochim)))))

## confirmed primers completed removed! 

## Original code: uniquesToFasta(getUniques(seqtab.nochim), fout="~/Dropbox/1_symITSdata/uniqueSeqs.fasta", ids=paste0("Seq", seq(length(getUniques(seqtab.nochim))))) 



getSequences(seqtab.nochim) 

#export seqs as fasta
uniquesToFasta(getUniques(seqtab.nochim), fout="./uniqueSeqs_16S_CBASS_18Oct2022.fasta", ids=paste0("Seq", seq(length(getUniques(seqtab.nochim)))))

## confirmed primers completed removed! 

## Original code: uniquesToFasta(getUniques(seqtab.nochim), fout="~/Dropbox/1_symITSdata/uniqueSeqs.fasta", ids=paste0("Seq", seq(length(getUniques(seqtab.nochim))))) 




taxa_16S <- assignTaxonomy(seqtab.nochim, "/home/vmglynn/projects/def-barrett/vmglynn/CBASS/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

taxa_16S <- addSpecies(taxa_16S, "/home/vmglynn/projects/def-barrett/vmglynn/CBASS/silva_species_assignment_v138.1.fa.gz")

##Check taxa db look okay

taxa.print <- taxa_16S # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

tax_silva_16s  = tax_table(taxa_16S)
write.csv(tax_silva_16s, file = "tax_silva_16s_CBASS_18Oct2022.csv")

tax_silva_otu  = otu_table(seqtab.nochim, taxa_are_rows = FALSE)
write.csv(tax_silva_otu, file = "tax_silva_otu_16S_CBASS_18Oct2022.csv")




## ----decipher phylogenetic tree-----------------------------------------------
# alignment using Decipher package note: may need to load
# This is for the length and abundance filtered samples.
##seqs_filt <- getSequences(st_len_abu_filt)
##names(seqs_filt) <- seqs_filt # This propagates to the tip labels of the tree
##alignment <- AlignSeqs(DNAStringSet(seqs_filt), anchor=NA)

## AlignSeqs : Align A Set Of Unaligned Sequences - Performs profile-to-profile alignment of multiple unaligned sequences following a guide tree.





## -----------------------------------------------------------------------------
## uniquesToFasta(getUniques(st_len_abu_filt),fout="./uniqueSeqs_16S_CBASS_19Sep2022.fasta")

#export seqs as fasta
uniquesToFasta(getUniques(st_len_abu_filt),fout= "./uniqueSeqs_asv_16s_CBASS_18Oct2022.fasta", ids=paste0("asv", seq(length(getUniques(st_len_abu_filt)))))


## -----------------------------------------------------------------------------
library(phyloseq); packageVersion("phyloseq")
library(ggplot2); packageVersion("ggplot2")
theme_set(theme_bw())


## ----handoff to Phyloseq------------------------------------------------------

#Import sample_data

samdf <- read.csv("/home/vmglynn/projects/def-barrett/vmglynn/CBASS/Sampledata_16S_CBASS_17Oct2022_MLG.csv", row.names = 1,  sep = ",")


### merge into ps object, save as Rmd files

ps <- phyloseq(otu_table(tax_silva_otu, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(tax_silva_16s))

saveRDS(ps, file = "symITSps_16S_CBASS_18Oct2022.rds")

readRDS("symITSps_16S_CBASS_18Oct2022.rds")