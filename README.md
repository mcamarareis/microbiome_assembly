# Pipeline used at the work: *"Exploring the phycosphere of Emiliania huxleyi: from bloom dynamics to microbiome assembly experiments"*

# Table of contents
* [Introduction](https://github.com/mcamarareis/microbiome_assembly/#Introduction)
* [Dataset characteristics](https://github.com/mcamarareis/microbiome_assembly/#Dataset-characteristics)
* [DADA2 pipeline](https://github.com/mcamarareis/microbiome_assembly/#DADA2-pipeline)
* [Demultiplexing](https://github.com/mcamarareis/microbiome_assembly/#Demultiplexing)
* [Primer trimming](https://github.com/mcamarareis/microbiome_assembly/#Primer-trimming)
* [DADA2: Filter and trim](https://github.com/mcamarareis/microbiome_assembly/#DADA2:-Filter-and-trim)
* [DADA2 learn errors and infer ASVs](https://github.com/mcamarareis/microbiome_assembly/#DADA2-learn-errors-and-infer-ASVs)
* [Merging all runs and running chimera detection](https://github.com/mcamarareis/microbiome_assembly/#Merging-all-runs-and-running-chimera-detection)
* [Taxonomic assignation IDTAXA and VSEARCH](https://github.com/mcamarareis/microbiome_assembly/#Taxonomic-assignation-IDTAXA-and-VSEARCH)
* [Filtering the table, comparing sequencing runs and generating final ASV table for analysis
](https://github.com/mcamarareis/microbiome_assembly/#Filtering-the-table,-comparing-sequencing-runs-and-generating-final-ASV-table-for-analysis)


## Introduction

This page details how the amplicon sequence variants (ASVs) table used in Câmara dos Reis et al., 2022 was built and downstream analysis. 
Briefly raw reads were demultiplexed and primers were trimmed using CUTADAPT.
These reads were analyzed using DADA2 to obtain the ASVs table. Downstream analysis were performed in the R environment.

![image](https://user-images.githubusercontent.com/54186891/153705647-5c934b14-d6d4-4202-a2ef-155a04f7d80e.png)

**Figure 1**: Bioinformatic processing to obtain the ASV table for microbiome analyses.


## Dataset characteristics

In total 126 samples were sequenced here.
- 30 were environmental samples from the TARA Breizh Bloom cruise 2019
- 96 were from the assembly experiment

This dataset was obtained from 3 library pools which were sequenced in three independent sequencing runs.
A particularity is that two of these libraries were sequenced together in two independent sequencing runs (consisting of technical replicates).
The third library was sequenced in another sequencing run without technical replicate.

The number of samples as well as the name of each sequencing run are summarized in the table:

![image](https://user-images.githubusercontent.com/54186891/153750517-e6601cf6-2a04-43d0-b0e3-ea083a0ed96d.png)


## DADA2 pipeline

### Demultiplexing

Demultiplexing and primer trimming have been performed using Cutadapt in separate steps.
Fastq files processed here, contain reads are in mixed orientation. 
So Cutadapt was used twice, first looking for tags in R1 and then looking for tags in the R2.

To do so, a text file in which are listed the samples and corresponding tags was used.

![image](https://user-images.githubusercontent.com/54186891/152359496-b25234b1-054d-4f65-adc9-45d70ad70c15.png)

Then the following script is applied to each line of this file. For each sample:

The results of demultiplexing are saved into ``demultiplexed/${SAMPLE}_${RUNID}_1_R1.fastq.gz``, ``demultiplexed/${SAMPLE}_${RUNID}_2_R1.fastq.gz``, ``demultiplexed/${SAMPLE}_${RUNID}_1_R2.fastq.gz`` and ``demultiplexed/${SAMPLE}_${RUNID}_2_R2.fastq.gz``. 

Notice that now the reads containing forward primers are the ones finished by ``R1``, those extracted from the first trimming run are named ``1_R1`` and those extracted from the second are named ``2_R1``. The reads containing reverse primers files are named accordingly to the forward reads (``1_R2`` and ``2_R2``).


```
#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -wd workdirectory
#$ -q long.q
#$ -N demultiplexing_200207
#$ -pe thread 16

N=16

export TMPDIR=scratch

source $CONDA3/activate cutadapt-2.8

#Demultiplexing run "200207_SN6662_A_L001_ADNU-41"

RUNID="200207_SN6662_A_L001_ADNU-41"
LOG="${RUNID}_demultiplexing.log"
FORWARD="200207_SN6662_A_L001_ADNU-41_AdapterTrimmed_R1.fastq.gz"
REVERSE="200207_SN6662_A_L001_ADNU-41_AdapterTrimmed_R2.fastq.gz"

echo 'sample in_reads out_reads w/adapters w/adapters2' > ${LOG}

demutrim() {
    TMP_TAG_F1=$(mktemp --tmpdir="/scratch")
    TMP_TAG_R1=$(mktemp --tmpdir="/scratch")
    TMP_TAG_F2=$(mktemp --tmpdir="/scratch")
    TMP_TAG_R2=$(mktemp --tmpdir="/scratch")
    TMP_LOG=$(mktemp --tmpdir="/scratch")
    
    cutadapt -g "^NNNN${TAG_SEQ}" --report=minimal --no-indels --discard-untrimmed -o ${TMP_TAG_F1} -p ${TMP_TAG_R1} ${FORWARD} ${REVERSE} 1> ${TMP_LOG}
    
    cutadapt -g "^NNNN${TAG_SEQ}" --report=minimal --no-indels --discard-untrimmed -o ${TMP_TAG_F2} -p ${TMP_TAG_R2} ${REVERSE} ${FORWARD} 1>> ${TMP_LOG}
   
    gawk -v a="${SPLE_NAME}" '{
        if(NR == 2)
            for (i = 1; i <= NF; i++)
                var[i] = $i
        else if(NR == 4)
            print a,$2,$7+var[7],$8+var[11],$11+var[8]
    }' ${TMP_LOG} >> ${LOG}
    
    cat ${TMP_TAG_F1} | gzip > "path/demultiplexed/${SPLE_NAME}-${RUNID}_1_R1.fastq.gz"
    cat ${TMP_TAG_R1} | gzip > "path/demultiplexed/${SPLE_NAME}-${RUNID}_1_R2.fastq.gz"
    cat ${TMP_TAG_F2} | gzip > "path/demultiplexed/${SPLE_NAME}-${RUNID}_2_R1.fastq.gz"
    cat ${TMP_TAG_R2} | gzip > "path/demultiplexed/${SPLE_NAME}-${RUNID}_2_R2.fastq.gz"
    
    rm -f "${TMP_TAG_F1}" "${TMP_TAG_F2}" "${TMP_TAG_R1}" "${TMP_TAG_R2}" "${TMP_LOG}"
}



while read SPLE_NAME TAG_SEQ; do
    ((i=i%N)); ((i++==0)) && wait
    demutrim &
done < "tags_POOL.txt"

wait
```

REMARK:

The demultiplexing run for submiting the reads to public databases were performed as mentioned before. However, all reads with the tag were sent to the ``R1`` file and all reads without the tags were sent to the ``R2`` file.
The results of demultiplexing alone are saved into ``${SAMPLE}_${RUNID}_R1.fastq.gz and ${SAMPLE}_${RUNID}_R2.fastq.gz``. **The reads are not reorientated. For further analysis primer trimming should be used to separate reads with forward primer from reads with the reverse primer**.



### Primer trimming

Now that the reads are in their good files, the primer trimming will act removing forward primers from ``R1`` reads and reverse primers from ``R2`` reads.

Here CUTADAPT was also used.

In general the same parameters were used for all the sequencing runs, except for the run 201218, for which we removed the flag *--no_indels* to allow more reads to be kept


```
#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -wd /path/demultiplexed/
#$ -q long.q
#$ -pe thread 16

N=16

export TMPDIR=/path/demultiplexed/
  
source $CONDA3/activate cutadapt-2.8

PRIMER_F="^GTGYCAGCMGCCGCGGTAA"
PRIMER_R="^CCGYCAATTYMTTTRAGTTT"


#####For run 1_1
RUNID="200207_SN6662_A_L001_ADNU-41"
LOG="primer_trim-${RUNID}_1.log"
MANIFEST="manifest_1.txt"


ls /path/folder_demultiplexed_files/ | sort | grep _1_R[12].fastq.gz | \
gawk '{
for (i =1; i <= NF; i++){
sample = $0
sub(/_200207_SN6662_A_L001_ADNU-41_1_R[12].+$/,"",sample)
if (/_R1\.fastq.gz/)
forward[sample]=$0
if (/_R2\.fastq.gz/)
reverse[sample]=$0
}
}
END {
for (i in forward){
print i " " "/path/demultiplexed/" forward[i] " " "/path/demultiplexed/" reverse[i]
}
}' > ${MANIFEST}

echo 'sample in_reads out_reads w/adapters w/adapters2' > ${LOG}

trimonly() {
  TMP_CUT_F1=$(mktemp --tmpdir="/scratch")
  TMP_CUT_R1=$(mktemp --tmpdir="/scratch")
  TMP_LOG=$(mktemp --tmpdir="/scratch")
  
  cutadapt -g "${PRIMER_F}" -G "${PRIMER_R}" --report=minimal --discard-untrimmed --minimum-length 75 --no-indels -o ${TMP_CUT_F1} -p ${TMP_CUT_R1} ${FORWARD} ${REVERSE} 1>> ${TMP_LOG}
  
  gawk -v a="${OUTPUT}" '{
  if(NR == 1)
  for (i = 1; i <= NF; i++)
  var[i] = $i
  else if(NR == 2)
  print a,$2,$7+var[7],$8+var[11],$11+var[8]
}' ${TMP_LOG} >> ${LOG}

  cat ${TMP_CUT_F1}  | gzip > "primers_trimmed/${OUTPUT}_1_R1_trimmed.fastq.gz"
  cat ${TMP_CUT_R1}  | gzip > "primers_trimmed/${OUTPUT}_1_R2_trimmed.fastq.gz"
  
  rm -f "${TMP_CUT_F1}" "${TMP_CUT_R1}" "${TMP_LOG}"
  }
while read OUTPUT FORWARD REVERSE; do
((i=i%N)); ((i++==0)) && wait
trimonly &    
  
  done < "${MANIFEST}"

#####For reads that came from second demultiplexing run

RUNID="200207_SN6662_A_L001_ADNU-41"
LOG="primer_trim-${RUNID}_2.log"
MANIFEST="manifest_2.txt"
echo 'sample in_reads out_reads w/adapters w/adapters2' > ${LOG}

ls /path/folder_demultiplexed_files/ | sort | grep _2_R[12].fastq.gz | \
gawk '{
for (i =1; i <= NF; i++){
sample = $0
sub(/_200207_SN6662_A_L001_ADNU-41_2_R[12].+$/,"",sample)
if (/_R1\.fastq.gz/)
forward[sample]=$0
if (/_R2\.fastq.gz/)
reverse[sample]=$0
}
}
END {
for (i in forward){
print i " " "/path/folder_demultiplexed_files/" forward[i] " " "/path/folder_demultiplexed_files/" reverse[i]
}
}' > ${MANIFEST}

trimonly() {
  TMP_CUT_F2=$(mktemp --tmpdir="/path/folder_demultiplexed_files/")
  TMP_CUT_R2=$(mktemp --tmpdir="/path/folder_demultiplexed_files/")
  TMP_LOG=$(mktemp --tmpdir="/path/folder_demultiplexed_files/")
  
  cutadapt -g "${PRIMER_F}" -G "${PRIMER_R}" --report=minimal --discard-untrimmed --minimum-length 75 --no-indels -o ${TMP_CUT_F2} -p ${TMP_CUT_R2} ${FORWARD} ${REVERSE} 1>> ${TMP_LOG}
  
  gawk -v a="${OUTPUT}" '{
  if(NR == 1)
  for (i = 1; i <= NF; i++)
  var[i] = $i
  else if(NR == 2)
  print a,$2,$7+var[7],$8+var[11],$11+var[8]
}' ${TMP_LOG} >> ${LOG}
    
  cat ${TMP_CUT_F2}  | gzip > "primers_trimmed/${OUTPUT}_2_R1_trimmed.fastq.gz"
  cat ${TMP_CUT_R2}  | gzip > "primers_trimmed/${OUTPUT}_2_R2_trimmed.fastq.gz"
  
  rm -f "${TMP_CUT_F2}" "${TMP_CUT_R2}" "${TMP_LOG}"
  }

while read OUTPUT FORWARD REVERSE; do
((i=i%N)); ((i++==0)) && wait
trimonly &    
  done < "${MANIFEST}"

wait
```

## DADA2: Filter and trim

For all the first steps of DADA2 pipeline reads from different sequencing runs and reads that were mixed-oriented were kept independent since they can have different error models. 
The function ``filterAndTrim`` was run using default parameters and truncation lenghs of ``truncLen=c(215, 190)`` for reads which came from the first demultiplexing run (``1_R1`` and ``1_R2``) and the parameters were inversed ``truncLen=c(190, 215)`` for the reads coming from the second demultiplexing run.

This is necessary because reads from the original ``R2`` had lower quality (``1_R2`` and ``2_R1``).


```
library(dada2); packageVersion("dada2")

# Input fastq directory
path <-"/primers_trimmed/"

TRUNCLEN_1=215
TRUNCLEN_2=190

#### Reads from first demultiplexing run #####

# List fastqs 

fastqFs1 <- sort(list.files(path, pattern="_1_R1_trimmed.fastq.gz",  full.names = TRUE))
fastqRs2 <- sort(list.files(path, pattern="_1_R2_trimmed.fastq.gz",  full.names = TRUE))

# Filtered fastqs name

filtFs1 <- gsub("primers_trimmed","quality_filtered",fastqFs1)
filtRs2 <- gsub("primers_trimmed","quality_filtered",fastqRs2)


# Filtering

out <- filterAndTrim(fastqFs1, filtFs1, fastqRs2, filtRs2, truncLen=c(TRUNCLEN_1,TRUNCLEN_2),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)


###### Reads coming from second demultiplexing run ######

# List fastqs 
fastqFs2 <- sort(list.files(path, pattern="_2_R1_trimmed.fastq.gz",  full.names = TRUE))
fastqRs1 <- sort(list.files(path, pattern="_2_R2_trimmed.fastq.gz",  full.names = TRUE))

# Filtered fastqs name
filtFs2 <- gsub("primers_trimmed","quality_filtered",fastqFs2)
filtRs1 <- gsub("primers_trimmed","quality_filtered",fastqRs1)


# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
out_2 <- filterAndTrim(fastqFs2, filtFs2, fastqRs1, filtRs1, truncLen=c(TRUNCLEN_2,TRUNCLEN_1),
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)
out
out_2
```

This command was used for each sequencing run.

## DADA2 learn errors and infer ASVs 

To lear the errors I increased the number of bases to include more samples to ``nbases = 2e+08`` and also used the flag ``randomize=TRUE``.
The DADA command was run in pooled mode ``pool=TRUE``

It's important to notice that reads obtained from the second demultiplexing run ``(_2)`` need to be merged in opposite position.

``mergers_2 <- mergePairs(dadaRs_2, filtRs_2, dadaFs_2, filtFs_2, verbose=TRUE) ##here I invert them``


```
library(dada2); packageVersion("dada2")

# File parsing RUNID reads from first demultiplexing run (_1)

filtpath <- "/quality_filtered"

filtFs <- list.files(filtpath, pattern="_1_R1_trimmed.fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpath, pattern="_1_R2_trimmed.fastq.gz", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) 
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) 

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names

##Learning errors##
###################
set.seed(100)

errF <- learnErrors(filtFs, nbases = 2e+08, randomize=TRUE, verbose=TRUE)
errR <- learnErrors(filtRs, nbases = 2e+08, randomize=TRUE, verbose=TRUE)

##ASVs inference##
##################

dadaFs <- dada(filtFs, err=errF, pool=TRUE)
dadaRs <- dada(filtRs, err=errR, pool=TRUE)

####Merging pairs###
####################

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

##exporting all intermediary files
saveRDS(errF, "errF_runID_1.rds")
saveRDS(errR, "errR_runID_1.rds")
saveRDS(mergers,"mergers_runID_1.rds")
saveRDS(dadaFs, "ddF_runID_1.rds")
saveRDS(dadaRs, "ddR_runID_1.rds")

# Construct sequence table 
seqtab_runID_1 <- makeSequenceTable(mergers)
saveRDS(seqtab_runID, "seqtab_runID_1.rds")


# File parsing RUNID reads from first demultiplexing run (_2)

filtpath <- "/quality_filtered"

filtFs_2 <- list.files(filtpath, pattern="_2_R1_trimmed.fastq.gz", full.names = TRUE)
filtRs_2 <- list.files(filtpath, pattern="_2_R2_trimmed.fastq.gz", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs_2), "_"), `[`, 1) 
sample.namesR <- sapply(strsplit(basename(filtRs_2), "_"), `[`, 1) 

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs_2) <- sample.names
names(filtRs_2) <- sample.names

set.seed(100)

##Learning errors##
###################
errF_2 <- learnErrors(filtFs_2, nbases = 2e+08, randomize=TRUE, verbose=TRUE)
errR_2 <- learnErrors(filtRs_2, nbases = 2e+08, randomize=TRUE, verbose=TRUE)

##ASVs inference##
##################

dadaFs_2 <- dada(filtFs_2, err=errF_2, pool=TRUE)
dadaRs_2 <- dada(filtRs_2, err=errR_2, pool=TRUE)

####Merging pairs###
####################

mergers_2<- mergePairs(dadaRs_2, filtRs_2, dadaFs_2, filtFs_2, verbose=TRUE) ##here I invert them


##exporting all intermediary files
saveRDS(errF_2, "errF_runID_2.rds")
saveRDS(errR_2, "errR_runID_2.rds")
saveRDS(mergers_2,"mergers_runID_2.rds")
saveRDS(dadaFs_2, "ddF_runID_2.rds")
saveRDS(dadaRs_2, "ddR_runID_2.rds")

# Construct sequence table 
seqtab_runID_2 <- makeSequenceTable(mergers_2)
saveRDS(seqtab_runID, "seqtab_runID_2.rds")


```

Then a step of reverse-complementing the sequences from the second demultiplexing run ``_2`` is added so we can merge these sequencing tables for each sequencing run.

```
library(dada2); packageVersion("dada2")
library(Biostrings); packageVersion("Biostrings") 

##MERGING SEQSTABS FROM RUNID###

seqtab_runID_1 <- readRDS("seqtab_runID_1.rds")
seqtab_runID_2 <- readRDS("seqtab_runID_2.rds")

##Reverse-complement seqs from ``_2``.

reverse_complement <- function(x) {
  toString(reverseComplement(DNAString(x)))
}

reverse_complemented_seq <- lapply(colnames(seqtab_runID_2), reverse_complement)
colnames(seqtab_runID_2) <- reverse_complemented_seq

st.all_runID <- mergeSequenceTables(seqtab_runID_1, seqtab_runID_2, repeats="sum")
saveRDS(st.all_runID, "seqtab_runID_merged.rds")

```

## Merging all runs and running chimera detection

Here the command `mergeSequenceTables` is used to merge all the sequencing runs together for chimera detection.
It is recommended to merge them to control the chimera detection across sequencing runs.

It is important to notice that technical replicates will be kept separate for evaluation of consistency among sequencing runs.

```
library(dada2); packageVersion("dada2")
library("data.table"); packageVersion("data.table")

st.all_runID1 <- readRDS("seqtab_runID_merged.rds.rds")
st.all_runID2 <- readRDS("seqtab_runID2_merged.rds.rds")
st.all_runID3 <- readRDS("seqtab_runID3_merged.rds.rds")


seqtab_all_runs <- mergeSequenceTables(st.all_runID1, st.all_runID2, st.all_runID3)
saveRDS(seqtab_all_runs, "seqtab_allruns.rds")

#Remove chimeras

seqtab_nochimera <- removeBimeraDenovo(seqtab_all_runs, method="pooled",minFoldParentOverAbundance=8, verbose=TRUE, multithread=FALSE)

saveRDS(seqtab_nochimera, "seqtab_allruns_nochimera.rds")

```

## Taxonomic assignation IDTAXA and VSEARCH

Here the taxonomy is assigned using the trainingset of the Silva database V138: SILVA_SSU_r138_2019.RData for IDTAXA.
For vsearch the sequences in SILVA database were trimmed using the same primers of our study using an overlap of 2/3 of each primers length and default errors rate.

```

# First we give our seq headers more manageable names (ASV_1, ASV_2...)

asv_seqs <- colnames(seqtab_nochimera)
asv_headers <- vector(dim(seqtab_nochimera)[2], mode="character")

for (i in 1:dim(seqtab_nochimera)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs_sequences.fa")

# count table:
asv_tab <- t(seqtab_nochimera)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_table.tsv", sep="\t", quote=F, col.names=NA)

# tax table:

library(DECIPHER); packageVersion("DECIPHER")

dna <- DNAStringSet(getSequences(seqtab_nochimera)) # Create a DNAStringSet from the ASVs
load("/path_for_traininset/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET

set.seed(100)

ids <- IdTaxa(dna, trainingSet, strand="top", threshold=50, processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest

# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab_nochimera)
saveRDS(taxid, "taxaid_final_50.rds")
write.table(taxid,  "ASVs_taxonomy_IDTAXA_all_ASVs_50.tsv", sep="\t", quote=F, col.names=NA)

```

For vsearch 

```
#' # Assigining taxonomy vsearch

fas <- "path/ASVs_restricted_nochimera_final_coccobloom.fa"
tax_coccobloom <- "path/tax_vsearch_80.tsv"
reference_515F_926R <- "path/silva_v138/SILVA_138_SSURef_NR99_tax_silva_515F_926R.fasta"

system(paste0("/Users/marianacamaradosreis/Downloads/vsearch-2.15.0-macos-x86_64/bin/vsearch ",
              "--usearch_global ",  fas, 
              " --db ",              reference_515F_926R,
              " --id 0.8 ",
              " --blast6out ",          tax_coccobloom))
```

The taxonomy of best-hits obtained by vsearch was used to plot the taxonomic diversity in the environmental and experiment samples. 


## Filtering the table, comparing sequencing runs and generating final ASV table for analysis

Once the ASV table is ready the next step is to inspect the results and filter non bacterial sequences (Chloroplasts, mithocondria and not assigned sequences) and to filter the table by abundance. 
This filtering is to avoid working with low abundance ASVs, which can be spurious sequences and are not relevant in the context of this work.


```
seqtab <- readRDS("seqtab_allruns.rds") #the ASV table
seqtab_nochimera <- readRDS("seqtab_allruns_nochimera.rds") #the ASV table after removing chimera

sum(seqtab_nochimera)/sum(seqtab)

[1] 0.9833528

# 98% were kept after chimera removal

dim(seqtab_nochimera)/dim(seqtab)

# About 70% of the ASVs were removed as chimera


tapply(colSums(seqtab_nochimera), nchar(colnames(seqtab_nochimera)), sum)


    215     254     273     296     302     351     353     354     360     364     365     **366     367 
    195       7       2      37       3       3     382     444      30      36     241  128644    1396 
    368     369     370     371     372     373     374     375     376     377     378     384     386 
 148527 1210800   59585    7988  224914 3471548 2050039   25582    9815   15845       2       3       3 
    390     391     392 
      2       5       2 
  
#Most of the sequences are distributed between 366 and 377 

```

The reads out of this range are either assigned to chloroplasts and mithocondria, either not assigned to any bacterial group (using 50% identity threshold of VSEARCH), or low abundance ASVs (which would be removed from the table anyway).
So I proceeded to remove these ASVs.

```
seqtab.filt <- seqtab_nochimera[, (seqlens <= 377) & (seqlens >= 366)]

1-sum(seqtab.filt)/sum(seqtab_nochimera)

#This removed 0.01% of the reads

```

The next step is to filter the table to remove Chloroplasts, Mithocondria, non-classified sequences at domain level and ASVs accounting to less than 0.001% of the reads.

First I import the taxonomy and the metadata tables.


```
tax_data <- fread("path/ASVs_taxonomy_IDTAXA_restricted_allASVs.tsv") %>%
  replace(., is.na(.), "unclassified") 
  
#Here empty cells need to be filled otherwise the row is removed


metadata <- fread("path/metadata.csv")


#Getting metadata of samples for which I have technical replicates and of samples for which I don't

metadata_togroup <- metadata[time_course != "Day 392",]  #samples that I'll group

metadata_last_sampling <- metadata[time_course == "Day 392",] #samples which were sequenced only once

```


Chloroplasts and Mitochondria sequences are removed. For now I did not remove ASVs not classified at the domain level by IDTAXA. I'll check them latter.



```

#Filtering out chloroplasts, mitochondia
 

tax_filtered <- tax_data %>% 
  filter(., order != 'Chloroplast' & family != 'Mitochondria')
  
dt_filt_chloro <- dt.filt[rownames(filt) %in% tax_filtered$ASVId,] %>%
  t(.)

#Removing ASV 303 which is a chloroplast that was not classified by idtaxa

dt_filt_chloro <- dt_filt_chloro[,colnames(dt_filt_chloro) != "ASV_303"] 

```


Then I filter the ASVs with less than 0.001% of the total of reads, which is 61 reads

```

round(sum(dt_filt_chloro)*0.00001)
#61

dt_filt <- dt_filt_chloro[,which(colSums(dt_filt_chloro) >= 61)] %>%
  round(.)
min(colSums(dt_filt))
dim(dt_filt)
sum(dt_filt)

```

Consistency of the technical replicates (samples which were sequences twice) is checked by procrustes analysis.


```
metadata$sample <- gsub(metadata$sample , pattern = "\\-", replacement = "_") 
metadata$treatments <- gsub(metadata$treatments , pattern = "\\-", replacement = "_") 

print(metadata)
dim(metadata)

metadata_run1 <- metadata[run == "run1",] 
metadata_run2 <- metadata[run == "run2",] 
dim(metadata_run1)

run1 <- dt_filt[rownames(dt_filt) %in% metadata_run1$id,] #secting only samples sequenced in the first seuqencing run
run1 <- run1[,which(colSums(run1) > 0)] #removing ASVs with 0 reads

rownames(run1) <- metadata_run1$sample #defining rownames as the name of the samples withough run identification


run2 <- dt_filt[rownames(dt_filt) %in% metadata_run2$id,] #secting only samples sequenced in the first seuqencing run
run2 <- run2[,which(colSums(run2) >0)] #removing ASVs with 0 reads
rownames(run2) <- metadata_run2$sample  #defining rownames as the name of the samples withough run identification

all(rownames(run1) == rownames(run2)) #checking if rownames are the same in both runs (it's important for procrustes)


run1_hell = decostand(run1, method = "hellinger") #hellinger trasnformation of the abundances
run2_hell = decostand(run2, method = "hellinger") #hellinger trasnformation of the abundances

pca.run1=rda(run1) #runing the pca and saving in an object
pca.run2=rda(run2) #runing the pca and saving in an object

prot <- protest(pca.run1, pca.run2, permutations = 1000)

summary(prot)
Call:
protest(X = pca.run1, Y = pca.run2, permutations = 1000) 

Procrustes Sum of Squares (m12 squared):        0.0005705 
Correlation in a symmetric Procrustes rotation: 0.9997 
Significance:  0.000999 

Permutation: free
Number of permutations: 1000

```

Since the technical replicates are higly correlated they will be merged by summing the numbers of reads of both.
Here only ASVs found in **both** replicates are kept.


```

dt_togroup <- dt_filt[rownames(dt_filt) %in% metadata_togroup$id,] %>%
  na_if(., "0") %>% #here I replace 0 abundances by NA, this guarantees that ASVs not present in both replicates won't be kept
  data.table(., keep.rownames = T)
  
all(dt_togroup$rn == metadata_togroup$id) #checking if samples are in the same order

dt_togroup <- right_join(dt_togroup, metadata_togroup[,.(id, sample)], by=c("rn"="id")) %>%
  select(., -rn) #what this command does is to replace unique rownames (with run identification) by a column with the name of the samples

dt_togroup$sample

dt_grouped  <- dt_togroup[, lapply(.SD, sum), by=sample] %>% #grouping the samples with same name by the sum of the number of reads
  replace(., is.na(.), 0) %>% #now NAs are replaced back by O
  arrange(., sample) %>% #sorting the samples
  column_to_rownames(., "sample") %>% #setting rownames
  .[, which(colSums(.) > 0)] #removing ASVs with 0 reads
  
```

After the replicates are merged I check if there are still ASVs not classified at the domain level which were kept. If yes, I run VSEARCH to get their taxonomic assignation. If they can't be classified at a 80% identity threshold using VSEARCH I remove them.



```

tax_filtered_grouped <- tax_filtered[tax_filtered$ASVId %in% colnames(dt_grouped),] 

unknown <- tax_filtered_grouped %>% filter(., domain == "unclassified")

write.fasta(names=as.list(unknown[,ASVId]),
            sequences=as.list(unknown[,V1]), 
            file.out="tmp/asvs_domain_unknown.fasta")

reference_515F_926R <- "path/silva_v138/SILVA_138_SSURef_NR99_tax_silva_515F_926R.fasta"
tax_unknown <- "tmp/asvs_domain_unknown.txt"
fas <- "tmp/asvs_domain_unknown.fasta"

system(paste0("/Users/marianacamaradosreis/Downloads/vsearch-2.15.0-macos-x86_64/bin/vsearch ",
              "--usearch_global ",  fas, 
              " --db ",              reference_515F_926R,
              " --id 0.80 ",
              " --blast6out ",          tax_unknown))

dim(unknown)

```

From the 15 ASVs remaining without domain classification by IDTAXA, only ASV_115, ASV_172, ASV_294 seem like real bacterial sequences (confirmed by VSEARCH 80%). 
The others were not classified. If I blast them they look like chimera - low query cover and %similarity. So they will be removed.

```

keep<- c("ASV_115", "ASV_172", "ASV_294")

unknown <- unknown[!unknown$ASVId %in% keep,] #hehre I remove the 3 ASVs that I want to keep from the unknown vector
dim(unknown) # 12

tax_filtered_grouped <- tax_filtered_grouped[!tax_filtered_grouped$ASVId %in% unknown$ASVId,] #then I remove the unknown ASVs from the taxonomy
dim(tax_filtered_grouped)
tax_filtered_grouped$ASVId #checking if the ASV_115, ASV_172 and ASV_294 were kept

dt_grouped <- dt_grouped[,colnames(dt_grouped) %in% tax_filtered_grouped$ASVId] 

#Now removing the unknown asvs from the ASV table

```

Checking for the samples which without technical replicates

```

##Selecting last sampling (which were not merged)

last_sampling <- dt_filt[rownames(dt_filt) %in% metadata_last_sampling$id,]
last_sampling <- last_sampling[,which(colSums(last_sampling) >0)]


tax_filtered_last_sampling<- tax_filtered[tax_filtered$ASVId %in% colnames(last_sampling),] 
dim(tax_filtered_last_sampling)

unknown_last <- tax_filtered_last_sampling %>% filter(., domain == "unclassified")

#no ASVs without classification were detected

```

The last step is build the ASV table with the grouped replicates and with the samples without technical replicates

```

#Merging the last sampling with the grouped table

last_sampling <- data.frame(t(last_sampling)) %>%
  rownames_to_column(., "ASVId")

dt_grouped <- data.frame(t(dt_grouped)) %>%
  rownames_to_column(., "ASVId")

all_samplings <- full_join(dt_grouped, last_sampling, by=("ASVId"))

last_sampling$ASVId %in% dt_grouped$ASVId #there was only one ASV which appeared only in last sampling

all_samplings[is.na(all_samplings)] <- 0 #replacing NA by 0

all_samplings <- all_samplings %>%
  column_to_rownames(., "ASVId")
  
sum(all_samplings)

sum(all_samplings)/sum(dt_filt_chloro)

```
With all the filters of abundance and non-classified ASVs I still kept 99.59% of the reads.
The last step of the pipeline is to select only the samples which will be used for the analysis.
Those are 
11 environmental samples: 0.2-3 µm for AM and PM DCM sites and outside of bloom DCM site.
96 culture samples: 8 time poits of the experiment

```

metadata <- metadata[size_fraction != "3.0-20" & treatments != "D1S1-surf",]
metadata <- metadata[size_fraction != ">20",]
dim(metadata) #107 samples

dt_final <- all_samplings[rownames(all_samplings) %in% metadata$sample,] 

write.table(dt_final, "path/coccobloom_nochloroplast_replicates_merged_filt_abund_sum_all_samplings.txt", sep="\t", row.names=TRUE)

```

The final ASV table contains 107 samples, 294 ASVs and 6,017,019. This table was used for the analysis of the paper.
