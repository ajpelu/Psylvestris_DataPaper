#Processing of Illumina (MiSeq reads). 
#ITS2 (fungal) reads

#Install all the packages needeed
#dada2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")

#phyloseq
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

#Biocgenerics
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BiocGenerics")

#Biostrings
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")

#ShortRead
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ShortRead")

#ggplot2 (for plotting)
install.packages("ggplot2")

#tidyverse
install.packages("tidyverse")

#Load corresponding libraries
library(dada2)
library(phyloseq)
library(BiocGenerics)
library(Biostrings)
library(ShortRead)
library(devtools)
library(ggplot2)
library(tidyverse)

path <- "C:/Users/radik/Desktop/Psylvestris_DataPper" #set the path where the fastq.gz files are
list.files(path)#check that all the needed files are in the path

#NOTE: if our scientific project includes different sequencing runs (for instance, if you have a lot of samples that do not fit well into just one run), 
#we should repeat steps 1) to 6) for each sequencing run, and then (step 7), merge all the sequence tables obtained in the step 6)

#1) ADAPTION OF THE FILE NAMES------------------------------------------------
#Sort Forward (F) and Reverse (R) files separatedly
#(F and R fastq filenames should have the following format: SAMPLENAME_R1_001.fastq.gz;SAMPLENAME_R2_001.fastq.gz)
#in this example, the names of F and R file are the following: NGS006-22-EXP4-ITS2-SAMPLENAME_yyy-R1_001.fastq.gz
fnFs = sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs = sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: xxx-SAMPLENAME_yyy_R1_001.fastq.gz
sample.names_raw = sapply(strsplit(basename(fnFs), "_"), `[`, 1)
  
sample.names = gsub("NGS006-22-EXP4-ITS2-", "", sample.names_raw) #replace the xxx part of the string, by nothing (or by the characters that allow you to get the sample names
#output in this example: SAMPLENAME
    
#2) CHECK THE QUALITY OF THE READS-----------------------------------------------------------------

#a) Check the number of reads
raw_reads_count = NULL

for (i in 1:length(fnFs)){
    raw_reads_count <- rbind(raw_reads_count, length(ShortRead::readFastq(fnFs[i])))
  } #it counts the number of reads in fnFs

rownames(raw_reads_count)= sample.names
colnames(raw_reads_count)= "Number of reads"
a=data.frame("_"=rownames(raw_reads_count),raw_reads_count)
colnames(a)= c("Sample", "Number_of_reads")

#Check which samples have the highlest and lowest amount of reads
cbind(row.names(raw_reads_count)[which.min(raw_reads_count)],min(raw_reads_count))
cbind(row.names(raw_reads_count)[which.max(raw_reads_count)],max(raw_reads_count))

write.table(data.frame("_"=rownames(raw_reads_count),raw_reads_count),file="NumberRawReads.txt", sep="\t",row.names =F)

#b) Check the length of the reads
reads=ShortRead::readFastq(fnFs) #it saves the reads into a new variable
uniques = unique(reads@quality@quality@ranges@width)#to get the length of the reads (bp)

counts= NULL#it counts the number of reads of each length
for (i in 1:length(uniques)) {
  counts<- rbind(counts,length(which(reads@quality@quality@ranges@width==uniques[i])))
  
}

#check the histogram
histogram =  cbind(uniques,counts)
colnames(histogram) = c("Seq.length", "counts")

write.table(histogram,"Sequences_length.txt", sep="\t",row.names =F)

head(histogram[order(histogram[,1],decreasing = TRUE),]) 

#plot the histogram
hist(reads@quality@quality@ranges@width, main="Forward length distribution", xlab="Sequence length", ylab="Raw reads")

#b) #c) Check the quality plots
plotQualityProfile(fnFs[4:5])#here you can plot whichever of the samples you want (in this example, samples number 4 and 5)
plotQualityProfile(fnRs[4:5])

#3)CUTADAPT, primer removal--------------------------------------------------
#define the sequence of the primers
FWD ="GTGARTCATCGAATCTTTG" #primer fITS7
REV = "TCCTCCGCTTATTGATATGC"  #primer ITS4

#Create a function to calculate all the orientations of F and R primers
allOrients =function(primer) {
  require(Biostrings)
  dna= DNAString(primer)  #BioStrings package works with strings but not with characters
  orients = c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString)) #change to string
}

FWD.orients = allOrients(FWD)#pass the function to primer F
REV.orients = allOrients(REV)#pass the function to primer R
FWD.orients #check the calculated orientations
REV.orients

#Count the number of times both primers (in all the possible orientations) appear in our dataset
primerHits = function(primer, fn) {
  nhits = vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

#check the number of reads in which the primers F and R (in all possible orientations) are found in sample number "5"
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[5]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[5]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[5]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[5]]))
#On the left, we have the presence of forward and reverse primers (FWD and REV) in forward and reverse reads (FWD.ForwardReads, FWD.ReverseReads, REV.ForwardReads and REV.ReverseReads). On the top of the table, we see the four possible orientations of the primers (Forward, Complement, Reverse and RevComp).

cutadapt = "/usr/local/bin/cutadapt" #path to cutadapt in your computer/machine

system2(cutadapt, args = c("--version")) # Run shell commands from R

#Create a directoy or folder where your reads will be saved
path.cut = file.path(path, "cutadapt") 

if(!dir.exists(path.cut)) dir.create(path.cut)

fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

#Produce arguments for cutadapt (visit the cutadapt website for more information)
FWD.RC <- dada2:::rc(FWD)  #calculate the reverse complementary sequence 
REV.RC <- dada2:::rc(REV)

#add the specific adapters to F and R reads, which are need for cutadapt
R1.flags <- paste0("-a", " ", "^",FWD,"...", REV.RC) 
R2.flags <- paste0("-A"," ","^", REV, "...", FWD.RC)

#Run cutadapt itself
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,"-m", 1, 
                           "--discard-untrimmed",
                             "-j",0,
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                           fnFs[i], fnRs[i],# input files
                             "--report=minimal")) #Report minimal reports a summary
  
}

#  -n 2: remove the primers
#   -m 1: remove empty reads

#cutadapt will remove the primers and also the sequences in which none of the primers are found

#after running cutadapt, check the number of primers found now in the output reads, in sample number "5"
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[5]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[5]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[5]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[5]]))

#be CAREFUL! it could be possible to find still some primers in your dataset. Do not worry if there are no many. Cutadapt and primerHits function
#do not work in the same way, so some primers can be found.

#4) FILTER AND TRIMMING STEPS----------------------------------------------------------------
#Create the folder "filtered" to save the trimmed reads
filtFs = file.path(path.cut, "filtered", basename(fnFs.cut)) 
filtRs = file.path(path.cut, "filtered", basename(fnRs.cut))

out = filterAndTrim(fnFs.cut, filtFs, fnRs.cut, filtRs, maxN = 0, 
                      maxEE = c(2, 3), truncQ = 2, minLen = 50, rm.phix = TRUE, 
                      compress = TRUE, multithread = TRUE) 
  #input reads F without primers, output reads F, input reads R without primers, output reads R,
  #maxN: number of allowed ambiguities (N) in the reads
  #maxEE: maximum expected errors in F and R reads; they should be selected according to reads' quality
  #trunqQ: the minimum quality score of each nucleotides; remove reads with at least one nucleotide with an associated quality score under the selected value
  #minLen: the minimum size of the reads; remove all the reads under this size (bp)

View(out)

#5) DADA2--------------------------------------------
#(visit DADA2 website for more information about the specific steps)
#Learn error rates................
errF =learnErrors(filtFs, multithread=T, verbose=1 )
errR = learnErrors(filtRs, multithread=T, verbose=1)

#View Error plots................
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#Sample inference.................
#it calculates the inferred composition of the samples. Removes all the possible sequence errors
#and keeps just the real sequences.
dadaFs = dada(filtFs, err=errF, multithread=TRUE)
dadaRs = dada(filtRs, err=errR, multithread=TRUE)
dadaFs[[5]]#check the inferred sample number 5. We can see the number of ASVs and unique sequences

names(dadaFs)=sample.names#giving the correct name of the samples
names(dadaRs)=sample.names

#Merging F and R reads.............................
mergers <=mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[5]]) #check the results for sample number 5

#Calculate the sequence tables......................
seqtab_ITS_primavera_Rizo = makeSequenceTable(mergers)
dim(seqtab_ITS_primavera_Rizo) #it shows the number of samples (including MOCK samples, negative control and so on) and the number of ASVs
saveRDS(seqtab_ITS_primavera_Rizo, file="C:/Users/radik/Desktop/Psylvestris_DataPper/seqtab_ITS_Primavera_rizo.rds")
#save the corresponding RDS file 

#IMPORTANT: if we have different sequencing runs included in the same scientific project, we should merge then
#at this point. For that purpose, we should merge the sequence tables coming from each sequencing runs

#6) MERGE SEQTABS OBTAINED FROM DIFFERENT SEQUENCING RUNS--------------------------------------------
#(skip this step if you included all the samples in one sequencing run. In that case, go to step7)

seqtab_pri_endo=readRDS("seqtab_ITS_Primavera_endo.rds")#load each sequence table
seqtab_pri_rizo=readRDS("seqtab_ITS_Primavera_rizo.rds")
seqtab_ve_endo=readRDS("seqtab_ITS_Verano_endo.rds")
seqtab_ve_endo_rizo=readRDS("seqtab_ITS_Verano_EndoRizo.rds")
seqtab_ve_rizo=readRDS("seqtab_ITS_Rizo_verano.rds")

#merge all the seqTabs
mergedSeqTab=mergeSequenceTables(seqtab_pri_endo,seqtab_pri_rizo,
                                 seqtab_ve_endo,seqtab_ve_endo_rizo,seqtab_ve_rizo1,
                                 repeats="error")
#7) CHIMERA REMOVAL-------------------------------------------------

#Eliminacion de quimeras..................
seqtab.nochim = removeBimeraDenovo(mergedSeqTab, method="consensus", multithread=TRUE, verbose=TRUE)
#in case you have just one seqtab, replace "mergedSeqTab" by "seqtab_ITS_primavera_Rizo"
dim(seqtab.nochim)#indicates the number of samples and ASVs, but not the number of sequences

#8)TAXONOMICAL CLASSIFICATION--------------------------------------------------------
taxa_UNITE = assignTaxonomy(seqtab.nochim, "/mnt/datos/databases/sh_general_release_dynamic_04.02.2020.fasta", multithread=TRUE)
#indicate the path where the database is located.

#9) OBTAINING THE ASV TABLE
#Format the data
ASV = seqtab.nochim
ASVt = t(ASV)

#Replace NA values by "unclassified" and remove the column "species" created during the classification
taxa_UNITE_na = apply(taxa_UNITE,2, tidyr::replace_na, "unclassified")[,-7]

#Rename the ASVs (if we have 100 ASVs, we will name them from ASV001 to ASV100)
number.digit = nchar(as.integer(nrow(ASVt)))
names =paste0("ASV%0", number.digit, "d") #As many 0 as digits
ASV_names<- sprintf(names, 1:nrow(ASVt))

#Join taxonomy and ASV table
ASV_table_classified_raw = cbind(as.data.frame(taxa_UNITE_na,stringsAsFactors = FALSE),as.data.frame(ASV_names, stringsAsFactors = FALSE),as.data.frame(ASVt,stringsAsFactors = FALSE))

#Add the sequence of each ASV to the ASV table in the correct position
ASV_seqs = rownames(ASV_table_classified_raw)
rownames(ASV_table_classified_raw) <- NULL
ASV_table_classified_raw = cbind(ASV_seqs, ASV_table_classified_raw)

#10) REMOVE THE MOCK COMMUNITY SAMPLES------------------------------------
#(we will use as trimming cut-off the cut-off value obtained for bacterial dataset)
#(in our case, the bacterial cut-off value was 0.001792%)

#obtain the number of sequences per ASV
ASV_sums = rowSums(ASV_table_classified_raw[,9:ncol(ASV_table_classified_raw)])

#Obtain the total number of sequences
sum.total=sum(ASV_sums)

#Calculate the number of sequences that suppose that percentage in the fungal dataset
nseq_cutoff=(0.001792/100)*sum.total

#Filtering according to the stablished cut-offf
ASV_filtered_MOCK= ASV_table_classified_raw[which(ASV_sums>nseq_cutoff),]

#Order the corresponding table
ASV_filtered_MOCK=ASV_filtered_MOCK[order(ASV_filtered_MOCK[["ASV_names"]]),]

#11) REMOVAL OF PLASTIDS, PLANT, CHLOROPLAST AND OTHER ARTIFACTS SEQUENCES
ASV_table_MOCK=ASV_filtered_MOCK[(which(ASV_filtered_MOCK$Genus!="Streptophyta"
                                         & ASV_filtered_MOCK$Genus!="Chlorophyta"
                                         & ASV_filtered_MOCK$Genus!="Bacillariophyta"
                                         & ASV_filtered_MOCK$Family!="Streptophyta"
                                         & ASV_filtered_MOCK$Family!="Chlorophyta"
                                         & ASV_filtered_MOCK$Family!="Bacillariophyta"
                                         & ASV_filtered_MOCK$Family!="Mitochondria"
                                         & ASV_filtered_MOCK$Class!="Chloroplast"
                                         & ASV_filtered_MOCK$Order!="Chloroplast"
                                         & ASV_filtered_MOCK$Kingdom!="Eukaryota"
                                         & ASV_filtered_MOCK$Kingdom!="unclassified")),]

ASV_final=ASV_table_MOCK[which(ASV_table_MOCK$Kingdom=="k__Fungi"),]#we also remove those ASV not classified at Kingdom level

#Save the definitive ASV Table 
write.table(ASV_final, file="ASV_final_ITS.txt", sep="\t")