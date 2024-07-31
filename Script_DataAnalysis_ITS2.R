#Script for the analysis of already processed ITS2 (fungal) reads obtained by Illumina (MiSeq)
#Install packages 
BiocManager::install("phyloseq")
install.packages("dplyr")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("vegan")
install.packages("purrr")

#load the libraries
library(phyloseq)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(vegan)
library(purrr)

#LOAD ALL THE DATA ----------------------------------------------------------
data=read.table("ASV_final_ITS.txt", header = T, sep="\t")
#we modify the taxonomy
data$Kingdom=gsub("k__", "",  data$Kingdom)
data$Phylum=gsub("p__", "",  data$Phylum)
data$Class=gsub("c__", "",  data$Class)
data$Order=gsub("o__", "",  data$Order)
data$Family=gsub("f__", "",  data$Family)
data$Genus=gsub("g__", "",  data$Genus)
colnames_data=colnames(data)

mt=read.table("metadatos.txt",sep="\t", header=T)#load the metadata describing the main characteristics of each sample
row.names(mt)=mt$Sample

tax =data[,2:8]#select the columns where the bacterial taxonomy is indicated (From Phylum level to the ASV)
ASV=data[,9:ncol(data)]#select the columns where the number of sequences of each ASV in each sample is indicated

dna=Biostrings::DNAStringSet(data$ASV_seqs) #load the sequence of each ASV
names(dna)= data$ASV_names

row.names(tax)=data$ASV_names
row.names(ASV)=data$ASV_names
identical(rownames(ASV), rownames(tax))#it should be TRUE

#2)CREATE A PHYLOSEW OBJECT-------------------------------------------------------------
phy_OTUtable=otu_table(ASV, taxa_are_rows = T)
phy_taxonomy=tax_table(as.matrix(tax))
phy_metadata=sample_data(mt) 

#Put everything into a phyloseq object
phy_data_total=phyloseq(phy_OTUtable,phy_taxonomy,phy_metadata)

#Check the number of samples included in each object
nsamples(phy_OTUtable); sum(sample_names(phy_metadata) %in% sample_names(phy_OTUtable)); nsamples(phy_metadata)

#3) FIRST TRIMMING OF THE DATASET-------------------------------------------------------------------
#(our dataset included samples obtained in spring and summer, samples coming from P.sylvestris, P.pinaster and P. sylvestris subsp. nevadensis)
#(here we only wanted the samples corresponding to P. sylvestris sampled in spring)
phy_data1=subset_samples(phy_data_total, Season != "Summer")#removal of summer samples
phy_data2=subset_samples(phy_data1, Species == "Psylvestris")#selecting the plant host

#check the total number of sequences
numsec=as.data.frame(colSums(otu_table(phy_data2)))
View(numsec) #it is convenient to order this table

write.table(data.frame(" "=rownames(numsec),numsec),file="NumberofWorkingSequences.txt", sep="\t",row.names =F)

#4) SUBSET SAMPLES-----------------------------------------------------------------------
#(we want to separate the data by mountainous area under study and by plant compartment)
Almijara=subset_samples(phy_data2, Site=="Almijara" )
Almijara=prune_taxa(taxa_sums(Almijara)>0.0, Almijara)

Almijara_rizo=subset_samples(Almijara, Compartment=="Rhizosphere" )
Almijara_rizo=prune_taxa(taxa_sums(Almijara_rizo)>0.0, Almijara_rizo)

Almijara_endo=subset_samples(Almijara, Compartment=="Endosphere" )
Almijara_endo=prune_taxa(taxa_sums(Almijara_endo)>0.0, Almijara_endo)

Baza=subset_samples(phy_data2, Site=="Baza")
Baza=prune_taxa(taxa_sums(Baza)>0.0, Baza)

Baza_rizo=subset_samples(Baza, Compartment=="Rhizosphere" )
Baza_rizo=prune_taxa(taxa_sums(Baza_rizo)>0.0, Baza_rizo)

Baza_endo=subset_samples(Baza, Compartment=="Endosphere" )
Baza_endo=prune_taxa(taxa_sums(Baza_endo)>0.0, Baza_endo)

Nevada=subset_samples(phy_data2, Site=="SierraNevada")
Nevada=prune_taxa(taxa_sums(Nevada)>0.0, Nevada)

Nevada_rizo=subset_samples(Nevada, Compartment=="Rhizosphere" )
Nevada_rizo=prune_taxa(taxa_sums(Nevada_rizo)>0.0, Nevada_rizo)

Nevada_endo=subset_samples(Nevada, Compartment=="Endosphere" )
Nevada_endo=prune_taxa(taxa_sums(Nevada_endo)>0.0, Nevada_endo)

rizo=subset_samples(phy_data2, Compartment=="Rhizosphere" )
rizo=prune_taxa(taxa_sums(rizo)>0.0, rizo)

endo=subset_samples(phy_data2, Compartment=="Endosphere" )
endo=prune_taxa(taxa_sums(endo)>0.0, endo)

#5) RAREFACTION CURVES----------------------------------------
#1.Almijara Rhizosphere..........................................
mt=as.data.frame(sample_data(Almijara_rizo))[order(as.character(rownames(as.data.frame(sample_data(Almijara_rizo)))),decreasing=F),]#save the metadata of this sub-dataset
otu = otu_table(Almijara_rizo)#save the ASV table of this sub-dataset
otu = as.data.frame(t(otu))
otu=otu[order(as.character(rownames(otu)),decreasing=FALSE),]

rownames(otu)==rownames(mt) 
rownames(otu)=paste0(rownames(otu),"/",mt$Condition)#prepare the data for a better plotting. Condition= "Symptomatic" and "Asymptomatic"
sample_names = rownames(otu)

out <- rarecurve(otu, step = 100, label = F) #rarefaction itself
rare <- lapply(out, function(x){  #create a table from the rarefaction
  b <- as.data.frame(x)
  b <- data.frame(ASV = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})

names(rare) <- sample_names
rare <- map_dfr(rare, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "Sample")

rare$Condition=rare$Sample  #data preparation for the plotting
rare$Condition=gsub(".*/", "", rare$Condition)#remove all the characters before "/"
rare$raw.read=as.numeric(rare$raw.read)
#Esto que acabamos de hacer solo vale para las muestras en concreto que tengo ahora

#Plotting......................
p<-ggplot(rare, aes(x=raw.read, y=ASV, colour=Condition, group=Sample)) +#coloring by condition of the trees
  theme_bw()+
  geom_point(aes(colour=Condition), size=0.85)+
  geom_line(aes(colour=Condition),linewidth=1.2)+
  scale_fill_manual(values = c("Symptomatic"="black","Asymptomatic"="#3399FF"))+#specific colors
  scale_color_manual(values = c("Symptomatic"="black", "Asymptomatic"="#3399FF"),
                     name="Condition",
                     breaks=c("Symptomatic", "Asymptomatic"),
                     labels=c("Symptomatic", "Asymptomatic"))+
  labs(title= "Sierra de Almijara (Rhizosphere)", x="Number of sequences", y="Number of ASV")+
  guides(alpha=FALSE)+
  theme(legend.key=element_blank(),
        legend.title.align = 0.85,
        legend.title = element_text(face="bold",size=14),
        axis.text = element_text(size=14),
        axis.title = element_text(size = 16),
        plot.title = element_text(hjust=0.5, face="bold", size=16),
        legend.text = element_text(size = 16))
x11()
p

#save the plot as .tiff
ggsave(filename = "Rarefacc_rizo_Almijara_SciData_Fungi.tiff", plot = p,device = tiff(),width = 18, height = 16, units = "cm", dpi = 800)
rare=NULL 

#2. Almijara Endosphere..........................................
mt=as.data.frame(sample_data(Almijara_endo))[order(as.character(rownames(as.data.frame(sample_data(Almijara_endo)))),decreasing=F),]
otu = otu_table(Almijara_endo)
otu = as.data.frame(t(otu))
otu=otu[order(as.character(rownames(otu)),decreasing=FALSE),]

rownames(otu)==rownames(mt) 
rownames(otu)=paste0(rownames(otu),"/",mt$Condition)
sample_names = rownames(otu)

out <- rarecurve(otu, step = 100, label = F) 
rare <- lapply(out, function(x){  
  b <- as.data.frame(x)
  b <- data.frame(ASV = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})

names(rare) <- sample_names
rare <- map_dfr(rare, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "Sample")

rare$Condition=rare$Sample 
rare$Condition=gsub(".*/", "", rare$Condition)
rare$raw.read=as.numeric(rare$raw.read)

p<-ggplot(rare, aes(x=raw.read, y=ASV, colour=Condition, group=Sample)) +
  theme_bw()+
  geom_point(aes(colour=Condition), size=0.85)+
  geom_line(aes(colour=Condition), size =1.2)+
  scale_fill_manual(values = c("Symptomatic"="black","Asymptomatic"="deepskyblue"))+
  scale_color_manual(values = c("Symptomatic"="black", "Asymptomatic"="deepskyblue"),
                     name="Condition",
                     breaks=c("Symptomatic", "Asymptomatic"),
                     labels=c("Symptomatic", "Asymptomatic"))+
  labs(title= "Sierra de Almijara (Endosphere)", x="Number of sequences", y="Number of ASV")+
  guides(alpha=FALSE)+
  theme(legend.key=element_blank(),
        legend.title.align = 0.85,
        legend.title = element_text(face="bold",size=14),
        axis.text = element_text(size=14),
        axis.title = element_text(size = 16),
        plot.title = element_text(hjust=0.5, face="bold", size=16),
        legend.text = element_text(size = 16))
x11()
p

ggsave(filename = "Rarefacc_Almijara_endo_SciData_fungi.tiff", plot = p,device = tiff(),width = 18, height = 16, units = "cm", dpi = 800)
rare=NULL 

#3.Baza. Rhizosphere..........................................
mt=as.data.frame(sample_data(Baza_rizo))[order(as.character(rownames(as.data.frame(sample_data(Baza_rizo)))),decreasing=F),]
otu = otu_table(Baza_rizo)
otu = as.data.frame(t(otu))
otu=otu[order(as.character(rownames(otu)),decreasing=FALSE),]

rownames(otu)==rownames(mt)
rownames(otu)=paste0(rownames(otu),"/",mt$Condition)
sample_names = rownames(otu)

out <- rarecurve(otu, step = 100, label = F) 
rare <- lapply(out, function(x){
  b <- as.data.frame(x)
  b <- data.frame(ASV = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})

names(rare) <- sample_names
rare <- map_dfr(rare, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "Sample")

rare$Condition=rare$Sample 
rare$Condition=gsub(".*/", "", rare$Condition)
rare$raw.read=as.numeric(rare$raw.read)


p<-ggplot(rare, aes(x=raw.read, y=ASV, colour=Condition, group=Sample)) +
  theme_bw()+
  geom_point(aes(colour=Condition), size=0.85)+
  geom_line(aes(colour=Condition),size=1.2)+
  scale_fill_manual(values = c("Symptomatic"="#66CC33","Asymptomatic"="#66CC33"))+
  scale_color_manual(values = c("Symptomatic"="#66CC33", "Asymptomatic"="#66CC33"),
                     name="Condition",
                     breaks=c("Symptomatic", "Asymptomatic"),
                     labels=c("Symptomatic", "Asymptomatic"))+
  labs(title= "Sierra de Baza (Rhizosphere)", x="Number of sequences", y="Number of ASV")+
  guides(alpha=FALSE)+
  theme(legend.key=element_blank(),
        legend.title.align = 0.85,
        legend.title = element_text(face="bold",size=14),
        axis.text = element_text(size=14),
        axis.title = element_text(size = 16),
        plot.title = element_text(hjust=0.5, face="bold", size=16),
        legend.text = element_text(size = 16))
x11()
p

ggsave(filename = "Rarefacc_Baza_Rizo_SciData_fungi.tiff", plot = p,device = tiff(),width = 18, height = 16, units = "cm", dpi = 800)
rare=NULL

#4.Baza Endosphere..........................................
mt=as.data.frame(sample_data(Baza_endo))[order(as.character(rownames(as.data.frame(sample_data(Baza_endo)))),decreasing=F),]
otu = otu_table(Baza_endo)
otu = as.data.frame(t(otu))
otu=otu[order(as.character(rownames(otu)),decreasing=FALSE),]

rownames(otu)==rownames(mt)
rownames(otu)=paste0(rownames(otu),"/",mt$Condition)
sample_names = rownames(otu)

out <- rarecurve(otu, step = 100, label = F) 
rare <- lapply(out, function(x){  
  b <- as.data.frame(x)
  b <- data.frame(ASV = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})

names(rare) <- sample_names
rare <- map_dfr(rare, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "Sample")

rare$Condition=rare$Sample 
rare$Condition=gsub(".*/", "", rare$Condition)
rare$raw.read=as.numeric(rare$raw.read)


p<-ggplot(rare, aes(x=raw.read, y=ASV, colour=Condition, group=Sample)) +
  theme_bw()+
  geom_point(aes(colour=Condition), size=0.85)+
  geom_line(aes(colour=Condition),size=1.2)+
  scale_fill_manual(values = c("Symptomatic"="green","Asymptomatic"="#3399FF"))+
  scale_color_manual(values = c("Symptomatic"="green", "Asymptomatic"="#3399FF"),
                     name="Condition",
                     breaks=c("Symptomatic", "Asymptomatic"),
                     labels=c("Symptomatic", "Asymptomatic"))+
  labs(title= "Sierra de Baza (Endosphere)", x="Number of sequences", y="Number of ASV")+
  guides(alpha=FALSE)+
  theme(legend.key=element_blank(),
        legend.title.align = 0.85,
        legend.title = element_text(face="bold",size=14),
        axis.text = element_text(size=14),
        axis.title = element_text(size = 16),
        plot.title = element_text(hjust=0.5, face="bold", size=16),
        legend.text = element_text(size = 16))
x11()#estas dos lineas correrlas juntas. En el servidor no funciona
p

ggsave(filename = "Rarefacc_Baza_endosphere_SciData_fungi.tiff", plot = p,device = tiff(),width = 18, height = 16, units = "cm", dpi = 800)
rare=NULL

#5.Sierra Nevada Rhizosphere......................................................
mt=as.data.frame(sample_data(Nevada_rizo))[order(as.character(rownames(as.data.frame(sample_data(Nevada_rizo)))),decreasing=F),]
otu = otu_table(Nevada_rizo)
otu = as.data.frame(t(otu))
otu=otu[order(as.character(rownames(otu)),decreasing=FALSE),]

rownames(otu)==rownames(mt) 
rownames(otu)=paste0(rownames(otu),"/",mt$Condition)
sample_names = rownames(otu)

out <- rarecurve(otu, step = 100, label = F) 
rare <- lapply(out, function(x){  
  b <- as.data.frame(x)
  b <- data.frame(ASV = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})

names(rare) <- sample_names
rare <- map_dfr(rare, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "Sample")

rare$Condition=rare$Sample 
rare$Condition=gsub(".*/", "", rare$Condition)
rare$raw.read=as.numeric(rare$raw.read)


p<-ggplot(rare, aes(x=raw.read, y=ASV, colour=Condition, group=Sample)) +
  theme_bw()+
  geom_point(aes(colour=Condition), size=0.85)+
  geom_line(aes(colour=Condition),size=1.2)+
  scale_fill_manual(values = c("Symptomatic"="black","Asymptomatic"="magenta"))+
  scale_color_manual(values = c("Symptomatic"="black", "Asymptomatic"="magenta"),
                     name="Condition",
                     breaks=c("Symptomatic", "Asymptomatic"),
                     labels=c("Symptomatic", "Asymptomatic"))+
  labs(title= "Sierra Nevada (Rhizosphere)", x="Number of sequences", y="Number of ASV")+
  guides(alpha=FALSE)+
  theme(legend.key=element_blank(),
        legend.title.align = 0.85,
        legend.title = element_text(face="bold",size=14),
        axis.text = element_text(size=14),
        axis.title = element_text(size = 16),
        plot.title = element_text(hjust=0.5, face="bold", size=16),
        legend.text = element_text(size = 16))
x11()#estas dos lineas correrlas juntas. En el servidor no funciona
p

ggsave(filename = "Rarefacc_SNevada_Rizo_SciData_fungi.tiff", plot = p,device = tiff(),width = 18, height = 16, units = "cm", dpi = 800)
rare=NULL

#6. Sierra Nevada Endosphere......................................................
mt=as.data.frame(sample_data(Nevada_endo))[order(as.character(rownames(as.data.frame(sample_data(Nevada_endo)))),decreasing=F),]
otu = otu_table(Nevada_endo)
otu = as.data.frame(t(otu))
otu=otu[order(as.character(rownames(otu)),decreasing=FALSE),]

rownames(otu)==rownames(mt) 
rownames(otu)=paste0(rownames(otu),"/",mt$Condition)
sample_names = rownames(otu)

out <- rarecurve(otu, step = 100, label = F) 
rare <- lapply(out, function(x){  
  b <- as.data.frame(x)
  b <- data.frame(ASV = b[,1], raw.read = rownames(b))
  b$raw.read <- as.numeric(gsub("N", "",  b$raw.read))
  return(b)
})

names(rare) <- sample_names
rare <- map_dfr(rare, function(x){
  z <- data.frame(x)
  return(z)
}, .id = "Sample")

rare$Condition=rare$Sample 
rare$Condition=gsub(".*/", "", rare$Condition)
rare$raw.read=as.numeric(rare$raw.read)

p<-ggplot(rare, aes(x=raw.read, y=ASV, colour=Condition, group=Sample)) +
  theme_bw()+
  geom_point(aes(colour=Condition), size=0.85)+
  geom_line(aes(colour=Condition),size=1.2)+
  scale_fill_manual(values = c("Symptomatic"="black","Asymptomatic"="pink"))+
  scale_color_manual(values = c("Symptomatic"="black", "Asymptomatic"="pink"),
                     name="Condition",
                     breaks=c("Symptomatic", "Asymptomatic"),
                     labels=c("Symptomatic", "Asymptomatic"))+
  labs(title= "Sierra Nevada (Endosphere)", x="Number of sequences", y="Number of ASV")+
  guides(alpha=FALSE)+
  theme(legend.key=element_blank(),
        legend.title.align = 0.85,
        legend.title = element_text(face="bold",size=14),
        axis.text = element_text(size=14),
        axis.title = element_text(size = 16),
        plot.title = element_text(hjust=0.5, face="bold", size=16),
        legend.text = element_text(size = 16))
x11()
p

ggsave(filename = "Rarefacc_SNevada_endo_SciData_fungi.tiff", plot = p,device = tiff(),width = 18, height = 16, units = "cm", dpi = 800)
rare=NULL