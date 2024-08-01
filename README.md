# Psylvestris_DataPaper
Processing of Illumina MiSeq sequencing reads (16S V3-V4 and ITS2) from the root endosphere and rhizosphere of decline-affected and asymptomatic Pinus sylvestris trees
# Welcome to SUMHAL Project 

The SUMHAL Project aims at deciphering the changes in the microbial communities associated to the roots of different *Pinus* species as a consequence of forest decline. Three different mountainous areas covered by *Pinus* spp. stands in the Southeast Spain were selected for this purpose.
Since different pine species in different conditions were selected, different sub-project were created, which resulted in different scientific publications. Below you will find specific information about the main project and sub-projects

# Experimental design

In the following table, all the information regarding the project is summarized 

|Species  |Season  |Plant compartment  |Area  | Condition  |Organism  |
|----------|----------|------------------------|------ |--------------|-
|*P. sylvestris* | Summer |Rhizosphere|Sierra Nevada|Symptomatic|Bacteria
|*P. pinaster* |Spring|Root endosphere|Sierra de Baza|Asymptomatic|Fungi|
|*P. sylvestris* subsp. *nevadensis*|||Sierra de Almijara, Tejeda y Alhama|

Here, "*Condition*" refers to the presence (or absence) of pine decline symptoms

12 tree replicates were selected in all cases, however, the development of the symptoms depended on the considered mountainous area and plant species. On the other hand, not all the plant host species were found in all the areas. Have a look at the following table to know the distribution of the samples:
|Area | Species | Condition (Acronym)|
|------|----------|--------|
|Sierra Nevada| *P. sylvestris*|Symptomatic (NSD)|
|||Aymptomatic (NSH)|
|| *P. pinaster*|Symptomatic (NPD)|
|||Aymptomatic (NPH)|
|| *P. sylvestris* subsp. *nevadensis*|Symptomatic (NAD)
|||Symptomatic (NAD)
|Sierra de Baza|*P. sylvestris*| Symptomatic (BSD)
||*P. pinaster*| Symptomatic (BPD)|
||*P. sylvestris* subsp. *nevadensis*|Asymptomatic (BAH)|
|Sierra de Almijara, Tejeda y Alhama|*P. sylvestris*| Asymptomatic (ASH)|
||*P. pinaster* | Asymptomatic (APH)|

Thus, we have 11 groups of samples x 2 Seasons x 2 Plant compartments x 2 Organisms x 12 trees per group of samples = 1056 samples.

# Sample sequencing scheme
Total DNA was extracted from rhizosphere and root endosphere samples according to Lasa et al., (2024; https://doi.org/10.1016/j.scitotenv.2024.171858). Due to the high number of samples, DNA had to be sequenced in different MiSeq runs. So, finally, we got the following sequencing scheme:

 - Bacteria + Fungal samples:
	 - Run 1: Samples from rhizosphere, spring 2022.
	 - Run 2: Samples from endosphere, spring 2022.
	 - Run 3: Samples from rhizosphere, summer 2021.
	 - Run 4: samples from endosphere, summer 2021.
	 - Run 5: some samples from endosphere and other samples from root endosphere, summer 2021.
It should be taken into account that in each sequencing run, 3 replicates of a MockCommunity were included.

The Data descriptor manuscript describes the data corresponding to the sub-project that included just ***P. sylvestris*** trees from **Sierra Nevada, Sierra de Baza and Sierra de Almijara, Tejeda y Alhama** (both symptomatic and decline-affected trees, according the Table above) obtained in **spring 2022**. In this project, both **fungal** and **bacterial** communities inhabiting the **rhizosphere** and **root endosphere** of the selected trees were considered.

# Sample sequencing scheme
According to DADA2 developers, when our scientific projects include 2 or more sequencing runs, they should be processed separatedly. Once the Sequence Tables are obtained, all of them have to be merged into a one merged table, and we can go on processing all the samples of the same scientific project from this unique merged table.
In our case, the 528 bacterial samples (+replicates of the MockCommunity) were thus separatedly processed according to the corresponding runs. Then, the 5 Sequence Table were joined into a just one table, all we followed processing that table. The same procedure was applied for fungal dataset.

Thus, all the samples of SUMHAL project were futher processed altogether, although the Data Descriptor **just** included samples corresponding to *P. sylvestris* trees (spring 2022) located in the three selected areas. Accordingly, the BioProject PRJNA993625 (NCBI SRA) included just the fastq files corresponding to this sub-project. In a near future, a second BioProject will be created again in NCBI SRA including the fastq files corresponding to another sub-project of the SUMHAL project.

# What will I find in this GitHub project?

**IMPORTANT**: This GitHub project is divided into 2 different branches: one related to Fungal dataset, and second related to Bacterial reads. In each brand, we can find two different typs of scripts:

 - "script_SequencingReadsProcessing_XXX.R": where all the steps from the fastq files to a filtered ASV table are explained.
 - "script_DataAnalysis_XXX.R": where all the steps from the previously obtained ASV table to the obtaining of the rarefaction curves are explained.
 
XXX: 16S or ITS2.
