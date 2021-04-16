# PERD
# Status

Active development

# Introduction

Regulatory elements are hypothezied to play an important role in individual responses to drug treatment. The relative little knowledge of the location and function of regulatory elements limits our understanding of drug-dependent regulatory elements. PERD is developed to predict the drug responsive regulatory element from their nearby gene and transcript factor transcriptome after drug treatment.

# Usage

1. Installation

   Prerequisites of PERD includes the following: 

   - R is properly installed; 

   - Rscript is available in your system path ($PATH);

   - git (2.21.1)

    Installation of PERD includes the following steps:

    - step 1: git clone https://github.com/wangyc82/PERD;

    - step 2: download the example data (example-Roadmap-data-netopen.RData,example-CMAP-data-PERD.RData) from https://onedrive.live.com/?id=20874D8228EBFF1E%21106&cid=20874D8228EBFF1E, and put it in PERD folder.

    Dependencies of DeepDRK includes the following: 

    - Readr1.3.1 and all its dependencies;

    - glmnet R package (Version 4.1-1 [glmnet_4.1-1.tar.gz](https://cran.r-project.org/web/packages/glmnet/index.html)) and its dependencies.


2. Preparation of the input files

PERD applies a computational model, NetOpen, to predict the openness value for regulatory elements before/after drug treatment, then treats the regulatory elements that display the significantly different openness value after drug treatment as the drug responsive elements. The training data for build the NetOpen model came from RNA-seq and DNase-seq data in ENCODE 167 cells, which can be download from [ENCODE-training-data](https://github.com/WeiqiangZhou/BIRD-data/releases/download/v3.0/BIRD_data_ENCODE.zip). Put the RNA_data and DNase_data in PERD folder. The nearby genes for enhancers came form GeneCard sub-package GeenHancer (GeneHancer_version_4_4), which can be find in PERD repository. The TFs that are located in the enhancer region came from ENCODE Chip-seq data and were summarized in [tfbs_info](https://onedrive.live.com/?id=20874D8228EBFF1E%21106&cid=20874D8228EBFF1E). Put all these files in PERD folder.
    
   - Preparing training data
    
    > DNase_167_cells <- readRDS("~/PERD/ENCODE-training-data/DNase_data_167_cells.rds")
    > RNA_167_cells <- readRDS("~/PERD/ENCODE-training-data/RNA_data_167_cells.rds")
    > library("OrganismDbi")
    > library("org.Hs.eg.db")
    > ensembleIDs<-sapply(1:length(strsplit(rownames(RNA_167_cells),"[.]")),function(x) strsplit(rownames(RNA_167_cells),"[.]")[[x]][1])
    > gene_symbol <- select(org.Hs.eg.db,keys = ensembleIDs,columns = "SYMBOL",keytype = "ENSEMBL")
    > RNAdata_trn<-RNA_167_cells
    > RNAdata.gene<-gene_symbol$SYMBOL[match(ensembleIDs,gene_symbol$ENSEMBL)]
    > rownames(RNAdata_trn)<-RNAdata.gene
    
    > library(readxl)
    > GeneHancer_version_4_4 <- read_excel("~/PERD/GeneHancer_version_4-4.xlsx")
    > A<-strsplit(GeneHancer_version_4_4$attributes,";")
    > enhancer_gene_list<-lapply(1:nrow(GeneHancer_version_4_4),function(x) sapply(seq(2,length(A[[x]]),2),function(y) substr(A[[x]][y],16,nchar(A[[x]][y]))))
    > names(enhancer_gene_list)<-paste(GeneHancer_version_4_4$chrom,paste(GeneHancer_version_4_4$start,GeneHancer_version_4_4$end,sep = "-"),sep = ":")

    > enhancer.withOpen.lab<-lapply(1:nrow(GeneHancer_version_4_4),function(x) which(DNase_167_cells$chromosome == GeneHancer_version_4_4$chrom[x] & DNase_167_cells$start>=GeneHancer_version_4_4$start[x] & DNase_167_cells$end<=GeneHancer_version_4_4$end[x]))
    > len<-sapply(1:length(enhancer.withOpen.lab),function(x) length(enhancer.withOpen.lab[[x]]))
    > enhancer.withOpen.lab1<-enhancer.withOpen.lab[which(len!=0)]
    #prepare the openness for enhancers in GeneHancer
    > DNase_data<-data.matrix(DNase_167_cells[,-c(1,2,3)])
    > rownames(DNase_data)<-paste(DNase_167_cells$chromosome,paste(DNase_167_cells$start,DNase_167_cells$end,sep = "-"),sep = ":")
    > enhancer.withOpen.openness<-NULL
    # calculate the openness for enhancers by find the locus with maximun openness
    for (i in 1:length(enhancer.withOpen.lab1)) {
    m<-DNase_data[enhancer.withOpen.lab1[[i]],]
    if (length(m)==ncol(DNase_data)) {enhancer.withOpen.openness<-rbind(enhancer.withOpen.openness,m)} else {enhancer.withOpen.openness<-rbind(enhancer.withOpen.openness,apply(m,2,max))}
    rm(m)
    #cat(i,"\n")
    }
    > rownames(enhancer.withOpen.openness)<-paste(GeneHancer_info$chrom,paste(GeneHancer_info$start,GeneHancer_info$end,sep = "-"),sep = ":")
    
    # prepare the target gene list
    > enhancer_gene_list_open<-enhancer_gene_list[rownames(enhancer.withOpen.openness)]
    # preparing the binding TF list
    # The example TFBS information cames from ENCODE chip-seq data, which can be obtained in PERD repository
    > library(readr)
    > tfbsInfo <- read_delim("~/PERD/tfbsInfo.txt",  "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
    > tfbs_info<-tfbsInfo[-c(1,2),]
    > TFgene<-sapply(1:length(strsplit(tfbs_info$description," ")),function(x) strsplit(tfbs_info$description," ")[[x]][1])
    > enhancer.withOpen.TFlab<-lapply(1:nrow(enhancer.withOpen.openness),function(x) which(tfbs_info$chrom==GeneHancer_info$chrom[x] & tfbs_info$start>=GeneHancer_info$start[x] & tfbs_info$end<=GeneHancer_info$end[x]))
    > len.t<-sapply(1:length(enhancer.withOpen.TFlab),function(x) length(enhancer.withOpen.TFlab[[x]]))
    > enhancer.withOpen.TFlab1<-enhancer.withOpen.TFlab[which(len.t!=0)]
    > GeneHancer_info1<-GeneHancer_info[which(len.t!=0),]
    > enhancer.withOpen.lab2<-enhancer.withOpen.lab1[which(len.t!=0)]
    > enhancer.withOpen.TG.list<-enhancer_gene_list[which(len.t!=0)]
    > enhancer.withOpen.openness1<-enhancer.withOpen.openness[which(len.t!=0),]
    > enhancer.withOpen.TFgene<-lapply(1:length(enhancer.withOpen.lab2),function(x) unique(TFgene[enhancer.withOpen.TFlab1[[x]]]))
    > enhancer.withOpen.TF.list<-lapply(1:length(enhancer.withOpen.lab2),function(x) intersect(enhancer.withOpen.TFgene[[x]],rownames(RNAdata)))


3. Running PERD

The main function of DeepDRK is DeepDRKpredictor.R. Get your input files prepared, and run it like this:

Usage example:
    
    # using NetOpen to predict the enhancers' openness value
    > load("~/PERD/example-Roadmap-data-netopen.RData") #load the training RData
    > source('~/PERD/netopen.R')
    > preEopen<-netopen(enhancer.withOpen.TG.list,enhancer.withOpen.TF.list,enhancer.withOpen.openness.train2,RNAdata.train,RNAdata.test)
    
    # using PERD to predict drug responsive enhancers
    > load("~/PERD/example-CMAP-data-PERD.RData") #load the training RData
    > source('~/PERD/perd.R')
    > diffEn<-perd(enhancerGeneList_withTGexp_test2,enhancerList_withTGexp_TFexp_test2,enhancer.withOpen.openness.train1,RNAdata,drug_with10_binary_GEmat,drug_with10_binary_instance_info)



# Contact

For technical issues please send an email to ycwang@nwipb.cas.cn.
