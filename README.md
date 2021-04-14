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

    - step 2: download the example data (example-Roadmap-data-netopen.RData,example-CMAP-data-PERD.RData) from https://wanglab.shinyapps.io/DeepDRK/, and put it in the PERD folder.

    Dependencies of DeepDRK includes the following: 

    - Readr1.3.1 and all its dependencies;

    - glmnet R package (Version 4.1-1 [glmnet_4.1-1.tar.gz](https://cran.r-project.org/web/packages/glmnet/index.html)) and its dependencies.


2. Preparation of the input files

PERD applies a computational model, NetOpen, to predict the openness value for regulatory elements before/after drug treatment, then treats the regulatory elements that display the significantly different openness value after drug treatment as the drug responsive elements. The training data for build the NetOpen model came from RNA-seq and DNase-seq data in ENCODE 167 cells, which can be download from [ENCODE-training-data](https://github.com/WeiqiangZhou/BIRD-data/releases/download/v3.0/BIRD_data_ENCODE.zip).

              A1BG A1CF  A2M …
    201T       0     0    0
    22RV1      0     1    0
    42-MG-BA   0     0    1
      .
      .
      .

Similarly, the input files that contain the copy number alteration data (CN.csv), the status of DNA methylation (methylation.csv), and the gene expression of the cancer cells (expression.csv) are also data matrices with a row representing a cancer cell line and a column representing a gene. The elements of these matrices are respectively integers for gene copy numbers, and float numbers for level of gene methylation and expression.

expression.csv

               A1BG   A1CF    A2M …
    201T      3.162   2.919   3.379
    22RV1     3.531   6.336   5.331
    42-MG-BA  6.002   3.137   3.237
      .
      .
      .

The input file chem.csv that describes the chemical properties of the cancer drugs is a matrix with each row representing one cancer drug and each column representing one feature to describe drug’s chemical properties. The descriptors of the chemical properties of a cancer drug were inferred from its chemical structure. Particularly, to describe a drug, such as Erlotinib, we will need to first download the sdf file of this drug from PubChem, and then upload the chemical structure into “StarVue” (StarVue-macinstall-1.4.dmg) software to extract the 2D Molecular Operating Environment (MOE)) descriptors, including physical properties, atom counts, and bond counts.

chem.csv

              PUBCHEM_MOLECULAR_WEIGHT   PUBCHEM_EXACT_MASS    PUBCHEM_CACTVS_TPSA …
    Erlotinib            393.4                   393.2                  74.4
    Rapamycin            917.2                   913.6                  195.0
    Sunitinib            398.5                   398.2                  77.2
       .
       .
       .


The input file DT.csv includes the known targeting proteins of the cancer drugs. Each row represents one cancer drug, and each column represents a target protein. “1” indicates a potential drug-gene interaction reported in DrugBank or KEGG.

DT.csv

                         EGFR                     KIT                   PDGRA …
    Erlotinib            393.4                   393.2                  74.4
    Rapamycin            917.2                   913.6                  195.0
    Sunitinib            398.5                   398.2                  77.2
       .
       .
       .

The example of all input files can be found in the “data” folder of the Github repository.

3. Running DeepDRK

The main function of DeepDRK is DeepDRKpredictor.R. Get your input files prepared, and run it like this:

Usage example:

    > cell_tst<-list()
    > library(readr)
    > mutation <- read_csv("~/DeepDRK/data/mutation.csv");A<-data.matrix(mutation[,-1]);rownames(A)<-mutation$X1;cell_tst[[1]]<-A
    > CN <- read_csv("~/DeepDRK/data/CN.csv");A<-data.matrix(CN[,-1]);rownames(A)<- CN$X1;cell_tst[[2]]<-A
    > Methy <- read_csv("~/DeepDRK/data/methylation.csv");A<-data.matrix(Methy[,-1]);rownames(A)<- Methy$X1;cell_tst[[3]]<-A
    > Exp <- read_csv("~/DeepDRK/data/expression.csv");A<-data.matrix(Exp[,-1]);rownames(A)<- Exp$X1;cell_tst[[4]]<-A
    > drug_tst<-list()
    > chem <- read_csv("~/DeepDRK/data/chem.csv");A<-data.matrix(chem[,-1]);rownames(A)<- chem$X1;drug_tst[[1]]<-A
    > DT <- read_csv("~/DeepDRK/data/DT.csv");A<-data.matrix(DT[,-1]);rownames(A)<- DT$X1;drug_tst[[2]]<-A
    > load("~/DeepDRK/combination_data.RData") #load the training RData
    > source('~/DeepDRK/DeepDRKpredictor.R')
    > predictions<-DeepDRKpredictor(cell_tst,drug_tst)
     Are you sure you want to shutdown the H2O instance running at http://localhost:54321/ (Y/N)? y
     TRUE

Moreover, DeepDRK could also handle task with missing features using the DeepDRKpredictor.e R function. Here is the example showing how to use it:

In case the mutation, methylation and target proteins are missing

    > cell_tst<-list()
    > library(readr)
    > CN <- read_csv("~/DeepDRK/data/CN.csv");A<-data.matrix(CN[,-1]);rownames(A)<- CN$X1;cell_tst[[2]]<-A
    > Exp <- read_csv("~/DeepDRK/data/expression.csv");A<-data.matrix(Exp[,-1]);rownames(A)<- Exp$X1;cell_tst[[4]]<-A
    > drug_tst<-list()
    > chem <- read_csv("~/DeepDRK/data/chem.csv");A<-data.matrix(chem[,-1]);rownames(A)<- chem$X1;drug_tst[[1]]<-A
    > drug_tst[[2]]<-matrix()
    > missCtype=c(1,3)
    > missDtype=2
    > load("~/DeepDRK/combination_data.RData") #load the training RData
    > source('~/DeepDRK/DeepDRKpredictor.e.R')
    > predictions<-DeepDRKpredictor.e(cell_tst,drug_tst,missCtype,missDtype)        
    Are you sure you want to shutdown the H2O instance running at http://localhost:54321/ (Y/N)? y
    TRUE


# Contact

For technical issues please send an email to ycwang@nwipb.cas.cn.
