###########################################################################
# This program identifies 2Mb hotspots on a chromosome based on QTL density
###########################################################################

###########################################################################
#### Step 1. Load Necessary R Functions and Libraries and Install QHOT
###########################################################################

## Load R functions and libraries and download required files, replace ?????? with specific directory
source('??????/WNV_rix_qtl_mapping_functions_publication.r')
install.packages("QHOT")
library(QHOT)

### All of the files below need to have destination directories added by replacing ??????. ####

# Download allsnps_genesonly_final_infected.rda
download.file(url = "https://figshare.com/ndownloader/files/51509954", destfile = "??????/allsnps_genesonly_final_infected.rda")

###########################################################################
#### Step 2. Create the QTL list to be used in the hotspot analysis
###########################################################################

## Load the master annotated variant file and derive the QTL information, all QTL are used
load("??????/allsnps_genesonly_final_infected.rda")

#QTL that match criteria for QTL hotspot detection
allqtl<-distinct(allsnps_genesonly_final_infected,finalpheno, time, tissue, panel, start.y, end.y, chr,label1,maxlodval.y,finalcat)

###########################################################################
#### Step 3. Run the hotpot analyses for each tissue/timepoint separately. 
#### Panels will be annotated on each analysis output.
###########################################################################

library(QHOT)

#Hotspots by time and tissue combined

#Brain
WNVQTLBraind7<- filter(allqtl, time=="D7" & tissue=='Brain' & !is.na(start.y)) %>%
  dplyr::mutate(.,trait=panel,chr="X",L=start.y,R=end.y) %>%
  dplyr::select(.,trait,chr,L,R)
DataCrop.WNV<-data.frame(chr="X",center=0,length=181)
wnvhs_braind7<-QHOT(WNVQTLBraind7, DataCrop.WNV, ScanStep=2, NH=3, NP=1000)
wnvhs_braind7df<-data.frame(wnvhs_braind7[[1]])

WNVQTLBraind12<- filter(allqtl, time=="D12" & tissue=='Brain' & !is.na(start.y)) %>%
  dplyr::mutate(.,trait=panel,chr="X",L=start.y,R=end.y) %>%
  dplyr::select(.,trait,chr,L,R)
DataCrop.WNV<-data.frame(chr="X",center=0,length=181)
wnvhs_braind12<-QHOT(WNVQTLBraind12, DataCrop.WNV, ScanStep=2, NH=3, NP=1000)
wnvhs_braind12df<-data.frame(wnvhs_braind12[[1]])

WNVQTLBraind21<- filter(allqtl, time=="D21" & tissue=='Brain' & !is.na(start.y)) %>%
  dplyr::mutate(.,trait=panel,chr="X",L=start.y,R=end.y) %>%
  dplyr::select(.,trait,chr,L,R)
DataCrop.WNV<-data.frame(chr="X",center=0,length=181)
wnvhs_braind21<-QHOT(WNVQTLBraind21, DataCrop.WNV, ScanStep=2, NH=3, NP=1000)
wnvhs_braind21df<-data.frame(wnvhs_braind21[[1]])

WNVQTLBraind28<- filter(allqtl, time=="D28" & tissue=='Brain' & !is.na(start.y)) %>%
  dplyr::mutate(.,trait=panel,chr="X",L=start.y,R=end.y) %>%
  dplyr::select(.,trait,chr,L,R)
DataCrop.WNV<-data.frame(chr="X",center=0,length=181)
wnvhs_braind28<-QHOT(WNVQTLBraind28, DataCrop.WNV, ScanStep=2, NH=3, NP=1000)
wnvhs_braind28df<-data.frame(wnvhs_braind28[[1]])

#spleen
WNVQTLSpleend7<- filter(allqtl, time=="D7" & tissue=='Spleen' & !is.na(start.y)) %>%
  dplyr::mutate(.,trait=panel,chr="X",L=start.y,R=end.y) %>%
  dplyr::select(.,trait,chr,L,R)
DataCrop.WNV<-data.frame(chr="X",center=0,length=181)
wnvhs_Spleend7<-QHOT(WNVQTLSpleend7, DataCrop.WNV, ScanStep=2, NH=3, NP=1000)
wnvhs_Spleend7df<-data.frame(wnvhs_Spleend7[[1]])

WNVQTLSpleend12<- filter(allqtl, time=="D12" & tissue=='Spleen' & !is.na(start.y)) %>%
  dplyr::mutate(.,trait=panel,chr="X",L=start.y,R=end.y) %>%
  dplyr::select(.,trait,chr,L,R)
DataCrop.WNV<-data.frame(chr="X",center=0,length=181)
wnvhs_Spleend12<-QHOT(WNVQTLSpleend12, DataCrop.WNV, ScanStep=2, NH=3, NP=1000)
wnvhs_Spleend12df<-data.frame(wnvhs_Spleend12[[1]])

WNVQTLSpleend21<- filter(allqtl, time=="D21" & tissue=='Spleen' & !is.na(start.y)) %>%
  dplyr::mutate(.,trait=panel,chr="X",L=start.y,R=end.y) %>%
  dplyr::select(.,trait,chr,L,R)
DataCrop.WNV<-data.frame(chr="X",center=0,length=181)
wnvhs_Spleend21<-QHOT(WNVQTLSpleend21, DataCrop.WNV, ScanStep=2, NH=3, NP=1000)
wnvhs_Spleend21df<-data.frame(wnvhs_Spleend21[[1]])

WNVQTLSpleend28<- filter(allqtl, time=="D28" & tissue=='Spleen' & !is.na(start.y)) %>%
  dplyr::mutate(.,trait=panel,chr="X",L=start.y,R=end.y) %>%
  dplyr::select(.,trait,chr,L,R)
DataCrop.WNV<-data.frame(chr="X",center=0,length=181)
wnvhs_Spleend28<-QHOT(WNVQTLSpleend28, DataCrop.WNV, ScanStep=2, NH=3, NP=1000)
wnvhs_Spleend28df<-data.frame(wnvhs_Spleend28[[1]])

 



 

 