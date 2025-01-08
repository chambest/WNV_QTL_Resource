###########################################################################
# This program identifies 2Mb hotspots on a chromosome based on QTL density
###########################################################################

## Load R functions and libraries
 
source('/Users/chambest/Documents/BACKUP/Shannons Work/Immunogenetics/Publications/QTL Paper/Nature Scientific Data/Workflows for Publication/Final Workflow/Code/WNV_rix_qtl_mapping_functions_publication.r')

load("/Users/chambest/Documents/BACKUP/Shannons Work/Immunogenetics/Publications/QTL Paper/Nature Scientific Data/Workflows for Publication/Final Workflow/Data/QTL/allsnps_genesonly_final_infected.rda")

#QTL that match criteria for QTL hotspot detection
allqtl<-distinct(allsnps_genesonly_final_infected,finalpheno, time, tissue, panel, start.y, end.y, chr,label1,maxlodval.y,finalcat)

####################################################################################################### 
# Statistical QTL hotspot detection, this example uses 356 QTL selected as described in the paper
#######################################################################################################

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

 



 

 