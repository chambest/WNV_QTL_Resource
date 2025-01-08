################################################################################################### 
# This program uses kmeans clustering to group phenotypes together for candidate gene analysis
###################################################################################################

## Load R functions and libraries
#setwd('/Users/chambest/Documents/OHSU/Shannons Work/Immunogenetics/Code')
source('/Users/chambest/Documents/BACKUP/Shannons Work/Immunogenetics/Publications/QTL Paper/Nature Scientific Data/Workflows for Publication/Final Workflow/Code/WNV_rix_qtl_mapping_functions_publication.r')

########################################################################### 
## Plot overlapping intervals for a panel, tissue, timepoint, chromosome of interest
###########################################################################
allqtls<-read.csv('/Users/chambest/Documents/BACKUP/Shannons Work/Immunogenetics/Publications/QTL Paper/Nature Scientific Data/Workflows for Publication/Final Workflow/Data/QTL/Final Annotated QTL.csv')
 
qtls<-filter(allqtls,Panel=='Treg' & Tissue=='Spleen' & Time=='D07' & Chromosome=='X')
minint=min(qtls$Start.Pos)
maxint=max(qtls$End.Pos)

set.seed(20000)
cluster<-qtls[qtls$Peak.LOD.Value>=6,]
 
clust<-ifelse(nrow(cluster)>8,8,ifelse(nrow(cluster)>4 & nrow(cluster)<8,5,2)) 

fit <- kmeans(cluster$Peak.LOD.Position, clust)  

# append cluster assignment
cluster<- data.frame(cluster$Phenotype, fit$cluster)

ints<-left_join(qtls,cluster,by=c("Phenotype"="cluster.Phenotype"))
save(ints,file=' ')
cc <- function(x){ifelse(x==1, "red1",ifelse(x==2,"red2",ifelse(x==3,"red3",ifelse(x==4,"red4",
                                                                                   ifelse(x==5,'red5',ifelse(x==6,"red6",ifelse(x==7,"red7",ifelse(x==8,"red8","green"))))))))}
ff <- function(x){ifelse(x < 6, "grey3","black")}

ints$pheno=paste0(ints$Phenotype,'(c=',ints$fit.cluster,')')

gg<-ggplot(ints[ints$Peak.LOD.Value>=6,]) +
  geom_segment( aes(x=reorder(pheno, Peak.LOD.Position, FUN=median), xend=pheno, y=Start.Pos, yend=End.Pos), color="grey") +
  geom_point( aes(x=pheno, y=Start.Pos, color=cc(fit.cluster)), size=3 ) +
  geom_point( aes(x=pheno, y=End.Pos, color=cc(fit.cluster)),size=3 ) +
  geom_point( aes(x=pheno, y=Peak.LOD.Position, color=ff(Peak.LOD.Value)),size=2) +
  coord_flip() +
  theme_classic() +  
  theme(
    legend.position = "none",
  ) +
  xlab("Pheno (LOD,Pos") +
  ylab("Position B38")  + ggtitle(paste0('Chromosome= ',chr))
print(gg)
 
