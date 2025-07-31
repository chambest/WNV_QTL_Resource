####################################################################### 
# This program runs the RWR for five seed genesets and annotates candidate genes
#######################################################################
setwd()
dir.create('DATA')
#######################################################################
#### Step 1. Load Necessary R Functions and Libraries
#######################################################################
 
## Load R functions and libraries and download required files 
source('WNV_rix_qtl_mapping_functions_publication.r')

# Download allsnps_genesonly_final_infected.rda
download.file(url = "https://figshare.com/ndownloader/files/51509954", destfile = "../DATA/allsnps_genesonly_final_infected.rda")

# Download 10090.protein.links.v11.0.txt
download.file(url = "https://figshare.com/ndownloader/files/51509903", destfile = "./DATA/10090.protein.links.v11.0.txt")

# Download PPI_MGIBatchReport_20210624_181229.xlsx
download.file(url = "https://figshare.com/ndownloader/files/51509900", destfile = "./DATA/PPI_MGIBatchReport_20210624_181229.xlsx")

# Download NIHMS154918-supplement-Table_1_genes.csv
download.file(url = "https://figshare.com/ndownloader/files/51509840", destfile = "./DATA/NIHMS154918-supplement-Table_1_genes.csv")

# Download NIHMS154918-supplement-Table_2_genes.csv
download.file(url = "https://figshare.com/ndownloader/files/51509825", destfile = "./DATA/NIHMS154918-supplement-Table_2_genes.csv")

# Download TCR_Signaling_Pathway.txt
download.file(url = "https://figshare.com/ndownloader/files/51509822", destfile = "./DATA/TCR_Signaling_Pathway.txt")

# Download Interferons.txt
download.file(url = "https://figshare.com/ndownloader/files/51509828", destfile = "./DATA/Interferons.txt")

# Download TNF_Family_Members.txt
download.file(url = "https://figshare.com/ndownloader/files/51509831", destfile = "./DATA/TNF_Family_Members.txt")

#######################################################################
#### Step 2. Load the annotated variant file and derive the list of candidate genes
#######################################################################

## Load the master annotated variant file
load("./DATA/allsnps_genesonly_final_infected.rda")

#### Select high impact genes from the snp file, as defined in the paper

Mb15<-dplyr::filter(allsnps_genesonly_final_infected,(end.y-start.y>0 & end.y-start.y<=15 & finalcat != "Cluster"))

highimpact<-filter(Mb15, (csq1=="High Impact" | csq1=='Moderate Impact' |
                            SIFT_PREDICTION=="DELETERIOUS" | SIFT_PREDICTION=="DELETERIOUS (*WARNING! Low confidence)") &
                     (finalcat=="5" | finalcat=="10")  )

highimpactgenes<-distinct(highimpact,Marker.Symbol, genome.coordinate.start,genome.coordinate.end,chr)
highimpactqtl<-distinct(highimpact,finalpheno, time, tissue, panel, start.y, end.y, chr,label1,maxlodval.y,finalcat)

#######################################################################
#### Step 3. Create the STRING PPI in iGRAPH and get human homologue genes
#######################################################################

#Read the base STRING interaction file and filter based on a physical interaction score>=700
string<-read.delim('./DATA/10090.protein.links.v11.0.txt',stringsAsFactors=F,header=T, sep=" ")
string_real<-string[string$combined_score>700,]

string_real$newprot1<-substr(string_real$protein1,str_locate(string_real$protein1,"10090")[,2]+2,length(string_real$protein1))
string_real$newprot2<-substr(string_real$protein2,str_locate(string_real$protein2,"10090")[,2]+2,length(string_real$protein2))

jax<-read.xls("./DATA/PPI_MGIBatchReport_20210624_181229.xlsx")
jax<-filter(jax,MGI.Gene.Marker.ID != 'No associated gene')


string_real<-inner_join(string_real,dplyr::select(jax,Input,Symbol),by=c("newprot1"="Input"))
string_real<-inner_join(string_real,dplyr::select(jax,Input,Symbol),by=c("newprot2"="Input"))

string_real$Symbol.x<-toupper(string_real$Symbol.x)
string_real$Symbol.y<-toupper(string_real$Symbol.y)

string_uniq1<-distinct(string_real,Symbol.x) %>%
  mutate(.,gene=Symbol.x) %>%
  dplyr::select(.,-Symbol.x)
string_uniq2<-distinct(string_real,Symbol.y) %>%
  mutate(.,gene=Symbol.y) %>%
  dplyr::select(.,-Symbol.y)

uniq<-distinct(rbind(string_uniq1,string_uniq2),gene)

#Create string PPI and perform various calculations
library(igraph)
string_ppi<-graph_from_data_frame(string_real[,c(6,7)],directed=FALSE)
clusters<-components(string_ppi)
clusters$csize
clusters$no

### Get human genes orthologues for the JAX mouse genes symbols
library('biomaRt')
ensembl <- useMart('ensembl', dataset="mmusculus_gene_ensembl")
 
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

genesppi = getLDS(attributes = c("mgi_symbol","ensembl_peptide_id"), filters = "mgi_symbol", values = jax$Symbol , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

#######################################################################
#### Step 4. Annotate candidate gene file with RWR scores for the 305 iRNA seed genes
#######################################################################

library(RandomWalkRestartMH)

# 305 iRNA seed genes for WNV infection involvement
# Annotates the high impact genes with the percentile for the RWR score for this seed set

seeds1<-read.csv("./DATA/NIHMS154918-supplement-Table_1_genes.csv",header=T)
seeds2<-read.csv("./DATA/NIHMS154918-supplement-Table_2_genes.csv",header=T)
seeds1$Symbol <- toupper(seeds1$Symbol)
seeds2$Symbol <- toupper(seeds2$Symbol)

seedstemp1<-inner_join(seeds1,distinct(genesppi,HGNC.symbol,MGI.symbol),by=c("Symbol"="HGNC.symbol"))  
seedstemp2<-inner_join(seeds2,distinct(genesppi,HGNC.symbol,MGI.symbol),by=c("Symbol"="HGNC.symbol"))  

seedstemp1$MGI.symbol<-toupper(seedstemp1$MGI.symbol)
seedstemp2$MGI.symbol<-toupper(seedstemp2$MGI.symbol)

SeedGene1<-toupper(inner_join(seedstemp1,uniq,by=c("MGI.symbol"="gene"))$MGI.symbol) 
SeedGene2<-toupper(inner_join(seedstemp2,uniq,by=c("MGI.symbol"="gene"))$MGI.symbol) 

SeedGene<-c(SeedGene1,SeedGene2)
SeedGene_305<-SeedGene

PPI_Monoplex_Obj <- create.multiplex(list(string_ppi),Layers_Name=c("PPI"))
PPI_Monoplex_Obj

AdjMatrix_PPI <- compute.adjacency.matrix(PPI_Monoplex_Obj)
AdjMatrixNorm_PPI <- normalize.multiplex.adjacency(AdjMatrix_PPI)

RWR_PPI_Results_305 <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI, PPI_Monoplex_Obj,SeedGene)

seeddf<-data.frame(NodeNames=SeedGene,Score=rep(1,times=length(SeedGene)))
df305<-rbind(RWR_PPI_Results_305$RWRM_Results,seeddf)

df305$Score=df305$Score+1.0e-15

highimpactgenes$Marker.Symbol<-toupper(highimpactgenes$Marker.Symbol)
highimpactgenes<-left_join(highimpactgenes,dplyr::select(df305,NodeNames,Score),by=c("Marker.Symbol"="NodeNames"))
 
highimpactgenes$pct305<-ecdf(df305$Score)(highimpactgenes$Score)  
df305$pct305<-ecdf(df305$Score)(df305$Score) 

highimpactgenes$Score305<-highimpactgenes$Score
highimpactgenes<-dplyr::select(highimpactgenes,-Score)

#######################################################################
#### Step 5. Annotate candidate gene file with RWR scores for GeneWeaver West Nile virus seed genes.
#######################################################################

# This is the geneweaver geneset for WNV
# Annotates the high impact genes with the percentile for the RWR score for this seed set
 
SeedGenedf <- data.frame(gene=c("BCL2L12",
              "CCR5",
              "CXCL10",
              "CXCR3",
              "DDX58",
              "EIF2AK2",
              "IFIT1",
              "IFIT2",
              "IFNAR1",
              "IFNB1",
              "IFNG",
              "IRF3",
              "IRF7",
              "MAVS",
              "MMP9",
              "MYD88",
              "OAS1",
              "PRF1",
              "RAG1",
              "TLR7",
              "TNF"))

seedstemp<-inner_join(SeedGenedf,distinct(genesppi,HGNC.symbol,MGI.symbol),by=c("gene"="HGNC.symbol"))
seedstemp$MGI.symbol<-toupper(seedstemp$MGI.symbol)
SeedGene<-toupper(inner_join(seedstemp,uniq,by=c("MGI.symbol"="gene"))$MGI.symbol)
SeedGene_GW<-SeedGene

RWR_PPI_Results_DW <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI, PPI_Monoplex_Obj,SeedGene)
 
seeddf<-data.frame(NodeNames=SeedGene,Score=rep(1,times=length(SeedGene)))
dfDW<-rbind(RWR_PPI_Results_DW$RWRM_Results,seeddf)

dfDW$Score=dfDW$Score+1.0e-15

highimpactgenes<-left_join(highimpactgenes,dplyr::select(dfDW,NodeNames,Score),by=c("Marker.Symbol"="NodeNames"))

highimpactgenes$pctDW<-ecdf(dfDW$Score)(highimpactgenes$Score)  
dfDW$pctDW<-ecdf(dfDW$Score)(dfDW$Score)  

highimpactgenes$ScoreGW<-highimpactgenes$Score
highimpactgenes<-dplyr::select(highimpactgenes,-Score)
 

#######################################################################
#### Step 6. Annotate candidate gene file with RWR scores for the t-cell 
#### receptor signaling pathway seed genes from Immport
#######################################################################

# Immport random walk for t-cell receptor signaling pathway genes
# Annotates the high impact genes with the percentile for the RWR score for this seed set

tcr<-read.table('./DATA/TCR_Signaling_Pathway.txt',header=T,sep="\t") 
tcrgenes<-dplyr::distinct(tcr,Symbol)

seedstemp<-inner_join(tcrgenes,distinct(genesppi,HGNC.symbol,MGI.symbol),by=c("Symbol"="HGNC.symbol"))
seedstemp$MGI.symbol<-toupper(seedstemp$MGI.symbol)
SeedGene<-toupper(inner_join(seedstemp,uniq,by=c("MGI.symbol"="gene"))$MGI.symbol)
SeedGene_tcr<-SeedGene

RWR_PPI_Results_tcr <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI, PPI_Monoplex_Obj,SeedGene)

seeddf<-data.frame(NodeNames=SeedGene,Score=rep(1,times=length(SeedGene)))
dftcr<-rbind(RWR_PPI_Results_tcr$RWRM_Results,seeddf)

dftcr$Score=dftcr$Score+1.0e-15

highimpactgenes<-left_join(highimpactgenes,dplyr::select(dftcr,NodeNames,Score),by=c("Marker.Symbol"="NodeNames"))

highimpactgenes$pcttcr<-ecdf(dftcr$Score)(highimpactgenes$Score)  
dftcr$pcttcr<-ecdf(dftcr$Score)(dftcr$Score) 

highimpactgenes$Scoretcr<-highimpactgenes$Score
highimpactgenes<-dplyr::select(highimpactgenes,-Score)

#######################################################################
#### Step 7. Annotate candidate gene file with RWR scores for the interferon seed genes from Immport.
#######################################################################

# Immport random walk for ifn genes
# Annotates the high impact genes with the percentile for the RWR score for this seed set

ifn<-read.table('./DATA/Interferons.txt',header=T,sep="\t") 
ifngenes<-dplyr::distinct(ifn,Symbol)

seedstemp<-inner_join(ifngenes,distinct(genesppi,HGNC.symbol,MGI.symbol),by=c("Symbol"="HGNC.symbol"))
seedstemp$MGI.symbol<-toupper(seedstemp$MGI.symbol)
SeedGene<-toupper(inner_join(seedstemp,uniq,by=c("MGI.symbol"="gene"))$MGI.symbol)
SeedGene_ifn<-SeedGene

RWR_PPI_Results_ifn <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI, PPI_Monoplex_Obj,SeedGene)

seeddf<-data.frame(NodeNames=SeedGene,Score=rep(1,times=length(SeedGene)))
dfifn<-rbind(RWR_PPI_Results_ifn$RWRM_Results,seeddf)

dfifn$Score=dfifn$Score+1.0e-15

highimpactgenes<-left_join(highimpactgenes,dplyr::select(dfifn,NodeNames,Score),by=c("Marker.Symbol"="NodeNames"))

highimpactgenes$pctifn<-ecdf(dfifn$Score)(highimpactgenes$Score)  
dfifn$pctifn<-ecdf(dfifn$Score)(dfifn$Score)  
 
highimpactgenes$Scoreifn<-highimpactgenes$Score
highimpactgenes<-dplyr::select(highimpactgenes,-Score)

#######################################################################
#### Step 8. Annotate candidate gene file with RWR scores for the tumor 
#### necrosis factor seed genes from Immport.
#######################################################################

# Immport random walk for tfn genes
# Annotates the high impact genes with the percentile for the RWR score for this seed set

tnf<-read.table('./DATA/TNF_Family_Members.txt',header=T,sep="\t") 
tnfgenes<-dplyr::distinct(tnf,Symbol)


seedstemp<-inner_join(tnfgenes,distinct(genesppi,HGNC.symbol,MGI.symbol),by=c("Symbol"="HGNC.symbol"))
seedstemp$MGI.symbol<-toupper(seedstemp$MGI.symbol)
SeedGene<-toupper(inner_join(seedstemp,uniq,by=c("MGI.symbol"="gene"))$MGI.symbol)
SeedGene_tnf<-SeedGene

RWR_PPI_Results_tnf <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI, PPI_Monoplex_Obj,SeedGene)

seeddf<-data.frame(NodeNames=SeedGene,Score=rep(1,times=length(SeedGene)))
dftnf<-rbind(RWR_PPI_Results_tnf$RWRM_Results,seeddf)

dftnf$Score=dftnf$Score+1.0e-15

highimpactgenes<-left_join(highimpactgenes,dplyr::select(dftnf,NodeNames,Score),by=c("Marker.Symbol"="NodeNames"))

highimpactgenes$pcttnf<-ecdf(dftnf$Score)(highimpactgenes$Score)  
dftnf$pcttnf<-ecdf(dftnf$Score)(dftnf$Score)  

highimpactgenes$Scoretnf<-highimpactgenes$Score
highimpactgenes<-dplyr::select(highimpactgenes,-Score)

 
 

 