## Functions to be used for the WNV QTL gene candidate processing

#The DOQTL package is the main package needed to run the QTL in the project
#library(htmltools)
## Need to install DOQTL from GITHUB
#devtools::install_github('dmgatti/DOQTL')
#library(matrixStats)
library(DOQTL)

#The following is used for checking overlap of QTL and gene coordinates
library(GenomicRanges)

#The following is used for any graphics produced in these functions
library(RColorBrewer)

#The following are used for the various data wrangling and graphics tasks done in these functions
library(plyr)  
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(gplots)
library(ggbeeswarm)
library(cluster)    
library(stringi) 
library(R.utils)
library(textclean)

showCols <- function(cl=colors(), bg = "grey",
                     cex = 0.75, rot = 30) {
  m <- ceiling(sqrt(n <-length(cl)))
  length(cl) <- m*m; cm <- matrix(cl, m)
  require("grid")
  grid.newpage(); vp <- viewport(w = .92, h = .92)
  grid.rect(gp=gpar(fill=bg))
  grid.text(cm, x = col(cm)/m, y = rev(row(cm))/m, rot = rot,
            vp=vp, gp=gpar(cex = cex, col = cm))
}
 
FoundProbAssocPlotsComplexShort_Publication<-function(phenocol,untransphenocol,sex,chr,start,end,qtl){
  print('in founder report loop')
  rm(tcoef, tmp, ma, sub_grp, g1,g2,coef,final_snps_phenos_all,tmpmrg)
  outputdir='../Data/Output/'

  #make parameters global for the boxplots and final candidate selection
  phenocol<<-phenocol
  untransphenocol<<-untransphenocol
  minmark<-get_min_marker(qtl, chr=chr)  
  
  #################################################################################################
  # Run assoc plot
  #################################################################################################
  #chr='X'
  
  print("assocplot")
  ma = assoc.map(pheno = pheno, 
                 pheno.col = phenocol, 
                 probs = model.probs, 
                 K = K,
                 addcovar = covar[, c('sex', "Oas1b_High", "Oas1b_Mod"), drop=F], 
                 snps = marker_pos, 
                 chr = chr,
                 start=start, 
                 end=end)
  tmp<-ma[ma$LOD>=3,]
  #tmp <<- print(assoc.plot(ma, thr = 3))
  
  ##############################################
  ##Identify Gene Candidates, Create Boxplot
  ############################################## 
  
  #Reformat the coeficient data
  colors = do.colors
  num.founders = nrow(colors)
  if (chr == "X") {
    lod = qtl$lod$X
    coef = qtl$coef$X
    if (sex == "F") {
      columns = match(paste("F", colors[, 1], sep = "."), 
                      colnames(coef))
    }
    else {
      columns = match(paste("F", colors[, 1], sep = "."), 
                      colnames(coef))
    }
    columns = columns[!is.na(columns)]
    coef = coef[, c(1, columns)]
    colnames(coef)[1] = "A"
    colnames(coef) = sub("^[MF]\\.", "", colnames(coef))
    coef = coef[rownames(coef) %in% lod[, 1], ]
    coef[, 2:ncol(coef)] = coef[, 2:ncol(coef)] + coef[,1]
    coef = coef - rowMeans(coef)
  }
  if (chr != "X") {
    lod = qtl$lod$A
    lod = lod[lod[, 2] == chr, ]
    intercept = qtl$coef$A[, 1]
    coef = qtl$coef$A[, (ncol(qtl$coef$A) - num.founders + 
                                   1):ncol(qtl$coef$A)]
    coef[, 1] = intercept
    colnames(coef)[1] = "A"
    coef = coef[rownames(coef) %in% lod[, 1], ]
    coef[, 2:ncol(coef)] = coef[, 2:ncol(coef)] + coef[,1]
    coef = coef - rowMeans(coef)
  }
  coef<-data.frame(coef)
  coef$marker<-rownames(coef)
  
  if(chr == "X"){coef<-inner_join(coef,marker_pos[marker_pos$position.B38.>=(start*1000000) & marker_pos$position.B38.<=(end*1000000) & marker_pos$chromosome=="X",],by="marker")}
  if(chr != "X"){coef<-inner_join(coef,marker_pos[marker_pos$position.B38.>=(start*1000000) & marker_pos$position.B38.<=(end*1000000) & marker_pos$chromosome==as.character(chr),],by="marker")}
  
  tcoef<-t(dplyr::select(coef,A,B,C,D,E,F,G,H))
  d <- dist(tcoef, method = "euclidean")
  hc1 <- hclust(d, method = "average" )
  sub_grp <-cutree(hc1, k = 2)
  
  print("cluster plot")         
  plot(hc1, cex = 0.6)
  rect.hclust(hc1, k = 2, border = 2:5)
 
  tcoef<<-tcoef
  sub_grp<<-sub_grp
  
  #Find founder variant patterns for LOD>3 snps that match the clustering
  
  tmp$posm<-tmp$pos*1000000
  
  ## This section retrieves SNP classiications for those returned by the ASSOC function
  
  # This is specifically for SNPs only on the X chromosone, the PATH variable needs to be set here
  Sys.setenv(PATH="?????")
  write.table(tmp[,c('chr','posm')],file="../Data/Output/tmp.txt",sep="\t",col.names=FALSE,quote=FALSE,row.names = FALSE)
  getsnps<-"tabix -f '../Data/Sanger/mgp.v5.merged.snps_all.dbSNP142_chrX.recode.vcf.gz' -R '../Data/Output/tmp.txt' | cat > '../Data/Output/tmp_matched.txt'"
  system(getsnps)
  
  c<-"X"
  region_variant <- read.delim(paste0('../Data/Output/tmp_matched.txt'),stringsAsFactors=F,header=F)  
  header<-paste0("tabix -H -f '../Data/Sanger/mgp.v5.merged.snps_all.dbSNP142_chrX.recode.vcf.gz'  > '",outputdir,"/header.out'")
  system(header)
  header<-read.delim(file.path(outputdir,paste0("header.out")),stringsAsFactors=F,skip=103)
    
  colnames(region_variant)<-colnames(header)  
  #reformat snp file
  short_genotype<-apply(region_variant[,10:ncol(region_variant)],2,setGenotypeShort)
  missense_indels<-data.frame(region_variant, short_genotype)
  colnames(missense_indels)<-gsub("_","/",colnames(missense_indels))
  f_interest<-c("CAST/EiJ","PWK/PhJ","WSB/EiJ","C57BL/6J","NZO/HlLtJ","129S1/SvImJ","A/J","NOD/ShiLtJ")
  f_interest[which(f_interest=="129S1/SvImJ")]<-"X129S1/SvImJ"
  f_idx<-unlist(lapply(f_interest,function(x,y) { return (grep(paste("^",x,sep=""),colnames(y)))},missense_indels))
  missense_indels_trunc<<-missense_indels[,c(1:9,f_idx)]
  
  final_snps_phenos_all<-left_join(tmp,missense_indels_trunc,by=c("chr"="X.CHROM","posm"="POS"))  
  final_snps_phenos_all$csq1='N'
  
  print("Processing snps")         
  variants<-c("CAST/EiJ.1","PWK/PhJ.1","WSB/EiJ.1","NZO/HlLtJ.1","X129S1/SvImJ.1","A/J.1","NOD/ShiLtJ.1")
  variants2<-c("CAST_EiJ","C57BL_6J","PWK_PhJ","WSB_EiJ","NZO_HlLtJ","129S1_SvImJ","A_J","NOD_ShiLtJ")
  
  if(nrow(final_snps_phenos_all)>0){
    for (i in 1:nrow(final_snps_phenos_all)){ 
      #i=1
      reffounders<- vector(mode="character", length=8)
      altfounders<- vector(mode="character", length=8)
      ridx=1
      aidx=0
      for (j in variants2){ 
        if(final_snps_phenos_all[i,j]==0){
          ridx=ridx+1
          if(j=="C57BL_6J"){reffounders[ridx]="B"}
          if(j=="A_J"){reffounders[ridx]="A"}
          if(j=="129S1_SvImJ"){reffounders[ridx]="C"}
          if(j=="NOD_ShiLtJ"){reffounders[ridx]="D"}
          if(j=="NZO_HlLtJ"){reffounders[ridx]="E"}
          if(j=="CAST_EiJ"){reffounders[ridx]="F"}
          if(j=="PWK_PhJ"){reffounders[ridx]="G"}
          if(j=="WSB_EiJ"){reffounders[ridx]="H"}
        }
        if(final_snps_phenos_all[i,j]==1){
          aidx=aidx+1
          if(j=="C57BL_6J"){altfounders[aidx]="B"}
          if(j=="A_J"){altfounders[aidx]="A"}
          if(j=="129S1_SvImJ"){altfounders[aidx]="C"}
          if(j=="NOD_ShiLtJ"){altfounders[aidx]="D"}
          if(j=="NZO_HlLtJ"){altfounders[aidx]="E"}
          if(j=="CAST_EiJ"){altfounders[aidx]="F"}
          if(j=="PWK_PhJ"){altfounders[aidx]="G"}
          if(j=="WSB_EiJ"){altfounders[aidx]="H"}
        }
      }
      #print("assigning alt and ref")
      cluster1=names(sub_grp)[sub_grp==1]
      cluster2=names(sub_grp)[sub_grp==2]
      altfounders=altfounders[altfounders !=""]
      alt<-altfounders[order(altfounders)]
      reffounders=reffounders[reffounders !=""]
      ref<-reffounders[order(reffounders)]
      
      #print("creating labels")
      label1=""
      label2=""
      if((all(alt==cluster1)==TRUE & all(ref==cluster2)==TRUE) | (all(alt==cluster2)==TRUE & all(ref==cluster1)==TRUE)) {
        label1=alt
        label2=ref
      }
      
      final_snps_phenos_all[i,"label1"]=paste0(label1,collapse="")
      final_snps_phenos_all[i,"label2"]=paste0(label2,collapse="") 
      final_snps_phenos_all[i,"n_grp1"]=0
      final_snps_phenos_all[i,"n_grp2"]=0
      
      high_impact<-c("transcript_ablation","splice_acceptor_variant","splice_donor_variant","stop_gained","frameshift_variant","stop_lost","transcript_amplification")
      moderate_impact<-c("inframe_insertion","inframe_deletion", "missense_variant","regulatory_region_ablation")
      csqs<-paste(c(high_impact,moderate_impact))
      found=0
      for(n in csqs){
        #print(n)
        if(grepl(n,final_snps_phenos_all$INFO[i],fixed=TRUE)==TRUE){
          found=found+1
          #print("found") 
          if(found==1){final_snps_phenos_all[i,"csq1"]=n} 
          if(found==2){final_snps_phenos_all[i,"csq2"]=n} 
          if(found==3){final_snps_phenos_all[i,"csq3"]=n}  
          if(found==4){final_snps_phenos_all[i,"csq4"]=n} 
          if(found>4){print('high csq')}
        }
      }
      if(final_snps_phenos_all[i,"csq1"]=='N'){final_snps_phenos_all[i,"csq1"]="LOW"}
       
      if(!(label1[1]=="")){  
        if(length(label1)==1){finalcrit=paste0("data.frame(model.probs[,,minmark[1,1]])[,'",label1,"']>=.4")}
        if(length(label1)>1){ 
          finalcrit=""
          for(n in label1){
            selcrit=""
            if(n != label1[length(label1)]){selcrit<-paste0("data.frame(model.probs[,,minmark[1,1]])[,'",n,"']>=.4 | ")}
            if(n == label1[length(label1)]){selcrit<-paste0("data.frame(model.probs[,,minmark[1,1]])[,'",n,"']>=.4")}
            finalcrit=paste0(finalcrit,selcrit)
          }   
        }
        g2code=paste0("g2<-data.frame(model.probs[,,minmark[1,1]])[!(",finalcrit,"),]")
        g1code=paste0("g1<-data.frame(model.probs[,,minmark[1,1]])[(",finalcrit,"),]")
        eval(parse(text=g1code))
        eval(parse(text=g2code))
        
        if(nrow(g1)>0 & nrow(g2)>0){ 
           
          g1$id=rownames(g1)
          g2$id=rownames(g2)
          nrow(g1)
          pheno_g1<-inner_join(pheno,g1,by=c("ID"="id")) 
          pheno_g2<-inner_join(pheno,g2,by=c("ID"="id"))
          
          pheno_g1$Pheno=pheno_g1[,phenocol]
          pheno_g2$Pheno=pheno_g2[,phenocol]
          
          pheno_g1$HapGroup=paste0(label1,collapse="")
          pheno_g2$HapGroup=paste0(label2,collapse="")
          final_snps_phenos_all[i,"n_grp1"]=nrow(pheno_g1)
          final_snps_phenos_all[i,"n_grp2"]=nrow(pheno_g2)
          pheno_g1<<-pheno_g1
          pheno_g2<<-pheno_g2
        }
      }
    }
  }
   
  print("finished snp processing")
  if(exists("g1") & exists("g2")){ 
    if(nrow(g1)>0 & nrow(g2)>0){
      for(i in 1:nrow(g1)){
        maxval=0
        Max_Founder_Prob=""
        for(n in c("A","B","C","D","E","F","G","H")){
          #print(g1[i,n])
          if(g1[i,n]>maxval){
            maxval=g1[i,n]
            Max_Founder_Prob=n
          }
        }
        # print(max)
        g1[i,'Max_Founder_Prob']=Max_Founder_Prob
      }
      for(i in 1:nrow(g2)){
        maxval=0
        Max_Founder_Prob=""
        for(n in c("A","B","C","D","E","F","G","H")){
          #print(g2[i,n])
          if(g2[i,n]>maxval){
            maxval=g2[i,n]
            Max_Founder_Prob=n
          }
        }
        # print(max)
        g2[i,'Max_Founder_Prob']=Max_Founder_Prob
      }
      pheno_g1<-inner_join(pheno_g1,g1,by=c("ID"="id")) 
      pheno_g2<-inner_join(pheno_g2,g2,by=c("ID"="id"))
      phenofinal<-rbind(dplyr::select(pheno_g1,HapGroup,Pheno,Max_Founder_Prob),dplyr::select(pheno_g2,HapGroup,Pheno,Max_Founder_Prob))
    }}
  
  # Load the protein coding gene coordinates
  annotationdir="../Data/Annotation"
  mgi_annotation<-read.delim(file.path(annotationdir,"mgi_annotation.rpt"),stringsAsFactors=F)
  uniprot_annotation<-read.delim(file.path(annotationdir,"MOUSE_10090_idmapping.dat"),stringsAsFactors=F,header=F)
  
  colnames(uniprot_annotation)<-c("Uniprot","ID_type","ID")
  uniprot_annotation_gene<-uniprot_annotation[which(uniprot_annotation[,"ID_type"]=="Gene_Name"),]
  uniprot_annotation_mgi<-uniprot_annotation[which(uniprot_annotation[,"ID_type"]=="MGI"),]
  uniprot_annotation_gene_mgi<-merge(uniprot_annotation_gene[,c("Uniprot","ID")],uniprot_annotation_mgi[,c("Uniprot","ID")],suffix=c("_GeneName","_MGI"),by="Uniprot")
  uniprot_annotation_gene_mgi<-distinct(uniprot_annotation_gene_mgi,ID_MGI, ID_GeneName)
  
  region_interest_mgi<-mgi_annotation[mgi_annotation$Feature.Type=="protein coding gene",]
  region_interest_mgi_uniprot<<-merge(region_interest_mgi,uniprot_annotation_gene_mgi,by.x="MGI.Accession.ID",by.y="ID_MGI",all.x=T)
  clean_genes<-distinct(region_interest_mgi_uniprot[!is.na(region_interest_mgi_uniprot$genome.coordinate.start) &
                                                      !is.na(region_interest_mgi_uniprot$genome.coordinate.start),],Marker.Symbol,.keep_all = TRUE)
  
  #Create iRanges objects to identify overlap between these snps and gene regions
  usnps<-distinct(data.frame(final_snps_phenos_all),posm,chr) 
  
  chrs<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X")
  
  allsnpsgenes<-data.frame(matrix(ncol = 8, nrow = 0))
  genecols<-c("genome.coordinate.start","genome.coordinate.end","Marker.Symbol","Marker.Name","ID_GeneName",
              "Feature.Type")
  colnames(allsnpsgenes)<-c("chr","posm", "genome.coordinate.start","genome.coordinate.end","Marker.Symbol","Marker.Name","ID_GeneName",
                            "Feature.Type")
  
  for(c in chrs){
    rm(temp,igenes,isnps,clean_genes_tmp, usnps_tmp,ol,mol)
    clean_genes_tmp<-clean_genes[clean_genes$Chr==c,]
    usnps_tmp<-usnps[usnps$chr==c,]
    igenes<-IRanges(start=clean_genes_tmp$genome.coordinate.start,end=clean_genes_tmp$genome.coordinate.end)
    isnps<-IRanges(start=usnps_tmp$posm, end=usnps_tmp$posm)
    ol<-findOverlaps(igenes, isnps)
    mol<-data.frame(as.matrix(ol))
    temp<-cbind(usnps_tmp[mol$subjectHits,], clean_genes_tmp[mol$queryHits,genecols])
    allsnpsgenes<-rbind(allsnpsgenes,temp)
    
  }
  
  tmpmrg <- left_join(data.frame(final_snps_phenos_all),data.frame(allsnpsgenes),by=c('posm','chr'))  %>%
            mutate(.,ingene=ifelse(is.na(Marker.Symbol),"N","Y")) 
  final_snps_phenos_all<-tmpmrg
  
  final_snps_phenos_all$Marker.Symbol[final_snps_phenos_all$ingene=="N"]=""
  final_snps_phenos_all$Feature.Type[final_snps_phenos_all$ingene=="N"]=""
  final_snps_phenos_all$genome.coordinate.end[final_snps_phenos_all$ingene=="N"]=""
  final_snps_phenos_all$genome.coordinate.start[final_snps_phenos_all$ingene=="N"]=""
  
  #nrow(distinct(final_snps_phenos_check,posm))
  
  rm(candrep)
  print("PC Plot")   
  print(fviz_cluster(list(data = tcoef, cluster = sub_grp)))
  
  genes=""
  geneshimod=""
  if(nrow(filter(final_snps_phenos_all,label1 != "")) >0){
    if(nrow(filter(final_snps_phenos_all,n_grp1>=5 & n_grp2>=5 & !is.na(Marker.Symbol) & Feature.Type=="protein coding gene")) >0){
      print("Gene candidates found")
      boxplottitle=paste0(phenoclass, ", Prob>=.4, ",untransphenocol)
      candrep<-distinct(filter(final_snps_phenos_all,n_grp1>=5 & n_grp2>=5 & !is.na(Marker.Symbol)),Marker.Symbol)
      candrephimod<-distinct(filter(final_snps_phenos_all,n_grp1>=5 & n_grp2>=5 & !is.na(Marker.Symbol) & csq1 != "LOW"),Marker.Symbol)
      genes=""
      geneshimod=""
      if(nrow(candrep)>0){for(i in 1:nrow(candrep)){genes=paste(genes,candrep[i,])}}
      if(nrow(candrephimod)>0){for(i in 1:nrow(candrephimod)){geneshimod=paste(geneshimod,candrephimod[i,])}}
      print("before boxplot")
      print(ggplot(phenofinal, aes(x=HapGroup, y=Pheno, fill=HapGroup)) +
              geom_boxplot() +
              geom_jitter(aes(color=Max_Founder_Prob), size=1.2, alpha=0.9) +
              
              theme(
                legend.position="bottom",
                plot.title = element_text(size=11),
                plot.subtitle = element_text(size=8) 
              ) + 
              ggtitle(boxplottitle, subtitle=paste0('Candidate genes= ',genes,"\n High/Mod Csq= ",geneshimod) ) +
              
              xlab(""))
     }}
  
  print("Assoc plot final") 
  ma=ma[!is.na(ma[,12]),]
  reportgenes=""
  if(exists("candrep")){reportgenes=candrep$Marker.Symbol} 
  #print(candrep$ID_GeneName) 
  tmp = assoc.plot(ma, thr = 3, highlight=reportgenes, highlight.col='red', show.sdps = T)
  save(final_snps_phenos_all, qtl, minmark, tcoef, sub_grp,file =paste0(outputdir,"/QTL_Mapping_FounderGrouping_",untransphenocol,".rda") )
}

setGenotypeShort<-function(x) {
  y<-rep("Alt",length(x))
  y[grep("0/0",x)]<-"Ref"
  y[grep("\\./\\.",x)]<-"No Call"
  return(y)
}

get_min_marker = function(doqtl, chr=NULL, window=0) {
  if (is.null(chr)) {
    lod = doqtl$lod$A
  } else {
    if (chr == 'X') {
      lod = doqtl$lod$X
    } else {
      lod = doqtl$lod$A[doqtl$lod$A[,2]==chr, ]
    }
  }
  min_idx = which.min(lod$p)
  lod[(min_idx-window):(min_idx+window),]
}
attr(get_min_marker, 'help') = "
This function returns the marker with the minimum p-value, either across the 
autosomes (if chr==NULL), or within a specific chromosome. Adjacent markers can 
also be returned.

Parameters:
doqtl: A QTL object returned by scanone()
chr: A chromosome (1-19, or 'X')
window: An integer specifying how many markers on either side of the minimum marker
should be returned.

Returns:
A dataframe containing the minimum marker (and adjacent markers).   
"

 
scanone.perm_v2 = function(pheno, pheno.col, probs, K, addcovar, snps, value='lod', B=100, parallel=F, ...) {
  
  ## Check that samples and matings match across pheno, K, and probs
  if (!identical(rownames(pheno), rownames(K)) | !identical(rownames(pheno), dimnames(probs)[[1]])) {
    stop("Sample names do not match!")
  }
  
  ## Get sample IDs
  sample_ids = rownames(pheno)
  
  ## Get matings and check that they match IDs
  if ('Mating' %in% colnames(pheno)) {
    sample_matings = pheno$Mating
    if (!identical(sample_matings, sapply(strsplit(sample_ids, "_"), function(x) {x[1]}))) {
      stop("Matings in pheno do not match sample IDs!")
    }
  } else {
    sample_matings = sapply(strsplit(sample_ids, "_"), function(x) {x[1]})
  }
  
  unique_matings = unique(sample_matings)
  matings_table = table(sample_matings)[unique_matings]
  
  max_vals = c()
  
  ## Assign matings as names for K and probability array so they can be subset easily during permutations
  dimnames(probs)[[1]] = sample_matings
  rownames(K) = sample_matings
  colnames(K) = sample_matings
  
  if (!parallel) {
    for (i in 1:B) {
      ## Print progress
      if (i==1) {
        cat('\n Permutation:', i)
        flush.console()
      } else {
        cat('\r','Permutation:', i)
        flush.console()
      }
      
      ## Permute population and create new IDs
      new_unique_matings = sample(unique_matings)
      new_matings = rep(new_unique_matings, matings_table)
      new_matings_table = table(new_matings)[new_unique_matings]
      new_sample_ids = unname(unlist(sapply(names(new_matings_table), function(x) {paste(x, 1:new_matings_table[x], sep='_')})))
      
      ## Permute probability array and assign new IDs
      new_probs = probs[new_matings, , ]
      dimnames(new_probs)[[1]] = new_sample_ids
      
      ## Permute K matrix and assign new IDs
      new_K = K[new_matings, new_matings]
      rownames(new_K) = new_sample_ids
      colnames(new_K) = new_sample_ids
      
      ## Assign new sample IDs to pheno dataframe
      new_pheno = pheno
      rownames(new_pheno) = new_sample_ids
      
      ## Assign new sample IDs to covariate dataframe
      new_addcovar = addcovar
      rownames(new_addcovar) = new_sample_ids
      
      ## Perform QTL scan
      msg = capture.output({qtl_perm = scanone(pheno=new_pheno, pheno.col=pheno.col, probs=new_probs, K=new_K, addcovar=new_addcovar, snps=snps)})
      
      if (value == "lod") {
        max_A = max(qtl_perm$lod$A[,'lod'])
        if ('X' %in% names(qtl_perm$lod)) {
          max_X = max(qtl_perm$lod$X[,'lod'])
        } else {
          max_X = NA
        }
        perm_max_val = max(max_A, max_X, na.rm=T)
      } else {
        max_A = min(qtl_perm$lod$A[,'p'])
        if ('X' %in% names(qtl_perm$lod)) {
          max_X = min(qtl_perm$lod$X[,'p'])
        } else {
          max_X = NA
        }
        perm_max_val = min(max_A, max_X, na.rm=T)
      }
      max_vals = c(max_vals, perm_max_val)
      gc(verbose=F)
    }
  } else {
    if (is.numeric(parallel)) {
      num_cores = parallel
    } else {
      num_cores = detectCores() - 1
    }
    registerDoParallel(num_cores)
    cat('\nNumber of parallel processes:', num_cores)
    max_vals = foreach(i=1:B, .combine=c, ...) %dopar% {
      gc(verbose=F)
      
      ## Print progress
      if (i==1) {
        cat('\n Permutation:', i)
        flush.console()
      } else if (i %% num_cores == 0) {
        cat('\r','Permutation:', i)
        flush.console()
      }
      
      ## Permute population and create new IDs
      new_unique_matings = sample(unique_matings)
      new_matings = rep(new_unique_matings, matings_table)
      new_matings_table = table(new_matings)[new_unique_matings]
      new_sample_ids = unname(unlist(sapply(names(new_matings_table), function(x) {paste(x, 1:new_matings_table[x], sep='_')})))
      
      ## Permute probability array and assign new IDs
      new_probs = probs[new_matings, , ]
      dimnames(new_probs)[[1]] = new_sample_ids
      
      ## Permute K matrix and assign new IDs
      new_K = K[new_matings, new_matings]
      rownames(new_K) = new_sample_ids
      colnames(new_K) = new_sample_ids
      
      ## Assign new sample IDs to pheno dataframe
      new_pheno = pheno
      rownames(new_pheno) = new_sample_ids
      
      ## Assign new sample IDs to covariate dataframe
      new_addcovar = addcovar
      rownames(new_addcovar) = new_sample_ids
      
      ## Perform QTL scan
      msg = capture.output({qtl_perm = scanone(pheno=new_pheno, pheno.col=pheno.col, probs=new_probs, K=new_K, addcovar=new_addcovar, snps=snps)})
      
      if (value == "lod") {
        max_A = max(qtl_perm$lod$A[,'lod'])
        if ('X' %in% names(qtl_perm$lod)) {
          max_X = max(qtl_perm$lod$X[,'lod'])
        } else {
          max_X = NA
        }
        perm_max_val = max(max_A, max_X, na.rm=T)
      } else {
        max_A = min(qtl_perm$lod$A[,'p'])
        if ('X' %in% names(qtl_perm$lod)) {
          max_X = min(qtl_perm$lod$X[,'p'])
        } else {
          max_X = NA
        }
        perm_max_val = min(max_A, max_X, na.rm=T)
      }
      perm_max_val
    }
    gc(verbose=F)
    stopImplicitCluster()
  }
  
  ## Return max values for all permutations
  cat(' ... Done. \n')
  return(max_vals)
}
attr(scanone.perm_v2, 'help') = "
This function permutes the genomes of the mapping population and returns the max LOD
(or min P-value) for each permutation.

Parameters:
pheno: the following 6 parameters are those passed to scanone()
pheno.col:
probs:
K:
addcovar:
snps:
value: the value to collect from each permutation ('lod' or 'p'; default='lod').
B: number of permutations to perform (default=100)
parallel: logical indicating whether parallel processes should be run (default=F), 
or numeric indicating number of cores to use. If set to TRUE, detectCores() will be used
to automatically start n-1 parallel processes.
...: optional parameters to pass to foreach

Returns:
A vector of the specified return values for all permutations
"


prob.plot <- function(pheno, pheno.col, probs, marker, qtl) {
  ## Sort the pheno dataframe
  pheno = pheno[order(pheno[, pheno.col]),]
  
  ## Remove NAs from pheno data
  keep = rownames(pheno[!is.na(pheno[,pheno.col]),])
  if (length(keep) < 3) {
    stop("Too many missing phenotype values!")
  }
  probmat = probs[keep, , marker]
  
  ## set margins
  old_mar = par("mar")
  par(mar=c(5,8,4,4))
  
  ## Create heatmap
  grays = gray(seq(from=0, to=0.94, by=0.01))
  oranges = brewer.pal(9, 'Oranges')[3:7]
  probmat[,8:1] -> probmat
  image(z=probmat, zlim=c(0, 1), col=rev(c(grays, oranges, '#FFFFFF')), axes=FALSE)
  #image(z=probmat, zlim=c(0, 1), col=rev(gray(seq(from=0, to=1, by=0.01))), axes=FALSE)
  founder_id = c("A/J","C57BL/6J","129S1/SvImJ","NOD/ShiLtJ","NZO/HlLtJ","CAST/EiJ","PWK/PhJ","WSB/EiJ")
  founder_id[8:1] -> founder_id2
  states <- paste0(colnames(probmat)," (",founder_id2,")")
  nk <- ncol(probmat)
  axis(2, at=seq(from=0, to=1, len=nk), labels=states, las=1, tck=0)
  
  ## Calculate z-scores
  y = pheno[keep, pheno.col]
  z <- (y - mean(y))/sd(y)
  
  axis(3, at=ecdf(z)(pretty(z)), labels=pretty(z))
  mtext("z-score scale", side=3, at=0.5, line=2)
  
  if(sum(qtl$lod$X[,"marker"] == marker) > 0){
    
    marker_pos = paste0(qtl$lod$X[qtl$lod$X[,"marker"] == marker,"chromosome"],":",round(qtl$lod$X[qtl$lod$X[,"marker"] == marker,"position.B38."],2))
    
    lod = round(qtl$lod$X[qtl$lod$X[,"marker"] == marker,"lod"],2)
    
  }else if(sum(qtl$lod$A[,"marker"] == marker) > 0){
    
    marker_pos = paste0(qtl$lod$A[qtl$lod$A[,"marker"] == marker,"chromosome"],":",round(qtl$lod$A[qtl$lod$A[,"marker"] == marker,"position.B38."],2))
    
    lod = round(qtl$lod$A[qtl$lod$A[,"marker"] == marker,"lod"],2)
    
  }else{
    stop("Missing Marker ID")
  }
  
  mtext(paste0("Marker: ", marker," at Chr", marker_pos,"Mb; LOD score: ",lod), side=1, at=0.5, line=3)
  
  xpos <- c(0, mean(y<mean(y)), 1)
  xtxt <- signif(c(min(y), mean(y), max(y)), 4)
  xtxt <- ifelse(abs(xtxt) < 1e-8, 0, xtxt)
  axis(1, at=xpos, labels=xtxt)
  box()
  par(mar=old_mar)
}
attr(prob.plot, 'help') = "
This function creates a probability plot for the mapping population.

Parameters:
pheno: The phenotype dataframe
pheno.col: The column in the pheno dataframe to use as the phenotype
probs: The 3D probability array
marker: The ID of the marker to plot
qtl: A doqtl object; results from scanone()
"


coefplot_v2 = function(doqtl, chr = 1, start=NULL, end=NULL, stat.name = "LOD", remove.outliers=NULL, conf.int = FALSE, legend = TRUE, colors = "DO", sex, ...) {
  
  old.par = par(no.readonly = TRUE)
  
  cross = attr(doqtl, "cross")
  if(is.null(cross)) {
    if(colors[1] == "DO") {    
      colors = do.colors
    } else if(colors[1] == "HS") {
      colors = hs.colors
    } # else if(colors[1] == "HS")
  } else {
    if(cross == "DO") {    
      colors = do.colors
    } else if(cross == "HS") {
      colors = hs.colors
    } # else if(cross == "HS")
  } # else
  
  num.founders = nrow(colors)
  call = match.call()
  
  # Keep only the founder coefficients from the coef.matrix.
  lod = NULL
  coef = NULL
  if(chr == "X") {
    if(missing(sex)) {
      stop("Sex (either M or F) must be specified on X chromosome.")
    } # if(missing(sex))
    lod  = doqtl$lod$X
    if (is.null(start) | is.null(end)) {
      lod = lod[lod[,2] == chr, ]
    } else {
      lod = lod[lod[,2] == chr & lod[,3] >= start & lod[,3] <= end,]
    }
    coef = doqtl$coef$X
    if(sex == "F") {
      columns = match(paste("F", colors[,1], sep = "."), colnames(coef))
    } else {
      columns = match(paste("M", colors[,1], sep = "."), colnames(coef))
    } # else
    columns = columns[!is.na(columns)]
    coef = coef[,c(1,columns)]
    colnames(coef)[1] = "A"
    colnames(coef) = sub("^[MF]\\.", "", colnames(coef))
    coef = coef[rownames(coef) %in% lod[,1],]
    
    ## Remove outliers
    if (!is.null(remove.outliers) & is.numeric(remove.outliers)) {
      coef[abs(coef) > remove.outliers | is.na(coef)] = 0
    }
    
    # Center the coefficient values.
    coef[,2:ncol(coef)] = coef[,2:ncol(coef)] + coef[,1]
    coef = coef - rowMeans(coef)
  } else {
    lod = doqtl$lod$A
    if (is.null(start) | is.null(end)) {
      lod = lod[lod[,2] == chr, ]
    } else {
      lod = lod[lod[,2] == chr & lod[,3] >= start & lod[,3] <= end,]
    }
    ## Remove markers with outlier effects
    #coef_subset = doqtl$coef$A[rownames(doqtl$coef$A) %in% lod[,1], 3:9]
    #markers_to_remove = apply(coef_subset, 1, function(x) {any(abs(x) > remove.outliers | is.na(x))})
    #lod = lod[!markers_to_remove,]
    
    intercept = doqtl$coef$A[,1]
    coef = doqtl$coef$A[,(ncol(doqtl$coef$A)-num.founders+1):ncol(doqtl$coef$A)]
    
    ## Remove outliers
    if (!is.null(remove.outliers) & is.numeric(remove.outliers)) {
      coef[abs(coef) > remove.outliers | is.na(coef)] = 0
    }
    
    coef[,1] = intercept
    colnames(coef)[1] = "A"
    coef = coef[rownames(coef) %in% lod[,1],]
    
    # Center the coefficient values.
    coef[,2:ncol(coef)] = coef[,2:ncol(coef)] + coef[,1]
    coef = coef - rowMeans(coef)
  } # else 
  # Verify that the SNP IDs in the lod & coef matrices match.
  if(!all(lod[,1] == rownames(coef))) {
    stop(paste("The SNP IDs in column 1 of the qtl data frame must match",
               "the SNP IDs in the rownames of the coef matrix."))
  } # if(!all(lod[,1] == rownames(coef)))
  # Verify that the coefficient column names are in the colors.
  if(!all(colnames(coef) %in% colors[,1])) {
    stop(paste("The founder names in the colnames of the coefficient matrix",
               "must be in column 1 of the colors matrix."))
  } # if(!all(colnames(coef) %in% colors[,1]))
  # Convert the chromosome locations to Mb.
  if(max(lod[,3], na.rm = TRUE) > 200) {
    lod[,3] = lod[,3] * 1e-6
  } # if(max(lod[,3], na.rm = TRUE) > 200)
  # Set the layout to plot the coefficients on top and the p-values on the 
  # bottom.
  layout(mat = matrix(1:2, 2, 1), heights = c(0.66666667, 0.3333333))
  par(font = 2, font.lab = 2, font.axis = 2, las = 1, plt =
        c(0.12, 0.99, 0, 0.85), xaxs = "i", lwd = 2)
  # Plot the coefficients.
  plot(lod[,3], coef[,colors[1,1]], type = "l", col = colors[1,3], lwd = 2,
       ylim = c(min(coef, na.rm = TRUE), max(coef * 2, na.rm = TRUE)), xlab = 
         paste("Chr", chr, "(Mb)"), ylab = "Founder Effects", axes = FALSE, ...)
  #abline(v = 0:20 * 10, col = "grey80")
  for(i in 1:nrow(colors)) {
    points(lod[,3], coef[,colors[i,1]], type = "l", col = colors[i,3],
           lwd = 2)
  } # for(i)
  # Draw a legend for the founder names and colors.
  if(legend) {
    legend.side = "topleft"
    if(which.max(lod[,7]) < nrow(lod) / 2) {
      legend.side = "topright"
    } # if(which.max(apply(coef, 1, max)) < nrow(lod) / 2)
    legend(legend.side, colors[,2], col = colors[,3], lty = 1, lwd = 2,
           x.intersp = 0.75, y.intersp = 0.75, bg = "white", cex = 0.8)
  } # if(legend)
  # Add the axis.
  axis(2)
  # Plot a rectangle around the plot.
  par(xpd = NA)
  usr = par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], lwd = 2)
  par(xpd = FALSE)
  # Plot the mapping statistic.
  par(plt = c(0.12, 0.99, 0.35, 1))
  # Single chromosome plot.
  plot(lod[,3], lod[,7], type = "l", lwd = 2, xlab = "",
       ylab = stat.name, ...)
  #abline(v = 0:20 * 10, col = "grey80")
  points(lod[,3], lod[,7], type = "l", lwd = 2)
  # Shade the confidence interval.
  if(conf.int) {
    interval = bayesint_v2(doqtl, chr = chr)
    usr = par("usr")
    rect(interval[1,3], usr[3], interval[3,3], usr[4], col = rgb(0,0,1,0.1), 
         border = NA)
  } # if(!is.na(conf.int))
  mtext(paste("Chr", chr, "(Mb)"), 1, 2)
  usr = par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], lwd = 2)
  par(old.par)
}
attr(coefplot_v2, 'help') = "
A slight modification of the coefplot() DOQTL function that ignores very large effects from rare haplotypes.

Parameters:
doqtl: A doqtl object; results from scanone()
chr: the chromosome to plot
start: the start position (Mb) of the interval to be plotted
end: the end position (Mb) of the interval to be plotted
stat.name: default is 'LOD'
remove.outliers: numeric; sets coefficients above this value to zero (default is NULL--do nothing)
conf.int: draw a box showing the QTL confidence interval (calculated using bayesint_v2())
legend: draw a legend (default is TRUE)
colors: default is 'DO'
sex: required only for plotting the X chromosome
...: additional arguments passed to plot()
"


bayesint_v2 = function (qtl, chr, prob = 0.95, expandtomarkers = TRUE, ignore_genetic_dist=TRUE) 
{
  if (missing(qtl)) {
    stop("bayesint: The qtl cannot be null. Please supply a QTL object.")
  }
  if (missing(chr)) {
    stop(paste("bayesint: The chromosome cannot be null."))
  }
  else if (!chr %in% c(1:19, "X")) {
    stop(paste("bayesint: The chromosome must be 1 to 19 or X."))
  }
  if (prob < 0 | prob > 1) {
    stop(paste("bayesint: The probability must between 0 and 1."))
  }
  old.warn = options("warn")$warn
  options(warn = -1)
  if (is.numeric(chr)) {
    qtl = qtl$lod$A
  }
  else {
    qtl = qtl$lod$X
  }
  options(warn = old.warn)
  qtl[, 1] = as.character(qtl[, 1])
  qtl[, 2] = as.character(qtl[, 2])
  qtl[, 3] = as.numeric(qtl[, 3])
  qtl[, 7] = as.numeric(qtl[, 7])
  qtl = qtl[qtl[, 2] == chr, ]
  pos = qtl[, 3]
  if (any(is.na(pos))) {
    remove = which(is.na(pos))
    qtl = qtl[-remove, ]
    pos = pos[-remove]
  }
  
  breaks = approx(x = pos, y = 10^qtl[, 7], xout = seq(pos[1], 
                                                       pos[length(pos)], length.out = 1e+05))
  widths = diff(breaks$x)
  heights = breaks$y[-1] + breaks$y[-length(breaks$y)]
  trapezoids = 0.5 * heights * widths
  trapezoids = trapezoids/sum(trapezoids)
  ord = order(breaks$y[-length(breaks$y)], decreasing = TRUE)
  wh = min(which(cumsum(trapezoids[ord]) >= prob))
  int = range(ord[1:wh])
  if (ignore_genetic_dist) {
    left.snp = c(NA, qtl[1, 2], breaks$x[int][1], NA, approx(qtl[, 3], qtl[, 
                                                                           5], breaks$x[int][1])$y, approx(qtl[, 3], qtl[, 6], breaks$x[int][1])$y, 
                 approx(qtl[, 3], qtl[, 7], breaks$x[int][1])$y)
    max.snp = qtl[which.max(qtl[, 7]), ]
    right.snp = c(NA, qtl[1, 2], breaks$x[int][2], NA, approx(qtl[, 3], qtl[, 
                                                                            5], breaks$x[int][2])$y, approx(qtl[, 3], qtl[, 6], breaks$x[int][2])$y, 
                  approx(qtl[, 3], qtl[, 7], breaks$x[int][2])$y)
    
    ## TESTING/DEBUGGING
    #print(dim(qtl))
    #print(which.max(qtl[, 7]))
  } else {
    left.snp = c(NA, qtl[1, 2], breaks$x[int][1], approx(qtl[, 
                                                             3], qtl[, 4], breaks$x[int][1])$y, approx(qtl[, 3], qtl[, 
                                                                                                                     5], breaks$x[int][1])$y, approx(qtl[, 3], qtl[, 6], breaks$x[int][1])$y, 
                 approx(qtl[, 3], qtl[, 7], breaks$x[int][1])$y)
    max.snp = qtl[which.max(qtl[, 7]), ]
    right.snp = c(NA, qtl[1, 2], breaks$x[int][2], approx(qtl[, 
                                                              3], qtl[, 4], breaks$x[int][2])$y, approx(qtl[, 3], qtl[, 
                                                                                                                      5], breaks$x[int][2])$y, approx(qtl[, 3], qtl[, 6], breaks$x[int][2])$y, 
                  approx(qtl[, 3], qtl[, 7], breaks$x[int][2])$y)
  }
  if (expandtomarkers) {
    left.snp = qtl[max(which(breaks$x[int][1] >= qtl[, 3])), 
                   ]
    max.snp = qtl[which.max(qtl[, 7]), ]
    right.snp = qtl[min(which(breaks$x[int][2] <= qtl[, 3])), 
                    ]
  }
  retval = rbind(left.snp, max.snp, right.snp)
  retval[, 3] = round(as.numeric(retval[, 3]), digits = 6)
  retval[, 4] = round(as.numeric(retval[, 4]), digits = 6)
  retval[, 5] = round(as.numeric(retval[, 5]), digits = 6)
  retval[, 6] = round(as.numeric(retval[, 6]), digits = 6)
  retval$lod = as.numeric(retval[, 7])
  return(retval)
}
attr(bayesint_v2, 'help') = "
A slight modification of the bayesint() DOQTL function that does not require genetic distances for input.

Parameters:
qtl: A doqtl object; results from scanone()
chr: The column in the pheno dataframe to use as the phenotype
prob: A probability for the desired confidence interval
expandtomarkers: logical (default is TRUE)
ignore_genetic_dist: logical (default is TRUE)
"

## Below is the same code from DOQTL, except code for automatically writing results to a text file is commented out
scanone.perm = function (pheno, pheno.col = 1, probs, K = K, addcovar, intcovar, snps, model = c("additive", "full"), path = ".", nperm = 1000, return.val = c("lod", "p"))
{
  return.val = match.arg(return.val)
  if (!missing(intcovar)) {
    stop("Interactive covariates not yet implemented")
  }
  if (missing(addcovar)) {
    stop(paste("'addcovar' is required and must contain a variable called 'sex' in",
               "order to map correctly on the X chromosome. This is required even if all",
               "of your samples have the same sex."))
  }
  if (length(grep("sex", colnames(addcovar), ignore.case = TRUE)) == 0) {
    stop(paste("'addcovar' is required and must contain a variable called 'sex' in",
               "order to map correctly on the X chromosome. This is required even if all",
               "of your samples have the same sex."))
  }
  if (is.null(rownames(addcovar))) {
    stop("rownames(addcovar) is null. The sample IDs must be in rownames(pheno).")
  }
  if (is.null(rownames(pheno))) {
    stop("rownames(pheno) is null. The sample IDs must be in rownames(pheno).")
  }
  if (is.character(pheno.col)) {
    pheno.col = match(pheno.col, colnames(pheno))
  }
  if (!file.exists(path)) {
    stop(paste("The path", path, "does not exist."))
  }
  probs = filter.geno.probs(probs)
  snps = snps[snps[, 1] %in% dimnames(probs)[[3]], ]
  probs = probs[, , match(snps[, 1], dimnames(probs)[[3]])]
  addcovar = as.matrix(addcovar)
  addcovar = addcovar[rowMeans(is.na(addcovar)) == 0, , drop = FALSE]
  samples = intersect(rownames(pheno), rownames(probs))
  samples = intersect(samples, rownames(addcovar))
  if (length(samples) == 0) {
    stop(paste("There are no matching samples in pheno, addcovar and probs.",
               "Please verify that the sample IDs are in rownames(pheno),",
               "rownames(addcovar) and rownames(probs) and that they all match."))
  }
  pheno = pheno[samples, , drop = FALSE]
  probs = probs[samples, , , drop = FALSE]
  addcovar = addcovar[samples, , drop = FALSE]
  if (is.null(colnames(addcovar))) {
    colnames(addcovar) = paste("addcovar", 1:ncol(addcovar), sep = ".")
  }
  num.auto = get.num.auto(snps)
  auto.snps = which(snps[, 2] %in% 1:num.auto)
  X.snps = which(snps[, 2] == "X")
  xprobs = array(0, c(nrow(probs), 2 * ncol(probs), length(X.snps)), dimnames = list(rownames(probs), paste(rep(c("F", "M"), each = 8), LETTERS[1:8], sep = "."), snps[X.snps, 1]))
  sex.col = grep("^sex$", colnames(addcovar), ignore.case = TRUE)
  sex = as.numeric(factor(addcovar[, sex.col])) - 1
  females = which(sex == 0)
  males = which(sex == 1)
  xprobs[females, 1:8, ] = probs[females, , X.snps]
  xprobs[males, 9:16, ] = probs[males, , X.snps]
  perms = array(0, c(nperm, 2, length(pheno.col)), dimnames = list(1:nperm, c("A", "X"), colnames(pheno)[pheno.col]))
  index = 1
  for (i in pheno.col) {
    print(colnames(pheno)[i])
    p = pheno[, i]
    names(p) = rownames(pheno)
    keep = which(!is.na(p))
    auto.perms = permutations.qtl.LRS(pheno = p[keep], probs = probs[keep, , auto.snps], snps = snps[auto.snps, ], addcovar = addcovar[keep, , drop = FALSE], nperm = nperm, return.val = return.val)
    X.perms = permutations.qtl.LRS(pheno = p[keep], probs = xprobs[keep, , ], snps = snps[X.snps, ], addcovar = addcovar[keep, , drop = FALSE], nperm = nperm, return.val = return.val)
    if (return.val == "lod") {
      auto.perms = auto.perms/(2 * log(10))
      X.perms = X.perms/(2 * log(10))
    }
    perms[, , index] = cbind(auto.perms, X.perms)
    #write.table(perms[, , index], file = paste(path, "/", colnames(pheno)[i], ".perms.txt", sep = ""), sep = "\t")
    index = index + 1
  }
  return(perms)
}


abind = function (..., along = N, rev.along = NULL, new.names = NULL,
                  force.array = TRUE, make.names = use.anon.names, use.anon.names = FALSE,
                  use.first.dimnames = FALSE, hier.names = FALSE, use.dnns = FALSE)
{
  if (is.character(hier.names))
    hier.names <- match.arg(hier.names, c("before", "after",
                                          "none"))
  else hier.names <- if (hier.names)
    "before"
  else "no"
  arg.list <- list(...)
  if (is.list(arg.list[[1]]) && !is.data.frame(arg.list[[1]])) {
    if (length(arg.list) != 1) 
      stop("can only supply one list-valued argument for ...")
    if (make.names) 
      stop("cannot have make.names=TRUE with a list argument")
    arg.list <- arg.list[[1]]
    have.list.arg <- TRUE
  }
  else {
    N <- max(1, sapply(list(...), function(x) length(dim(x))))
    have.list.arg <- FALSE
  }
  if (any(discard <- sapply(arg.list, is.null))) 
    arg.list <- arg.list[!discard]
  if (length(arg.list) == 0) 
    return(NULL)
  N <- max(1, sapply(arg.list, function(x) length(dim(x))))
  if (!is.null(rev.along)) 
    along <- N + 1 - rev.along
  if (along < 1 || along > N || (along > floor(along) && along < 
                                 ceiling(along))) {
    N <- N + 1
    along <- max(1, min(N + 1, ceiling(along)))
  }
  if (length(along) > 1 || along < 1 || along > N + 1) 
    stop(paste("\"along\" must specify one dimension of the array,", 
               "or interpolate between two dimensions of the array", 
               sep = "\n"))
  if (!force.array && N == 2) {
    if (!have.list.arg) {
      if (along == 2) 
        return(cbind(...))
      if (along == 1) 
        return(rbind(...))
    }
    else {
      if (along == 2) 
        return(do.call("cbind", arg.list))
      if (along == 1) 
        return(do.call("rbind", arg.list))
    }
  }
  if (along > N || along < 0) 
    stop("along must be between 0 and ", N)
  pre <- seq(from = 1, len = along - 1)
  post <- seq(to = N - 1, len = N - along)
  perm <- c(seq(len = N)[-along], along)
  arg.names <- names(arg.list)
  if (is.null(arg.names)) 
    arg.names <- rep("", length(arg.list))
  if (is.character(new.names)) {
    arg.names[seq(along = new.names)[nchar(new.names) > 0]] <- new.names[nchar(new.names) > 
                                                                           0]
    new.names <- NULL
  }
  if (any(arg.names == "")) {
    if (make.names) {
      dot.args <- match.call(expand.dots = FALSE)$...
      if (is.call(dot.args) && identical(dot.args[[1]], 
                                         as.name("list"))) 
        dot.args <- dot.args[-1]
      arg.alt.names <- arg.names
      for (i in seq(along = arg.names)) {
        if (arg.alt.names[i] == "") {
          if (object.size(dot.args[[i]]) < 1000) {
            arg.alt.names[i] <- paste(deparse(dot.args[[i]], 
                                              40), collapse = ";")
          }
          else {
            arg.alt.names[i] <- paste("X", i, sep = "")
          }
          arg.names[i] <- arg.alt.names[i]
        }
      }
    }
    else {
      arg.alt.names <- arg.names
      arg.alt.names[arg.names == ""] <- paste("X", seq(along = arg.names), 
                                              sep = "")[arg.names == ""]
    }
  }
  else {
    arg.alt.names <- arg.names
  }
  use.along.names <- any(arg.names != "")
  names(arg.list) <- arg.names
  arg.dimnames <- matrix(vector("list", N * length(arg.names)), 
                         nrow = N, ncol = length(arg.names))
  dimnames(arg.dimnames) <- list(NULL, arg.names)
  arg.dnns <- matrix(vector("list", N * length(arg.names)), 
                     nrow = N, ncol = length(arg.names))
  dimnames(arg.dnns) <- list(NULL, arg.names)
  dimnames.new <- vector("list", N)
  arg.dim <- matrix(integer(1), nrow = N, ncol = length(arg.names))
  for (i in seq(len = length(arg.list))) {
    m <- arg.list[[i]]
    m.changed <- FALSE
    if (is.data.frame(m)) {
      m <- as.matrix(m)
      m.changed <- TRUE
    }
    else if (!is.array(m) && !is.null(m)) {
      if (!is.atomic(m)) 
        stop("arg '", arg.alt.names[i], "' is non-atomic")
      dn <- names(m)
      m <- as.array(m)
      if (length(dim(m)) == 1 && !is.null(dn)) 
        dimnames(m) <- list(dn)
      m.changed <- TRUE
    }
    new.dim <- dim(m)
    if (length(new.dim) == N) {
      if (!is.null(dimnames(m))) {
        arg.dimnames[, i] <- dimnames(m)
        if (use.dnns && !is.null(names(dimnames(m)))) 
          arg.dnns[, i] <- as.list(names(dimnames(m)))
      }
      arg.dim[, i] <- new.dim
    }
    else if (length(new.dim) == N - 1) {
      if (!is.null(dimnames(m))) {
        arg.dimnames[-along, i] <- dimnames(m)
        if (use.dnns && !is.null(names(dimnames(m)))) 
          arg.dnns[-along, i] <- as.list(names(dimnames(m)))
        dimnames(m) <- NULL
      }
      arg.dim[, i] <- c(new.dim[pre], 1, new.dim[post])
      if (any(perm != seq(along = perm))) {
        dim(m) <- c(new.dim[pre], 1, new.dim[post])
        m.changed <- TRUE
      }
    }
    else {
      stop("'", arg.alt.names[i], "' does not fit: should have `length(dim())'=", 
           N, " or ", N - 1)
    }
    if (any(perm != seq(along = perm))) 
      arg.list[[i]] <- aperm(m, perm)
    else if (m.changed) 
      arg.list[[i]] <- m
  }
  conform.dim <- arg.dim[, 1]
  for (i in seq(len = ncol(arg.dim))) {
    if (any((conform.dim != arg.dim[, i])[-along])) {
      stop("arg '", arg.alt.names[i], "' has dims=", paste(arg.dim[, 
                                                                   i], collapse = ", "), "; but need dims=", paste(replace(conform.dim, 
                                                                                                                           along, "X"), collapse = ", "))
    }
  }
  if (N > 1) 
    for (dd in seq(len = N)[-along]) {
      for (i in (if (use.first.dimnames) 
        seq(along = arg.names)
        else rev(seq(along = arg.names)))) {
        if (length(arg.dimnames[[dd, i]]) > 0) {
          dimnames.new[[dd]] <- arg.dimnames[[dd, i]]
          if (use.dnns && !is.null(arg.dnns[[dd, i]])) 
            names(dimnames.new)[dd] <- arg.dnns[[dd, 
                                                 i]]
          break
        }
      }
    }
  for (i in seq(len = length(arg.names))) {
    if (arg.dim[along, i] > 0) {
      dnm.along <- arg.dimnames[[along, i]]
      if (length(dnm.along) == arg.dim[along, i]) {
        use.along.names <- TRUE
        if (hier.names == "before" && arg.names[i] != 
            "") 
          dnm.along <- paste(arg.names[i], dnm.along, 
                             sep = ".")
        else if (hier.names == "after" && arg.names[i] != 
                 "") 
          dnm.along <- paste(dnm.along, arg.names[i], 
                             sep = ".")
      }
      else {
        if (arg.dim[along, i] == 1) 
          dnm.along <- arg.names[i]
        else if (arg.names[i] == "") 
          dnm.along <- rep("", arg.dim[along, i])
        else dnm.along <- paste(arg.names[i], seq(length = arg.dim[along, 
                                                                   i]), sep = "")
      }
      dimnames.new[[along]] <- c(dimnames.new[[along]], 
                                 dnm.along)
    }
    if (use.dnns) {
      dnn <- unlist(arg.dnns[along, ])
      if (length(dnn)) {
        if (!use.first.dimnames) 
          dnn <- rev(dnn)
        names(dimnames.new)[along] <- dnn[1]
      }
    }
  }
  if (!use.along.names) 
    dimnames.new[along] <- list(NULL)
  out <- array(unlist(arg.list, use.names = FALSE), dim = c(arg.dim[-along, 
                                                                    1], sum(arg.dim[along, ])), dimnames = dimnames.new[perm])
  if (any(order(perm) != seq(along = perm))) 
    out <- aperm(out, order(perm))
  if (!is.null(new.names) && is.list(new.names)) {
    for (dd in seq(len = N)) {
      if (!is.null(new.names[[dd]])) {
        if (length(new.names[[dd]]) == dim(out)[dd]) 
          dimnames(out)[[dd]] <- new.names[[dd]]
        else if (length(new.names[[dd]])) 
          warning(paste("Component ", dd, " of new.names ignored: has length ", 
                        length(new.names[[dd]]), ", should be ", 
                        dim(out)[dd], sep = ""))
      }
      if (use.dnns && !is.null(names(new.names)) && names(new.names)[dd] != 
          "") 
        names(dimnames(out))[dd] <- names(new.names)[dd]
    }
  }
  if (use.dnns && !is.null(names(dimnames(out))) && any(i <- is.na(names(dimnames(out))))) 
    names(dimnames(out))[i] <- ""
  out
}