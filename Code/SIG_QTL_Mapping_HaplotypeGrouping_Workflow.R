
####################################################################################
#### Step 1. Load Necessary R Functions and Libraries, Download All Required Files
####################################################################################

## Load R functions and libraries
source('../Code/WNV_rix_qtl_mapping_functions_publication.r')

### All of the files below need to have destination directories added by replacing ??????. ####

# Download Oas1b_status_recoded.txt
download.file(url = "https://figshare.com/ndownloader/files/51510053", destfile = "??????/Data/Mapping/Oas1b_status_recoded.txt")

# Download PhenotypesWNV Spleen Treg D7.csv. Both the file name and the directory destination need to be 
# added for this since there are 30 phenotype files
download.file(url = "https://figshare.com/ndownloader/files/51509957", destfile = "??????/Data/Phenotypes/PhenotypesWNV Spleen Treg D7.csv")
 
# Download rix_universal_model_prob_males_27-Jun-2016.rda
download.file(url = "https://figshare.com/ndownloader/files/51510062", destfile = "??????/Data/Mapping//rix_universal_model_prob_males_27-Jun-2016.rda")

# Download CC001-Uncb38V01.csv
download.file(url = "https://figshare.com/ndownloader/files/51510059", destfile = "??????/Data/Mapping//CC001-Uncb38V01.csv")

# Download mgp.v5.merged.snps_all.dbSNP142_chrX.recode.vcf.gz
download.file(url = "https://figshare.com/ndownloader/files/51509939", destfile = "??????/Data/Sanger/mgp.v5.merged.snps_all.dbSNP142_chrX.recode.vcf.gz")

# Download mgp.v5.merged.snps_all.dbSNP142_chrX.recode.vcf.gz.tbi
download.file(url = "https://figshare.com/ndownloader/files/51509936", destfile = "??????/Data/Sanger/mgp.v5.merged.snps_all.dbSNP142_chrX.recode.vcf.gz.tbi")

# Download mgi_annotation.rpt
download.file(url = "https://figshare.com/ndownloader/files/51510080", destfile = "??????/Data/Annotation/mgi_annotation.rpt")

# Download MOUSE_10090_idmapping.dat
download.file(url = "https://figshare.com/ndownloader/files/51510083", destfile = "??????/Data/Annotation/MOUSE_10090_idmapping.dat")

####################################################################################
#### Step 2. Create Phenotype Dataframe
####################################################################################

## Read phenotype data, also need to make sure these directories are created and fully notated here
cleaned_data_dir = '?????/Data/Phenotypes'
pheno_dir = '?????/Data/Phenotypes'
mapping_dir = '?????/Data/Mapping'
data_dir = '?????/Data/Output'

## For this example we're working with Spleen Treg Day 7 phenotypes
pheno = read.csv(file.path(cleaned_data_dir, 'PhenotypesWNV Spleen Treg D7.csv'), header=T, as.is=T)
dim(pheno)

## The phenotype will be D7 WNV in spleen, treg_CD4pos_Foxp3neg_box_cox_7_WNV_Spleen
pheno = pheno[with(pheno, Virus=='WNV' & Timepoint=='7' & Tissue=='spleen'), ]
dim(pheno)

## View the phenotype distribution accross the mapping population
library(lattice)
dotplot(reorder(pheno[,'UW_Line'], pheno[,'treg_CD4pos_Foxp3neg_box_cox_7_WNV_Spleen'], mean, na.rm=T) ~ 
        pheno[,'treg_CD4pos_Foxp3neg_box_cox_7_WNV_Spleen'], 
        panel = function(x,y,...) {panel.dotplot(x,y,...); panel.abline(v=0, col.line="red")}, 
        pch=19, ylab='UW Line', xlab="treg_CD4pos_Foxp3neg_box_cox_7_WNV_Spleen")

####################################################################################
#### Step 3. Update Phenotypes and Create Covariate Dataframe
####################################################################################

## Sort pheno dataframe and set rownames
pheno = pheno[with(pheno, order(Mating, RIX_ID)),]
rownames(pheno) = pheno$ID

## Add sex column
pheno$sex = 'M'

## Create covariate dataframe (must include sex)
covar = data.frame(sex = as.numeric(pheno$sex == 'M'))
#covar = data.frame(sex = pheno$sex)
rownames(covar) = pheno$ID

## Get IDs and Matings for each sample
samples = pheno$ID
matings = unlist(lapply(strsplit(samples, '_'), function(x) {x[1]}))
matings = unique(matings)

## Read strain ID mapping file (with Oas1b status for WNV mapping analyses)
## Note: use the Mx1 file (linked at top) for Flu analyses
strain_map = read.delim(file.path(mapping_dir, 'Oas1b_status_recoded.txt'), header=T, as.is=T, sep='\t')
head(strain_map)

## Check that matings match the strain mapping file
setdiff(matings, strain_map$Mating)

## Update covar with Oas1b status (only for WNV analyses)
covar$Mating = pheno$Mating
covar$Oas1b = sapply(covar$Mating, function(x) {if (x %in% strain_map$Mating) strain_map$Oas1b[strain_map$Mating == x] else NA})
covar$Oas1b_High = sapply(covar$Mating, function(x) {if (x %in% strain_map$Mating) strain_map$Oas1b_High[strain_map$Mating == x] else NA})
covar$Oas1b_Mod = sapply(covar$Mating, function(x) {if (x %in% strain_map$Mating) strain_map$Oas1b_Mod[strain_map$Mating == x] else NA})
covar$Oas1b_Low = sapply(covar$Mating, function(x) {if (x %in% strain_map$Mating) strain_map$Oas1b_Low[strain_map$Mating == x] else NA})
pheno$Oas1b = covar$Oas1b
pheno$Oas1b_High = covar$Oas1b_High
pheno$Oas1b_Mod = covar$Oas1b_Mod
pheno$Oas1b_Low = covar$Oas1b_Low

####################################################################################
#### Step 4. Construct 3D Array of Probabilities
####################################################################################

## Load universal model probabilities (loads a model.probs object containing all RIX lines)
load(file.path(mapping_dir, 'rix_universal_model_prob_males_27-Jun-2016.rda'))

## Create model.probs array specific to the mapping population
model.probs = model.probs[pheno$Mating, , ]
dimnames(model.probs)[[1]] = pheno$ID

## Check model.probs object
dim(model.probs)
names(dimnames(model.probs))
dim(model.probs)[1] == dim(pheno)[1]
model.probs[1,,1:5]

## Fix very small genotype probabilities
model.probs[model.probs < 0.005] = 1e-20

## Check model.probs object
model.probs[1,,1:5]

####################################################################################
#### Step 5. Calculate Kinship Matrix
####################################################################################

## Create kinship probability matrix
K = kinship.probs(model.probs)

## Garbage collection might help free memory
gc()

## Check kinship matrix
K[1:5, 1:5]

####################################################################################
#### Step 6. Perform QTL Scan
####################################################################################

## First get marker positions
marker_pos = read.csv(file.path(mapping_dir, 'CC001-Uncb38V01.csv'), as.is=T)
marker_pos = marker_pos[,1:3]
marker_pos$position_cM = NA
head(marker_pos)

## Run QTL scan
qtl = scanone(pheno=pheno, pheno.col='treg_CD4pos_Foxp3neg_box_cox_7_WNV_Spleen', probs=model.probs, K=K, 
              addcovar=covar[, c("Oas1b_High", "Oas1b_Mod",'sex'), drop=F], snps=marker_pos)

## Run permutations to calculate significance threshold
## Note: nperm=10 is for demonstration purposes only. Normally, at least 1000 permutations should be done.
perms = scanone.perm(pheno=pheno, pheno.col='treg_CD4pos_Foxp3neg_box_cox_7_WNV_Spleen', probs=model.probs, 
                     addcovar=covar[, c('sex', "Oas1b_High", "Oas1b_Mod"), drop=F], snps=marker_pos, nperm = 10)

## Get the significance thresholds
sig.thr = get.sig.thr(perms[,,'treg_CD4pos_Foxp3neg_CXCR3pos_biexp_21_WNV_Spleen'], alpha=0.05)
sugg.thr = get.sig.thr(perms[,,'treg_CD4pos_Foxp3neg_CXCR3pos_biexp_21_WNV_Spleen'], alpha=0.1)
med.thr = get.sig.thr(perms[,,'treg_CD4pos_Foxp3neg_CXCR3pos_biexp_21_WNV_Spleen'], alpha=0.5)

## Plot QTL results
plot(qtl, sig.thr = sig.thr, main = 'Day 7 Spleen Treg CD4pos Foxp3neg')

## Add lines for suggestive peaks (alpha = 0.1)
lines(x = c(0, 2462), y = rep(sugg.thr[1,"A"], 2), lwd = 2, lty=2, col = 'red')
lines(x = c(2462, 2633), y = rep(sugg.thr[1,"X"], 2), lwd = 2, lty=2, col = 'red')

####################################################################################
####Step 7. Identify QTL Intervals
####################################################################################

## Identify the Bayes Credible Interval for the peak on chromosome X
interval = bayesint_v2(qtl, chr = 'X')
interval

####################################################################################
#### Step 8. Create a Coefficient Plot and Probability Plot for Significant Peaks
####################################################################################

## Create founder effects (coefficient) plot for peak interval
coefplot_v2(qtl, chr='X', start=interval[1,3], end=interval[3,3], sex='M', main='treg_CD4pos_Foxp3neg_box_cox_7_WNV_Spleen', cex.main=.7)
 
## Get marker at peak
get_min_marker(qtl, chr='X')

## Create probability plot
phenoclass="Treg Spleen Day 7" #used in report title names
prob.plot(pheno=pheno, pheno.col='treg_CD4pos_Foxp3neg_box_cox_7_WNV_Spleen', probs=model.probs, marker='UNC31155388', qtl=qtl)

####################################################################################
#Step 9. Run haplotype groupings, identify variants of interest in QTL interval, and identify gene candidates
####################################################################################

## Create the association plots and identify gene candidates in the interval
## Directories in the following function will need to be changed to map the directories in the above file downloads
phenocol<-'treg_CD4pos_Foxp3neg_box_cox_7_WNV_Spleen'
phenocolend="WNV_Spleen"
untransphenocol<-substr(phenocol,1,str_locate(phenocol,phenocolend)[,1]-2)
print("inner function")
start=interval[1,3]
end=interval[3,3]
sex='M'
outputdir=data_dir
chr='X'
# The following function will need to have the system PATH variable set to run the 'tabix' command. Tabix will also need to be installed to 
# run this R function

FoundProbAssocPlotsComplexShort_Publication(phenocol=phenocol,untransphenocol=untransphenocol,sex=sex,chr=chr,start=start,end=end,qtl=qtl) 




