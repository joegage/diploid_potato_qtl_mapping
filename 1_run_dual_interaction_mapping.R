##--------------------------##
## QTL mapping interactions ##
##--------------------------##




rm(list=ls())
library(doParallel)
setwd("/Volumes/WD elements 25a3/218_data/w4m6_yield/")
index <- read.table("/Volumes/WD elements 25a3/218_data/potato_references/DM_v404/potato_dm_v404_all_pm_un.chr1-12_len.bed")

## Change this to 0 if you don't want to run permutations:
## Perms are not particularly fast to run, but the multithreading helps
nPerm=2
nCore=2
sigThreshold=0.05
registerDoParallel(cores=nCore)

## Read files
bedMat = read.table("W4M6_haplotype_bins_split.bed", header=TRUE, sep="\t", stringsAsFactors=FALSE)
pheno <- read.csv("yield_phenological_traits_allGenos.csv", stringsAsFactors=FALSE, header=TRUE, row.names=1)
pheno <- t(pheno)

## Order phenotypes
genoOrder <- grep("X[0-9]{1,3}$", colnames(bedMat), value=TRUE)
pheno <- pheno[genoOrder, ]

## Pull in functions and make output dir
source("map_dual_interaction.R")
source("runPermutations.R")
system("mkdir -p ./output/")

##--------------------------##
## begin iteration by trait ##
##--------------------------##

## set plot layout
layout(matrix(c(1:72), nrow=6, byrow=T), widths=index$V2)
par(mar=c(1,0,0.5,0), oma=c(3,3,1,1))

## loop over phenotypes
for(i in 1:ncol(pheno)){
    
        ## trait name
        trait <- colnames(pheno)[i]
        print(paste0("Starting trait: ", trait))
        b <- bedMat[,1:3]
        fileBase <- paste0("output/dual_parent_mapping_", trait)
        
        ## Get phenotype
        y <- pheno[,trait]
        
        ## Do mapping
        pValues <- dual_mapping(y, bedMat)
        colnames(pValues) <- c("W4", "M6", "Interaction", "Overall")
        pvals <- as.data.frame(pValues)
        b$W4 <- pvals$W4
        b$M6 <- pvals$M6
        b$Interaction <- pvals$Interaction
        b$Overall <- pvals$Overall
        
        ## Do permutations if desired
        if(nPerm > 0){
                perms <- foreach(i=1:nPerm, .combine=rbind) %dopar% runPerm(y, bedMat, i)
                
                ## Write file
                colnames(perms) <- colnames(pValues)
                permFile <- paste0(fileBase, "_perms.txt")
                write.table(perms, permFile, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
        }
        
        ## iterate over each chromosome
        for (j in 1:12){
                if(j < 10){
                        name <- paste("chr0",j,sep="")
                }else{
                        name <- paste("chr",j,sep="")
                }
                subb <- subset(b, b$chr==name)
                
                # convert to Mb 
                subb$start <- subb$start/1000000
                yRange <- range(-log10(b[,4:7]), na.rm=TRUE) + c(0,1)
                cols <- c("darkorange", "darkorchid1", "grey70", "dodgerblue")
                
                ## pvalue profiles
                chrlength <- subset(index, index$V1==name)
                
                ## new plots
                plot.new()
                plot.window(ylim=yRange, xlim=c(0,max(chrlength$V2)))
                
                ## plot -log10 significance
                lines(subb$start,-log10(subb$W4), col=cols[1])
                lines(subb$start, -log10(subb$M6), type="l", col=cols[2])
                lines(subb$start, -log10(subb$Interaction), type="l", col=cols[3])
                lines(subb$start, -log10(subb$Overall), type="l", col=cols[4])
                
                ## Annotate plots in specific locations
                if(j == 1){
                        axis(2)
                        title(ylab=expression(-log[10](p-value)), main=trait)
                        if(i==6){
                            axis(1)
                        }
                } else if(i == 6){
                        axis(1)
                }
                
                ## Sig thresholds from permutation, if any
                if(nPerm > 0){
                        thresholds <- apply(perms, 2, function(x) -log10(quantile(x, sigThreshold)))
                        abline(h=thresholds, col=cols, lty=2)
                }
                
                
        }
        ## Write file
        out <- cbind(bedMat[,1:4], pValues)
        pFile <- paste0(fileBase, "_pVals.txt")
        write.table(out, pFile, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
}

print("All Done.")
