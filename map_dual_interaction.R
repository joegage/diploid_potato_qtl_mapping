## Testing lines:
## bedMat = read.table("W4M6_haplotype_bins.bed", header=FALSE, sep="\t", stringsAsFactors=FALSE)
## nind <- ( ncol(bedMat) - 3 ) / 2
## y <- rnorm(nind)

make_model <- function(y, x1, x2){
    ## y is vector of phenotypes
    ## x1 is vector of parent1 calls, as a factor
    ## x2 is vector of parent2 calls, as a factor    

    pX1 <- NA
    pX2 <- NA
    pInteraction <- NA
    overall <- NA
    #fvals <- NA
    
    if(nlevels(x1) == 2 & nlevels(x2) == 2){
        reg <- lm(y~(factor(x1)*factor(x2)))
        pX1 <- summary(reg)$coefficients[2,4]
        pX2 <- summary(reg)$coefficients[3,4]
        pInteraction <- summary(reg)$coefficients[4,4]
        fval = summary(reg)$fstatistic
       # fvals <- fval[1]
        overall = pf(fval[1], fval[2], fval[3], lower.tail=FALSE)
    } else if(nlevels(x1) == 2 & nlevels(x2) == 1){
        reg <- lm(y~(factor(x1)))
        pX1 <- summary(reg)$coefficients[2,4]
        fval = summary(reg)$fstatistic
        #fvals <- fval[1]
        overall = pf(fval[1], fval[2], fval[3], lower.tail=FALSE)
    } else if(nlevels(x1) == 1 & nlevels(x2) == 2){
        reg <- lm(y~(factor(x2)))
        pX2 <- summary(reg)$coefficients[2,4]
        fval = summary(reg)$fstatistic
        #fvals <- fval[1]
        overall = pf(fval[1], fval[2], fval[3], lower.tail=FALSE)
    }
    return(c(pX1, pX2, pInteraction, overall))
}


dual_mapping <- function(y, bedMat){

    info <- bedMat[,1:3]
    bedMat <- as.matrix(bedMat[,-c(1:3)])

    genoParentA <- bedMat[,seq(1, ncol(bedMat), 2)]
    genoParentB <- bedMat[,seq(2, ncol(bedMat), 2)]

    ## Testing lines
    ## i <- 10
    ## x1 <- factor(genoParentA[i,])
    ## x2 <- factor(genoParentB[i,])

    testHaplos <- function(j){
        x1 <- factor(genoParentA[j,])
        x2 <- factor(genoParentB[j,])
        return(make_model(y, x1, x2))
    }
    
    pVals <- sapply(1:nrow(bedMat), function(j) testHaplos(j))
    pVals <- t(pVals)
}

