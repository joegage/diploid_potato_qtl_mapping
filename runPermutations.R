runPerm <- function(y, bedMat, i){
            print(paste0("Permutation number: ", i))
            yShuf <- sample(y, length(y))
            pValPerm <- dual_mapping(yShuf, bedMat)
            return(apply(pValPerm, 2, min, na.rm=TRUE))
        }
