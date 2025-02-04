# perform the estimation of the mutations associated to specific mutational signatures
associationEstimation <- function( alterations, signatures ) {
    
    # initialize alpha with the mutations data, 
    # alpha is an N x K matrix, where rows (N) represent the samples and columns (K) represent the K observed mutations
    alpha <- as.matrix(alterations)

    # initialize beta with an empty matrix, 
    # beta is a K x M matrix, where rows represent the associations for the K mutations to specific mutational signatures, 
    # and columns represent the M considered mutational signatures
    beta <- matrix(NA, nrow = ncol(alpha), ncol = ncol(signatures))
    rownames(beta) <- colnames(alpha)
    colnames(beta) <- colnames(signatures)

    # initialize intercept estimates
    intercept <- rep(NA,ncol(beta))
    names(intercept) <- colnames(beta)

    # perform the estimation by regularized regression
    for (i in seq_len(ncol(beta))) {
        x <- as.numeric(signatures[,i]) # the mutations are the esponse (binary) variables
        beta[,i] <- tryCatch({ # try to perform fit
            fit <- cv.glmnet(x = alpha, y = x, type.measure = "default", family = "gaussian", 
                        nfolds = 10, nlambda = 100, maxit = 1e+05, alpha = 1)
            fit <- as.numeric(coef(fit,s=fit$lambda.min))
            intercept[i] <- fit[1]
            res <- fit[-1]
            res
        }, error = function(err) { # if the fit fails, discard the current gene
            message("Fit failed for signature ",colnames(beta)[i],": discarding the signature.","\n")
            return(NA)
        })
    }

    # remove any signature that cannot be associated to genes
    invalid_genes <- which(is.na(beta),arr.ind=TRUE)
    if (nrow(invalid_genes)>0) {
        invalid_genes <- unique(invalid_genes[,"col"])
        alterations <- alterations[,-invalid_genes,drop=FALSE]
        beta <- beta[,-invalid_genes,drop=FALSE]
        intercept <- intercept[-invalid_genes]
    }

    # compute the fold changes of each mutation to be associated to a specific signature
    fold_changes <- beta
    for(i in 1:nrow(fold_changes)) {
        for(j in 1:ncol(fold_changes)) {
            estimate_0 <- as.numeric(intercept[j])
            estimate_var <- (estimate_0+as.numeric(fold_changes[i,j])) # estimation assuming the considered mutation is present
            min_val <- min(estimate_0,estimate_var)
            if(min_val < 0) { # signatures can cause mutations not remove them
                estimate_0 <- estimate_0 + abs(min_val)
                estimate_var <- estimate_var + abs(min_val)
            }
            fold_changes[i,j] <- ((estimate_var+1)/(estimate_0+1))
        }
    }

    # return the assigned signatures
    results <- list(alterations = alterations, signatures = signatures, alpha = alpha, beta = beta, intercept = intercept, fold_changes = fold_changes)
    return(results)

}
