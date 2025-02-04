# perform the estimation of the mutational signatures associated to mutations in specific genes
associationEstimation <- function( alterations, signatures ) {
    
    # initialize alpha with the signatures esposure data, 
    # alpha is an N x K matrix, where rows (N) represent the samples and columns (K) represent the K mutational signatures
    alpha <- as.matrix(signatures)

    # initialize beta with an empty matrix, 
    # beta is a K x M matrix, where rows represent the associations for the K signatures to mutations in specific genes, 
    # and columns represent the M considered genes (usually a selected set of genes of interest)
    beta <- matrix(NA, nrow = ncol(alpha), ncol = ncol(alterations))
    rownames(beta) <- colnames(alpha)
    colnames(beta) <- colnames(alterations)

    # initialize intercept estimates
    intercept <- rep(NA,ncol(beta))
    names(intercept) <- colnames(beta)

    # perform the estimation by regularized logistic regression
    for (i in seq_len(ncol(beta))) {
        x <- as.numeric(alterations[,i]) # the mutations are the response (binary) variables
        beta[,i] <- tryCatch({ # try to perform fit
            fit <- cv.glmnet(x = alpha, y = x, type.measure = "default", family = "binomial", 
                        nfolds = 10, nlambda = 100, maxit = 1e+05, alpha = 1)
            fit <- as.numeric(coef(fit,s=fit$lambda.min))
            intercept[i] <- fit[1]
            res <- fit[-1]
            res
        }, error = function(err) { # if the fit fails, discard the current gene
            message("Fit failed for gene ",colnames(beta)[i],": discarding the gene.","\n")
            return(NA)
        })
    }

    # remove any gene that cannot be associated to signatures
    invalid_genes <- which(is.na(beta),arr.ind=TRUE)
    if (nrow(invalid_genes)>0) {
        invalid_genes <- unique(invalid_genes[,"col"])
        alterations <- alterations[,-invalid_genes,drop=FALSE]
        beta <- beta[,-invalid_genes,drop=FALSE]
        intercept <- intercept[-invalid_genes]
    }

    # compute the probability of each signature to be associated to a specific mutation
    probabilities <- beta
    fold_changes <- beta
    for(i in 1:nrow(probabilities)) {
        for(j in 1:ncol(probabilities)) {
            beta_0 <- as.numeric(intercept[j])
            beta_var <- as.numeric(probabilities[i,j])
            t_0 <- beta_0
            t_var <- (beta_0+beta_var)
            prob_0 <- (1/(1+((exp(1)^(-t_0)))))
            prob_var <- (1/(1+((exp(1)^(-t_var)))))
            probabilities[i,j] <- prob_var
            fold_changes[i,j] <- (prob_var/prob_0)
        }
    }

    # return the assigned signatures
    results <- list(alterations = alterations, alpha = alpha, beta = beta, intercept = intercept, probabilities = probabilities, fold_changes = fold_changes)
    return(results)

}
