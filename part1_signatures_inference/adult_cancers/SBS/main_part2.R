# load the required libraries and sources
library("glmnet")
library("lsa")
library("parallel")
source("R/mutational.signatures.discovery.R")
source("R/utils.R")

# load data and results
load("final_data/trinucleotides_counts.RData")
load("final_results/goodness_fit.RData")
load("final_results/signatures_assignment.RData")

# define signatures incremental ordering
beta = signatures_assignment$beta
normalized_alpha = (signatures_assignment$alpha/rowSums(signatures_assignment$alpha))
signatures_ordering = unique(names(sort(colSums(normalized_alpha,na.rm=TRUE),decreasing=TRUE)))
save(signatures_ordering,file="results/signatures_ordering.RData")

# perform incremental signatures assignment
cat(0.00,"\n")
for(i in 1:length(signatures_ordering)) {
    set.seed(i)
    curr_beta = beta[signatures_ordering[1:i],,drop=FALSE]
    incremental_assignment = signaturesAssignment(x=trinucleotides_counts,beta=curr_beta)
    save(incremental_assignment,file=paste0("results/incremental_assignment_",i,".RData"))
    cat(i/length(signatures_ordering),"\n")
}
