# load the required libraries and sources
library("glmnet")
library("lsa")
library("parallel")
source("R/mutational.signatures.significance.R")
source("R/utils.R")

# load data and results
load("final_data/analysis_data.RData")
counts = analysis_data[["Pan-Cancer"]]
load("final_results/signatures_assignment.RData")

# perform signatures significance estimation
set.seed(12345)
signatures_significance = signaturesSignificance(x=counts,beta=signatures_assignment$beta,cosine_thr=0.95,min_contribution=0.05,pvalue_thr=0.05,nboot=100,num_processes=100,verbose=TRUE)

# save the results
save(signatures_significance,file="final_results/signatures_significance.RData")
