# load the required libraries
library("lsa")

# load the data
load(file="final_results/signatures_significance.RData")

# compute signatures-based distance
alpha = signatures_significance$alpha
alpha = alpha[sort(unique(names(which(signatures_significance$goodness_fit>0.90)))),]
cos_similarity = cosine(t(alpha))
invalid = which(is.na(cos_similarity))
if(length(invalid)>0) {
    cos_similarity[invalid] = 0
}
signatures_distance = (1-cos_similarity)
invalid = which(signatures_distance<0)
if(length(invalid)>0) {
    signatures_distance[invalid] = 0
}
invalid = which(signatures_distance>1)
if(length(invalid)>0) {
    signatures_distance[invalid] = 1
}

# save the results
save(signatures_distance,file="results/signatures_distance.RData")
