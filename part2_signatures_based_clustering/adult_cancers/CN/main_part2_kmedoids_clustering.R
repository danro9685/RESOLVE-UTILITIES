# load the required libraries
library("cluster")
library("parallel")

# load the results
load(file="results/signatures_distance.RData")
sig_dist = as.dist(signatures_distance)

# setting up parallel execution
configurations = 1:25
parallel = makeCluster(length(configurations), outfile = "")
res_clusterEvalQ = clusterEvalQ(parallel,library("cluster", warn.conflicts = FALSE,
    quietly = TRUE, verbose = FALSE))
clusterExport(parallel, varlist = c("sig_dist"), envir = environment())
clusterSetRNGStream(parallel, iseed = 12345)

# perform K-Medoids clustering
inference_results = parLapply(parallel, configurations, function(i) {

    # perform the analysis
    set.seed(i)
    res_clustering = pam( x = sig_dist,
                          k = i,
                          diss = TRUE,
                          medoids = "random",
                          nstart = 100,
                          cluster.only = TRUE,
                          do.swap = TRUE,
                          keep.diss = FALSE,
                          keep.data = FALSE)

    # save the results
    save(res_clustering,file=paste0("results/res_clustering_",i,".RData"))
    return(NA)

})

# close parallel
stopCluster(parallel)
