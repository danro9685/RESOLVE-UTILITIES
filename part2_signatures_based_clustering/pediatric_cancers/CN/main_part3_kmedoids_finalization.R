# load the required libraries
library("cluster")
library("ggplot2")

# load the results
load(file="results/signatures_distance.RData")
sig_dist = as.dist(signatures_distance)

# set the seed
set.seed(12345)

# save K-Medoids clustering results
all_res_clustering = list()
for(i in 1:25) {
    load(file=paste0("results/res_clustering_",i,".RData"))
    all_res_clustering[[i]] = res_clustering
}
res_clustering = all_res_clustering

# compute silhouette and make plot
res_silhouette = list()
avg = 0
for(i in 2:25) {
    res_silhouette[[i]] = silhouette(x = res_clustering[[i]], dist = sig_dist)
    avg = c(avg,summary(res_silhouette[[i]])[["avg.width"]])
}
dev.new()
pdf(file="plots/kmedoids_avg_width.pdf", width = 10, height = 5)
RANGE = 1:25
VALUE = avg
ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("Average Silhouette Width") + xlab("Number of Clusters") + ylab("Average Silhouette Width (100 restarts)")
dev.off()

# save results
kmedoids_all_clustering = res_clustering
kmedoids_all_silhouette = res_silhouette
save(kmedoids_all_clustering,file="results/kmedoids_all_clustering.RData")
save(kmedoids_all_silhouette,file="results/kmedoids_all_silhouette.RData")
