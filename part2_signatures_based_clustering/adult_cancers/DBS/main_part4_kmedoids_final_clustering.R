# load the required libraries
library("cluster")
library("factoextra")
library("FactoMineR")
library("ggpubr")
library("NbClust")
library("RColorBrewer")

# load the results
load(file="final_results/signatures_significance.RData")
alpha = signatures_significance$alpha
alpha = alpha[sort(unique(names(which(signatures_significance$goodness_fit>0.90)))),]
normalized_alpha = (alpha/rowSums(alpha))
load(file="results/signatures_distance.RData")
sig_dist = as.dist(signatures_distance)
load(file="results/kmedoids_all_clustering.RData")
load(file="results/kmedoids_all_silhouette.RData")

# evaluate silhouette for each clustering result
avg = 0
for(i in 2:25) {
    res_silhouette = silhouette(x = kmedoids_all_clustering[[i]], dist = sig_dist)
    avg = c(avg,summary(res_silhouette)[["avg.width"]])
}

# save the final clustering results
best_res = 10
kmedoids_clustering = kmedoids_all_clustering[[best_res]]

# make silhouette plot for the best clustering results
kmedoids_silhouette = silhouette(x = kmedoids_clustering, dist = sig_dist)
dev.new()
pdf(file="final_results/kmedoids_silhouette.pdf", width = 10, height = 5)
fviz_silhouette(kmedoids_silhouette)
dev.off()

# plot signatures per cluster
plot_data_sigs = NULL
plot_data_sigs_name = NULL
plot_data_sigs_red = NULL
plot_data_clusts = NULL
for(i in names(kmedoids_clustering)) {
    curr_sig = as.numeric(normalized_alpha[i,])
    curr_sig_name = names(normalized_alpha[i,])
    curr_clust = rep(kmedoids_clustering[[i]],length(curr_sig))
    plot_data_sigs = c(plot_data_sigs,curr_sig)
    plot_data_sigs_name = c(plot_data_sigs_name,curr_sig_name)
    plot_data_sigs_red = c(plot_data_sigs_red,paste0("D",1:length(curr_sig)))
    plot_data_clusts = c(plot_data_clusts,curr_clust)
}
plot_data = data.frame(VALUE = plot_data_sigs, SIGNATURE = factor(plot_data_sigs_name,levels=unique(plot_data_sigs_name)), CLUSTER = plot_data_clusts, SIGNATURE_RED = plot_data_sigs_red)
sig_sign = NULL
for(i in sort(unique(plot_data$CLUSTER))) {
    curr_sig_sign = NULL
    for(j in sort(unique(plot_data$SIGNATURE))) {
        curr_sig_sign = c(curr_sig_sign,mean(plot_data$VALUE[which(plot_data$CLUSTER==i&plot_data$SIGNATURE==j)],na.rm=TRUE))
    }
    names(curr_sig_sign) = sort(unique(plot_data$SIGNATURE))
    sig_sign = c(sig_sign,names(curr_sig_sign)[which(curr_sig_sign>0.05)])
}
sig_sign = sort(unique(sig_sign))
plot_data = plot_data[which(plot_data$SIGNATURE%in%sig_sign),]
dev.new()
pdf(file="final_results/signatures_clusters.pdf", width = 15, height = 5)
plot_data_tmp = plot_data[which(plot_data$CLUSTER==1),]
curr_sil = round(summary(kmedoids_silhouette)[["clus.avg.widths"]][[1]],2)
curr_siz = table(kmedoids_clustering)[[1]]
ggboxplot(plot_data_tmp,x="SIGNATURE_RED",y="VALUE",fill="SIGNATURE",palette=colorRampPalette(brewer.pal(name="Paired",n=11))(length(sig_sign)),title=paste0("Signatures Impact - Cluster 1 (Samples: ",curr_siz,", Silhouette: ",curr_sil,")")) + xlab("Signatures") + ylab("Percentage of Mutations") + guides(fill=guide_legend(title=""))
plot_data_tmp = plot_data[which(plot_data$CLUSTER==2),]
curr_sil = round(summary(kmedoids_silhouette)[["clus.avg.widths"]][[2]],2)
curr_siz = table(kmedoids_clustering)[[2]]
ggboxplot(plot_data_tmp,x="SIGNATURE_RED",y="VALUE",fill="SIGNATURE",palette=colorRampPalette(brewer.pal(name="Paired",n=11))(length(sig_sign)),title=paste0("Signatures Impact - Cluster 2 (Samples: ",curr_siz,", Silhouette: ",curr_sil,")")) + xlab("Signatures") + ylab("Percentage of Mutations") + guides(fill=guide_legend(title=""))
plot_data_tmp = plot_data[which(plot_data$CLUSTER==3),]
curr_sil = round(summary(kmedoids_silhouette)[["clus.avg.widths"]][[3]],2)
curr_siz = table(kmedoids_clustering)[[3]]
ggboxplot(plot_data_tmp,x="SIGNATURE_RED",y="VALUE",fill="SIGNATURE",palette=colorRampPalette(brewer.pal(name="Paired",n=11))(length(sig_sign)),title=paste0("Signatures Impact - Cluster 3 (Samples: ",curr_siz,", Silhouette: ",curr_sil,")")) + xlab("Signatures") + ylab("Percentage of Mutations") + guides(fill=guide_legend(title=""))
plot_data_tmp = plot_data[which(plot_data$CLUSTER==4),]
curr_sil = round(summary(kmedoids_silhouette)[["clus.avg.widths"]][[4]],2)
curr_siz = table(kmedoids_clustering)[[4]]
ggboxplot(plot_data_tmp,x="SIGNATURE_RED",y="VALUE",fill="SIGNATURE",palette=colorRampPalette(brewer.pal(name="Paired",n=11))(length(sig_sign)),title=paste0("Signatures Impact - Cluster 4 (Samples: ",curr_siz,", Silhouette: ",curr_sil,")")) + xlab("Signatures") + ylab("Percentage of Mutations") + guides(fill=guide_legend(title=""))
plot_data_tmp = plot_data[which(plot_data$CLUSTER==5),]
curr_sil = round(summary(kmedoids_silhouette)[["clus.avg.widths"]][[5]],2)
curr_siz = table(kmedoids_clustering)[[5]]
ggboxplot(plot_data_tmp,x="SIGNATURE_RED",y="VALUE",fill="SIGNATURE",palette=colorRampPalette(brewer.pal(name="Paired",n=11))(length(sig_sign)),title=paste0("Signatures Impact - Cluster 5 (Samples: ",curr_siz,", Silhouette: ",curr_sil,")")) + xlab("Signatures") + ylab("Percentage of Mutations") + guides(fill=guide_legend(title=""))
plot_data_tmp = plot_data[which(plot_data$CLUSTER==6),]
curr_sil = round(summary(kmedoids_silhouette)[["clus.avg.widths"]][[6]],2)
curr_siz = table(kmedoids_clustering)[[6]]
ggboxplot(plot_data_tmp,x="SIGNATURE_RED",y="VALUE",fill="SIGNATURE",palette=colorRampPalette(brewer.pal(name="Paired",n=11))(length(sig_sign)),title=paste0("Signatures Impact - Cluster 6 (Samples: ",curr_siz,", Silhouette: ",curr_sil,")")) + xlab("Signatures") + ylab("Percentage of Mutations") + guides(fill=guide_legend(title=""))
plot_data_tmp = plot_data[which(plot_data$CLUSTER==7),]
curr_sil = round(summary(kmedoids_silhouette)[["clus.avg.widths"]][[7]],2)
curr_siz = table(kmedoids_clustering)[[7]]
ggboxplot(plot_data_tmp,x="SIGNATURE_RED",y="VALUE",fill="SIGNATURE",palette=colorRampPalette(brewer.pal(name="Paired",n=11))(length(sig_sign)),title=paste0("Signatures Impact - Cluster 7 (Samples: ",curr_siz,", Silhouette: ",curr_sil,")")) + xlab("Signatures") + ylab("Percentage of Mutations") + guides(fill=guide_legend(title=""))
plot_data_tmp = plot_data[which(plot_data$CLUSTER==8),]
curr_sil = round(summary(kmedoids_silhouette)[["clus.avg.widths"]][[8]],2)
curr_siz = table(kmedoids_clustering)[[8]]
ggboxplot(plot_data_tmp,x="SIGNATURE_RED",y="VALUE",fill="SIGNATURE",palette=colorRampPalette(brewer.pal(name="Paired",n=11))(length(sig_sign)),title=paste0("Signatures Impact - Cluster 8 (Samples: ",curr_siz,", Silhouette: ",curr_sil,")")) + xlab("Signatures") + ylab("Percentage of Mutations") + guides(fill=guide_legend(title=""))
plot_data_tmp = plot_data[which(plot_data$CLUSTER==9),]
curr_sil = round(summary(kmedoids_silhouette)[["clus.avg.widths"]][[9]],2)
curr_siz = table(kmedoids_clustering)[[9]]
ggboxplot(plot_data_tmp,x="SIGNATURE_RED",y="VALUE",fill="SIGNATURE",palette=colorRampPalette(brewer.pal(name="Paired",n=11))(length(sig_sign)),title=paste0("Signatures Impact - Cluster 9 (Samples: ",curr_siz,", Silhouette: ",curr_sil,")")) + xlab("Signatures") + ylab("Percentage of Mutations") + guides(fill=guide_legend(title=""))
plot_data_tmp = plot_data[which(plot_data$CLUSTER==10),]
curr_sil = round(summary(kmedoids_silhouette)[["clus.avg.widths"]][[10]],2)
curr_siz = table(kmedoids_clustering)[[10]]
ggboxplot(plot_data_tmp,x="SIGNATURE_RED",y="VALUE",fill="SIGNATURE",palette=colorRampPalette(brewer.pal(name="Paired",n=11))(length(sig_sign)),title=paste0("Signatures Impact - Cluster 10 (Samples: ",curr_siz,", Silhouette: ",curr_sil,")")) + xlab("Signatures") + ylab("Percentage of Mutations") + guides(fill=guide_legend(title=""))
dev.off()

# save results
save(kmedoids_clustering,file="final_results/kmedoids_clustering.RData")
save(kmedoids_silhouette,file="final_results/kmedoids_silhouette.RData")
