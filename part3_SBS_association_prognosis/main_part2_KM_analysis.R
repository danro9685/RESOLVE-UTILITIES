# load the required libraries
library("dynamicTreeCut")
library("glmnet")
library("RColorBrewer")
library("survival")
library("survminer")

# load results for the analysis considering normalized exposures
dev.new(width=15,height=10)
load(file="results/inference_results_normalized_exposures.RData")
inference_results = inference_results_normalized_exposures
inference_results_normalized_exposures = NULL
clustering_results = list()

# perform the analysis for the Esophagus analysis
curr_configuration = "Esophagus"
print(curr_configuration)
print(inference_results$significant_features[[curr_configuration]])
set.seed(12345)
data = inference_results$analysis_data[[curr_configuration]]
analyzed_samples = rownames(data)
fit = inference_results$cv.fit[[curr_configuration]]
x = as.numeric(coef(fit,s=fit$lambda.min))
y = data[,3:ncol(data)]
beta = as.numeric(x%*%t(y))
dist_beta = dist(beta)
dendogram = hclust(dist_beta)
Clustering = cutreeDynamic(dendro=dendogram,minClusterSize=50,distM=as.matrix(dist_beta),method="hybrid")
data$Clusters = as.numeric(Clustering)
data$Clusters[which(as.numeric(Clustering)==1)] = 2
data$Clusters[which(as.numeric(Clustering)==2)] = 1
data$Clusters[which(as.numeric(Clustering)==3)] = 3
valid_samples = which(data$Clusters!=0)
if(length(valid_samples)>=10) {
    data = data[valid_samples,]
    rownames(data) = 1:nrow(data)
    Cluster = paste0("C",data$Clusters)
    names(Cluster) = analyzed_samples[valid_samples]
    data$Clusters = NULL
    clustering_results[[curr_configuration]] = Cluster
    print(ggsurvplot(survfit(Surv(as.numeric(data$Times),as.numeric(data$Status))~Cluster),
            data = data,
            xlab = "Months",
            ylab = "Overall survival probability",
            palette = "Dark2",
            mark.time = TRUE,
            pval = TRUE,
            risk.table = TRUE,
            ggtheme = theme_bw(),
            title = paste0(curr_configuration," (Kaplan-Meier)"),
            font.main = 18,
            font.x = 18,
            font.y = 18,
            font.caption = 18,
            font.legend = 18,
            font.tickslab = 18))
}

# perform the analysis for the Pancreas analysis
curr_configuration = "Pancreas"
print(curr_configuration)
print(inference_results$significant_features[[curr_configuration]])
set.seed(12345)
data = inference_results$analysis_data[[curr_configuration]]
analyzed_samples = rownames(data)
fit = inference_results$cv.fit[[curr_configuration]]
x = as.numeric(coef(fit,s=fit$lambda.min))
y = data[,3:ncol(data)]
beta = as.numeric(x%*%t(y))
dist_beta = dist(beta)
dendogram = hclust(dist_beta)
Clustering = cutreeDynamic(dendro=dendogram,minClusterSize=50,distM=as.matrix(dist_beta),method="hybrid")
data$Clusters = as.numeric(Clustering)
data$Clusters[which(as.numeric(Clustering)==1)] = 2
data$Clusters[which(as.numeric(Clustering)==2)] = 3
data$Clusters[which(as.numeric(Clustering)==3)] = 1
valid_samples = which(data$Clusters!=0)
if(length(valid_samples)>=10) {
    data = data[valid_samples,]
    rownames(data) = 1:nrow(data)
    Cluster = paste0("C",data$Clusters)
    names(Cluster) = analyzed_samples[valid_samples]
    data$Clusters = NULL
    clustering_results[[curr_configuration]] = Cluster
    print(ggsurvplot(survfit(Surv(as.numeric(data$Times),as.numeric(data$Status))~Cluster),
            data = data,
            xlab = "Months",
            ylab = "Overall survival probability",
            palette = "Dark2",
            mark.time = TRUE,
            pval = TRUE,
            risk.table = TRUE,
            ggtheme = theme_bw(),
            title = paste0(curr_configuration," (Kaplan-Meier)"),
            font.main = 18,
            font.x = 18,
            font.y = 18,
            font.caption = 18,
            font.legend = 18,
            font.tickslab = 18))
}

# perform the analysis for the Prostate analysis
curr_configuration = "Prostate"
print(curr_configuration)
print(inference_results$significant_features[[curr_configuration]])
set.seed(12345)
data = inference_results$analysis_data[[curr_configuration]]
analyzed_samples = rownames(data)
fit = inference_results$cv.fit[[curr_configuration]]
x = as.numeric(coef(fit,s=fit$lambda.min))
y = data[,3:ncol(data)]
beta = as.numeric(x%*%t(y))
dist_beta = dist(beta)
dendogram = hclust(dist_beta)
Clustering = cutreeDynamic(dendro=dendogram,minClusterSize=50,distM=as.matrix(dist_beta),method="hybrid")
data$Clusters = as.numeric(Clustering)
data$Clusters[which(as.numeric(Clustering)==1)] = 1
data$Clusters[which(as.numeric(Clustering)==2)] = 2
valid_samples = which(data$Clusters!=0)
if(length(valid_samples)>=10) {
    data = data[valid_samples,]
    rownames(data) = 1:nrow(data)
    Cluster = paste0("C",data$Clusters)
    names(Cluster) = analyzed_samples[valid_samples]
    data$Clusters = NULL
    clustering_results[[curr_configuration]] = Cluster
    print(ggsurvplot(survfit(Surv(as.numeric(data$Times),as.numeric(data$Status))~Cluster),
            data = data,
            xlab = "Months",
            ylab = "Overall survival probability",
            palette = "Dark2",
            mark.time = TRUE,
            pval = TRUE,
            risk.table = TRUE,
            ggtheme = theme_bw(),
            title = paste0(curr_configuration," (Kaplan-Meier)"),
            font.main = 18,
            font.x = 18,
            font.y = 18,
            font.caption = 18,
            font.legend = 18,
            font.tickslab = 18))
}

# perform the analysis for the Skin analysis
curr_configuration = "Skin"
print(curr_configuration)
print(inference_results$significant_features[[curr_configuration]])
set.seed(12345)
data = inference_results$analysis_data[[curr_configuration]]
analyzed_samples = rownames(data)
fit = inference_results$cv.fit[[curr_configuration]]
x = as.numeric(coef(fit,s=fit$lambda.min))
y = data[,3:ncol(data)]
beta = as.numeric(x%*%t(y))
dist_beta = dist(beta)
dendogram = hclust(dist_beta)
Clustering = cutreeDynamic(dendro=dendogram,minClusterSize=50,distM=as.matrix(dist_beta),method="hybrid")
data$Clusters = as.numeric(Clustering)
data$Clusters[which(as.numeric(Clustering)==1)] = 1
data$Clusters[which(as.numeric(Clustering)==2)] = 2
data$Clusters[which(as.numeric(Clustering)==0)] = 3
valid_samples = which(data$Clusters!=0)
if(length(valid_samples)>=10) {
    data = data[valid_samples,]
    rownames(data) = 1:nrow(data)
    Cluster = paste0("C",data$Clusters)
    names(Cluster) = analyzed_samples[valid_samples]
    data$Clusters = NULL
    clustering_results[[curr_configuration]] = Cluster
    print(ggsurvplot(survfit(Surv(as.numeric(data$Times),as.numeric(data$Status))~Cluster),
            data = data,
            xlab = "Months",
            ylab = "Overall survival probability",
            palette = "Dark2",
            mark.time = TRUE,
            pval = TRUE,
            risk.table = TRUE,
            ggtheme = theme_bw(),
            title = paste0(curr_configuration," (Kaplan-Meier)"),
            font.main = 18,
            font.x = 18,
            font.y = 18,
            font.caption = 18,
            font.legend = 18,
            font.tickslab = 18))
}

# save results
clustering_results_normalized_exposures = clustering_results
save(clustering_results_normalized_exposures,file="results/clustering_results_normalized_exposures.RData")
