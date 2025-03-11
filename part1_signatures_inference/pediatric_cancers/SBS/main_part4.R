# load the required libraries and sources
library("ggpubr")
library("lsa")

# load the data
load("final_data/trinucleotides_counts.RData")
load("catalogues_assignments/COSMIC_goodness_fit.RData")
load("catalogues_assignments/COSMIC_assignment.RData")
load("catalogues_assignments/reference_signatures_goodness_fit.RData")
load("catalogues_assignments/reference_signatures_assignment.RData")

# consider each ingremental solution
SIGNATURE = NULL
PREDICTIONS = NULL
UNEXPLAINED_MUTATIONS = NULL
UNEXPLAINED_MUTATIONS_NORM = NULL
load("results/signatures_ordering.RData")
for(i in 1:length(signatures_ordering)) {
    load(paste0("results/incremental_assignment_",i,".RData"))
    curr_signature = rep(signatures_ordering[i],nrow(trinucleotides_counts))
    curr_predictions = ((incremental_assignment$alpha%*%incremental_assignment$beta)+incremental_assignment$unexplained_mutations)
    curr_cos = NULL
    for(j in rownames(trinucleotides_counts)) {
        curr_cos = c(curr_cos,as.numeric(cosine(curr_predictions[j,],trinucleotides_counts[j,])))
    }
    curr_unexplained_mutations = as.numeric(incremental_assignment$unexplained_mutations[rownames(trinucleotides_counts)])
    curr_unexplained_mutations_norm = as.numeric(curr_unexplained_mutations/rowSums(trinucleotides_counts))
    SIGNATURE = c(SIGNATURE,curr_signature)
    PREDICTIONS = c(PREDICTIONS,curr_cos)
    UNEXPLAINED_MUTATIONS = c(UNEXPLAINED_MUTATIONS,curr_unexplained_mutations)
    UNEXPLAINED_MUTATIONS_NORM = c(UNEXPLAINED_MUTATIONS_NORM,curr_unexplained_mutations_norm)
}
plot_data = data.frame(SIGNATURES=SIGNATURE,PREDICTIONS=PREDICTIONS,UNEXPLAINED_MUTATIONS=UNEXPLAINED_MUTATIONS,UNEXPLAINED_MUTATIONS_NORM=UNEXPLAINED_MUTATIONS_NORM,check.names=FALSE,stringsAsFactors=FALSE)

# make plots
dev.new()
pdf(file="plots/incremental_predictions_boxplot.pdf", width = 15, height = 10)
ggboxplot(plot_data,x="SIGNATURES",y="PREDICTIONS",fill="darkorange2",title="Predictions (Cosine Similarity)") + xlab("Signatures") + ylab("Predictions") + guides(fill=guide_legend(title="")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
dev.new()
pdf(file="plots/incremental_predictions_violin.pdf", width = 15, height = 10)
ggviolin(plot_data,x="SIGNATURES",y="PREDICTIONS",fill="darkorange2",title="Predictions (Cosine Similarity)") + xlab("Signatures") + ylab("Predictions") + guides(fill=guide_legend(title="")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
dev.new()
pdf(file="plots/incremental_num_mutations_boxplot.pdf", width = 15, height = 10)
ggboxplot(plot_data,x="SIGNATURES",y="UNEXPLAINED_MUTATIONS",fill="darkorange2",title="Unexplained Mutations (Raw Counts)") + xlab("Signatures") + ylab("Number of unexplained mutations") + guides(fill=guide_legend(title="")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
dev.new()
pdf(file="plots/incremental_num_mutations_violin.pdf", width = 15, height = 10)
ggviolin(plot_data,x="SIGNATURES",y="UNEXPLAINED_MUTATIONS",fill="darkorange2",title="Unexplained Mutations (Raw Counts)") + xlab("Signatures") + ylab("Number of unexplained mutations") + guides(fill=guide_legend(title="")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
dev.new()
pdf(file="plots/incremental_proportion_mutations_boxplot.pdf", width = 15, height = 10)
ggboxplot(plot_data,x="SIGNATURES",y="UNEXPLAINED_MUTATIONS_NORM",fill="darkorange2",title="Unexplained Mutations (Normalized)") + xlab("Signatures") + ylab("Percentage  of unexplained mutations") + guides(fill=guide_legend(title="")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
dev.new()
pdf(file="plots/incremental_proportion_mutations_violin.pdf", width = 15, height = 10)
ggviolin(plot_data,x="SIGNATURES",y="UNEXPLAINED_MUTATIONS_NORM",fill="darkorange2",title="Unexplained Mutations (Normalized)") + xlab("Signatures") + ylab("Percentage  of unexplained mutations") + guides(fill=guide_legend(title="")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# build line plots data
LINE_SIGNATURE = NULL
LINE_PREDICTIONS = NULL
LINE_UNEXPLAINED_MUTATIONS = NULL
LINE_UNEXPLAINED_MUTATIONS_NORM = NULL
for(i in unique(plot_data$SIGNATURES)) {
    LINE_SIGNATURE = c(LINE_SIGNATURE,i)
    LINE_PREDICTIONS = c(LINE_PREDICTIONS,(length(which(PREDICTIONS[which(SIGNATURE==i)]>0.95))/nrow(trinucleotides_counts)))
    LINE_UNEXPLAINED_MUTATIONS = c(LINE_UNEXPLAINED_MUTATIONS,median(UNEXPLAINED_MUTATIONS[which(SIGNATURE==i)],na.rm=TRUE))
    LINE_UNEXPLAINED_MUTATIONS_NORM = c(LINE_UNEXPLAINED_MUTATIONS_NORM,mean(UNEXPLAINED_MUTATIONS_NORM[which(SIGNATURE==i)],na.rm=TRUE))
}

# make plots
VALUE = LINE_PREDICTIONS
RANGE = LINE_SIGNATURE
COSMIC = (length(which(COSMIC_goodness_fit>0.95))/nrow(trinucleotides_counts))
reference = (length(which(reference_signatures_goodness_fit>0.95))/nrow(trinucleotides_counts))
plot_data = data.frame(VALUE=VALUE,RANGE=RANGE,check.names=FALSE,stringsAsFactors=FALSE)
dev.off()
dev.new()
pdf(file="plots/incremental_predictions_lineplot.pdf", width = 15, height = 10)
ggline(plot_data,x="RANGE",y="VALUE",fill="darkorange2",title="Goodness of predictions (% Cosine Similarity >0.95)") + xlab("Signatures") + ylab("Goodness of predictions") + guides(fill=guide_legend(title="")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_hline(yintercept=reference,colour="red",linetype="longdash") + geom_hline(yintercept=COSMIC,colour="blue",linetype="longdash") + ylim(0,1)
dev.off()
VALUE = LINE_UNEXPLAINED_MUTATIONS
RANGE = LINE_SIGNATURE
COSMIC = median(as.numeric(COSMIC_assignment$unexplained_mutations),na.rm=TRUE)
reference = median(as.numeric(reference_signatures_assignment$unexplained_mutations),na.rm=TRUE)
plot_data = data.frame(VALUE=VALUE,RANGE=RANGE,check.names=FALSE,stringsAsFactors=FALSE)
dev.off()
dev.new()
pdf(file="plots/incremental_num_mutations_lineplot.pdf", width = 15, height = 10)
ggline(plot_data,x="RANGE",y="VALUE",fill="darkorange2",title="Unexplained Mutations (Median Raw Counts)") + xlab("Signatures") + ylab("Median number of unexplained mutations") + guides(fill=guide_legend(title="")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_hline(yintercept=reference,colour="red",linetype="longdash") + geom_hline(yintercept=COSMIC,colour="blue",linetype="longdash")
dev.off()
VALUE = LINE_UNEXPLAINED_MUTATIONS_NORM
RANGE = LINE_SIGNATURE
COSMIC = as.numeric(COSMIC_assignment$unexplained_mutations[rownames(trinucleotides_counts)])
COSMIC = mean(as.numeric(COSMIC/rowSums(trinucleotides_counts)),na.rm=TRUE)
reference = as.numeric(reference_signatures_assignment$unexplained_mutations[rownames(trinucleotides_counts)])
reference = mean(as.numeric(reference/rowSums(trinucleotides_counts)),na.rm=TRUE)
plot_data = data.frame(VALUE=VALUE,RANGE=RANGE,check.names=FALSE,stringsAsFactors=FALSE)
dev.off()
dev.new()
pdf(file="plots/incremental_proportion_mutations_lineplot.pdf", width = 15, height = 10)
ggline(plot_data,x="RANGE",y="VALUE",fill="darkorange2",title="Unexplained Mutations (Mean Proportions)") + xlab("Signatures") + ylab("Mean percentage  of unexplained mutations") + guides(fill=guide_legend(title="")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + geom_hline(yintercept=reference,colour="red",linetype="longdash") + geom_hline(yintercept=COSMIC,colour="blue",linetype="longdash")
dev.off()
