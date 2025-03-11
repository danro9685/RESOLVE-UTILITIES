# load the required libraries and sources
library("ggpubr")

# load the data
load("final_data/cna_counts.RData")
counts = cna_counts
load("final_results/goodness_fit.RData")
load("final_results/signatures_assignment.RData")
load("final_results/signatures_significance.RData")
load("catalogues_assignments/COSMIC_goodness_fit.RData")
load("catalogues_assignments/COSMIC_assignment.RData")

# make plots data
IDs = sort(unique(rownames(counts)))
SAMPLES = c(IDs,IDs,IDs,IDs)
CATALOGUE = NULL
PREDICTIONS = NULL
UNEXPLAINED_MUTATIONS = NULL
UNEXPLAINED_MUTATIONS_NORM = NULL
CATALOGUE = c(CATALOGUE,rep("RESOLVE",length(IDs)))
PREDICTIONS = c(PREDICTIONS,as.numeric(goodness_fit[IDs]))
unexplained_mut = as.numeric(signatures_assignment$unexplained_mutations[IDs])
unexplained_mut_norm = as.numeric(unexplained_mut/rowSums(counts)[IDs])
UNEXPLAINED_MUTATIONS = c(UNEXPLAINED_MUTATIONS,unexplained_mut)
UNEXPLAINED_MUTATIONS_NORM = c(UNEXPLAINED_MUTATIONS_NORM,unexplained_mut_norm)
CATALOGUE = c(CATALOGUE,rep("RESOLVE (Significant Exposure)",length(IDs)))
curr_predictions = as.numeric(signatures_significance$goodness_fit[IDs])
curr_predictions[which(is.na(curr_predictions))] = 0
PREDICTIONS = c(PREDICTIONS,curr_predictions)
unexplained_mut = as.numeric(signatures_significance$unexplained_mutations[IDs])
unexplained_mut_norm = as.numeric(unexplained_mut/rowSums(counts)[IDs])
UNEXPLAINED_MUTATIONS = c(UNEXPLAINED_MUTATIONS,unexplained_mut)
UNEXPLAINED_MUTATIONS_NORM = c(UNEXPLAINED_MUTATIONS_NORM,unexplained_mut_norm)
CATALOGUE = c(CATALOGUE,rep("COSMIC Signatures",length(IDs)))
PREDICTIONS = c(PREDICTIONS,as.numeric(COSMIC_goodness_fit[IDs]))
unexplained_mut = as.numeric(COSMIC_assignment$unexplained_mutations[IDs])
unexplained_mut_norm = as.numeric(unexplained_mut/rowSums(counts)[IDs])
UNEXPLAINED_MUTATIONS = c(UNEXPLAINED_MUTATIONS,unexplained_mut)
UNEXPLAINED_MUTATIONS_NORM = c(UNEXPLAINED_MUTATIONS_NORM,unexplained_mut_norm)
plot_data = data.frame(CATALOGUES=CATALOGUE,PREDICTIONS=PREDICTIONS,UNEXPLAINED_MUTATIONS=UNEXPLAINED_MUTATIONS,UNEXPLAINED_MUTATIONS_NORM=UNEXPLAINED_MUTATIONS_NORM,check.names=FALSE,stringsAsFactors=FALSE)

# make plots
dev.new()
pdf(file="plots/catalogues_predictions_boxplot.pdf", width = 10, height = 5)
ggboxplot(plot_data,x="CATALOGUES",y="PREDICTIONS",fill="CATALOGUES",palette="Dark2",title="Predictions (Cosine Similarity)") + xlab("Catalogues") + ylab("Predictions") + guides(fill=guide_legend(title=""))
dev.off()
dev.new()
pdf(file="plots/catalogues_predictions_violin.pdf", width = 10, height = 5)
ggviolin(plot_data,x="CATALOGUES",y="PREDICTIONS",fill="CATALOGUES",palette="Dark2",title="Predictions (Cosine Similarity)") + xlab("Catalogues") + ylab("Predictions") + guides(fill=guide_legend(title=""))
dev.off()
dev.new()
pdf(file="plots/catalogues_num_mutations_boxplot.pdf", width = 10, height = 5)
ggboxplot(plot_data,x="CATALOGUES",y="UNEXPLAINED_MUTATIONS",fill="CATALOGUES",palette="Dark2",title="Unexplained Mutations (Raw Counts)") + xlab("Catalogues") + ylab("Number of unexplained mutations") + guides(fill=guide_legend(title=""))
dev.off()
dev.new()
pdf(file="plots/catalogues_num_mutations_violin.pdf", width = 10, height = 5)
ggviolin(plot_data,x="CATALOGUES",y="UNEXPLAINED_MUTATIONS",fill="CATALOGUES",palette="Dark2",title="Unexplained Mutations (Raw Counts)") + xlab("Catalogues") + ylab("Number of unexplained mutations") + guides(fill=guide_legend(title=""))
dev.off()
dev.new()
pdf(file="plots/catalogues_proportion_mutations_boxplot.pdf", width = 10, height = 5)
ggboxplot(plot_data,x="CATALOGUES",y="UNEXPLAINED_MUTATIONS_NORM",fill="CATALOGUES",palette="Dark2",title="Unexplained Mutations (Normalized)") + xlab("Catalogues") + ylab("Percentage  of unexplained mutations") + guides(fill=guide_legend(title=""))
dev.off()
dev.new()
pdf(file="plots/catalogues_proportion_mutations_violin.pdf", width = 10, height = 5)
ggviolin(plot_data,x="CATALOGUES",y="UNEXPLAINED_MUTATIONS_NORM",fill="CATALOGUES",palette="Dark2",title="Unexplained Mutations (Normalized)") + xlab("Catalogues") + ylab("Percentage  of unexplained mutations") + guides(fill=guide_legend(title=""))
dev.off()
