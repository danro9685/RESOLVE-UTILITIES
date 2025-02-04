# load the required libraries and sources
library("glmnet")
library("lsa")
library("parallel")
source("R/association.estimation_signatures_from_mutations.R")

# load and process the data
load(file="final_data/clinical_data.RData")
load(file="final_data/cosine_similarity.RData")
load(file="final_data/signatures_assignment.RData")
load(file="final_data/variants_matrix.RData")

# consider only samples with good goodness of fit
clinical_data = clinical_data[which(clinical_data$PATIENT%in%sort(unique(names(which(cosine_similarity>0.95))))),]

# process the mutations data
variants_matrix_driver = variants_matrix[which(variants_matrix$PATIENT%in%clinical_data$PATIENT&variants_matrix$CONSEQUENCE!="synonymous"),]
rownames(variants_matrix_driver) = 1:nrow(variants_matrix_driver)
clinical_data_driver = clinical_data[which(clinical_data$PATIENT%in%variants_matrix_driver$PATIENT),]
rownames(clinical_data_driver) = 1:nrow(clinical_data_driver)
samples = sort(unique(clinical_data_driver$PATIENT))
genes = sort(unique(variants_matrix_driver$GENE))
mutations = matrix(0, nrow = length(samples), ncol = length(genes))
rownames(mutations) = samples
colnames(mutations) = genes
for(i in 1:nrow(variants_matrix_driver)) {
    mutations[variants_matrix_driver$PATIENT[i],variants_matrix_driver$GENE[i]] = 1
}
mutations_driver = mutations[,sort(unique(names(which(colSums(mutations)>=10))))]

# process the signatures data
alpha = signatures_assignment$alpha
colnames(alpha) = c("SBS1","SBS2","SBS3","SBS4","SBS5","SBS7a","SBS7b","SBS8","SBS9","SBS10a","SBS10d","SBS11","SBS13","SBS14","SBS15","SBS17","SBS18","SBS19","SBS20","SBS22","SBS23","SBS26","SBS28","SBS31","SBS32","SBS44","SBS88","SBS92","SBS97")
raw_alpha_driver = alpha[clinical_data_driver$PATIENT,]
normalized_alpha_driver = raw_alpha_driver/rowSums(raw_alpha_driver)

# perform the analysis for functional mutations
analysis_data = list()
analysis_results = list()
for(i in c("Breast","Central Nervous System","Esophagus","Hematopoietic and Lymphoid Tissues","Liver","Pancreas","Prostate","Skin")) {
    cat(i,"\n")
    analysis_data[[i]] = list()
    curr_clinical_data = clinical_data_driver[which(clinical_data_driver$TISSUE==i),,drop=FALSE]
    rownames(curr_clinical_data) = 1:nrow(curr_clinical_data)
    curr_mutations = mutations_driver[curr_clinical_data$PATIENT,,drop=FALSE]
    valid_genes = sort(colSums(curr_mutations)[which(colSums(curr_mutations)>=3)],decreasing=TRUE)
    if(length(valid_genes)>100) {
        valid_genes = valid_genes[1:100]
    }
    valid_genes = sort(unique(names(valid_genes)))
    curr_mutations = curr_mutations[,valid_genes,drop=FALSE]
    curr_variants_matrix = variants_matrix_driver[which(variants_matrix_driver$PATIENT%in%clinical_data$PATIENT&variants_matrix_driver$GENE%in%colnames(curr_mutations)),,drop=FALSE]
    rownames(curr_variants_matrix) = 1:nrow(curr_variants_matrix)
    curr_raw_alpha = raw_alpha_driver[rownames(curr_mutations),,drop=FALSE]
    curr_normalized_alpha = normalized_alpha_driver[rownames(curr_mutations),,drop=FALSE]
    analysis_data[[i]][["clinical_data"]] = curr_clinical_data
    analysis_data[[i]][["mutations"]] = curr_mutations
    analysis_data[[i]][["variants_matrix"]] = curr_variants_matrix
    analysis_data[[i]][["raw_alpha"]] = curr_raw_alpha
    analysis_data[[i]][["normalized_alpha"]] = curr_normalized_alpha
    set.seed(12345)
    analysis_results[[i]] = associationEstimation(alterations = curr_mutations, signatures = curr_normalized_alpha)
}
analysis_data_driver = analysis_data
analysis_results_driver = analysis_results

# save the results
save(analysis_data_driver,file="results/signatures_from_mutations/analysis_data_driver.RData")
save(analysis_results_driver,file="results/signatures_from_mutations/analysis_results_driver.RData")
