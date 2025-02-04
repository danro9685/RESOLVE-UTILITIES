# load the required libraries and sources
library("glmnet")
library("parallel")
library("survival")

# load the data
load(file="final_data/clinical_data.RData")
load(file="final_data/cosine_similarity.RData")
load(file="final_data/signatures_assignment.RData")

# consider only samples with high goodness of fit
clinical_data = clinical_data[which(clinical_data$PATIENT%in%sort(unique(names(which(cosine_similarity>0.95))))),]

# process the survival data
valid_samples = which(!is.na(clinical_data$AGE)&!is.na(clinical_data$OS_STATUS)&!is.na(clinical_data$OS_YEARS))
clinical_data = clinical_data[valid_samples,]
clinical_data = clinical_data[which(clinical_data$AGE>=18),]
survival_data = clinical_data[,c("PATIENT","TISSUE","AGE","OS_YEARS","OS_STATUS")]
survival_data[,"OS_STATUS"][which(survival_data[,"OS_STATUS"]=="DECEASED")] = 1
survival_data[,"OS_STATUS"][which(survival_data[,"OS_STATUS"]=="ALIVE")] = 0
survival_data[,"OS_STATUS"][which(as.numeric(survival_data[,"OS_YEARS"])<0.08|as.numeric(survival_data[,"AGE"])>80)] = NA
survival_data[,"OS_YEARS"][which(as.numeric(survival_data[,"OS_YEARS"])<0.08|as.numeric(survival_data[,"AGE"])>80)] = NA
survival_data[,"OS_STATUS"][which(as.numeric(survival_data[,"OS_YEARS"])>5)] = 0
survival_data[,"OS_YEARS"][which(as.numeric(survival_data[,"OS_YEARS"])>5)] = 5
clinical_data = survival_data
clinical_data$AGE = NULL
invalid_samples = unique(which(is.na(clinical_data),arr.ind=TRUE)[,"row"])
if(length(invalid_samples)>0) {
    clinical_data = clinical_data[-invalid_samples,]
}
alpha = signatures_assignment$alpha
colnames(alpha) = c("SBS1","SBS2","SBS3","SBS4","SBS5","SBS7a","SBS7b","SBS8","SBS9","SBS10a","SBS10d","SBS11","SBS13","SBS14","SBS15","SBS17","SBS18","SBS19","SBS20","SBS22","SBS23","SBS26","SBS28","SBS31","SBS32","SBS44","SBS88","SBS92","SBS97")
raw_alpha = alpha[clinical_data$PATIENT,]
normalized_alpha = raw_alpha/rowSums(raw_alpha)

# perform the association analysis
analysis_data_normalized_exposures = list()
cv.fit_normalized_exposures = list()
significant_features_normalized_exposures = list()
set.seed(12345)
analysis_data = cbind(clinical_data,normalized_alpha)
rownames(analysis_data) = analysis_data$PATIENT
analysis_data$PATIENT = NULL
analysis_data$TISSUE = NULL
colnames(analysis_data)[1:2] = c("Times","Status")
analysis_data = data.frame(analysis_data)
x = as.matrix(analysis_data[,3:ncol(analysis_data)])
y = Surv(as.numeric(analysis_data$Times),as.numeric(analysis_data$Status))
cv.fit = cv.glmnet(x,y,family="cox",maxit=1000000000)
coeff_fit = coef(cv.fit,s=cv.fit$lambda.min)[colnames(analysis_data)[3:ncol(analysis_data)],1]
coeff_values = as.numeric(coeff_fit)
names(coeff_values) = names(coeff_fit)
significant_features = list()
significant_features[["Positive"]] = coeff_values[which(coeff_values>0)]
significant_features[["Negative"]] = coeff_values[which(coeff_values<0)]
analysis_data_normalized_exposures[["Pan-cancer"]] = analysis_data
cv.fit_normalized_exposures[["Pan-cancer"]] = cv.fit
significant_features_normalized_exposures[["Pan-cancer"]] = significant_features
analysis_data = NULL
cv.fit = NULL
significant_features = NULL
for(ct in names(which(table(clinical_data$TISSUE)>=10))) {
    set.seed(12345)
    valid_samples = clinical_data$PATIENT[which(clinical_data$TISSUE==ct)]
    analysis_data = cbind(clinical_data[which(clinical_data$PATIENT%in%valid_samples),],normalized_alpha[valid_samples,])
    rownames(analysis_data) = analysis_data$PATIENT
    analysis_data$PATIENT = NULL
    analysis_data$TISSUE = NULL
    colnames(analysis_data)[1:2] = c("Times","Status")
    analysis_data = data.frame(analysis_data)
    x = as.matrix(analysis_data[,3:ncol(analysis_data)])
    y = Surv(as.numeric(analysis_data$Times),as.numeric(analysis_data$Status))
    cv.fit = cv.glmnet(x,y,family="cox",maxit=1000000000)
    coeff_fit = coef(cv.fit,s=cv.fit$lambda.min)[colnames(analysis_data)[3:ncol(analysis_data)],1]
    coeff_values = as.numeric(coeff_fit)
    names(coeff_values) = names(coeff_fit)
    significant_features = list()
    significant_features[["Positive"]] = coeff_values[which(coeff_values>0)]
    significant_features[["Negative"]] = coeff_values[which(coeff_values<0)]
    analysis_data_normalized_exposures[[ct]] = analysis_data
    cv.fit_normalized_exposures[[ct]] = cv.fit
    significant_features_normalized_exposures[[ct]] = significant_features
    analysis_data = NULL
    cv.fit = NULL
    significant_features = NULL
}

# save the results
inference_results_normalized_exposures = list(analysis_data=analysis_data_normalized_exposures,cv.fit=cv.fit_normalized_exposures,significant_features=significant_features_normalized_exposures)
save(inference_results_normalized_exposures,file="results/inference_results_normalized_exposures.RData")
