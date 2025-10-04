# ULSL: Unified Latent and Similarity Learning for Robust Multi-omics Cancer Subtype Identification
## 
we can install and run _ULSL_ as follows:
```R
	library("devtools")
	devtools::install_github('codelzy-01/ULSL-1')
	# load ULSL library
	library("ULSL")
```
### Environment:

[](https://github.com/yushanqiu/MDICC#environment)

Anaconda, R>= 4.5.0;

### Connect Anaconda environment.
```R
library(reticulate)
use_python("C:\Users\anaconda3\python.exe", required = TRUE)
reticulate::py_install(c("numpy", "scipy"))
```
#  Example
This is a basic example which shows you how to solve a common problem:
```R
source("Data/Data processing function.R")

library("ULSL")
library(survival)
library(survminer)
library(parallel)


res <- load_tcga_data("AML",
                      data.dir = "D:/Data/tcga/AML",
                      clinical.dir = "D:/Data/tcga",
                      survival.dir = "D:/Data/tcga")
res[[1]]
omics.list = log.and.normalize(res$omics.data,normalize = c(TRUE,T,T),
                               log.transform = c(T,FALSE,T))
X <- omics.list
set.seed(10)

measures <- ULSL(X,alpha0 =2 , beta =0.1 , d =40 , maxIters =200 ,gt = NULL, c = NULL)

cl <- measures$cl

###**Survival Analysis

exp_common <- res$omics.data[[1]]
survival_common <- res$survival
survival_common$PatientID <- gsub("-", ".", survival_common$PatientID)
names(cl) <- colnames(exp_common)
df_cl <- data.frame(PatientID = names(cl), cluster = as.factor(cl))
survival2 <- merge(survival_common, df_cl, by = "PatientID")
surv_obj <- Surv(time = survival2$Survival, event = survival2$Death)


cluster_counts <- table(survival2$cluster)
cluster_labels <- paste0("Cluster ", names(cluster_counts), " (n=", cluster_counts, ")")


fit <- survfit(surv_obj ~ cluster, data = survival2)

survresult  <- survdiff(Surv(Survival, Death) ~ cluster, data = survival2)
p.val <- 1 - pchisq(survresult$chisq, df = length(survresult$n) - 1)


p <- ggsurvplot(
  fit,
  data = survival2,
  pval = TRUE,
  risk.table = FALSE,
  legend.title = "",
  legend.labs = cluster_labels,
  ggtheme = theme_bw(base_size = 14),
  legend = "top"
)


p$plot <- p$plot +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(color = "black", fill = "white")
  ) +
  guides(color = guide_legend(title = NULL)) +
  labs(x = "Time (days)") +
  ggtitle(cancer)


###Clinical Enrichment Analysis

clustering <- cl
clinical.params <-  clinical_common


clinical.metadata = list(gender='DISCRETE', age_at_initial_pathologic_diagnosis='NUMERIC',
                         pathologic_M='DISCRETE', pathologic_N='DISCRETE', pathologic_T='DISCRETE', pathologic_stage='DISCRETE')
pvalues = c()

params.being.tested = c()

for (clinical.param in names(clinical.metadata)) {

  if (!(clinical.param %in% colnames(clinical.params))) {

  next
  }

  clinical.values = clinical.params[names(clustering),clinical.param]
  is.discrete.param = clinical.metadata[clinical.param] == 'DISCRETE'
  is.numeric.param = clinical.metadata[clinical.param] == 'NUMERIC'
  stopifnot(is.discrete.param | is.numeric.param)



  if (is.numeric.param) {
    numeric.entries = !is.na(as.numeric(clinical.values))
    if (2 * sum(numeric.entries) < length(clinical.values)) {

   next
    }
  } else {
    not.na.entries = !is.na(clinical.values)
    should.skip = F
    if (2 * sum(not.na.entries) < length(clinical.values)) {
      should.skip = T
    } else if (length(table(clinical.values[not.na.entries])) == 1) {
      should.skip = T
    }
    if (should.skip) {

   next
    }
  }

  params.being.tested = c(params.being.tested, clinical.param)

  if (is.discrete.param) {

  pvalue = get.empirical.clinical(clustering[!is.na(clinical.values)], clinical.values[!is.na(clinical.values)], T)

  } else if (is.numeric.param) {

   pvalue = get.empirical.clinical(clustering[numeric.entries], as.numeric(clinical.values[numeric.entries]), F)
  }

  pvalues = c(pvalues, pvalue)

}
names(pvalues) = params.being.tested
pvalues
```
