load_tcga_data <- function(cancer,
                           data.dir,
                           clinical.dir,
                           survival.dir,
                           only.primary = FALSE) {
  # ---- cancer-specific settings ----
  if (cancer=="BIC") {
    cancer2 <- "breast"; nn <- 12
  } else if (cancer=="COAD") {
    cancer2 <- "colon"; nn <- 15
  } else if (cancer=="AML") {
    cancer2 <- "aml"; nn <- 12
  } else if (cancer=="GBM") {
    cancer2 <- "gbm"; nn <- 12
  } else if (cancer=="KIRC") {
    cancer2 <- "kidney"; nn <- 15
  } else if (cancer=="LIHC") {
    cancer2 <- "liver"; nn <- 15; only.primary <- TRUE
  } else if (cancer=="OV") {
    cancer2 <- "ovarian"; nn <- 15
  } else if (cancer=="SARC") {
    cancer2 <- "sarcoma"; nn <- 15
  } else if (cancer=="SKCM") {
    cancer2 <- "melanoma"; nn <- 15
  } else if (cancer=="LUSC") {
    cancer2 <- "lung"; nn <- 12
  } else {
    stop("Unsupported cancer type.")
  }

  # ---- load data ----
  clinical  <- read.table(file.path(clinical.dir, "clinical",cancer2), sep = "\t", header = TRUE)
  exp       <- read.table(file.path(data.dir, "exp"), sep = " ", header = TRUE)
  methy     <- read.table(file.path(data.dir, "methy"), sep = " ", header = TRUE)
  mirna     <- read.table(file.path(data.dir, "mirna"), sep = " ", header = TRUE)
  survival  <- read.table(file.path(survival.dir, cancer,"survival"), sep = "\t", header = TRUE)

  # ---- filter tumor samples ----
  exp   <- filter.non.tumor.samples(exp, only.primary = only.primary)
  methy <- filter.non.tumor.samples(methy, only.primary = only.primary)
  mirna <- filter.non.tumor.samples(mirna, only.primary = only.primary)

  # ---- special fix for KIRC ----
  if (cancer=="KIRC") {
    colnames(clinical)[1] <- "sampleID"
  }

  # ---- harmonize IDs ----
  clinical$sampleID <- toupper(substr(clinical$sampleID, 1, nn))
  colnames(exp)     <- toupper(substr(colnames(exp), 1, nn))
  colnames(methy)   <- toupper(substr(colnames(methy), 1, nn))
  colnames(mirna)   <- toupper(substr(colnames(mirna), 1, nn))
  survival$PatientID <- toupper(substr(survival$PatientID, 1, nn))

  survival <- na.omit(survival)

  # ---- replace "-" with "." ----
  clinical_ids <- gsub("-", ".", clinical$sampleID)
  survival_ids <- gsub("-", ".", survival$PatientID)
  exp_ids      <- colnames(exp)
  methy_ids    <- colnames(methy)
  mirna_ids    <- colnames(mirna)

  # ---- find common IDs ----
  common_ids <- Reduce(intersect, list(clinical_ids, survival_ids, exp_ids, methy_ids, mirna_ids))

  # ---- subset matched samples ----
  clinical_common <- clinical[clinical_ids %in% common_ids, ]
  clinical_common <- clinical_common[!duplicated(clinical_common$sampleID), ]
  clinical_common$sampleID <- gsub("-", ".", clinical_common$sampleID)
  rownames(clinical_common) <- clinical_common$sampleID

  survival_common <- survival[survival_ids %in% common_ids, ]
  survival_common <- survival_common[!duplicated(survival_common$PatientID), ]

  exp_common   <- as.matrix(exp[, common_ids])
  methy_common <- as.matrix(methy[, common_ids])
  mirna_common <- as.matrix(mirna[, common_ids])

  # ---- return ----
  return(list(
    omics.data = list(exp_common, methy_common, mirna_common),
    clinical   = clinical_common,
    survival   = survival_common,
    common.ids = common_ids
  ))
}



plot_survival_by_cluster <- function(exp_common,
                                     survival_common,
                                     cl,
                                     cancer = "Cancer") {


  library(survival)
  library(survminer)


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

  return(list(plot = p, p.value = p.val, survfit = fit, merged.data = survival2))
}


filter.non.tumor.samples <- function(raw.datum, only.primary=TRUE) {
  # 01 is primary, 06 is metastatic, 03 is blood derived cancer
  if (!only.primary)
    return(raw.datum[,substring(colnames(raw.datum), 14, 15) %in% c('01', '03', '06')])
  else
    return(raw.datum[,substring(colnames(raw.datum), 14, 15) %in% c('01')])
}


normalize.matrix <- function(data.matrix) {
  temp = data.matrix - rowMeans(data.matrix)
  should.keep = (apply(temp, 1, sd) != 0)
  return ((temp / apply(temp, 1, sd))[should.keep, ])
}
