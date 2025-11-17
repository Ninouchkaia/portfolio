### To do different gene expression by limma package

DGEs <- function(data, Phenotype, lFC, FDR) {
  require(limma)
  require(dplyr)
  # require(dendextend)


  message("[===========================]")
  message("[<<<<<<< DGEs START >>>>>>>>>]")
  message("[<<<< Pairwise analysis >>>>>]")
  message("------------------------------")

  # This function creates the pairs for the pairwise matrices
  design.pairs <- function(levels) {
    n <- length(levels)
    design <- matrix(0, n, choose(n, 2))
    rownames(design) <- levels
    colnames(design) <- 1:choose(n, 2)
    k <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        k <- k + 1
        design[i, k] <- 1
        design[j, k] <- -1
        colnames(design)[k] <- paste(levels[i], "-", levels[j], sep = "")
      }
    }
    design
  }

  # This function creates the pairs for the pairwise matrices

  design <- model.matrix(~ 0 + Phenotype)
  contr.matrix <- design.pairs(levels(factor(Phenotype)))
  colnames(design) <- rownames(contr.matrix)

  # Removing heteroscedascity from data
  v <- voom(log2(data + 1), design, plot = F)

  # Fitting linear models for comparisons of interest
  Fit <- lmFit(v, design) %>%
    contrasts.fit(., contr.matrix) %>%
    eBayes(.) %>%
    treat(., lfc = lFC)

  # this code includes list of DEGs for pairwise comparisons
  FitList <- list()
  for (i in 1:ncol(contr.matrix)) {
    FitList[[i]] <- topTreat( Fit[ which( decideTests( Fit)[, i] != 0), ], coef = i, adjust.method = "BH", number = nrow(v$E) ) %>%
      mutate(ID = rownames(.)) %>%
      filter(., adj.P.Val < FDR)
  }

  names(FitList) <- colnames(contr.matrix)
  message("thank you for waiting")
  return(FitList)
}

### To do signature ###

Signfeature <- function(data_TPM, Phenotype, FitList, FCcuoff, MaxDMRs) {
  require(dplyr)
  require(matrixStats)
  require(limma)

  message("[===========================]")
  message("[<<<< Signature START >>>>>]")
  message("-----------------------------")

  # design <- model.matrix(~0 + Phenotype)
  # v_expr <- voom(data, design, plot = FALSE)

  # using positive FC (=over-expressed in cell type of interest)
  # FitList is DGEs which calculate by DGEs function

  #this code explained how to choose the genes for final signature gens
  A1 <- list()
  for (i in 1:length(FitList)) {
    A1[[i]] <- filter(FitList[[i]], abs(logFC) > FCcuoff) %>% arrange(., desc(logFC)) # %>% top_n(ceiling(n()*0.5), wt = logFC)
    if (nrow(A1[[i]]) > MaxDMRs) {
      A1[[i]] <- A1[[i]][1:MaxDMRs, ]
    }
  }


  sig <- lapply(A1, function(x) dplyr::select(x, ID))
  sig <- do.call(rbind, sig)
  sig <- filter(sig, !duplicated(ID))
  data1 <- data_TPM[rownames(data_TPM) %in% sig$ID, ]

  # Print number of selected probes (signature)
  nrow(data1)

  result <- getMedVal(data1, Phenotype)
  message("[===========================]")
  message("[<<<< Signature END >>>>>]")
  message("-----------------------------")
  return(result)
}

### Function to get median value of each gene from an input dataframe ###

getMedVal <- function(data, Phenotype) {
  library(matrixStats)

  Trans <- data.frame(t(data))
  Mt.Split <- split(Trans, Phenotype)
  Mt.Split <- lapply(Mt.Split, function(x) colMedians(data.matrix(x)))
  Mt.Split <- do.call(cbind, Mt.Split)
  rownames(Mt.Split) <- rownames(data)
  return(Mt.Split)
}

### To choose suitable threshold of LogFC (save as PDF file)

DGEs.QC <- function(FitList, FDR) {
  require(ggpubr)

  p <- list()

  for (i in 1:length(FitList)) {
    p[[i]] <- ggplot(FitList[[i]][FitList[[i]]$adj.P.Val < FDR, ], aes(x = logFC)) +
      geom_density(colour = "red") +
      labs(title = names(FitList)[i], x = "logFC", y = "Density") +
      theme(plot.title = element_text(hjust = 0.5))
  }
  return(do.call(ggarrange, p))
}


getExpMatrix <- function(geneList, expression) {
  df_out <- expression[geneList, ]
  return(df_out)
}

fpkm2tpm <- function(fpkm) {
  tpm <- exp(log(fpkm) - log(sum(fpkm, na.rm = T)) + log(1e6))
  tpm[is.na(tpm)] <- 0
  return(tpm)
}

TPM_normalization <- function(data) {

  # TPM normalization

  TPM_data <- t(t(data) * 1e6 / apply(data, 2, sum))
  return(TPM_data)
}
