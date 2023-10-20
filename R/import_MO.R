#' Import, format and option to preprocess omics data for integration
#'
#' @param rnaseq_counts Matrix of RNAseq counts
#' @param metab_peaks Matrix of metabolite peaks
#' @param rrbs_mvals Matrix of RRBS m-values
#' @param mirna_counts Matrix of microRNA counts
#' @param pirna_counts Matrix of piRNA counts
#' @param meta Metadata, must have sample ids in first column matching those in omics data and treatment col called TRT
#' @param norm Normalization, True or False
#' @param batch Omics layers with batches corresponding to "batch_omic" columns in metadata, column labels as follows: batch_rna, batch_metab, batch_rrbs
#'
#' @return formatted data_list for integration
#' @export
#'
import_MO <- function(rnaseq_counts = NULL,
                        metab_peaks = NULL,
                        rrbs_mvals = NULL,
                        mirna_counts = NULL,
                        pirna_counts = NULL,
                        meta = NULL,
                        norm = TRUE,
                        batch = NULL) {

    # stop if metadata not provided
    if(is.null(meta)) stop("No metadata provided")

    #stop if data are not provided
    data_list <- list(rnaseq_counts, metab_peaks, rrbs_mvals, mirna_counts, pirna_counts)
    names(data_list) <- c("rnaseq_counts", "metab_peaks", "rrbs_mvals", "mirna_counts", "pirna_counts")
    if(is.null(rnaseq_counts) & is.null(metab_peaks) & is.null(rrbs_mvals) & is.null(mirna_counts)) {
      stop("No omics layers provided")
    }

    #check sample names match rows of metadata
    for(j in names(data_list)) {
      if(!is.null(data_list[[j]])) {
        a <- dim(data_list[[j]]) == dim(meta)[1]
        if (a[1] == TRUE){

          #rownames of data should match metadata
          for(i in 1:length(rownames(data_list[[j]]))) {
            if (rownames(data_list[[j]])[i] != rownames(meta)[i])stop(j, " sample ", i, " name does not match metadata")
          }
          data_list[[j]] <- as.data.frame(t(data_list[[j]]))
        }
        if (a[2] == TRUE){

          #colnames of data should match metadata
          for(i in 1:length(colnames(data_list[[j]]))) {
            if (colnames(data_list[[j]])[i] != rownames(meta)[i])stop(j, " sample ", i, " name does not match metadata")
          }}
      }}

    if(norm == TRUE){
      cdir <- paste("MOnorm", Sys.time(), sep = "_")
      cdir <- gsub(":", ".", cdir)
      dir.create(cdir)

    ## filter, normalize, scale
    #rnaseq
    if(!is.null(data_list$rnaseq_counts)) {
      print("normalizing RNAseq counts with edgeR TMM")
      grDevices::pdf(file = paste0(cdir, "/", "rnaseq_norm_filter.pdf"))
      graphics::par(mfrow = c(1, 2))
      # sub-setting %20 of rows for plots
      select_rows <- sample(nrow(data_list$rnaseq_counts), length(rownames(data_list$rnaseq_counts)) * .2)
      suppressWarnings(graphics::boxplot(log2(data_list$rnaseq_counts[select_rows, ] + 1),
                               main = paste("log2 raw counts\n n=", length(data_list$rnaseq_counts[, 1])), las = 2))
      y <- edgeR::DGEList(counts = data_list$rnaseq_counts, group = meta$TRT)
      keep <- edgeR::filterByExpr(y)
      y <- y[keep,,keep.lib.sizes = FALSE]
      y <- edgeR::calcNormFactors(y)
      print("RNAseq viz")
      data_list$rnaseq_counts <- edgeR::cpm(y, log = TRUE, normalized.lib.sizes = TRUE) #trying log transformed data
      if("rnaseq_counts" %in% batch){
        print("RNAseq batch correction with ComBat")
        data_list$rnaseq_counts <- sva::ComBat(data_list$rnaseq_counts, meta$batch_rna, mod = as.factor(meta$TRT))

      }
      #sub-setting %20 of rows for plots
      select_rows <- sample(nrow(data_list$rnaseq_counts), length(rownames(data_list$rnaseq_counts)) * .2)
      suppressWarnings(graphics::boxplot(data_list$rnaseq_counts[select_rows, ],
                               main = paste("log2 filter & TMM norm counts\n n=",
                               length(data_list$rnaseq_counts[, 1])), las = 2))
      grDevices::dev.off()
      #Normalization independent visualizations
      #PCA
      grDevices::pdf(file = paste0(cdir, "/", "rnaseq_norm_PCA.pdf"))
      TRTColor <- as.numeric(y$samples$group) + .01
      #TRT can only be 5 here
      color_options<-c("red3", "coral","green3", "purple", "black")
      for(i in 1:length(unique(TRTColor))) {
        TRTColor <- gsub(unique(TRTColor)[i], color_options[i], TRTColor)
      }
      pr <- stats::prcomp(t(data_list$rnaseq_counts), center = TRUE)
      proportion_variance <- pr$sdev^2 / sum(pr$sdev^2) * 100
      par(mar=c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
      plot(pr$x[, 1], pr$x[, 2], col = TRTColor, main = 'RNAseq norm principal components 1 & 2', pch = 16, cex = 2,
           xlab = paste0("Principal Component 1 %", round(proportion_variance[1], 2)),
           ylab = paste0("Principal Component 2 %", round(proportion_variance[2], 2)))
      legend('topright', inset = c(-0.3, 0),legend = unique(y$samples$group), fill = unique(TRTColor),
             bg = "transparent", bty = "n", title = "Treatment")
      plot(pr$x[, 3], pr$x[, 4], col = TRTColor, main = 'RNAseq norm principal components 3 & 4', pch = 16, cex = 2,
          xlab = paste0("Principal Component 3 %", round(proportion_variance[3], 2)),
          ylab = paste0("Principal Component 4 %", round(proportion_variance[4], 2)))
      legend('topright', inset = c(-0.3, 0),legend = unique(y$samples$group), fill = unique(TRTColor),
             bg = "transparent", bty = "n", title = "Treatment")
      par(mar = c(5.1, 4.1, 4.1, 2.1))
      #scree plot
      bp <- barplot(proportion_variance[1:20], ylab = "Proportion variance", xlab = "PCs",
                    names.arg = c(1:20), col = rep("black", 20), main = "RNAseq scree plot")
      text(bp, proportion_variance[1:20] - 0.5, labels = round(proportion_variance[1:20],
                                                            digits = 2), col = "white", cex = 0.5)
      grDevices::dev.off()
      #MDS
      grDevices::pdf(file = paste0(cdir, "/", "rnaseq_norm_PCoA_MDS.pdf"))
      par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
      limma::plotMDS(y, col = TRTColor, pch = 16, gene.selection = "pairwise", top = 500)
      legend("topright", inset = c(-0.25, 0), legend = unique(y$samples$group), fill = unique(TRTColor),
             bg = "transparent", bty = "n", title = "Treatment", cex = 0.5)
      grDevices::dev.off()
    }

    #metab
    if(!is.null(data_list$metab_peaks)) {
      print("normalizing metabolite sum peak area with mean centering and pareto scaling")
      data_list$metab_peaks <- data_list$metab_peaks[rowSums(data_list$metab_peaks) > 0, ]
      grDevices::pdf(file = paste0(cdir, "/" , "metab_impute_norm_filter.pdf"))
      graphics::par(mfrow = c(1, 2))
      graphics::boxplot(data_list$metab_peaks, main = "Raw values", las = 2)
      #add zero impute step
      data_list$metab_peaks <- sweep(data_list$metab_peaks * 1000, 1, rowSums(data_list$metab_peaks), FUN = "/")
      data_list$metab_peaks <- as.data.frame(IMIFA::pareto_scale(data_list$metab_peaks), centering = FALSE)
      if("metab_peaks" %in% batch){
        print("Metabolite batch correction with ComBat")
        data_list$rnaseq_counts <- sva::ComBat(data_list$metab_peaks, meta$batch_metab, mod = as.factor(meta$TRT))

      }
      graphics::boxplot(data_list$metab_peaks, main = "Sum peak norm\n & Pareto scale counts", las = 2)
      grDevices::dev.off()
    }

    #rrbs
    if(!is.null(data_list$rrbs_mvals)) {
      if("rrbs_mvals" %in% batch){
        print("RRBS batch correction with ComBat")
        data_list$rrbs_mvals <- sva::ComBat(data_list$rnaseq_counts, meta$batch_rrbs, mod = as.factor(meta$TRT))

      }
      print("RRBS viz")
      grDevices::pdf(file = paste0(cdir, "/", "rrbs_boxplot.pdf"))
      #sub-setting %20 of rows for plots
      select_rows <- sample(nrow(data_list$rrbs_mvals), length(rownames(data_list$rrbs_mvals)) * .2)
      suppressWarnings(graphics::boxplot(data_list$rrbs_mvals[select_rows, ],
                               main = paste("rrbs features=",
                               length(data_list$rrbs_mvals[, 1])), las = 2))
      grDevices::dev.off()
      #MDS
      grDevices::pdf(file = paste0(cdir, "/", "rrbs_PCoA_MDS.pdf"))
      par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
      limma::plotMDS(data_list$rrbs_mvals, col = TRTColor, pch = 16, gene.selection = "pairwise", top = 500)
      legend('topright', inset = c(-0.25, 0), legend = unique(y$samples$group), fill = unique(TRTColor),
             bg = "transparent", bty = "n", title = "Treatment", cex = .5)
      grDevices::dev.off()
    }
    #mirna
    #pirna
    assign("data_list", data_list, envir = .GlobalEnv)
    assign("meta", meta, envir = .GlobalEnv)
    }
    if(norm == FALSE){
      assign("data_list", data_list, envir = .GlobalEnv)
      assign("meta", meta, envir = .GlobalEnv)
    }

  }
