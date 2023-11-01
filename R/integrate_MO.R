#' Multi-omics integration for toxicology data
#'
#' @param int_method sPLS-DA, MOFA2, WGCNA, SNF, or iPCA
#' @param RRBS_feature_map dataframe of RRBS features to genes, must match format of Unique_Features_to_Genes.csv from formatRRBS output
#'
#' @return output
#' @importFrom WGCNA cor
#' @export
integrate_MO <- function(int_method = c("sPLS-DA", "MOFA", "WGCNA", "SNF", "iPCA"), RRBS_feature_map = NULL){

  #stop if metadata not provided
  if (is.null(meta)) stop("No metadata provided, run import function")

  #stop if data are not provided
  if (is.null(data_list)) stop("No data provided, run import function")

  cdir <- paste("integrateMO", int_method, Sys.time(), sep = "_")
  cdir <- gsub(":", ".", cdir)
  dir.create(cdir)

  # Integrate

  int_method <- match.arg(int_method)
  cat("Using",int_method,"in as integration method\n")
  X <- list()
  for (i in names(data_list[lengths(data_list) != 0])){
    X[[i]] <- as.data.frame(t(data_list[[i]]))
  }

  ## Prepare Unique_Features_to_Genes for RRBS

  if (!is.null(RRBS_feature_map)){
  UF2G <- utils::read.csv(RRBS_feature_map)
  UF2G$id <- paste(UF2G$chrom, UF2G$loc, sep = "-")
  }

  ## mixOmics

  ### limit to top variable features by MAD
  if (int_method == "sPLS-DA"){
    for (i in names(X)){
      if (!is.null(X[[i]])){
        if (ncol(X[[i]]) > 10000){
          X[[i]] <- X[[i]][, order(apply(X[[i]], 2, mad), decreasing = TRUE)[1:10000]]
          cat("Limiting ", i,"to top 10,000 most variable features using median absolute deviation\n")
        }else{
          cat("Using all", ncol(X[[i]]), "rnaseq_counts features\n")
        }
      }
    }
    Y <- as.factor(meta$TRT)
    # doing "data driven" approach to setting design matrix
    test_combo <- gtools::combinations(n = length(X), r = 2, v = 1:length(X))
    test_vec <- vector()
    for (i in nrow(test_combo)){
      res1 <- mixOmics::pls(X[[test_combo[i, 1]]], X[[test_combo[i ,2]]], ncomp = 1)
      cor1 <- stats::cor(res1$variates$X, res1$variates$Y)
      test_vec <- append(test_vec, cor1)
    }
    test_vec <- mean(test_vec)
    design <- matrix(test_vec, ncol = length(X), nrow = length(X),
                     dimnames = list(names(X), names(X)))
    diag(design) <- 0
    diablo.mo <- mixOmics::block.plsda(X, Y, ncomp = length(unique(meta$TRT)) + 1, design = design)
    perf.diablo.mo <- mixOmics::perf(diablo.mo, validation = "Mfold", folds = min(table(meta$TRT)), nrepeat = 10)
    ncomp <- min(perf.diablo.mo$choice.ncomp$WeightedVote[2,])
    # dist_name <- colnames(perf.diablo.mo$choice.ncomp$WeightedVote[,perf.diablo.mo$choice.ncomp$WeightedVote[2,]==ncomp])[1]
    dist_name <- colnames(perf.diablo.mo$error.rate[[1]])[which(perf.diablo.mo$error.rate[[1]] == min(perf.diablo.mo$error.rate[[1]]),
                                                                  arr.ind = TRUE)[[2]]]
    svglite::svglite(file = paste0(cdir, "/", "mixOmics_classification_error_rate.svg"))
    plot(perf.diablo.mo)
    graphics::mtext(paste0("comp = ", ncomp, " dist = ", dist_name), side = 3)
    grDevices::dev.off()
    test.keepX <- list()
    for (i in names(X)){
      test.keepX[[i]] <- c(5:9, seq(10, 25, 5))
    }
    for (i in dist_name){
      tune.diablo.mo <- mixOmics::tune.block.splsda(X, Y, ncomp = ncomp,
                                            test.keepX = test.keepX, design = design,
                                            validation = "Mfold", folds = 5, nrepeat = 10,
                                            # use two CPUs for faster computation
                                            BPPARAM = BiocParallel::SnowParam(workers = parallel::detectCores()),
                                            dist = i) # should this error rate match multiomics_pdf?
    }
    list.keepX <- tune.diablo.mo$choice.keepX
    diablo.mo <- mixOmics::block.splsda(X, Y, ncomp = ncomp,
                                          keepX = list.keepX)
    for (i in 1:ncomp){
      values<-vector()
      for (j in 1:length(X)){
        values<-rbind(values, as.data.frame(mixOmics::selectVar(diablo.mo, comp = i)[[j]][2]))
      }
      if (exists("UF2G")){
        values[, 3] <- UF2G[match(values[, 1], UF2Gg$id), 4]
      }
      utils::write.csv(values, paste0(cdir, "/", "splsda_variables_comp", i, ".csv"))
    }
    if (ncomp > 1){
      svglite::svglite(file = paste0(cdir, "/", "mixOmics_plotVar1-2.svg"))
      mixOmics::plotVar(diablo.mo, var.names = c(TRUE), style = "graphics", legend = TRUE,
                        title = "DIABLO components 1 - 2")
      grDevices::dev.off()
      if (ncomp == 4){
        svglite::svglite(file = paste0(cdir, "/", "mixOmics_plotVar3-4.svg"))
        mixOmics::plotVar(diablo.mo, var.names = c(TRUE), style = "graphics", legend = TRUE,
                          comp = c(3, 4), title = "DIABLO components 3 - 4")
        grDevices::dev.off()
      }
      svglite::svglite(file = paste0(cdir, "/", "mixOmics_plotIndiv1-2.svg"))
      mixOmics::plotIndiv(diablo.mo, ind.names = FALSE, legend = TRUE, comp = c(1, 2),
                          title = "DIABLO components 1 - 2", block = "weighted.average")
      grDevices::dev.off()
      if (ncomp == 4){
        svglite::svglite(file = paste0(cdir, "/", "mixOmics_plotIndiv3-4.svg"))
        mixOmics::plotIndiv(diablo.mo, ind.names = FALSE, legend = TRUE, comp = c(3, 4),
                            title = "DIABLO components 3 - 4", block = "weighted.average")
        grDevices::dev.off()
      }
      for (i in 1:ncomp){
        svglite::svglite(file = paste0(cdir, "/", "mixOmics_cim", i, ".svg"))
        mixOmics::cimDiablo(diablo.mo, comp = i, margin=c(11, 15), legend.position = "topright",
                            trim = FALSE, size.legend = 0.7)
        grDevices::dev.off()
      }
      for (i in names(X)){
        svglite::svglite(file = paste0(cdir, "/", "mixOmics_auroc_", i, ".svg"))
        auc.diablo.mo <- mixOmics::auroc(diablo.mo, roc.block = i, roc.comp = ncomp,
                                         print = FALSE)
        grDevices::dev.off()
      }
    }
  }

  # WGCNA - counts need to be vst or log2
  if (int_method == "WGCNA"){
    if (!is.null(X$rrbs_mvals)){
      X$rrbs_mvals <- X$rrbs_mvals[, order(apply(X$rrbs_mvals, 2, mad), decreasing = TRUE)[1:(.25 * ncol(X$rrbs_mvals))]]
      cat("Limiting rrbs_mvals to top", ncol(X$rrbs_mvals), "most variable using median absolute deviation\n")
      }
    lowcount_omics_MEs <- list()
    block_MEs <- list()
    TRT_number <- meta$TRT
    for (i in 1:length(unique(meta$TRT))){
      TRT_number <- gsub(unique(meta$TRT)[i], i, TRT_number)
      }
    TRT_number <- as.numeric(TRT_number)
    traitColors <- WGCNA::numbers2colors(TRT_number, signed = FALSE)
    # rnaseq counts are already logged
    for (i in names(X)){
      svglite::svglite(file = paste0(cdir, "/", "hclust_", i, "_sampleTree.svg"))
      sampleTree <- fastcluster::hclust(stats::dist(X[[i]]), method = "average")
      #cluster with metadata
      WGCNA::plotDendroAndColors(fastcluster::hclust(stats::dist(X[[i]]), method = "single"), traitColors,
                                 groupLabels = colnames(meta)[2],
                                 main = paste("Single clustering dendro", i))
      WGCNA::plotDendroAndColors(fastcluster::hclust(stats::dist(X[[i]]), method = "average"), traitColors,
                                 groupLabels = colnames(meta)[2],
                                 main = paste("Average clustering dendro", i))
      WGCNA::plotDendroAndColors(fastcluster::hclust(stats::dist(X[[i]]), method = "complete"), traitColors,
                                 groupLabels = colnames(meta)[2],
                                 main = paste("Complete clustering dendro", i))
      WGCNA::plotDendroAndColors(fastcluster::hclust(stats::dist(X[[i]]), method = "median"), traitColors,
                                 groupLabels = colnames(meta)[2],
                                 main = paste("Median clustering dendro", i))
      WGCNA::plotDendroAndColors(fastcluster::hclust(stats::dist(X[[i]]), method = "mcquitty"), traitColors,
                                 groupLabels = colnames(meta)[2],
                                 main = paste("Mcquitty clustering dendro", i))
      WGCNA::plotDendroAndColors(fastcluster::hclust(stats::dist(X[[i]]), method = "ward.D"), traitColors,
                                 groupLabels = colnames(meta)[2],
                                 main = paste("ward.D clustering dendro", i))
      grDevices::dev.off()
      # Module detection one step - for low feature count omic layers but if very low dont need to cluster
      if (length(X[[i]]) > 100 && length(X[[i]]) < 5000){
        powers = c(c(1:10), seq(from = 12, to = 40, by = 2))
        #allowWGCNAThreads()
        sft <- WGCNA::pickSoftThreshold(X[[i]], powerVector = powers, verbose = 0)
        power_from_sft <- sft$fitIndices[sft$fitIndices$SFT.R.sq == max(sft$fitIndices$SFT.R.sq[1:10]), 1]
        color_forplot <- rep("black", length(sft$fitIndices$Power))
        color_forplot[power_from_sft] <- "red"
        svglite::svglite(file = paste(cdir, "/", "soft threasholding", i))
        graphics::par(mfrow = c(1, 2))
        plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
             xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n",
             main = paste("Scale independence"))
        graphics::text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
             labels = powers, cex = 0.9, col = color_forplot)
        plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
             xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
             main = paste("Mean connectivity"))
        graphics::text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = 0.9, col = color_forplot)
        grDevices::dev.off()
        adjacency <- WGCNA::adjacency(X[[i]], power = power_from_sft)
        TOM <- WGCNA::TOMsimilarity(adjacency)
        dissTOM <- 1 - TOM
        # Call the hierarchical clustering function
        geneTree <- fastcluster::hclust(stats::as.dist(dissTOM), method = "average")
        # Module size too large will result in one module, make it dependent on features
        minModuleSize <- floor(length(X[[i]]) / 20)
        # Module identification using dynamic tree cut:
        dynamicMods <- dynamicTreeCut::cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                            deepSplit = 2, pamRespectsDendro = FALSE,
                                            minClusterSize = minModuleSize)
        dynamicColors <- WGCNA::labels2colors(dynamicMods)
        # Calculate eigengenes
        MEList <- WGCNA::moduleEigengenes(X[[i]], colors = dynamicColors)
        MEs <- MEList$eigengenes
        # Calculate dissimilarity of module eigengenes
        MEDiss <- 1 - WGCNA::cor(MEs)
        # Cluster module eigengenes
        METree <- fastcluster::hclust(stats::as.dist(MEDiss), method = "average")
        # Plot the result
        svglite::svglite(file = paste0(cdir, "/", "module_eigengene_clustering_cutoff_"), i)
        plot(METree, main = paste("Clustering of", i, "module eigengenes"),
             xlab = "", sub = "")
        MEDissThres <- 0.25
        # Plot the cut line into the dendrogram
        graphics::abline(h = MEDissThres, col = "red")
        grDevices::dev.off()
        # Call an automatic merging function
        merge <- WGCNA::mergeCloseModules(X[[i]], dynamicColors, cutHeight = MEDissThres, verbose = 0)
        # The merged module colors
        mergedColors <- merge$colors
        # Eigengenes of the new merged modules:
        mergedMEs <- merge$newMEs
        svglite::svglite(file = paste0(cdir, "/", "dendrogram_", i))
        WGCNA::plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                                    c("Dynamic Tree Cut", "Merged dynamic"),
                                    dendroLabels = FALSE, hang = 0.03,
                                    addGuide = TRUE, guideHang = 0.05)
        grDevices::dev.off()
        lowcount_omics_MEs[[i]] <- WGCNA::moduleEigengenes(X[[i]], mergedColors)$eigengenes
        utils::write.csv(lowcount_omics_MEs[[i]], paste0(cdir, "/", "Module_Eigengens_", i, ".csv"))
      }else{
        # Module detection one step - for high feature count omic layers
        # Choose a set of soft-thresholding powers
        powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
        # Call the network topology analysis function
        sft <- WGCNA::pickSoftThreshold(X[[i]], powerVector = powers, verbose = 0)
        power_from_sft <- sft$fitIndices[sft$fitIndices$SFT.R.sq == max(sft$fitIndices$SFT.R.sq[1:10]), 1]
        color_forplot <- rep("black", length(sft$fitIndices$Power))
        color_forplot[power_from_sft] <- "red"
        svglite::svglite(file = paste0(cdir, "/", "soft_threasholding_", i, ".svg"))
        graphics::par(mfrow = c(1, 2))
        # Plot the results:
        # Scale-free topology fit index as a function of the soft-thresholding power
        plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
             xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
             main = paste("Scale independence"))
        graphics::text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
             labels = powers, cex = .9, col = color_forplot)
        # Mean connectivity as a function of the soft-thresholding power
        plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
             xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
             main = paste("Mean connectivity"))
        graphics::text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = 0.9, col = color_forplot)
        grDevices::dev.off()
        bwnet <- WGCNA::blockwiseModules(X[[i]], maxBlockSize = 5000,
                                         power = power_from_sft, TOMType = "unsigned",
                                         minModuleSize = 30,
                                         reassignThreshold = 0, mergeCutHeight = 0.25,
                                         numericLabels = TRUE, saveTOMs = FALSE,
                                         saveTOMFileBase = paste0(i, "_TOM"),
                                         verbose = 0)

        bwLabels <- bwnet$colors
        bwModuleColors <- WGCNA::labels2colors(bwLabels)
        svglite::svglite(file = paste0(cdir, "/", i, "_dendro.svg"))
        if (length(bwnet$dendrograms) > 1){
          for (j in 1:length(bwnet$dendrograms)){
            WGCNA::plotDendroAndColors(bwnet$dendrograms[[j]], bwModuleColors[bwnet$blockGenes[[j]]],
                                       "Module colors", main = paste(i,"Dendro & module colors in block", j),
                                       dendroLabels = FALSE, hang = 0.03,
                                       addGuide = TRUE, guideHang = 0.05)
          }
        }else{
          graphics::par(mfrow = c(1, 1))
          WGCNA::plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
                                     "Module colors", main = "Dendrogram and module colors in block 1",
                                     dendroLabels = FALSE, hang = 0.03,
                                     addGuide = TRUE, guideHang = 0.05)
        }
        grDevices::dev.off()
        block_MEs[[i]] <- WGCNA::moduleEigengenes(X[[i]], bwModuleColors)$eigengenes
        utils::write.csv(block_MEs[[i]], paste0(cdir, "/", "Module_Eigengens_", i, ".csv"))

      }
    }
    # This all needs to be modified for having +2 omic layers
    all_MEs <- c(block_MEs, lowcount_omics_MEs)
    MEs_combo <- gtools::combinations(n = length(all_MEs), r = 2, v = 1:length(X))
    for (v in 1:nrow(MEs_combo)){
      nGenes <- ncol(X[[MEs_combo[v, 1]]])
      nSamples <- nrow(X[[MEs_combo[v, 1]]])
      moduleTraitCor <- WGCNA::cor(WGCNA::orderMEs(all_MEs[[MEs_combo[v, 1]]]),
                                   WGCNA::orderMEs(all_MEs[[MEs_combo[v, 2]]]), method = "pearson")
      moduleTraitPvalue <- WGCNA::corPvalueStudent(moduleTraitCor, nSamples)
      textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")")
      dim(textMatrix) <- dim(moduleTraitCor)
      pdf_name <- paste(names(X)[MEs_combo[v, 1]], names(X)[MEs_combo[v, 2]], "module_relationship.svg", sep = "_")
      svglite::svglite(file = paste0(cdir, "/", pdf_name))
      # Display the correlation values within a heatmap plot
      graphics::par(mar = c(7, 7, 2, 2))
      WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                            xLabels = names(WGCNA::orderMEs(all_MEs[[MEs_combo[v, 2]]])),
                            yLabels = names(WGCNA::orderMEs(all_MEs[[MEs_combo[v, 1]]])),
                            xSymbols = names(WGCNA::orderMEs(all_MEs[[MEs_combo[v, 2]]])),
                            ySymbols = names(WGCNA::orderMEs(all_MEs[[MEs_combo[v, 1]]])),
                            colorLabels = FALSE,
                            colors = WGCNA::blueWhiteRed(50),
                            textMatrix = textMatrix,
                            setStdMargins = FALSE,
                            cex.text = 0.5,
                            zlim = c(-1, 1),
                            main = paste0("y =", names(X)[MEs_combo[v, 1]], "x =", names(X)[MEs_combo[v, 2]], "/nModule-module relationships"))
      grDevices::dev.off()

      # Write correlation and p-values of module to module relationships
      utils::write.csv(moduleTraitCor, paste0(cdir, "/", "Module-Module_cor.csv"))
      utils::write.csv(moduleTraitPvalue, paste0(cdir, "/", "Module-Module_pval.csv"))

      # Test extract module membership
      geneModuleMembership <- as.data.frame(WGCNA::cor(X[[MEs_combo[v, 1]]], WGCNA::orderMEs(all_MEs[[MEs_combo[v, 1]]]), method = "pearson"))
      if (names(X)[MEs_combo[v, 1]]=="rrbs_mvals"){
        if (exists("UF2G")){
          rownames(geneModuleMembership) <- paste(rownames(geneModuleMembership), UF2G[match(rownames(geneModuleMembership), UF2Gg$id), 4])
        }
      }
      utils::write.csv(geneModuleMembership, paste0(cdir, "/", "Module_membership_", names(X)[MEs_combo[v, 1]], ".csv"))
      geneModuleMembership2 <- as.data.frame(WGCNA::cor(X[[MEs_combo[v, 2]]], WGCNA::orderMEs(all_MEs[[MEs_combo[v, 2]]]), method = "pearson"))
      utils::write.csv(geneModuleMembership2, paste0(cdir, "/", "Module_membership_", names(X)[MEs_combo[v, 2]], ".csv"))
      # Testing module membership viz
      svglite::svglite(file = paste0(cdir, "/", paste(names(X)[MEs_combo[v, 1]], names(X)[MEs_combo[v, 2]], "module_membership_subset.svg", sep = "_")))
      graphics::par(cex.main = 1)
      df <- unique(c(unlist(lapply(geneModuleMembership, function(x) utils::tail(rownames(geneModuleMembership)[order(x, decreasing = TRUE)], n = 3))),
                     unlist(lapply(geneModuleMembership, function(x) utils::tail(rownames(geneModuleMembership)[order(x, decreasing = FALSE)], n = 3)))))
      stats::heatmap(as.matrix(geneModuleMembership[rownames(geneModuleMembership) %in% df, ]), main = paste("module membership of extreme", names(X)[MEs_combo[v, 1]]),
              margins = c(8, 8))
      df <- unique(c(unlist(lapply(geneModuleMembership2, function(x) utils::tail(rownames(geneModuleMembership2)[order(x, decreasing = TRUE)], n = 3))),
                     unlist(lapply(geneModuleMembership2, function(x) utils::tail(rownames(geneModuleMembership2)[order(x, decreasing = FALSE)], n = 3)))))
      stats::heatmap(as.matrix(geneModuleMembership2[rownames(geneModuleMembership2) %in% df, ]), main = paste("module membership of extreme", names(X)[MEs_combo[v, 2]]),
              margins = c(8, 8))
      grDevices::dev.off()
      }
  }

  # SNF https://github.com/cran/SNFtool
  if (int_method == "SNF"){
    # for (i in names(X)){
    #    X[[i]]<-as.data.frame(t(X[[i]]))
    #  }
    if (!is.null(X$rrbs_mvals)){
      X$rrbs_mvals <- X$rrbs_mvals[, order(apply(X$rrbs_mvals, 2, mad), decreasing = TRUE)[1:(.25 * ncol(X$rrbs_mvals))]]
      cat("Limiting rrbs_mvals to top", ncol(X$rrbs_mvals), "most variable using median absolute deviation\n")}
    if (length(meta$sample) < 20){
      K <- length(meta$sample) - 1
    }else{K <- 20}
    #number of nearest neighbors, must be >1, usually 10-20
    alpha <- 0.5 ##hyperparameter, usually (0.3~0.8)
    t <- 20 ##Number of Iterations, usually (10~50)
    clusterNum <- length(unique(meta$TRT))
    W_temp <- list()
    for (i in names(X)){
      distance <- SNFtool::dist2(as.matrix((X[[i]])), as.matrix((X[[i]])))
      W_temp[[i]] <- SNFtool::affinityMatrix(distance, K, alpha)
    }
    W <- SNFtool::SNF(W_temp, K = K, t = t)
    group <- SNFtool::spectralClustering(W, clusterNum)
    TRT_number <- meta$TRT
    for (i in 1:length(unique(meta$TRT))){
      j <- unique(group)[i]
      TRT_number <- gsub(unique(meta$TRT)[i], j, TRT_number)
    }
    svglite::svglite(file = paste0(cdir, "/" , "SNF_cluster_heatmap.svg"))
    SNFtool::displayClustersWithHeatmap(W, group, ColSideColors = cbind(as.character(TRT_number), as.character(group))) #hack this to be better https://rdrr.io/cran/SNFtool/src/R/displayClustersWithHeatmap.R
    graphics::legend("topleft", legend = unique(group), fill = unique(group), cex = 0.8)
    grDevices::dev.off()
    #diag(W)=0 #which of these is correct?
    #diag(W)=max(W) #is correct? https://rdrr.io/bioc/CancerSubtypes/src/R/ClusteringMethod.R
    distanceMatrix <- W
    attr(distanceMatrix, "class") <- "Similarity"
    result <- list(group = group, distanceMatrix = distanceMatrix, originalResult = group)
    #rownames(meta)==rownames(result$distanceMatrix)
    truelabel <- TRT_number
    SNFNMI <- SNFtool::calNMI(group, truelabel)
    net <- W * 100
    diag(net) <- 0
    net[net < min(apply(net, 2, max))] <- 0 #whatever number the dark red color is in heatplot
    net <- network::network(net, directed = FALSE, ignore.eval = FALSE, names.eval = "weights")
    network::set.vertex.attribute(net, "group", group)
    # net %v% "group" <- group
    # vertex names
    network::network.vertex.names(net) <- rownames(meta)
    svglite::svglite(file = paste0(cdir, "/" , "network_map.svg"))
    #normalized mutual information close to 1, good, close to zero not similar
    #plot works but is not saving correctly
    GGally::ggnet2(net, label = rownames(meta), color = "group", label.size = 2, label.color = "yellow", edge.size = "weights") +
    ggplot2::ggtitle(paste0("NMI=", round(SNFtool::calNMI(group, truelabel), 2)))
    grDevices::dev.off()

  }
  if (int_method == "MOFA2"){

  }
}
