#' Multi-omics integration for toxicology data
#'
#' @param int_method sPLS-DA, MOFA2, WGCNA, or SNF
#'
#' @return output
#' @export
integrate_MO <- function(int_method = c("sPLS-DA", "MOFA", "WGCNA", "SNF", "iPCA")){

  #stop if metadata not provided
  if(is.null(meta)) stop("No metadata provided, run import function")

  #stop if data are not provided
  if(is.null(data_list)) stop("No data provided, run import function")

  cdir <- paste("FKAmoli", int_method, Sys.time(), sep = "_")
  cdir <- gsub(":", ".", cdir)
  dir.create(cdir)

  ##Integrate
  int_method <- match.arg(int_method)
  cat("Using",int_method,"in as integration method\n")
  X <- list()
  for (i in names(data_list[lengths(data_list) != 0])){
    X[[i]] <- as.data.frame(t(data_list[[i]]))
  }

  #mixOmics
  if(int_method == "sPLS-DA"){
    if(!is.null(X$rnaseq_counts)){
      if(ncol(X$rnaseq_counts) > 10000){
        X$rnaseq_counts <- X$rnaseq_counts[, order(apply(X$rnaseq_counts, 2, mad), decreasing = TRUE)[1:10000]]
        cat("Limiting rnaseq_counts to top 10,000 most variable features using MAD\n")
      }else{
        cat("Using all", ncol(X$rnaseq_counts), "rnaseq_counts features\n")
        }
    }
    if(!is.null(X$rrbs_counts)){
      if(ncol(X$rrbs_counts) > 10000){
        X$rrbs_counts <- X$rrbs_counts[, order(apply(X$rrbs_counts, 2, mad), decreasing = TRUE)[1:10000]]
        cat("Limiting rrbs_counts to top 10,000 most variable features using MAD\n")
      }else{
        cat("Using all", ncol(X$rrbs_counts), "rrbs_counts features\n")}
    }
    Y <- as.factor(meta$TRT)
    #doing "data driven" approach to setting design matrix
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
    diablo.tcga <- mixOmics::block.plsda(X, Y, ncomp = length(unique(meta$TRT)) + 1, design = design)
    perf.diablo.tcga <- mixOmics::perf(diablo.tcga, validation = "Mfold", folds = min(table(meta$TRT)), nrepeat = 10)
    ncomp <- min(perf.diablo.tcga$choice.ncomp$WeightedVote[2,])
    #dist_name <- colnames(perf.diablo.tcga$choice.ncomp$WeightedVote[,perf.diablo.tcga$choice.ncomp$WeightedVote[2,]==ncomp])[1]
    dist_name <- colnames(perf.diablo.tcga$error.rate[[1]])[which(perf.diablo.tcga$error.rate[[1]] == min(perf.diablo.tcga$error.rate[[1]]), arr.ind = TRUE)[[2]]]
    grDevices::pdf(file = paste0(cdir, "/", "mixOmics_classification_error_rate.pdf"))
    plot(perf.diablo.tcga)
    graphics::mtext(paste0("comp=", ncomp, " dist=", dist_name), side = 3)
    grDevices::dev.off()
    test.keepX <- list()
    for( i in names(X)){
      test.keepX[[i]] <- c(5:9, seq(10, 25, 5))
    }
    for (i in dist_name){
      tune.diablo.tcga <- mixOmics::tune.block.splsda(X, Y, ncomp = ncomp,
                                            test.keepX = test.keepX, design = design,
                                            validation = "Mfold", folds = 5, nrepeat = 1,
                                            # use two CPUs for faster computation
                                            #BPPARAM = BiocParallel::SnowParam(workers = 2),
                                            dist = i) # should this error rate match multiomics_pdf?
    }
    list.keepX <- tune.diablo.tcga$choice.keepX
    diablo.tcga <- mixOmics::block.splsda(X, Y, ncomp = ncomp,
                                          keepX = list.keepX)
    for ( i in 1:ncomp){
      values<-vector()
      for ( j in 1:length(X)){
        values<-rbind(values, as.data.frame(mixOmics::selectVar(diablo.tcga, comp = i)[[j]][2]))
      }
      utils::write.table(values, paste(cdir, "/", "splsda_variables_comp", i, ".csv"))
    }
    if(ncomp > 1){
      grDevices::pdf(file = paste0(cdir, "/", "mixOmics_plotVar.pdf"))
      mixOmics::plotVar(diablo.tcga, var.names = c(TRUE), style = "graphics", legend = TRUE,
                        title = "DIABLO comp 1 - 2")
      grDevices::dev.off()
      grDevices::pdf(file = paste0(cdir, "/", "mixOmics_plotIndiv.pdf"))
      mixOmics::plotIndiv(diablo.tcga, ind.names = FALSE, legend = TRUE, comp=c(1, 2),
                          title = "DIABLO comp 1 - 2", block = "weighted.average")
      grDevices::dev.off()
      for(i in names(X)){
        grDevices::pdf(file = paste0(cdir, "/", "mixOmics_auroc_", i, ".pdf"))
        auc.diablo.tcga <- mixOmics::auroc(diablo.tcga, roc.block = i, roc.comp = ncomp,
                                           print = FALSE)
        grDevices::dev.off()
      }
    }
  }

  #WGCNA - counts need to be vst or log2
  if(int_method == "WGCNA"){
    if(!is.null(X$rrbs_counts)){
      X$rrbs_counts <- X$rrbs_counts[, order(apply(X$rrbs_counts, 2, mad), decreasing = TRUE)[1:(.25 * ncol(X$rrbs_counts))]]
      cat("Limiting rrbs_counts to top", ncol(X$rrbs_counts), "most variable using MAD\n")}
    lowcount_omics_MEs <- list()
    block_MEs <- list()
    TRT_number <- meta$TRT
    for ( i in 1:length(unique(meta$TRT))){
      TRT_number <- gsub(unique(meta$TRT)[i], i, TRT_number)
    }
    TRT_number <- as.numeric(TRT_number)
    traitColors <- WGCNA::numbers2colors(TRT_number, signed = FALSE)
    #rnaseq counts are already logged
    for (i in names(X)){
      grDevices::pdf(file = paste0(cdir, "/", "WGCNA_hclust_", i, "_sampleTree.pdf"))
      sampleTree <- fastcluster::hclust(stats::dist(X[[i]]), method = "average")
      #cluster with metadata
      WGCNA::plotDendroAndColors(sampleTree, traitColors,
                                 groupLabels = colnames(meta)[2], #want this to say TRT, col depends on format
                                 main = paste("Sample dendrogram and trait heatmap", i)) #could remove? or could use to relate two omic layers
      grDevices::dev.off()
      #module detection one step - for low feature count omic layers
      if(length(X[[i]]) > 100 && length(X[[i]]) < 5000){ #if less than 100 features should not cluster?
        powers = c(c(1:10), seq(from = 12, to = 40, by = 2))
        #allowWGCNAThreads()
        sft <- WGCNA::pickSoftThreshold(X[[i]], powerVector = powers, verbose = 0)
        power_from_sft <- sft$fitIndices[sft$fitIndices$SFT.R.sq == max(sft$fitIndices$SFT.R.sq[1:10]), 1]
        color_forplot <- rep("black", length(sft$fitIndices$Power))
        color_forplot[power_from_sft] <- "red"
        grDevices::pdf(file = paste(cdir, "/", i, "soft threasholding"))
        graphics::par(mfrow = c(1, 2))
        plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
             xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n",
             main = paste("Scale independence"));
        graphics::text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
             labels = powers, cex = .9, col = color_forplot)
        plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
             xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
             main = paste("Mean connectivity"))
        graphics::text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = .9, col = color_forplot)
        grDevices::dev.off()
        adjacency <- WGCNA::adjacency(X$metab_peaks, power = power_from_sft)
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
        grDevices::pdf(file = paste0(cdir, "/", i, "module_eigengene_clustering_cutoff"))
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
        grDevices::pdf(file = paste0(cdir, "/", i, "_dendrogram"))
        WGCNA:: plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                                    c("Dynamic Tree Cut", "Merged dynamic"),
                                    dendroLabels = FALSE, hang = 0.03,
                                    addGuide = TRUE, guideHang = 0.05)
        grDevices::dev.off()
        lowcount_omics_MEs[[i]] <- WGCNA::moduleEigengenes(X[[i]], mergedColors)$eigengenes
      }else{
        #module detection one step - for high feature count omic layers
        #X[[i]]<-X[[i]][,order(apply(X[[i]],2,mad),decreasing=T)[1:10000]]#could remove
        # Choose a set of soft-thresholding powers
        powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
        # Call the network topology analysis function
        sft <- WGCNA::pickSoftThreshold(X[[i]], powerVector = powers, verbose = 0)
        power_from_sft <- sft$fitIndices[sft$fitIndices$SFT.R.sq == max(sft$fitIndices$SFT.R.sq[1:10]), 1]
        color_forplot <- rep("black", length(sft$fitIndices$Power))
        color_forplot[power_from_sft] <- "red"
        grDevices::pdf(file = paste0(cdir, "/", i, "_soft_threasholding.pdf"))
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
        graphics::text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = .9, col = color_forplot)
        grDevices::dev.off()
        #could do saveTOMd=FALSE
        bwnet <- WGCNA::blockwiseModules(X[[i]], maxBlockSize = 5000,
                                         power = power_from_sft, TOMType = "unsigned", minModuleSize = 30,
                                         reassignThreshold = 0, mergeCutHeight = 0.25,
                                         numericLabels = TRUE, saveTOMs = FALSE,
                                         saveTOMFileBase = paste0(i, "_TOM"),
                                         verbose = 0)

        bwLabels <- bwnet$colors
        bwModuleColors <- WGCNA::labels2colors(bwLabels)
        grDevices::pdf(file = paste0(cdir, "/", i, "_dendro.pdf"))
        if(length(bwnet$dendrograms) > 1){
          #par does not work
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
      }
    }
    #this all needs to be modified for having +2 omic layers
    all_MEs <- c(block_MEs, lowcount_omics_MEs)
    MEs_combo <- gtools::combinations(n = length(all_MEs), r = 2, v = 1:length(X))
    for (v in 1:nrow(MEs_combo)){
      nGenes <- ncol(X[[MEs_combo[v, 1]]])
      nSamples <- nrow(X[[MEs_combo[v, 1]]])
      moduleTraitCor <- WGCNA::cor(WGCNA::orderMEs(all_MEs[[MEs_combo[v, 1]]]), WGCNA::orderMEs(all_MEs[[MEs_combo[v, 2]]]), use = "pearson")
      moduleTraitPvalue <- WGCNA::corPvalueStudent(moduleTraitCor, nSamples)
      textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")")
      dim(textMatrix) <- dim(moduleTraitCor)
      pdf_name <- paste(names(X)[MEs_combo[v, 1]], names(X)[MEs_combo[v, 2]], "module_relationship.pdf", sep = "_")
      grDevices::pdf(file = paste0(cdir, "/" ,pdf_name))
      # Display the correlation values within a heatmap plot
      graphics::par(mar = c(4, 7, 2, 2))
      WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                            xLabels = names(WGCNA::orderMEs(all_MEs[[MEs_combo[v, 2]]])),
                            yLabels = names(WGCNA::orderMEs(all_MEs[[MEs_combo[v, 1]]])),
                            ySymbols = names(WGCNA::orderMEs(all_MEs[[MEs_combo[v, 1]]])),
                            colorLabels = FALSE,
                            colors = WGCNA::greenWhiteRed(50),
                            textMatrix = textMatrix,
                            setStdMargins = FALSE,
                            cex.text = 0.5,
                            zlim = c(-1, 1),
                            main = paste0(names(X)[MEs_combo[v, 1]], names(X)[MEs_combo[v, 2]],"Module-trait relationships"))
      grDevices::dev.off() #need to have more metabolites or not cluster metabolites, and need to reduce number of rna features
    }
  }

  #SNF https://github.com/cran/SNFtool
  if(int_method == "SNF"){
    # for (i in names(X)){
    #    X[[i]]<-as.data.frame(t(X[[i]]))
    #  }
    if(!is.null(X$rrbs_counts)){
      X$rrbs_counts <- X$rrbs_counts[, order(apply(X$rrbs_counts, 2, mad), decreasing = TRUE)[1:(.25 * ncol(X$rrbs_counts))]]
      cat("Limiting rrbs_counts to top", ncol(X$rrbs_counts), "most variable using MAD\n")}
    if(length(meta$sample) < 20){
      K <- length(meta$sample) - 1
    }else{K <- 20}
    #number of nearest neighbors, must be >1, usually 10-20
    alpha <- 0.5 ##hyperparameter, usually (0.3~0.8)
    t <- 20 ##Number of Iterations, usually (10~50)
    clusterNum <- length(unique(meta$TRT))
    W_temp <- list()
    for(i in names(X)){
      distance <- SNFtool::dist2(as.matrix((X[[i]])), as.matrix((X[[i]])))
      W_temp[[i]] <- SNFtool::affinityMatrix(distance, K, alpha)
    }
    W <- SNFtool::SNF(W_temp, K = K, t = t)
    group <- SNFtool::spectralClustering(W, clusterNum)
    TRT_number <- meta$TRT
    for ( i in 1:length(unique(meta$TRT))){
      j <- unique(group)[i]
      TRT_number <- gsub(unique(meta$TRT)[i], j, TRT_number)
    }
    grDevices::pdf(file = paste0(cdir, "/" , "SNF_cluster_heatmap.pdf"))
    SNFtool::displayClustersWithHeatmap(W, group, ColSideColors = cbind(as.character(TRT_number), as.character(group))) #hack this to be better https://rdrr.io/cran/SNFtool/src/R/displayClustersWithHeatmap.R
    graphics::legend("topleft", legend = unique(group), fill = unique(group), cex = .8)
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
    grDevices::pdf(file = paste0(cdir, "/" , "network_map.pdf"))
    #normalized mutual information close to 1, good, close to zero not similar
    #plot works but is not saving correctly
    GGally::ggnet2(net, label = rownames(meta), color = "group", label.size = 2, label.color = "yellow", edge.size = "weights") +
    ggplot2::ggtitle(paste0("NMI=", round(SNFtool::calNMI(group, truelabel), 2)))
    grDevices::dev.off()

  }
  if(int_method == "MOFA2"){

  }
}
