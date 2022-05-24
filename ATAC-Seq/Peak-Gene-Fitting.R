getDfFit <-
  function (vecGeneIDs = NULL,
            scaNTopIDs = NULL,
            objectImpulseDE2,
            boolCaseCtrl,
            dirOut = NULL,
            strFileName = "ImpulseDE2_Trajectories.pdf",
            boolMultiplePlotsPerPage = TRUE,
            boolSimplePlot = FALSE,
            vecRefPval = NULL,
            strNameRefMethod = NULL)
  {
    library(ImpulseDE2)
    dfAnnot <- get_dfAnnotationProc(obj = objectImpulseDE2)
    scaNPlotsPerPage <- 4
    cbPalette <- c(
      "#000000",
      "#E69F00",
      "#56B4E9",
      "#009E73",
      "#F0E442",
      "#0072B2",
      "#D55E00",
      "#CC79A7"
    )
    if (is.null(get_lsModelFits(objectImpulseDE2))) {
      error(
        paste0(
          "objectImpulseDE2 does not contain model fits. ",
          "Run ImpulseDE2_main first."
        )
      )
    }
    if (is.null(vecGeneIDs) & is.null(scaNTopIDs)) {
      stop("Supply either vecGeneIDs or scaNTopIDs.")
    }
    if (!is.null(vecGeneIDs) & !is.null(scaNTopIDs)) {
      stop("Only one of the two: vecGeneIDs or scaNTopIDs.")
    }
    if (is.null(get_lsModelFits(objectImpulseDE2)$IdxGroups$case$lsvecBatchUnique) &
        !is.null(boolSimplePlot)) {
      #print("Setting boolSimplePlot=TRUE as no batch structure was found.")
      boolSimplePlot <- TRUE
    }
    if (!is.null(dirOut) && !file.exists(dirOut)) {
      stop("Output directory dirOut not available.")
    }
    if (!is.null(vecRefPval) && (names(vecRefPval) != vecGeneIDs)) {
      stop("Names of vecRefPval have to be IDs from vecGeneIDs.")
    }
    if (!boolSimplePlot & boolCaseCtrl & boolMultiplePlotsPerPage) {
      warning(
        paste0(
          "Plots are likely overloaded. ",
          "Consider switching to boolSimplePlot=TRUE ",
          "or boolMultiplePlotsPerPage=FALSE."
        )
      )
    }
    if (is.null(vecGeneIDs))
      vecGeneIDs <-
      objectImpulseDE2$dfImpulseDE2Results[with(objectImpulseDE2$dfImpulseDE2Results,
                                                order(padj)),]$Gene[1:scaNTopIDs]
    scaNIDs <- length(vecGeneIDs)
    lsgplotsID <- list()
    for (id in vecGeneIDs) {
      vecTimePointsFit <-
        seq(
          min(get_dfAnnotationProc(obj = objectImpulseDE2)$Time),
          max(get_dfAnnotationProc(obj = objectImpulseDE2)$Time),
          length.out = 100
        )
      vecCaseImpulseParam <-
        get_lsModelFits(obj = objectImpulseDE2)$case[[id]]$lsImpulseFit$vecImpulseParam
      vecCaseImpulseValue <-
        ImpulseDE2:::evalImpulse_comp(vecImpulseParam = vecCaseImpulseParam,
                         vecTimepoints = vecTimePointsFit)
      if (!is.null(get_vecConfounders(obj = objectImpulseDE2))) {
        vecCaseBatchFactors <-
          get_lsModelFits(obj = objectImpulseDE2)$case[[id]]$lsImpulseFit$lsvecBatchFactors[[1]]
        vecBatchLabelsCase <-
          get_lsModelFits(obj = objectImpulseDE2)$IdxGroups$case$lsvecBatchUnique[[1]]
      }
      else {
        vecCaseBatchFactors <- 1
        vecBatchLabelsCase <- " "
      }
      vecValueToPlotCase <- do.call(c,
                                    lapply(vecCaseBatchFactors,
                                           function(f)
                                             vecCaseImpulseValue * f))
      if (boolCaseCtrl) {
        vecSamples <- dfAnnot$Sample
        dfRaw <-
          data.frame(
            normCounts = get_matCountDataProc(obj = objectImpulseDE2)[id,
                                                                      vecSamples] /
              get_vecSizeFactors(objectImpulseDE2)[vecSamples],
            time = dfAnnot[vecSamples,]$Time,
            Batch = dfAnnot[vecSamples,
                            get_vecConfounders(obj = objectImpulseDE2)[1]],
            Condition = dfAnnot[vecSamples,]$Condition
          )
        vecControlImpulseParam <-
          get_lsModelFits(obj = objectImpulseDE2)$control[[id]]$lsImpulseFit$vecImpulseParam
        vecControlImpulseValue <-
          ImpulseDE2:::evalImpulse_comp(vecImpulseParam = vecControlImpulseParam,
                           vecTimepoints = vecTimePointsFit)
        vecCombinedImpulseParam <-
          get_lsModelFits(obj = objectImpulseDE2)$combined[[id]]$lsImpulseFit$vecImpulseParam
        vecCombinedImpulseValue <-
          ImpulseDE2:::evalImpulse_comp(vecImpulseParam = vecCombinedImpulseParam,
                           vecTimepoints = vecTimePointsFit)
        vecControlBatchFactors <-
          get_lsModelFits(obj = objectImpulseDE2)$control[[id]]$lsImpulseFit$lsvecBatchFactors[[1]]
        vecCombinedBatchFactors <-
          get_lsModelFits(obj = objectImpulseDE2)$combined[[id]]$lsImpulseFit$lsvecBatchFactors[[1]]
        vecBatchLabelsCtrl <-
          get_lsModelFits(obj = objectImpulseDE2)$IdxGroups$control$lsvecBatchUnique[[1]]
        vecBatchLabelsComb <-
          get_lsModelFits(obj = objectImpulseDE2)$IdxGroups$combined$lsvecBatchUnique[[1]]
        vecValueToPlotCtrl <-
          do.call(c,
                  lapply(vecControlBatchFactors,
                         function(f)
                           vecControlImpulseValue * f))
        vecValueToPlotComb <-
          do.call(c,
                  lapply(vecCombinedBatchFactors,
                         function(f)
                           vecCombinedImpulseValue * f))
        if (boolSimplePlot) {
          dfFit <- data.frame(
            time =
              rep(vecTimePointsFit,
                  3),
            value = c(
              vecCaseImpulseValue,
              vecControlImpulseValue,
              vecCombinedImpulseValue
            ),
            Condition = c(
              rep("case",
                  length(vecTimePointsFit)),
              rep("control", length(vecTimePointsFit)),
              rep("combined", length(vecTimePointsFit))
            )
          )
          if (!is.null(vecBatchLabelsComb)) {
            gplotGene <- ggplot() + geom_point(data = dfRaw,
                                               aes(
                                                 x = time,
                                                 y = normCounts,
                                                 colour = Condition,
                                                 shape = Batch
                                               )) + geom_line(data = dfFit,
                                                              aes(
                                                                x = time,
                                                                y = value,
                                                                colour = Condition
                                                              ))
          }
          else {
            gplotGene <- ggplot() + geom_point(data = dfRaw,
                                               aes(
                                                 x = time,
                                                 y = normCounts,
                                                 colour = Condition
                                               )) +
              geom_line(data = dfFit, aes(
                x = time,
                y = value,
                colour = Condition
              ))
          }
        }
        else {
          dfFit <- data.frame(
            time = c(
              rep(vecTimePointsFit,
                  length(vecCaseBatchFactors)),
              rep(vecTimePointsFit,
                  length(vecControlBatchFactors)),
              rep(vecTimePointsFit,
                  length(vecCombinedBatchFactors))
            ),
            value = c(
              vecValueToPlotCase,
              vecValueToPlotCtrl,
              vecValueToPlotComb
            ),
            Condition = c(
              rep(rep(
                "case",
                length(vecTimePointsFit)
              ), length(vecCaseBatchFactors)),
              rep(
                rep("control", length(vecTimePointsFit)),
                length(vecControlBatchFactors)
              ),
              rep(
                rep("combined",
                    length(vecTimePointsFit)),
                length(vecCombinedBatchFactors)
              )
            ),
            BatchFit = c(
              sapply(vecBatchLabelsCase, function(label)
                rep(label,
                    length(
                      vecTimePointsFit
                    ))),
              sapply(vecBatchLabelsCtrl,
                     function(label)
                       rep(label, length(
                         vecTimePointsFit
                       ))),
              sapply(vecBatchLabelsComb, function(label)
                rep(label,
                    length(
                      vecTimePointsFit
                    )))
            )
          )
          gplotGene <- ggplot() + geom_point(data = dfRaw,
                                             aes(
                                               x = time,
                                               y = normCounts,
                                               colour = Condition,
                                               shape = Batch
                                             )) + geom_line(data = dfFit,
                                                            aes(
                                                              x = time,
                                                              y = value,
                                                              colour = Condition,
                                                              linetype = BatchFit
                                                            ))
        }
      }
      else {
        vecSamples <- dfAnnot[dfAnnot$Condition == "case",]$Sample
        if (!is.null(get_vecConfounders(obj = objectImpulseDE2))) {
          vecBatches <-
            dfAnnot[vecSamples, get_vecConfounders(obj = objectImpulseDE2)[1]]
        }
        else {
          vecBatches <- rep(" ", length(vecSamples))
        }
        dfRaw <-
          data.frame(
            normCounts = get_matCountDataProc(obj = objectImpulseDE2)[id,
                                                                      vecSamples] /
              get_vecSizeFactors(obj = objectImpulseDE2)[vecSamples],
            time = dfAnnot[vecSamples,]$Time,
            Batch = vecBatches,
            condition = dfAnnot[vecSamples,]$Condition
          )
        if (boolSimplePlot) {
          dfFit <- data.frame(
            time = vecTimePointsFit,
            value = vecCaseImpulseValue,
            Batch = rep(vecBatchLabelsCase[1],
                        length(vecTimePointsFit))
          )
          if (!is.null(vecBatchLabelsCase)) {
            gplotGene <- ggplot() + geom_point(data = dfRaw,
                                               aes(
                                                 x = time,
                                                 y = normCounts,
                                                 shape = Batch
                                               )) +
              geom_line(data = dfFit, aes(
                x = time,
                y = value,
                linetype = Batch
              ))
          }
          else {
            gplotGene <- ggplot() + geom_point(data = dfRaw,
                                               aes(x = time, y = normCounts)) + geom_line(data = dfFit,
                                                                                          aes(x = time, y = value))
          }
        }
        else {
          dfFit <- data.frame(
            time = rep(vecTimePointsFit,
                       length(vecCaseBatchFactors)),
            value = vecValueToPlotCase,
            condition = rep(rep(
              "case", length(vecTimePointsFit)
            ),
            length(vecCaseBatchFactors)),
            BatchFit = as.vector(sapply(vecBatchLabelsCase,
                                        function(label)
                                          rep(label, length(vecTimePointsFit))))
          )
          gplotGene <- ggplot() + geom_point(data = dfRaw,
                                             aes(
                                               x = time,
                                               y = normCounts,
                                               colour = Batch,
                                               shape = Batch
                                             )) + geom_line(data = dfFit,
                                                            aes(
                                                              x = time,
                                                              y = value,
                                                              colour = BatchFit,
                                                              linetype = BatchFit
                                                            ))
        }
      }
      gplotID <-
        gplotGene + scale_colour_manual(values = cbPalette) +
        xlab("time [hours]") + ylab("Read counts")
      if (!is.null(strNameRefMethod)) {
        gplotID <- gplotID + labs(title = paste0(
          id,
          ": ImpulseDE2 [",
          round(
            log(objectImpulseDE2$dfImpulseDE2Results[id,]$padj) / log(10),
            2
          ),
          "] ",
          strNameRefMethod,
          " [",
          round(log(vecRefPval[id]) / log(10), 2),
          "]"
        ))
      }
      else {
        gplotID <- gplotID + labs(title = paste0(id, ": log10 q = ",
                                                 round(
                                                   log(objectImpulseDE2$dfImpulseDE2Results[id,]$padj) / log(10),
                                                   2
                                                 )))
      }
      lsgplotsID[[length(lsgplotsID) + 1]] <- gplotID
    }
    if (!is.null(dirOut)) {
      dirFileOut <- paste0(dirOut, strFileName)
      print(paste0("Creating ", dirFileOut))
      graphics.off()
      if (boolMultiplePlotsPerPage) {
        pdf(dirFileOut)
        scaNPages <- scaNIDs %/% scaNPlotsPerPage
        if (scaNIDs %% scaNPlotsPerPage == 0)
          scaNPages <- scaNPages - 1
        for (p in seq(0, scaNPages)) {
          if (p < scaNIDs %/% scaNPlotsPerPage) {
            vecidxPlots <- seq((p * scaNPlotsPerPage +
                                  1), ((p + 1) * (scaNPlotsPerPage)))
          }
          else {
            vecidxPlots <- seq((p * scaNPlotsPerPage +
                                  1), scaNIDs)
          }
          print(
            plot_grid(
              plotlist = lsgplotsID[vecidxPlots],
              align = "h",
              nrow = scaNPlotsPerPage / 2,
              ncol = 2,
              rel_widths = c(1, 1),
              rel_heights = c(1, 1,
                              1)
            )
          )
        }
      }
      else {
        pdf(dirFileOut)
        for (p in seq(1, scaNIDs)) {
          print(lsgplotsID[[p]])
        }
      }
      dev.off()
      graphics.off()
    }
    return(dfFit)
  }

compareGeneToPeaks <-
  function(gene = "") {
    # Get the RNA-Seq fit
    gene_fit <-
      getDfFit(
        vecGeneIDs = gene,
        objectImpulseDE2 = impulse_obj,
        boolCaseCtrl = T,
        boolSimplePlot = T
      )
    gene_fit <- gene_fit[gene_fit$Condition == "case", 2]
    
    # TODO
    # We need to switch from a single peak of interest to peak-gene interactions. We now have the same PeakID with different genes.
    peaks_of_interest <-
      closest_5_genes[closest_5_genes$geneId %in% gene, "PeakID"]
    
    # Make sure it's in the ATAC-Seq obj.
    peaks_of_interest <-
      peaks_of_interest[peaks_of_interest %in% impulse_obj2@vecAllIDs]
    
    # If NONE of the peaks are in the object, return NA.
    
    if (!length(peaks_of_interest)) {
      return(data.frame(
        spearman = NA,
        distanceToTSS = NA,
        gene = gene,
        PeakID = NA
      ))
    }
    
    # Get the ATAC-Seq fit
    peak_res <- as.data.frame(matrix(nrow = length(peaks_of_interest), ncol = 4))
    colnames(peak_res) <- c("spearman","distanceToTSS", "gene", "PeakID")
    peak_count <- 1
    
    for (peak in peaks_of_interest) {
      peak_fit <-
        getDfFit(
          vecGeneIDs = peak,
          objectImpulseDE2 = impulse_obj2,
          boolCaseCtrl = F,
          boolSimplePlot = T
        )
      peak_fit <- peak_fit[,2]
      #spearman correlation
      peak_res[peak_count, 1] <- cor(gene_fit, peak_fit, method = "spearman")
      peak_res[peak_count, 2] <-
        closest_5_genes[closest_5_genes$PeakID == peak &
                          closest_5_genes$geneId == gene, "distanceToTSS"]
      peak_res[peak_count, 3] <- gene
      peak_res[peak_count, 4] <- peak
      
      peak_count <- peak_count +1
      
    }
    return(peak_res)
  }

RunMultipleGenes <- function(genelist = "") {
  n_missing_rna <-
    length(genelist[!genelist %in% impulse_obj@vecAllIDs])
  print(paste0("Removing ", n_missing_rna, ", genes not in RNA."))
  genelist <- genelist[genelist %in% impulse_obj@vecAllIDs]
  n_missing_peaks <-
    length(genelist[!genelist %in% closest_5_genes$geneId])
  print(paste0("Removing ", n_missing_peaks, " genes, no associated peaks."))
  genelist <- genelist[genelist %in% closest_5_genes$geneId]
  # Split into 5000 item chunks if necessary.
  if (length(genelist) > 5000) {
    message("More than 5,000 genes detected. Chunking...")
    chunks <- split(genelist, ceiling(seq_along(genelist)/5000))
    chunk_res <- list()
    for (chunk in chunks){
      chunk_res[chunk] <- bplapply(chunk, FUN = compareGeneToPeaks)
    }
    
    results.df <- do.call("rbind", chunk_res)
    results.df$Symbol <- lookup(results.df$gene, key.match = as.data.frame(GenesUniverse[, c(1, 4)]))
  }
  
  else{
    temp <- bplapply(genelist, FUN = compareGeneToPeaks)
    print("bplApply Done")
    results.df <- do.call("rbind", temp)
    results.df$Symbol <- lookup(results.df$gene, key.match = as.data.frame(GenesUniverse[, c(1, 4)]))
  }
  
  return(results.df)
}

get_conservation_2 <- function(x){
  split <- str_split(x, pattern = "-", simplify = T)
  if (split[1] %in% c("MT")){
    return(NA)
  }
  else{
    gr1 <- GRanges(seqnames = as(split[1], "Rle"), IRanges(start =seq(as.numeric(split[2]), as.numeric(split[3])), width=1))
    return(mean((subsetByOverlaps(x = Window25bp.gr, ranges = gr1))@elementMetadata$score))
  }
  
}