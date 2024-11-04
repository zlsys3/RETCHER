#' Run Subclone Clustering
#'
#' @description Run subclone clustering from a list of input data.
#'
#' @param inputList A list of data.frames, where each data.frame represents a
#' sample, and each sample has 8 columns of required data.
#' @param calcCCFMethod Methods for calculating Cancer Cell Fraction (CCF) are
#' supported ("distance", "threshold"). The difference between them lies in the
#' method of calculating mutation multiplicity, with the default being "distance".
#' @param maxCCFValue The maximum CCF value of SNV. SNV larger than this value
#' will be filtered. The default is 1.2.
#' @param sampleNames A vector of sample names, the length of which should be the
#' same as the number of samples, or the names of the elements of the "inputList"
#' if not specified.
#' @param minimumDepth Minimum sequencing depth below which SNVs will be filtered,
#' default is 20.
#' @param maximumClusters Maximum number of clusters, default is 10.
#' @param useSexChrs Boolean parameter to specify whether (TRUE) or not (FALSE)
#' mutations on the sex chromosomes are used in the clustering process. default is TRUE.
#' @param bootstrapIter The number of bootstrap iterations, default is 2000.
#' @param bootstrapConf Confidence interval for bootstrap, default is 0.95.
#' @param usePercentage Whether to multiply CCF by 100 to express percentage, default is FALSE.
#' @param LimitCCF Whether to limit CCF to 0 to 1, default is TRUE.
#'
#' @return Returns a list containing a data.frame of the clustering results,
#' a data.frame of the confidence intervals for each subclone cluster, a list
#' of CCF and mutation multiplicities for each SNV for each sample, and a
#' SciClone clustering results object.
#'
#' @export
#'
#' @importFrom utils capture.output
#' @importFrom reshape2 dcast
runTumorCluster <-
  function(inputList,
           calcCCFMethod = "distance",
           maxCCFValue = 1.2,
           sampleNames = NULL,
           minimumDepth = 20,
           maximumClusters = 10,
           useSexChrs = TRUE,
           bootstrapIter = 2000,
           bootstrapConf = 0.95,
           usePercentage = FALSE,
           LimitCCF = TRUE) {
    if (is.null(sampleNames)) {
      sampleNames <- names(inputList)
    }

    # Filter the snv of ccf > maxCCFValue and convert it to the input format of sciclone
    x <- convert2SciCloneInput(
      inputList,
      calcCCFMethod = calcCCFMethod,
      maxCCFValue = maxCCFValue,
      LimitCCF = LimitCCF
    )

    all_snv_list <- x$all_snv
    all_anno_list <- x$all_anno

    if (!usePercentage) {
      for (i in seq_along(all_snv_list)) {
        all_snv_list[[i]]$ccf <- all_snv_list[[i]]$ccf / 100
      }
    }

    if (nrow(all_snv_list[[1]]) < 10) {
      cat("\nToo few shared mutations for subclonal clustering.
          Next infer the sample tree.\n")
      return(list(mutMulti = x$inputList))
    }

    suppressMessages(library(sciClone))

    cat("\n**** Running SciClone ****\n")
    cat("\nA total of", nrow(all_snv_list[[1]]), "shared genes\n\n")
    start_time <- Sys.time()
    out <- capture.output(
      sc <- sciClone(
        vafs = all_snv_list,
        annotation = all_anno_list,
        sampleNames = sampleNames,
        minimumDepth = minimumDepth,
        maximumClusters = maximumClusters,
        useSexChrs = useSexChrs,
        clusterParams = NULL,
        clusterMethod = "bmm"
      )
    )
    end_time <- Sys.time()

    a <- gsub("\\[1\\] \\\"", "", out[length(out) - 1])
    print(gsub("\\\"", "", a))

    a <- gsub("\\[1\\] \\\"", "", out[length(out)])
    print(gsub("\\\"", "", a))

    cat("\nRun Time:", as.double(difftime(end_time, start_time, units = "secs")), "secs\n")

    # Extract relevant information from sc clustering results
    ccf_mat <- extractFromSC(sc = sc,usePercentage = usePercentage)

    # Calculate confidence intervals using bootstrap
    ci <- calcCI(
      ccfMat = ccf_mat,
      bootstrapIter = bootstrapIter,
      bootstrapConf = bootstrapConf
    )

    cat("The mean CCF for each cluster for each sample:\n")
    print(dcast(ci, cluster ~ sample, value.var = "mean"))

    return(list(
      mat = ccf_mat,
      ci = ci,
      mutMulti = x$inputList,
      sc = sc
    ))
  }


#' Convert to SciClone Input Format
#' @param inputList Same as 'inputList' parameter of runTumorCluster()
#' @param calcCCFMethod calcCCFMethod.
#' @param maxCCFValue maxCCFValue
#' @importFrom dplyr bind_rows filter distinct across
#' @importFrom tidyselect all_of
#' @return Returns a list containing the snv list and anno list needed for
#' sciclone clustering as well as the original input list with calculated CCF
#' and mutation multiplicity.
convert2SciCloneInput <-
  function(inputList,
           calcCCFMethod = "distance",
           maxCCFValue = 1.2,
           LimitCCF = TRUE) {
    if (any(is.na(unlist(lapply(inputList, function(x) {
      x[, 6]
    }))))) {
      stop("Major cannot be equal to NA")
    }
    if (any(unlist(lapply(inputList, function(x) {
      x[, 6]
    })) == 0)) {
      stop("Major cannot be equal to 0")
    }

    if (calcCCFMethod == "distance") {
      inputList <-
        lapply(inputList, function(x)
          calcCCFDistance(df = x, LimitCCF = LimitCCF))
    } else if (calcCCFMethod == "threshold") {
      inputList <-
        lapply(inputList, function(x)
          calcCCFThreshold(df = x, LimitCCF = LimitCCF))
    } else {
      cat(
        "\nParameter 'calcCCFMethod' only supports 'distance' and 'threshold',",
        "use the default parameter 'distance' now \n"
      )
      inputList <-
        lapply(inputList, function(x)
          calcCCFDistance(df = x, LimitCCF = LimitCCF))
    }

    # Display filtered snv
    filterSNV <- data.frame()
    for (i in seq_along(inputList)) {
      filterSNV <-
        bind_rows(filterSNV, filter(inputList[[i]], ccf > maxCCFValue))
    }
    if (nrow(filterSNV) == 0) {
      cat("No snv is filtered\n")
    } else {
      filterSNV <- filterSNV[!duplicated(filterSNV[, 2]), ]
      cat(nrow(filterSNV), "snv were filtered\n")
      print(filterSNV)
    }

    inputList <- lapply(inputList, function(df) {
      filter(df, ccf <= maxCCFValue)
    })

    # Intersection of gene and pos for all samples
    gene_pos <- Reduce(intersect, lapply(inputList, function(df) {
      paste(df[, 3], df[, 2], sep = "_")
    }))

    all_inputList <- lapply(inputList, function(df) {
      df[paste(df[, 3], df[, 2], sep = "_") %in% gene_pos, ]
    })

    all_inputList <- lapply(all_inputList, function(df) {
      distinct(df, across(all_of(names(df)[1:3])), .keep_all = TRUE)
    })

    # List of snv data.frame for 'chr','pos','ref','alt','ccf'
    all_snv_list <- lapply(all_inputList, function(df) {
      df$ccf <- round(df$ccf * 100, 2)
      df <- subset(df, select = c(1, 2, 4, 5, ccf))
      return(df)
    })

    # List of anno data.frame for 'chr','pos','gene'
    all_anno_list <- lapply(all_inputList, function(df) {
      df <- subset(df, select = c(1, 2, 3))
      return(df)
    })

    return(list(
      all_snv = all_snv_list,
      all_anno = all_anno_list,
      inputList = inputList
    ))
  }


# Calculate mutation multiplicity based on distance and then calculate CCF
calcCCFDistance <- function(df, LimitCCF = TRUE) {
  # First calculate the ccf of all snvs with m=1
  x <-
    df[(df[, 6] == 1 &
          df[, 7] == 1) | (df[, 6] == 1 & df[, 7] == 0), ]
  x$ccf <-
    (x[, 5] / (x[, 5] + x[, 4]) * ((x[, 6] + x[, 7]) * x[, 8] + 2 * (1 - x[, 8]))) / x[, 8]

  if (length(x$ccf) == 0) {
    cat(
      "Unable to use 'distance' to calculate CCF, will automatically use the 'threshold' to continue calculation\n"
    )
    df <- calcCCFThreshold(df)
    return(df)
  }

  x$mutMulti <- 1
  x <- x[, c(2, 9, 10)]
  df <- merge(df, x, by = names(x)[1], all.x = TRUE)
  df <- df[, c(2, 1, 3:ncol(df))]

  # For m there are multiple choices, choose the m that minimizes the distance
  # between the snv and its closest n snvs
  if (nrow(df) <= 200) {
    n <- 5
  } else {
    n <- round(nrow(df) * 0.025)
  }

  if (table(!is.na(df$ccf))["TRUE"] < n) {
    cat(
      "Unable to use 'distance' to calculate CCF, will automatically use the 'threshold' to continue calculation\n"
    )
    df <- calcCCFThreshold(df)
    return(df)
  }

  for (i in which(is.na(df$ccf))) {
    min.distance <- Inf
    mutMulti <- 1
    for (m in 1:df[, 6][i]) {
      tmp.ccf <-
        (df[, 5][i] / (df[, 5][i] + df[, 4][i]) *
           ((df[, 6][i] + df[, 7][i]) *
              df[, 8][i] + 2 * (1 - df[, 8][i]))) / df[, 8][i] / m

      df$ccf[i] <- tmp.ccf
      distance <-
        sum(sort(abs(df$ccf[i] - df$ccf[!is.na(df$ccf)]))[1:(n + 1)])
      if (distance < min.distance) {
        min.distance <- distance
        mutMulti <- m
      }
    }
    df$ccf[i] <- (df[, 5][i] / (df[, 5][i] + df[, 4][i]) *
                    ((df[, 6][i] + df[, 7][i]) *
                       df[, 8][i] + 2 *
                       (1 - df[, 8][i]))) / df[, 8][i] / mutMulti
    df$mutMulti[i] <- mutMulti
  }

  if (LimitCCF) {
    df$ccf <- ifelse(df$ccf > 1, 1, df$ccf)
  }

  return(df)
}


# Calculate mutation multiplicity based on threshold and then calculate CCF
calcCCFThreshold <- function(df, LimitCCF = TRUE) {
  m <- pmax(1, round(((df[, 5] / (
    df[, 5] + df[, 4]
  )) *
    (((df[, 6] + df[, 7]) * df[, 8]) + 2 * (1 - df[, 8])
    )) / df[, 8]))
  ccf <- ((df[, 5] / (df[, 5] + df[, 4])) *
            (((df[, 6] + df[, 7]) * df[, 8]) + 2 * (1 - df[, 8]))) / (df[, 8] * m)

  df$ccf <- ccf
  df$mutMulti <- m
  df <- df[order(df[, 2]), ]

  if (LimitCCF) {
    df$ccf <- ifelse(df$ccf > 1, 1, df$ccf)
  }

  return(df)
}


#' Extract the Required Columns From the SC
#' @param sc cluster result of sciclone.
#' @param usePercentage Whether to multiply CCF by 100 to express percentage, default is FALSE.
#' @importFrom dplyr select
#' @importFrom stats na.omit
extractFromSC <- function(sc, usePercentage = FALSE) {
  ccf_mat <- sc@vafs.merged |>
    select(
      chr,
      st,
      names(sc@vafs.merged)[length((sc@vafs.merged))],
      names(sc@vafs.merged)[grepl("*.var", names(sc@vafs.merged))],
      names(sc@vafs.merged)[grepl("*.ref", names(sc@vafs.merged))],
      names(sc@vafs.merged)[grepl("*.vaf", names(sc@vafs.merged))],
      cluster
    )

  colnames(ccf_mat) <- sub("vaf", "ccf", colnames(ccf_mat))

  if (!usePercentage) {
    ccf_mat[, c(names(ccf_mat)[grepl("*.ccf", names(ccf_mat))])] <-
      ccf_mat[, c(names(ccf_mat)[grepl("*.ccf", names(ccf_mat))])]
  } else {
    ccf_mat[, c(names(ccf_mat)[grepl("*.ccf", names(ccf_mat))])] <-
      ccf_mat[, c(names(ccf_mat)[grepl("*.ccf", names(ccf_mat))])] / 100
  }

  ccf_mat[, c(names(ccf_mat)[grepl("*.ref", names(ccf_mat))])] <-
    ccf_mat[, c(names(ccf_mat)[grepl("*.ref", names(ccf_mat))])] + ccf_mat[, c(names(ccf_mat)[grepl("*.var", names(ccf_mat))])]

  colnames(ccf_mat) <- sub("ref", "depth", colnames(ccf_mat))

  colnames(ccf_mat) <-
    c("chr", "pos", "gene", colnames(ccf_mat)[-c(1, 2, 3)])

  ccf_mat <- na.omit(ccf_mat)

  ccf_mat <- ccf_mat[ccf_mat$cluster != 0, ]

  return(ccf_mat)
}


#' Calculate The CI of The CCF of The Cluster
#' @param ccfMat ccf matrix.
#' @param bootstrapIter bootstrapIter.
#' @param bootstrapConf bootstrapConf.
#' @importFrom boot boot boot.ci
calcCI <- function(ccfMat,
                   bootstrapIter = 2000,
                   bootstrapConf = 0.95) {
  intervals <-
    data.frame(
      cluster = numeric(),
      sample = character(),
      lower = numeric(),
      mean = numeric(),
      upper = numeric()
    )

  for (i in 1:max(ccfMat$cluster)) {
    for (j in names(ccfMat)[grepl("*.ccf", names(ccfMat))]) {
      df <- ccfMat[ccfMat$cluster == i, j]

      if (all(df == df[1]) | any(is.na(df))) {
        ci_lower <- df[1]
        mean <- df[1]
        ci_upper <- df[1]
      } else {
        x <- boot(df, mean_fun, R = bootstrapIter)
        ci_lower <-
          boot.ci(x,
                  bootstrapConf = bootstrapConf,
                  type = "perc")$percent[4]
        mean <- mean(df)
        ci_upper <-
          boot.ci(x,
                  bootstrapConf = bootstrapConf,
                  type = "perc")$percent[5]
      }

      intervals <-
        rbind(
          intervals,
          data.frame(
            cluster = i,
            sample = j,
            lower = ci_lower,
            mean = mean,
            upper = ci_upper
          )
        )
    }
  }
  return(intervals)
}


# bootstrap mean function
mean_fun <- function(data, indices) {
  df <- data[indices]
  return(mean(df))
}


#' Clustering Post-Processing
#'
#' @description Post-processing of the clustering results included
#' removing outlier points; removing points with too few mutations in the cluster;
#' merging clusters that were not statistically different;
#' and merging irrational subclone clusters into clonal clusters.
#'
#' @param clusterResult A list returned by runTumorCluster().
#' @param removeOutliers Boolean value, if TRUE then outliers are removed based.
#' on quartile distance, if FALSE then they are not removed.
#' @param minSnvNum Minimum number of SNVs in the cluster, clusters with SNVs
#' less than this value will be deleted, default is 5, if you don't want to
#' delete the cluster you can set it to 0.
#' @param mergeUndiffClust Boolean value, if TRUE then clusters that are not
#' statistically different are merged according to the wilcox test,
#' if FALSE then they are not merged.
#' @param mergePvalue The p-value for the wilcox test, valid only if
#' mergeUndiffClust=TRUE. Default is 0.05.
#' @param bootstrapIter The number of bootstrap iterations, default is 2000.
#' @param bootstrapConf Confidence interval for bootstrap, default is 0.95.
#'
#' @return Returns a list containing one more data.frame, used to infer the
#' evolutionary tree, than the input list (clusterResult). The clustering
#' results and confidence intervals in the list have all been processed.
#'
#' @export
#'
#' @importFrom dplyr group_by filter n select
#' @importFrom reshape2 dcast
postProcess <- function(clusterResult,
                        removeOutliers = TRUE,
                        minSnvNum = 5,
                        mergeUndiffClust = TRUE,
                        mergePvalue = 0.05,
                        bootstrapIter = 2000,
                        bootstrapConf = 0.95) {
  removedCluster <- clusterResult

  # Remove outliers
  if (removeOutliers) {
    clusterResult <- removeOut(clusterResult)
  }

  # Delete clusters with a number of snvs less than minSnvNum.
  pre_remove_cluster <- sort(unique(clusterResult$mat$cluster))

  clusterResult$mat <- as.data.frame(clusterResult$mat %>%
                                       group_by(cluster) %>%
                                       filter(n() >= minSnvNum))

  after_remove_cluster <- sort(unique(clusterResult$mat$cluster))

  removed_cluster <-
    setdiff(pre_remove_cluster, after_remove_cluster)

  if (length(removed_cluster) == 0) {
    cat("\nAll cluster's snv number >",
        minSnvNum,
        ", no clusters were deleted\n")
  } else {
    cat(
      "\nA total of",
      length(removed_cluster),
      "clusters's snv number <",
      minSnvNum,
      ", the deleted clusters are:",
      removed_cluster,
      "\n"
    )
    print(removedCluster$mat[removedCluster$mat$cluster == removed_cluster, ])
  }

  # Merge clusters that have no statistical difference with Wilcox test
  if (mergeUndiffClust) {
    clusterResult <- mergeCluster(clusterResult, pvalue = mergePvalue)
  }

  # Regenerate the cluster number
  clusterResult$mat <- reClusterID(clusterResult$mat)

  # If there is a subclonal cluster whose ccf is greater than the ccf of the
  # clonal cluster, that subclonal cluster is merged into the clonal cluster.
  clusterResult <- checkCloneCCF(clusterResult)

  # Recalculate CI
  ci <-
    calcCI(clusterResult$mat,
           bootstrapIter = bootstrapIter,
           bootstrapConf = bootstrapConf)
  clusterResult$ci <- ci

  cat("The mean CCF of each cluster for each sample after post-processing:\n")
  print(dcast(ci, cluster ~ sample, value.var = "mean"))

  # Extracting the input data needed to infer evolutionary trees
  inferTree_input <- clusterResult$mat |>
    select(chr,
           pos,
           gene,
           names(clusterResult$mat)[grepl("*.ccf", names(clusterResult$mat))],
           cluster)

  clusterResult$inferTree <- inferTree_input

  return(clusterResult)
}


#' Remove Outliers
#' @param clusterResult clusterResult.
#' @importFrom stats quantile
#' @importFrom dplyr bind_rows
removeOut <- function(clusterResult) {
  cat("\n")
  snames <-
    names(clusterResult$mat[, grepl("*.ccf", names(clusterResult$mat))])
  removedf <- data.frame()
  # Calculate the interquartile distance of each cluster for each sample
  for (name in snames) {
    interval <-
      data.frame(cluster = numeric(),
                 lower = double(),
                 upper = double())

    for (j in 1:max(clusterResult$mat$cluster)) {
      qnt <-
        quantile(clusterResult$mat[clusterResult$mat$cluster == j, name], c(0.25, 0.75))

      h <- as.numeric(qnt[2] - qnt[1]) * 1.5

      interval <-
        rbind(interval,
              data.frame(
                cluster = j,
                lower = qnt[1] - h,
                upper = qnt[2] + h
              ))
    }

    rownames(interval) <-
      1:max(clusterResult$mat$cluster)

    # Find the coordinates of the outliers
    clusterResult$mat <-
      merge(clusterResult$mat, interval, by = "cluster")

    outlier_indices <-
      which(
        clusterResult$mat[, name] < clusterResult$mat$lower |
          clusterResult$mat[, name] > clusterResult$mat$upper
      )

    # Save deleted snv
    removedf <-
      bind_rows(removedf, clusterResult$mat[outlier_indices, ])
    removedf <- removedf[, -c(ncol(removedf) - 1, ncol(removedf))]

    if (length(outlier_indices) != 0) {
      clusterResult$mat <-
        clusterResult$mat[-outlier_indices, ]
    }

    clusterResult$mat <-
      clusterResult$mat[, -c(ncol(clusterResult$mat) -
                               1, ncol(clusterResult$mat))]

    clusterResult$mat <-
      clusterResult$mat[, c(2:ncol(clusterResult$mat), 1)]
  }

  if (nrow(removedf) == 0) {
    cat("No outliers\n")
  } else {
    cat(nrow(removedf), "outliers have been removed\n")
    print(removedf)
  }

  return(clusterResult)
}


#' Merge Cluster
#' @param clusterResult clusterResult.
#' @param pvalue pvalue.
#' @importFrom stats wilcox.test
mergeCluster <- function(clusterResult, pvalue = 0.05) {
  cat("\n")
  label <- 0
  name <-
    names(clusterResult$mat[, grepl("*.ccf", names(clusterResult$mat))])
  for (i in unique(clusterResult$mat$cluster)) {
    for (j in unique(clusterResult$mat$cluster)) {
      if (j > i) {
        if (length(unlist(clusterResult$mat[clusterResult$mat$cluster == i, name])) != 0 &
            length(unlist(clusterResult$mat[clusterResult$mat$cluster == j, name])) != 0) {
          wilcox <-
            wilcox.test(unlist(clusterResult$mat[clusterResult$mat$cluster == i, c(name)]),
                        unlist(clusterResult$mat[clusterResult$mat$cluster == j, c(name)]),
                        exact = F)

          if (wilcox$p.value >= pvalue) {
            clusterResult$mat[clusterResult$mat$cluster == j, ]$cluster <-
              i
            cat(
              "Cluster",
              i,
              "and Cluster",
              j,
              "have been merged, they are not statistical difference\n"
            )
            label <- label + 1
          }
        }
      }
    }
  }

  if (label == 0) {
    cat("All clusters are statistically different\n")
  }

  return(clusterResult)
}


#' Checking for Clone Cluster Problems
#'
#' @description Check for clonal cluster problems sample by sample, merging
#' clusters where the ccf of subclonal clusters is greater than the ccf of
#' clonal clusters
#' @param clusterResult clusterResult.
#' @importFrom reshape2 dcast
checkCloneCCF <- function(clusterResult) {
  cat("\n")
  while (T) {
    # Calculate the average ccf for each cluster for each sample
    sample_mean <-
      dcast(calcCI(clusterResult$mat, 100, 0.95),
            cluster ~ sample,
            value.var = "mean")

    c1 <- as.data.frame(sample_mean[sample_mean$cluster == 1, -1])
    lable <- 0

    for (i in 2:ncol(sample_mean)) {
      larger <- which(sample_mean[, i] > c1[[1, i - 1]])
      if (length(larger) > 0) {
        clusterResult$mat[clusterResult$mat$cluster == larger, ]$cluster <-
          1
        # Regenerate the cluster number
        clusterResult$mat <- reClusterID(clusterResult$mat)
        cat(
          "CCF of cluster",
          larger,
          "of sample",
          i - 1,
          "is greater than CCF of cluster 1, so cluster",
          larger,
          "is merged into cluster 1\n"
        )
        break
      } else {
        lable <- lable + 1
      }
    }

    if (lable == ncol(sample_mean) - 1) {
      break
    }
  }
  return(clusterResult)
}


# Regenerate the Cluster Number
reClusterID <- function(df) {
  oldClusterID <- sort(unique(df$cluster))
  newClusterID <- 1:length(sort(unique(df$cluster)))
  mapping <-
    data.frame(oldClusterID = oldClusterID, newClusterID = newClusterID)
  df$cluster <-
    mapping$newClusterID[match(df$cluster, mapping$oldClusterID)]
  return(df)
}



#' Run Samples Grouping
#'
#' @description Group SNVs according to their presence or absence in the sample
#' and remove groups with too few mutations in the group.
#'
#' @param ccfInputList The 'mutMulti' list in the list returned by the runTumorClutter()
#' or postProcess() function, both functions return the same 'mutMulti'.
#' @param sampleName A vector of sample names, the length of which should be the
#' same as the number of samples, or the names of the elements of the "ccfInputList"
#' if not specified.
#' @param minPresenceCCF Minimum CCF in the sample to mark SNV presence,
#' below which SNV will be marked as absent, default 0.05.
#' @param minSnvNum Minimum number of SNVs in the group, groups with SNVs less
#' than this value will be deleted, default is 5, if you don't want to delete
#' the group you can set it to 0.
#'
#' @return Returns a list containing a grouped data.frame and an average CCF
#' data.frame for each group.
#'
#' @export
#'
#' @importFrom dplyr rowwise filter group_by mutate across summarise ungroup n
#' @importFrom tidyselect all_of
runSampleGroup <-
  function(ccfInputList,
           sampleName = NULL,
           minPresenceCCF = 0.05,
           minSnvNum = 5) {
    if (class(ccfInputList) != "list" | length(ccfInputList) == 1) {
      cat("\nOnly one sample, cannot infer tree\n")
      return(NULL)
    }

    if (is.null(sampleName)) {
      sampleName <- paste0(names(ccfInputList), ".ccf")
    } else if (length(sampleName) != length(ccfInputList)) {
      stop(
        "\nThe number of sampleName (",
        length(sampleName),
        ") does not equal the number of ccfInputList samples (",
        length(ccfInputList),
        ")\n"
      )
    } else {
      sampleNacme <- paste0(sampleName, ".ccf")
    }

    # Extract ccf column
    ccfInputList <- lapply(ccfInputList, function(x) {
      x[, c(1:3, 9)]
    })

    ccfInputList <- lapply(seq_along(ccfInputList), function(i) {
      df <- ccfInputList[[i]]
      names(df)[4] <- paste0(names(df)[4], ".", i)
      df
    })

    # Multiple sample ccfs are organized into a data.frame
    mergeName <- names(ccfInputList[[1]])[1:3]
    group_df <- Reduce(function(x, y) {
      merge(x, y, by = mergeName, all = T)
    }, ccfInputList)
    colnames(group_df) <- c(mergeName, sampleName)
    group_df[is.na(group_df)] <- 0

    # Add binary grouping tags
    group_df <- group_df %>%
      rowwise() %>%
      mutate(label = paste0(ifelse(
        across(all_of(sampleName), ~ .x < minPresenceCCF), 0, 1
      ), collapse = "")) %>%
      ungroup()

    remove_index <- grepl("^0+$", group_df$label)
    if (any(remove_index)) {
      cat(
        "\nCCF is less than minPresenceCCF (",
        minPresenceCCF,
        "), the following SNVs are marked as absent\n"
      )
      print(as.data.frame(group_df[remove_index, ]))
      group_df <- group_df[!remove_index, ]
    }

    # Delete groups with fewer mutations than minSnvNum (but not root and leaf nodes)
    pre_remove_label <- unique(group_df$label)

    root <-
      pre_remove_label[sapply(pre_remove_label, function(x) {
        all(strsplit(x, "")[[1]] == "1")
      })]

    leaf <-
      pre_remove_label[sapply(pre_remove_label, function(x) {
        sum(strsplit(x, "")[[1]] == "1") == 1
      })]

    no_filter_label <- group_df[group_df$label %in% c(root, leaf), ]

    group_df <- group_df[!group_df$label %in% c(root, leaf), ]

    group_df <- as.data.frame(group_df %>%
                                group_by(label) |>
                                filter(n() > minSnvNum))

    group_df <- rbind(no_filter_label, group_df)

    after_remove_label <- unique(group_df$label)

    removed_label <-
      setdiff(pre_remove_label, after_remove_label)

    if (length(removed_label) != 0) {
      cat(
        "\nGroup",
        removed_label,
        "has been removed because its snv number is less than",
        minSnvNum,
        "\n"
      )
    }

    group_df$mean_ccf <-
      rowSums(group_df[, -c(1:3, ncol(group_df))]) / length(sampleName)
    mean_ccf <- group_df %>%
      group_by(label) %>%
      summarise(meanCCF = mean(mean_ccf))
    group_df <- group_df[, -ncol(group_df)]

    return(list(group = group_df, meanCCF = mean_ccf))
  }


#' Saving the Results of Tumor Subclone Clustering
#'
#' @description Save subclone clustering results, confidence intervals for
#' subclone clusters, mutation multiplicity and CCF values for each sample.
#'
#' @param clusterResult A list returned by runTumorCluster() or postProcess().
#' @param saveDir Save directory, automatically save all files to this directory,
#' and automatically create it if it does not exist.
#'
#' @export
#'
#' @importFrom utils write.table
saveTumorClusterRes <- function(clusterResult, saveDir) {
  if (!dir.exists(saveDir)) {
    dir.create(saveDir, recursive = T)
  }

  write.table(
    clusterResult$mat,
    paste0(saveDir, "/cluster_res.tsv"),
    sep = "\t",
    row.names = F
  )

  write.table(
    clusterResult$ci,
    paste0(saveDir, "/cluster_ci.tsv"),
    sep = "\t",
    row.names = F
  )

  for (i in names(clusterResult$mutMulti)) {
    df <- clusterResult$mutMulti[i]
    write.table(
      df,
      paste0(saveDir, "/", i, "_mutMulti.tsv"),
      sep = "\t",
      row.names = F,
      col.names = colnames(df[[1]])
    )
  }
}


#' Saving the Results of Samples Grouping
#'
#' @description Save sample grouping data.frame and average CCF per group data.frame.
#'
#' @param sampleGroup A list returned by runSampleGroup().
#' @param saveDir Save directory, automatically save all files to this directory,
#' and automatically create it if it does not exist.
#'
#' @export
#'
#' @importFrom utils write.table
saveSampleGroupRes <- function(sampleGroup, saveDir) {
  if (!dir.exists(saveDir)) {
    dir.create(saveDir, recursive = T)
  }

  write.table(
    sampleGroup$group,
    paste0(saveDir, "/sampleGroup.tsv"),
    sep = "\t",
    quote = T,
    row.names = F
  )

  write.table(
    sampleGroup$mean_ccf,
    paste0(saveDir, "/sampleGroup_meanCCF.tsv"),
    sep = "\t",
    quote = T,
    row.names = F
  )
}
